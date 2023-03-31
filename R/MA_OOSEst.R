## Trying Methods for OOS Estimation in Random Effects Meta-Analysis ##

library(tidyverse)
library(lme4)
library(multcomp)
library(MASS)
#library(pimeta)
#library(merTools)

source("R/MDD_Generation_OOSEst.R")

#simulate data ####
sim_dat <- gen_mdd(K=6, n_mean=200, n_sd=0, eps_study_m=1, eps_study_tau=0.05, 
                   distribution="same", target_dist="same")
train_dat <- sim_dat[["train_dat"]]
target_dat <- sim_dat[["target_dat"]]

covars <- c("sex", "smstat", "weight", "age", "madrs")


#model ####
mod <- lmer(Y ~ W + age + W:age +
              (W + age + W:age | S), data=train_dat)
sum <- summary(mod)

#explore var/covar
#coef(summary(mod)) #fixed effects and SEs
#fc <- vcov(mod) #covariance matrix of fixed effects
#print(VarCorr(mod), comp=c("Variance","Std.Dev")) #random effects
#rc <- as.data.frame(VarCorr(mod)) #variance and covariances of random effects
#sigma(mod) #residual sd

#example CATE PI ####
# example: let age be 40 years old
# x <- 40
# meantheta <- fixef(mod)["W"] + fixef(mod)["W:age"]*x
# vartheta <- fc["W","W"] + x^2*fc["W:age","W:age"] +
#   2*x*fc["W","W:age"] + rc[which(rc$var1=="W" & is.na(rc$var2)==T),"vcov"] +
#   x^2*rc[which(rc$var1=="W:age" & is.na(rc$var2)==T),"vcov"] +
#   2*x*rc[which(rc$var1=="W" & rc$var2=="W:age"),"vcov"]


#prediction interval by hand ####
#get variance components of model
matrix_var <- function(mod) {

  #calculate variances
  fc <- vcov(mod) #covariance matrix of fixed effects
  rc <- as.data.frame(VarCorr(mod)) #variance and covariances of random effects
  
  #matrix calculation
  beta <- matrix(fixef(mod)[c("W","W:age")], nrow=2) #beta-hat
  var_beta <- fc[c("W","W:age"),c("W","W:age")] #Var(beta-hat)
  rand <- rc[c(which(rc$var1=="W" & is.na(rc$var2)==T),
               which(rc$var1=="W" & rc$var2=="W:age"),
               which(rc$var1=="W" & rc$var2=="W:age"),
               which(rc$var1=="W:age" & is.na(rc$var2)==T)), "vcov"] #can I clean this up??
  var_rand <- matrix(rand, nrow=2, dimnames=list(c("W","W:age"),c("W","W:age"))) #Var(ranef) 
  
  return(list(beta=beta, var_beta=var_beta, var_rand=var_rand))
}

#add pis to dataset
manual_pi <- function(df, mod, K) {
  
  #get variances
  res <- matrix_var(mod)
  beta <- res$beta
  var_beta <- res$var_beta
  var_rand <- res$var_rand
  
  #calculate theta-hats
  X <- df %>%
    mutate(trt = 1) %>%
    dplyr::select(trt, age) %>%
    as.matrix() #X (same as Z in this model)
  mean_theta <- X %*% beta %>% c() #theta-hat
  vcov_theta <- X %*% var_beta %*% t(X) + X %*% var_rand %*% t(X)
  var_theta <- diag(vcov_theta) #Var(theta-hat)
  
  #prediction interval
  pred_lower <- mean_theta + qt(.025, K-2)*sqrt(var_theta)
  pred_upper <- mean_theta - qt(.025, K-2)*sqrt(var_theta)
  
  df <- df %>%
    mutate(lower = pred_lower,
           mean = mean_theta,
           upper = pred_upper)
  
  return(df)
}


#confidence interval by glht ####
#get ci based on x=age (as a character)
age_ci <- function(x, mod) {
  
  ci <- glht(mod, linfct=paste0("W +",x,"*W:age = 0"))
  mean <- coef(ci) %>% as.numeric()
  sd <- sqrt(vcov(ci)[1,1])
  lower <- mean - 1.96*sd %>% as.numeric()
  upper <- mean + 1.96*sd %>% as.numeric()
  
  return(c(lower=lower, mean=mean, upper=upper))
} 

#add cis to dataset
glht_ci <- function(df, mod) {
  
  ages <- as.character(df$age)
  cis <- map_dfr(.x=ages, .f=age_ci, mod=mod)
  df <- df %>%
    bind_cols(cis)
  
  return(df)
} 


#prediction interval by bootstrap ####
#create intervals from bootstrap
sample_cate <- function(x, boot_fix, boot_rand) {
  
  #get cate|x according to each coefficient
  cates <- (boot_fix$W + boot_rand$W) + 
    (boot_fix$W.age + boot_rand$W.age)*x
  
  #calculate interval over all iterations
  mean <- mean(cates)
  sd <- sd(cates)
  lower <- mean - 1.96*sd
  upper <- mean + 1.96*sd

  return(c(lower=lower, mean=mean, upper=upper))
}

#bootstrap and add intervals to data
boot_pi <- function(df, mod) {
  
  #get var-covar of fixed and random effects from model
  res <- matrix_var(mod) 
  
  #randomly sample fixed and random coefficients
  boot_fix <- mvrnorm(n=1000, mu=c(res$beta), Sigma=res$var_beta) %>%
    data.frame()
  boot_rand <- mvrnorm(n=1000, mu=c(0,0), Sigma=res$var_rand) %>%
    data.frame() #assume ranefs have mean 0
  
  #apply coefficients to estimate cate for all ages
  intervals <- map_dfr(.x=df$age, .f=sample_cate,
                             boot_fix=boot_fix, boot_rand=boot_rand)
  df <- df %>%
    bind_cols(intervals)
  
  return(df)
}


#check results for all methods ####
assess_interval <- function(train_dat, target_dat) {
  
  #calculate mse
  train_mse <- mean((train_dat$mean - train_dat$tau)^2)
  target_mse <- mean((target_dat$mean - target_dat$tau)^2)
  
  #calculate coverage
  train_coverage <- sum(train_dat$tau >= train_dat$lower & train_dat$tau <= train_dat$upper)/nrow(train_dat)
  target_coverage <- sum(target_dat$tau >= target_dat$lower & target_dat$tau <= target_dat$upper)/nrow(target_dat)
  
  #calculate length
  train_length <- mean(train_dat$upper - train_dat$lower)
  target_length <- mean(target_dat$upper - target_dat$lower)
  
  return(c(train_mse = train_mse, target_mse = target_mse,
           train_coverage = train_coverage, target_coverage = target_coverage,
           train_length = train_length, target_length = target_length))
  
}


#overall function ####
compare_oos <- function(N=100, K=6, n_mean=200, n_sd=0, eps_study_m=0.05, eps_study_tau=3, 
                        eps_study_age=0.05, distribution="same", target_dist="same", eps_target=0) {
  
  
  ## Simulate training and target (OOS) data
  sim_dat <- gen_mdd(K, n_mean, n_sd, eps_study_m, eps_study_tau, 
                     eps_study_age, distribution, target_dist, eps_target)
  train_dat <- sim_dat[["train_dat"]]
  target_dat <- sim_dat[["target_dat"]]
  
  covars <- c("sex", "smstat", "weight", "age", "madrs")
  
  
  ## Fit mixed effects model
  mod <- lmer(Y ~ age + madrs + sex + W + W:age +
                  (W + W:age | S), data=train_dat)
  sum <- summary(mod)
  
  
  ## Calculate mean and CIs for individual and assess accuracy
  #confidence interval
  glht_train <- glht_ci(train_dat, mod)
  glht_target <- glht_ci(target_dat, mod)
  glht_res <- assess_interval(glht_train, glht_target)
    
  #manual PI
  manual_train <- manual_pi(train_dat, mod)
  manual_target <- manual_pi(target_dat, mod)
  manual_res <- assess_interval(manual_train, manual_target)

  #bootstrap PI
  boot_train <- boot_pi(train_dat, mod)
  boot_target <- boot_pi(target_dat, mod)
  boot_res <- assess_interval(boot_train, boot_target)
  
  
  ## Save results
  return(list(sum=sum, glht_res=glht_res, manual_res=manual_res, boot_res=boot_res,
              N=N, K=K, n_mean=n_mean, n_sd=n_sd,
              eps_study_m=eps_study_m, eps_study_tau=eps_study_tau, 
              distribution=distribution, target_dist=target_dist))
}
