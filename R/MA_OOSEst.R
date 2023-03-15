## Trying Methods for OOS Estimation in Random Effects Meta-Analysis ##

library(tidyverse)
library(rsample)
library(grf)
library(fastDummies)
library(nnet)
library(lme4)
library(pimeta)
library(multcomp)
library(merTools)
library(marginaleffects)

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
} #get variance components of model

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
} #add pis to dataset


#confidence interval by glht ####
age_ci <- function(x, mod) {
  
  ci <- glht(mod, linfct=paste0("W +",x,"*W:age = 0"))
  mean <- coef(ci) %>% as.numeric()
  sd <- sqrt(vcov(ci)[1,1])
  lower <- mean - 1.96*sd %>% as.numeric()
  upper <- mean + 1.96*sd %>% as.numeric()
  
  return(c(lower=lower, mean=mean, upper=upper))
} #get ci based on x=age (as a character)

glht_ci <- function(df, mod) {
  
  ages <- as.character(df$age)
  cis <- map_dfr(.x=ages, .f=age_ci, mod=mod)
  df <- df %>%
    bind_cols(cis)
  
  return(df)
} #add cis to dataset


#prediction interval by bootstrap ####



#check results for all methods ####
assess_interval <- function(mod, train_dat, target_dat, method) {
  
  #add intervals to both datasets
  if (method == "glht") {
    train_dat <- glht_ci(train_dat, mod)
    target_dat <- glht_ci(target_dat, mod)
    
  } else if (method == "manual") {
    train_dat <- manual_pi(train_dat, mod)
    target_dat <- manual_pi(target_dat, mod)
    
  }
  
  #calculate mse
  train_mse <- mean((train_dat$mean - train_dat$tau)^2)
  target_mse <- mean((target_dat$mean - target_dat$tau)^2)
  
  #calculate coverage
  train_coverage <- sum(train_dat$tau >= train_dat$lower & train_dat$tau <= train_dat$upper)/nrow(train_dat)
  target_coverage <- sum(target_dat$tau >= target_dat$lower & target_dat$tau <= target_dat$upper)/nrow(target_dat)
  
  return(c(train_mse = train_mse, target_mse = target_mse,
           train_coverage = train_coverage, target_coverage = target_coverage))
  
}





#overall function ####
compare_oos <- function(N=100, K=6, n_mean=200, n_sd=0, eps_study_m=0.05, 
                        eps_study_tau=0.01, distribution="same", target_dist="same") {
  
  
  ## Simulate training and target (OOS) data
  sim_dat <- gen_mdd(K, n_mean, n_sd, eps_study_m, eps_study_tau, distribution, target_dist)
  train_dat <- sim_dat[["train_dat"]]
  target_dat <- sim_dat[["target_dat"]]
  
  covars <- c("sex", "smstat", "weight", "age", "madrs")
  
  ## Fit mixed effects model
  mod <- lmer(Y ~ W + age + W:age +
                       (W + age + W:age | S), data=train_dat)
  sum <- summary(mod)
  
  #mse for training data
  
  train_mse <- mean(( - tau_true)^2)
  #CI coverage for training data
  train_coverage <- sum(ifelse(tau_true >= lower & tau_true <= upper, 1, 0))/nrow(train_dat)
  train_res <- c(train_mse=train_mse, train_coverage=train_coverage)
  
  
  ## Calculate mean and CIs for each target individual according to each imputation method
  #bootstrap PI
  
  #manual calculation
  
  #confidence interval
  
  
  
  ## Save results
  return(list(
              N=N, K=K, n_mean=n_mean, n_sd=n_sd,
              eps_study_m=eps_study_m, eps_study_tau=eps_study_tau, 
              distribution=distribution, target_dist=target_dist))
}
