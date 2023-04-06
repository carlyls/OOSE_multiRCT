## Trying Methods for OOS Estimation in Random Effects Meta-Analysis ##

library(tidyverse)
library(lme4)
library(rsample)
library(multcomp)
library(MASS)
#library(pimeta)
#library(merTools)

source("R/MDD_Generation_OOSEst.R")

#prediction interval by hand ####
#get variance components of model
matrix_var <- function(mod, rand_int=T) {

  #calculate variances
  fc <- vcov(mod) #covariance matrix of fixed effects
  rc <- as.data.frame(VarCorr(mod)) #variance and covariances of random effects
  
  #fixed effects
  beta <- matrix(fixef(mod)[c("W","W:age")], nrow=2) #beta-hat
  var_beta <- fc[c("W","W:age"),c("W","W:age")] #Var(beta-hat)
  
  #random effects
  if (rand_int == T) {
    u <- ranef(mod)$S[c("W","W:age")] %>% t() %>% matrix() #u-hat
    rand <- rc[c(which(rc$var1=="W" & is.na(rc$var2)==T),
                 which(rc$var1=="W" & rc$var2=="W:age"),
                 which(rc$var1=="W" & rc$var2=="W:age"),
                 which(rc$var1=="W:age" & is.na(rc$var2)==T)), "vcov"] 
    var_rand <- matrix(rand, nrow=2, dimnames=list(c("W","W:age"),c("W","W:age"))) #Var(ranef) 
  } else {
    u <- ranef(mod)$S[c("W")] %>% t() %>% matrix() #u-hat
    rand <- rc[c(which(rc$var1=="W" & is.na(rc$var2)==T)), "vcov"] 
    var_rand <- matrix(rand, dimnames=list(c("W"),c("W"))) #Var(ranef) 
  }
  
  return(list(beta=beta, var_beta=var_beta, u=u, var_rand=var_rand))
}

#add pis to dataset
manual_pi <- function(df, mod, K, rand_int=T) {
  
  #get variances
  res <- matrix_var(mod, rand_int)
  beta <- res$beta
  var_beta <- res$var_beta
  u <- res$u
  var_rand <- res$var_rand
  
  #calculate theta-hats
  X <- df %>%
    mutate(trt = 1) %>%
    dplyr::select(trt, age) %>%
    as.matrix() #X
  
  if ("S" %in% names(df)) { #training data structure is different
    Zlist <- list()
    var_rand_train <- list()
    for (s in 1:K) {
      X_S <- X[df$S==s,]
      if (rand_int == T) { Z_S <- X_S } else { Z_S <- X_S[,1] }
      Zlist[[s]] <- Z_S
      var_rand_train[[s]] <- var_rand
    }
    Z <- do.call("bdiag", Zlist) %>% as.matrix()
    var_rand <- do.call("bdiag", var_rand_train) %>% as.matrix()
    
    mean_theta <- X %*% beta + Z %*% u %>% c()

  } else { #target data doesn't require block diagonal matrices
    if (rand_int == T) { Z <- X } else { Z <- X[,1] }
    
    mean_theta <- X %*% beta %>% c() #theta-hat
  }
  
  vcov_theta <- X %*% var_beta %*% t(X) + Z %*% var_rand %*% t(Z)
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
sample_cate <- function(x, boot_fix, boot_rand, rand_int=T) {
  
  #get cate|x according to each coefficient
  if (rand_int == T) {
    cates <- (boot_fix$W + boot_rand$W) + 
      (boot_fix$W.age + boot_rand$W.age)*x
  } else {
    cates <- (boot_fix$W + boot_rand$W) + 
      (boot_fix$W.age)*x
  }
  
  #calculate interval over all iterations
  mean <- mean(cates)
  sd <- sd(cates)
  lower <- mean - 1.96*sd
  upper <- mean + 1.96*sd

  return(c(lower=lower, mean=mean, upper=upper))
}

#bootstrap and add intervals to data
boot_pi <- function(df, mod, rand_int=T) {
  
  #get var-covar of fixed and random effects from model
  res <- matrix_var(mod, rand_int) 
  
  #randomly sample fixed and random coefficients
  boot_fix <- mvrnorm(n=1000, mu=c(res$beta), Sigma=res$var_beta) %>%
    data.frame()
  boot_rand <- mvrnorm(n=1000, mu=rep(0, nrow(res$var_rand)), Sigma=res$var_rand) %>%
    data.frame() #assume ranefs have mean 0
  
  #apply coefficients to estimate cate for all ages
  intervals <- map_dfr(.x=df$age, .f=sample_cate,
                       boot_fix=boot_fix, boot_rand=boot_rand, rand_int=rand_int)
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
  
  #calculate signficance
  train_significance <- sum(sign(train_dat$lower) == sign(train_dat$upper))/nrow(train_dat)
  target_significance <- sum(sign(target_dat$lower) == sign(target_dat$upper))/nrow(target_dat)
  
  return(c(train_mse = train_mse, target_mse = target_mse,
           train_coverage = train_coverage, target_coverage = target_coverage,
           train_length = train_length, target_length = target_length,
           train_significance = train_significance, target_significance = target_significance))
  
}


#overall function ####
compare_oos <- function(K=10, n_mean=200, n_sd=0, n_target=100, eps_study_m=0.05, eps_study_tau=0.05, 
                        eps_study_age=0.05, distribution="same", target_dist="same") {
  
  
  ## Simulate training and target (OOS) data
  sim_dat <- gen_mdd(K, n_mean, n_sd, n_target, eps_study_m, eps_study_tau, 
                     eps_study_age, distribution, target_dist)
  train_dat <- sim_dat[["train_dat"]]
  target_dat <- sim_dat[["target_dat"]]

  
  ## Fit mixed effects models
  #correct
  mod <- lmer(Y ~  madrs + sex + W*age +
                  (W + W:age | S), data=train_dat,
              control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000)))
  sum <- summary(mod)
  
  #incorrect
  mod_wrong <- lmer(Y ~  madrs + sex + W*age +
                (W | S), data=train_dat,
              control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000)))
  sum_wrong <- summary(mod_wrong)
  
  
  ## Calculate mean and CIs for individuals and assess accuracy
  #confidence interval
  glht_train <- glht_ci(train_dat, mod)
  glht_target <- glht_ci(target_dat, mod)
  glht_res <- assess_interval(glht_train, glht_target)
  
  #confidence interval - incorrectly specified
  glht_train_wrong <- glht_ci(train_dat, mod_wrong)
  glht_target_wrong <- glht_ci(target_dat, mod_wrong)
  glht_res_wrong <- assess_interval(glht_train_wrong, glht_target_wrong)
    
  #manual PI
  manual_train <- manual_pi(train_dat, mod, K, rand_int=T)
  manual_target <- manual_pi(target_dat, mod, K, rand_int=T)
  manual_res <- assess_interval(manual_train, manual_target)
  
  #manual PI - incorrectly specified
  manual_train_wrong <- manual_pi(train_dat, mod_wrong, K, rand_int=F)
  manual_target_wrong <- manual_pi(target_dat, mod_wrong, K, rand_int=F)
  manual_res_wrong <- assess_interval(manual_train_wrong, manual_target_wrong)

  #bootstrap PI
  boot_train <- boot_pi(train_dat, mod, rand_int=T)
  boot_target <- boot_pi(target_dat, mod, rand_int=T)
  boot_res <- assess_interval(boot_train, boot_target)
  
  #bootstrap PI - incorrectly specified
  boot_train_wrong <- boot_pi(train_dat, mod_wrong, rand_int=F)
  boot_target_wrong <- boot_pi(target_dat, mod_wrong, rand_int=F)
  boot_res_wrong <- assess_interval(boot_train_wrong, boot_target_wrong)
  
  
  ## Save results
  #data frame of parameters
  params <- data.frame(K=K, n_mean=n_mean, n_sd=n_sd, n_target=n_target, eps_study_m=eps_study_m, 
                       eps_study_tau=eps_study_tau, eps_study_age=eps_study_age,
                       distribution=distribution, target_dist=target_dist)
  
  #data frame of results
  all_res <- cbind(glht_res, glht_res_wrong, manual_res, 
               manual_res_wrong, boot_res, boot_res_wrong) %>%
    data.frame() %>%
    rownames_to_column("Metric") %>%
    cbind(params)
  
  return(list(sum=sum, sum_wrong=sum_wrong, all_res=all_res))
}
