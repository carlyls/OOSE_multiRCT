### Random Effects Meta-Analysis Prediction Intervals ###

library(tidyverse)
library(lme4)
library(rsample)
library(multcomp)
library(MASS)
library(grf)


#### OPTION 1: Prediction interval by hand ####

#get variance components of model
matrix_var <- function(mod) {
  
  #calculate variances
  fc <- vcov(mod) #covariance matrix of fixed effects
  rc <- Matrix::bdiag(VarCorr(mod)) #variance and covariances of random effects
  
  #fixed effects
  beta <- fixef(mod)[grep("W", names(fixef(mod)))] %>% matrix() #beta-hat
  var_beta <- fc[grep("W", rownames(fc)),
                 grep("W", colnames(fc))] #Var(beta-hat)
  
  #random effects
  u <- ranef(mod)$S[grep("W", colnames(ranef(mod)$S))] %>% t() %>% matrix()
  var_rand <- rc[grep("W", rownames(rc)),
                 grep("W", colnames(rc))]
  v_uhat <- attr(ranef(mod, condVar = TRUE)[[1]], "postVar")
  
  return(list(beta=beta, var_beta=var_beta, u=u, var_rand=var_rand, v_uhat=v_uhat))
}

#add pis to dataset
manual_pi_train <- function(df, res, K, covars_fix, covars_rand) {
  
  #get variances
  #res <- matrix_var(mod)
  beta <- res$beta
  var_beta <- res$var_beta
  u <- res$u
  #var_rand <- res$var_rand
  v_uhat <- res$v_uhat
  
  #calculate theta-hats
  X <- df %>%
    dplyr::select(W, all_of(covars_fix)) %>%
    mutate(W = 1) %>%
    as.matrix()
  Z <- df %>%
    dplyr::select(W, all_of(covars_rand)) %>%
    mutate(W = 1) %>%
    as.matrix()
  
  Zlist <- list()
  for (s in 1:K) {
    Z_S <- Z[df$S==s,]
    Zlist[[s]] <- Z_S
  }
  Z <- do.call("bdiag", Zlist) %>% as.matrix()
  
  var_uhat_train <- list()
  for (i in 1:K) {
    var_uhat_train[[i]] <- v_uhat[,,i][2:ncol(v_uhat[,,i]),2:ncol(v_uhat[,,i])]
  }
  var_uhat <- do.call("bdiag", var_uhat_train) %>% as.matrix()

  mean_theta <- X %*% beta + Z %*% u %>% c()
  
  vcov_theta <- X %*% var_beta %*% t(X) + Z %*% var_uhat %*% t(Z)
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

manual_pi_target <- function(df, res, K, covars_fix, covars_rand) {
  
  #get variances
  #res <- matrix_var(mod)
  beta <- res$beta
  var_beta <- res$var_beta
  u <- res$u
  var_rand <- res$var_rand
  #v_uhat <- res$v_uhat
  
  #calculate theta-hats
  X <- df %>%
    dplyr::select(W, all_of(covars_fix)) %>%
    mutate(W = 1) %>%
    as.matrix()
  Z <- df %>%
    dplyr::select(W, all_of(covars_rand)) %>%
    mutate(W = 1) %>%
    as.matrix()
  
  mean_theta <- X %*% beta %>% c() #theta-hat
  
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





##### Currently unused methods for main function

#confidence interval by glht ####
## still needs to be updated to account for other moderators ##
#get ci based on x=age (as a character)
# age_ci <- function(x, mod) {
# 
# ci <- glht(mod, linfct=paste0("W +",x,"*W:age = 0"))
# mean <- coef(ci) %>% as.numeric()
# sd <- sqrt(vcov(ci)[1,1])
# lower <- mean - 1.96*sd %>% as.numeric()
# upper <- mean + 1.96*sd %>% as.numeric()
# 
# return(c(lower=lower, mean=mean, upper=upper))
# } 

#add cis to dataset
# glht_ci <- function(df, mod) {
# 
# ages <- as.character(df$age)
# cis <- map_dfr(.x=ages, .f=age_ci, mod=mod)
# df <- df %>%
#   bind_cols(cis)
# 
# return(df)
# } 

#prediction interval by bootstrap ####
#create intervals from bootstrap
# sample_cate <- function(fix, rand, boot_fix, boot_rand) {
# 
# #get cate|x according to each coefficient
# cates <- boot_fix %*% t(fix) + boot_rand %*% t(rand)
# 
# #calculate interval over all iterations
# mean <- mean(cates)
# sd <- sd(cates)
# lower <- mean - 1.96*sd
# upper <- mean + 1.96*sd
# 
# return(c(lower=lower, mean=mean, upper=upper))
# }

#bootstrap and add intervals to data
# boot_pi <- function(df, mod, covars_fix, covars_rand) {
# 
# #get var-covar of fixed and random effects from model
# res <- matrix_var(mod) 
# 
# #randomly sample fixed and random coefficients
# boot_fix <- mvrnorm(n=1000, mu=c(res$beta), Sigma=res$var_beta)
# boot_rand <- mvrnorm(n=1000, mu=rep(0, 1+length(covars_rand)), Sigma=res$var_rand) #assume ranefs have mean 0
# 
# fix <- dplyr::select(df, W, all_of(covars_fix)) %>%
#   mutate(W=1)
# rand <- dplyr::select(df, W, all_of(covars_rand)) %>%
#   mutate(W=1)
# 
# #apply coefficients to estimate cate for all ages
# intervals <- c()
# for (i in 1:nrow(fix)) {
#   intervals <- bind_rows(intervals, sample_cate(fix[i,], rand[i,], boot_fix, boot_rand))
# }
# 
# df <- df %>%
#   bind_cols(intervals)
# 
# return(df)
# }
