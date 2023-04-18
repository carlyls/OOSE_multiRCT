### Simulation Approach for OOS Estimation ###
### Based on MDD Data ###

library(tidyverse)
library(rsample)
library(lme4)
library(multcomp)
library(MASS)
#library(pimeta)
#library(merTools)

#interior function
sample_dist <- function(K, k, n, Sigma, eps_study_m, eps_study_tau, eps_study_inter, distribution) {
  
  #define mu based on distribution input
  if (distribution == "same") {
    mu <- c(age=44.8971, sex=0.6784, smstat=0.3043, weight=79.0253, madrs=31.4088)

    } else if (distribution == "varying_madrs") {
      mu <- c(age=44.8971, sex=0.6784, smstat=0.3043, weight=79.0253, madrs=31.4088-k*1.5)
    
    } else if (distribution == "halfdiff_madrsage") {
      if (k%%2 == 0 ) {
        mu <- c(age=50, sex=0.6784, smstat=0.3043, weight=79.0253, madrs=40)
      } else {
        mu <- c(age=44.8971, sex=0.6784, smstat=0.3043, weight=79.0253, madrs=31.4088)
      }
    
    } else if (distribution == "separate_age") {
      ages <- seq(30, 55, length.out=K)
      mu <- c(age=ages[k], sex=0.6784, smstat=0.3043, weight=79.0253, madrs=31.4088)
      
    }
  
  #define random slopes for moderators
  eps_inter <- matrix(nrow=n, ncol=length(covars_rand))
  colnames(eps_inter) <- paste0("eps_",covars_rand)
  for (i in 1:length(covars_rand)) {
    eps_inter[,i] <- rnorm(n=1, mean=0, sd=eps_study_inter[i])
  }
  
  #simulate data
  dat <- MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma) %>%
    as.data.frame() %>%
    mutate(S = k,
           W = rbinom(n=n, size=1, prob=.5),
           sex = ifelse(sex > 1-0.6784, 1, 0),
           smstat = ifelse(smstat > 1-0.3043, 1, 0),
           eps = rnorm(n=n, mean=0, sd=.05),
           eps_m = rnorm(n=1, mean=0, sd=eps_study_m),
           eps_tau = rnorm(n=1, mean=0, sd=eps_study_tau)) %>%
    bind_cols(eps_inter)
  
  return(dat)
}

#main function
gen_mdd <- function (K=10, n_mean=500, n_sd=0, n_target=100, covars_fix="age", covars_rand="age",
                     eps_study_m=0.05, eps_study_tau=0.05, eps_study_inter=0.05,
                     distribution="same", target_dist="same") {
  
  #training data
  train_dat <- data.frame()
  n_study <- floor(rnorm(K, mean=n_mean, sd=n_sd))
  
  #define covariance matrix
  Sigma <- data.frame(age=c(165.6471, 0.2448, -0.5180, 1.6408, -0.9666),
                      sex=c(0.2448, 0.2183, -0.0218, -1.9030, 0.1380),
                      smstat=c(-0.5180, -0.0218, 0.2118, -0.1429, 0.1155),
                      weight=c(1.6408, -1.9030, -0.1428, 452.6100, -7.6864),
                      madrs=c(-0.9666, 0.1380, 0.1155, -7.6864, 17.5343),
                      row.names=c("age","sex","smstat","weight","madrs"))
  
  for (k in 1:K) {
    n <- n_study[k]
    
    #sample
    dat <- sample_dist(K, k, n, Sigma, eps_study_m, eps_study_tau, eps_study_inter, distribution)
    train_dat <- bind_rows(train_dat, dat)
    
  }
  
  #target data
  if (target_dist == "same") {
    target_dat <- train_dat[sample(nrow(train_dat), n_target),] %>%
      dplyr::select(-S, -contains("eps_"))
    
  } else if (target_dist == "upweight") {
    train_weight <- train_dat %>%
      mutate(study_weight = ifelse(S %in% c(3, 5), 3, 1))
    target_dat <- train_weight[sample(nrow(train_weight), n_target, prob=train_weight$study_weight),] %>%
      dplyr::select(-study_weight, -S, -contains("eps_"))
    
  } else if (target_dist == "different") {
    target_mean <- c(age=30, sex=0.6784, smstat=0.3043, weight=79.0253, madrs=25)
    target_dat <- MASS::mvrnorm(n=n_target, mu=target_mean, Sigma=Sigma) %>%
      as.data.frame() %>%
      mutate(sex = ifelse(sex > 1-0.6784, 1, 0),
             smstat = ifelse(smstat > 1-0.3043, 1, 0),
             eps = rnorm(n=n_target, mean=0, sd=.05),
             W = rbinom(n=n_target, size=1, prob=.5))
    
  }
  
  #define random slopes for moderators in target sample
  eps_inter_target <- matrix(nrow=n_target, ncol=length(covars_rand))
  colnames(eps_inter_target) <- paste0("eps_",covars_rand)
  for (i in 1:length(covars_rand)) {
    eps_inter_target[,i] <- rnorm(n=n_target, mean=0, sd=eps_study_inter[i])
  }
  
  #add error terms to target sample
  target_dat <- target_dat %>%
    mutate(eps_m = rnorm(n=n_target, mean=0, sd=eps_study_m),
           eps_tau = rnorm(n=n_target, mean=0, sd=eps_study_tau)) %>%
    bind_cols(eps_inter_target)
  
  #standardize variables
  #add m and tau
  train_dat <- train_dat %>% 
    mutate(age = (age - mean(age))/sd(age),
           madrs = (madrs - mean(madrs))/sd(madrs),
           weight = (weight - mean(weight))/sd(weight))
  target_dat <- target_dat %>% 
    mutate(age = (age - mean(age))/sd(age),
           madrs = (madrs - mean(madrs))/sd(madrs),
           weight = (weight - mean(weight))/sd(weight))
  
  if (length(covars_fix) == 1 & length(covars_rand) == 1) {
    train_dat <- train_dat %>% 
      mutate(m = (-17.40 + eps_m) - 0.13*age - 2.05*madrs - 0.11*sex,
             tau = (2.505 + eps_tau) + (0.82 + eps_age)*age)
    target_dat <- target_dat %>% 
      mutate(m = (-17.40 + eps_m) - 0.13*age - 2.05*madrs - 0.11*sex,
             tau = (2.505 + eps_tau) + (0.82 + eps_age)*age)
  } else if (covars_fix == c("age", "madrs") & covars_rand == c("age", "madrs")) {
    ## FIGURE THIS OUT
  } else if (covars_fix == c("age", "sex") & covars_rand == c("age", "sex")) {
    ## FIGURE THIS OUT
  } else if (covars_fix == c("age", "sex") & covars_rand == c("age")) {
    ## FIGURE THIS OUT
  }
  
  #outcome Y
  train_dat <- train_dat %>%
    mutate(Y = m + W*tau + eps,
           S = factor(S)) %>%
    dplyr::select(S, W, sex, smstat, weight, age, madrs, Y, tau)
  
  target_dat <- target_dat %>%
    mutate(Y = m + W*tau + eps) %>%
    dplyr::select(W, sex, smstat, weight, age, madrs, Y, tau)
  
  return(list(train_dat=train_dat, target_dat=target_dat))
}
