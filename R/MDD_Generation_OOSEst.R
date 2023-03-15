### Simulation Approach for OOS Estimation ###
### Based on MDD Data ###

library(tidyverse)
library(rsample)
library(truncnorm)

#interior function
sample_dist <- function(k, n, Sigma, eps_study_m, eps_study_tau, distribution) {
  
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
      ages <- c(30, 35, 40, 45, 50, 55)
      mu <- c(age=ages[k], sex=0.6784, smstat=0.3043, weight=79.0253, madrs=31.4088)
      
    }
  
  #simulate data
  dat <- MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma) %>%
    as.data.frame() %>%
    mutate(sex = ifelse(sex > 1-0.6784, 1, 0),
           smstat = ifelse(smstat > 1-0.3043, 1, 0),
           eps = rnorm(n=n, mean=0, sd=.05),
           W = rbinom(n=n, size=1, prob=.5),
           eps_m = rnorm(n=1, mean=0, sd=eps_study_m),
           eps_tau = rnorm(n=1, mean=0, sd=eps_study_tau),
           S = k)
  
  return(dat)
}

#main function
gen_mdd <- function (K=6, n_mean=200, n_sd=0, eps_study_m=0.05, eps_study_tau=0.01, 
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
    dat <- sample_dist(k, n, Sigma, eps_study_m, eps_study_tau, distribution)
    
    train_dat <- bind_rows(train_dat, dat)
  }
  
  #target data
  if (target_dist == "same") {
    target_dat <- train_dat[sample(nrow(train_dat), 100),] %>%
      dplyr::select(-S, -eps_m, -eps_tau)
    
  } else if (target_dist == "upweight") {
    train_weight <- train_dat %>%
      mutate(study_weight = ifelse(S %in% c(3, 5), 3, 1))
    target_dat <- train_weight[sample(nrow(train_weight), 100, prob=train_weight$study_weight),] %>%
      dplyr::select(-study_weight, -S, -eps_m, -eps_tau)
    
  } else if (target_dist == "different") {
    target_mean <- c(age=30, sex=0.6784, smstat=0.3043, weight=79.0253, madrs=25)
    target_dat <- MASS::mvrnorm(n=100, mu=target_mean, Sigma=Sigma) %>%
      as.data.frame() %>%
      mutate(sex = ifelse(sex > 1-0.6784, 1, 0),
             smstat = ifelse(smstat > 1-0.3043, 1, 0),
             eps = rnorm(n=100, mean=0, sd=.05),
             W = rbinom(n=100, size=1, prob=.5))
    
  }
  
  #m and tau
  train_dat <- train_dat %>% 
    mutate(m = -0.02*age - 0.7*madrs - 0.15*sex + eps_m,
           tau = -8.5 + 0.07*age + 0.20*madrs + eps_tau)
  target_dat <- target_dat %>% 
    mutate(m = -0.02*age - 0.7*madrs - 0.15*sex,
           tau = -8.5 + 0.07*age + 0.20*madrs)
  
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
