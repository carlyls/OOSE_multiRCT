### Bootstrapping Confidence Intervals for Out of Sample CATE Estimation ###

library(tidyverse)
library(rsample)
library(grf)
library(fastDummies)
library(nnet)

#source("Comparing_methods_functions.R")
source("MDD_Simulation_OOSEst.R")


#### OPTION 1: COMPLETELY RANDOM ####

impute_rand <- function(N, test_dat, tau_forest) {
  
  #assign study
  new_dat <- test_dat %>%
    slice(rep(1:n(), each=N)) %>%
    mutate(S = sample(1:K, nrow(test_dat)*N, replace = T)) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  new_feat <- new_dat %>%
    select(-c(W, tau, Y))
  
  #predict CATE
  new_dat$tau_hat <- predict(tau_forest, newdata = new_feat)$predictions
  
  #create confidence intervals
  cis <- new_dat %>%
    group_by(sex, smstat, weight, age, madrs, tau) %>%
    summarise(mean = mean(tau_hat),
              sd = sd(tau_hat)) %>%
    mutate(lower = mean + qt(.025, df=N-1)*sd,
           upper = mean + qt(.975, df=N-1)*sd)
  
  #calculate accuracy
  mse <- mean((cis$mean - cis$tau)^2)
  ci_coverage <- sum(ifelse(cis$tau >= cis$lower & cis$tau <= cis$upper, 1, 0))/nrow(cis)
  ci_length <- mean(cis$upper - cis$lower)
  
  return(list(mse=mse, ci_coverage=ci_coverage, ci_length=ci_length, cis=cis))
}

#### OPTION 2: STUDY MEMBERSHIP MODEL ####

impute_mem <- function(N, train_dat, test_dat, tau_forest) {
  
  #create membership model
  mem_mod <- multinom(S ~ sex + smstat + weight + age + madrs + Y, data=train_dat)
  #summary(mem_mod)
  #round(fitted(mem_mod), 2) #looks at probabilities in each class
  #preds <- predict(mem_mod, newdata = train_dat, "class")
  #tab <- table(train_dat$S, preds); #round((sum(diag(tab))/sum(tab))*100,2)
  
  #define probabilities
  mem_probs <- predict(mem_mod, newdata = test_dat, type = "probs")
  S_mem <- c()
  for (i in 1:nrow(mem_probs)) {
    S_mem <- c(S_mem, sample(1:K, N, replace=T, prob=mem_probs[i,]))
  }
  
  #assign study
  new_mem <- test_dat %>%
    slice(rep(1:n(), each=N)) %>%
    mutate(S = S_mem) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  new_feat <- new_mem %>%
    select(-c(W, tau, Y))
  
  #predict CATE
  new_mem$tau_hat <- predict(tau_forest, newdata = new_feat)$predictions
  
  #create confidence intervals
  cis <- new_mem %>%
    group_by(sex, smstat, weight, age, madrs, tau) %>%
    summarise(mean = mean(tau_hat),
              sd = sd(tau_hat)) %>%
    mutate(lower = mean + qt(.025, df=N-1)*sd,
           upper = mean + qt(.975, df=N-1)*sd)
  
  #calculate accuracy
  mse <- mean((cis$mean - cis$tau)^2)
  ci_coverage <- sum(ifelse(cis$tau >= cis$lower & cis$tau <= cis$upper, 1, 0))/nrow(cis)
  ci_length <- mean(cis$upper - cis$lower)
  
  return(list(mse=mse, ci_coverage=ci_coverage, ci_length=ci_length, cis=cis))
}

#### OPTION 3: WITHIN-FOREST ####

impute_default <- function(K, test_dat, tau_forest) {
  
  #default method: https://grf-labs.github.io/grf/REFERENCE.html#missing-values
  #assign study
  #we don't need to replicate because we will get the same prediction each time
  new_default <- test_dat %>%
    mutate(S = factor(NA, levels=1:K)) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T, ignore_na = T)
  new_feat <- new_default %>%
    select(-c(W, tau, Y))
  
  #predict CATE
  cate_default <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
  new_default$mean <- cate_default$predictions 
  new_default$sd <- sqrt(cate_default$variance.estimates)
  
  #create confidence intervals
  cis <- new_default %>%
    mutate(lower = mean + qt(.025, df=N-1)*sd,
           upper = mean + qt(.975, df=N-1)*sd)
  
  #calculate accuracy
  mse <- mean((cis$mean - cis$tau)^2)
  ci_coverage <- sum(ifelse(cis$tau >= cis$lower & cis$tau <= cis$upper, 1, 0))/nrow(cis)
  ci_length <- mean(cis$upper - cis$lower)
  
  return(list(mse=mse, ci_coverage=ci_coverage, ci_length=ci_length, cis=cis))
  
}

#### OPTION 4: WITHIN-FOREST RANDOM SAMPLING ####


#### OPTION 5: GRADE OF MEMBERSHIP MODEL ####



##################################

## FUNCTION

compare_oos <- function(N=100, K=6, n_mean=200, n_sd=0, eps_study_m=0.05, 
                        eps_study_tau=0.01, distribution="same", test_dist="same") {
  
  
  ## Simulate training and testing (OOS) data
  sim_dat <- gen_mdd(K, n_mean, n_sd, eps_study_m, eps_study_tau, distribution, test_dist)
  train_dat <- sim_dat[["train_dat"]]
  test_dat <- sim_dat[["test_dat"]]
  
  covars <- c("sex", "smstat", "weight", "age", "madrs")
  feat <- select(train_dat, c(S,all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  tau_true <- train_dat$tau
  
  
  ## Fit models (for now, causal forest with pooling with trial indicator)
  tau_forest <- causal_forest(X=feat, Y=train_dat$Y, W=train_dat$W, 
                              num.threads=3, honesty=F, num.trees=1000)
  tau_hat <- predict(tau_forest, estimate.variance=T)
  #mse for training data
  train_mse <- mean((tau_hat$predictions - tau_true)^2)
  #CI coverage for training data
  sd <- sqrt(tau_hat$variance.estimates)
  lower <- tau_hat$predictions + qt(.025, df=nrow(train_dat)-1)*sd
  upper <- tau_hat$predictions + qt(.975, df=nrow(train_dat)-1)*sd
  train_coverage <- sum(ifelse(tau_true >= lower & tau_true <= upper, 1, 0))/nrow(train_dat)
  train_res <- c(train_mse=train_mse, train_coverage=train_coverage)
  
  
  ## Calculate mean and CIs for each test individual according to each imputation method
  #random
  res_rand <- impute_rand(N, test_dat, tau_forest)
  #study membership model
  res_mem <- impute_mem(N, train_dat, test_dat, tau_forest)
  #within-forest default
  res_default <- impute_default(K, test_dat, tau_forest)
  
  
  ## Save results
  return(list(train_res=train_res, res_rand=res_rand, 
              res_mem=res_mem, res_default=res_default,
              N=N, K=K, n_mean=n_mean, n_sd=n_sd,
              scenario=scenario, distribution=distribution,
              test_dat=test_dat, test_scenario=test_scenario))
}

