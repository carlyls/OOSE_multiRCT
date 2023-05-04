### Bootstrapping Confidence Intervals for Out of Sample CATE Estimation ###

library(tidyverse)
library(rsample)
library(grf)
library(fastDummies)
library(nnet)

#source("Comparing_methods_functions.R")
source("R/MDD_Generation_OOSEst.R")


#### OPTION 1: COMPLETELY RANDOM ####

impute_rand <- function(N, target_dat, tau_forest) {
  
  #assign study
  new_dat <- target_dat %>%
    slice(rep(1:n(), each=N)) %>%
    mutate(S = sample(1:K, nrow(target_dat)*N, replace = T),
           tau_hat = NA) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  new_feat <- new_dat %>%
    dplyr::select(-c(W, tau, Y, tau_hat))
  
  #predict CATE
  target_pred <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
  for (i in 1:nrow(target_pred)) {
    new_dat$tau_hat[i] <- rnorm(1, mean=target_pred$predictions[i],
                                sd = sqrt(target_pred$variance.estimates[i]))
  }
  
  #create confidence intervals
  cis <- new_dat %>%
    group_by(sex, smstat, weight, age, madrs, tau) %>%
    summarise(mean = mean(tau_hat),
              sd = sd(tau_hat)) %>%
    mutate(lower = mean + qt(.025, df=N-1)*sd,
           upper = mean + qt(.975, df=N-1)*sd)
  
  return(cis)
}


#### OPTION 2: STUDY MEMBERSHIP MODEL ####

impute_mem <- function(N, train_dat, target_dat, tau_forest) {
  
  #create membership model
  mem_mod <- multinom(S ~ sex + smstat + weight + age + madrs + Y, data=train_dat)
  #summary(mem_mod)
  #round(fitted(mem_mod), 2) #looks at probabilities in each class
  #preds <- predict(mem_mod, newdata = train_dat, "class")
  #tab <- table(train_dat$S, preds); #round((sum(diag(tab))/sum(tab))*100,2)
  
  #define probabilities
  mem_probs <- predict(mem_mod, newdata = target_dat, type = "probs")
  S_mem <- c()
  for (i in 1:nrow(mem_probs)) {
    S_mem <- c(S_mem, sample(1:K, N, replace=T, prob=mem_probs[i,]))
  }
  
  #assign study
  new_mem <- target_dat %>%
    slice(rep(1:n(), each=N)) %>%
    mutate(S = S_mem,
           tau_hat = NA) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  new_feat <- new_mem %>%
    dplyr::select(-c(W, tau, Y, tau_hat))
  
  #predict CATE
  target_pred <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
  for (i in 1:nrow(target_pred)) {
    new_mem$tau_hat[i] <- rnorm(1, mean=target_pred$predictions[i],
                                sd = sqrt(target_pred$variance.estimates[i]))
  }
  
  #create confidence intervals
  cis <- new_mem %>%
    group_by(sex, smstat, weight, age, madrs, tau) %>%
    summarise(mean = mean(tau_hat),
              sd = sd(tau_hat)) %>%
    mutate(lower = mean + qt(.025, df=N-1)*sd,
           upper = mean + qt(.975, df=N-1)*sd)
  
  return(cis)
}


#### OPTION 3: WITHIN-FOREST ####

impute_default <- function(K, target_dat, tau_forest) {
  
  #default method: https://grf-labs.github.io/grf/REFERENCE.html#missing-values
  #assign study
  #we don't need to replicate because we will get the same prediction each time
  new_default <- target_dat %>%
    mutate(S = factor(NA, levels=1:K)) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T, ignore_na = T)
  new_feat <- new_default %>%
    dplyr::select(-c(W, tau, Y))
  
  #predict CATE
  cate_default <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
  new_default$mean <- cate_default$predictions 
  new_default$sd <- sqrt(cate_default$variance.estimates)
  
  #create confidence intervals
  cis <- new_default %>%
    mutate(lower = mean + qt(.025, df=N-1)*sd,
           upper = mean + qt(.975, df=N-1)*sd)
  
  return(cis)
}

