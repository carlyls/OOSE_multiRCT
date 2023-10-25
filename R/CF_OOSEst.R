### Bootstrapping Confidence Intervals for Out of Sample CATE Estimation ###

# library(tidyverse)
# library(rsample)
# library(grf)
# library(fastDummies)
# library(nnet)


#### TRAINING DATA ####

cf_pi_train <- function(df, tau_hat) {
  df <- df %>%
    mutate(mean = tau_hat$predictions,
           sd = sqrt(tau_hat$variance.estimates),
           lower = mean - 1.96*sd,
           upper = mean + 1.96*sd)
  return(df)
}


#### TARGET DATA ####

# Mimic meta-analytic prediction intervals
cf_pi_target <- function(K, target_dat, tau_forest, covars) {
  
  #assign study
  new_dat <- target_dat %>%
    slice(rep(1:n(), each=K)) %>%
    mutate(S = rep(1:K, nrow(target_dat)))
  new_feat <- new_dat %>%
    dplyr::select(c(S, all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  
  #predict CATE
  target_pred <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
  new_dat <- new_dat %>%
    mutate(tau_hat = target_pred$predictions,
           var = target_pred$variance.estimates)
  
  #for each X profile, calculate avg tau_hat, avg within-study var, variance of tau_hats
  cis <- new_dat %>%
    group_by(W, sex, smstat, weight, age, age2, madrs, Y, tau) %>%
    summarise(mean = mean(tau_hat),
              var_within = mean(var),
              var_btwn = var(tau_hat),
              n_K = n()) %>%
    mutate(var_tot = var_within + var_btwn,
           sd = sqrt(var_tot),
           lower = mean - qt(.975, K-2)*sd,
           upper = mean + qt(.975, K-2)*sd)
  
  #reorder to match target data
  cis_ord <- target_dat %>%
    left_join(cis, by = c("W","sex","smstat","weight","age",
                          "age2","madrs","Y","tau"))
  
  return(cis_ord)
}


# Random imputation
impute_rand <- function(N, K, target_dat, tau_forest, covars) {

  #assign study
  new_dat <- target_dat %>%
    slice(rep(1:n(), each=N)) %>%
    mutate(S = sample(1:K, nrow(target_dat)*N, replace = T),
           tau_hat = NA)
  new_feat <- new_dat %>%
    dplyr::select(c(S, all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)

  #predict CATE
  target_pred <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
  for (i in 1:nrow(target_pred)) {
    new_dat$tau_hat[i] <- rnorm(1, mean=target_pred$predictions[i],
                                sd = sqrt(target_pred$variance.estimates[i]))
  }

  #create confidence intervals
  cis <- new_dat %>%
    group_by(sex, smstat, weight, age, age2, madrs, tau) %>%
    summarise(mean = mean(tau_hat),
              sd = sd(tau_hat),
              q2.5 = quantile(tau_hat, .025),
              q97.5 = quantile(tau_hat, .975)) %>%
    mutate(lower = mean - qt(.975, K-2)*sd,
           upper = mean + qt(.975, K-2)*sd)

  #reorder to match target data
  cis_ord <- target_dat %>%
    left_join(cis, by = c("sex","smstat","weight","age",
                          "age2","madrs","tau"))

  return(cis_ord)
}





#### OPTION 2: STUDY MEMBERSHIP MODEL ####

# impute_mem <- function(N, train_dat, target_dat, tau_forest, covars) {
#   
#   #create membership model
#   mem_mod <- multinom(S ~ sex + smstat + weight + age + madrs + Y, data=train_dat)
#   #summary(mem_mod)
#   #round(fitted(mem_mod), 2) #looks at probabilities in each class
#   #preds <- predict(mem_mod, newdata = train_dat, "class")
#   #tab <- table(train_dat$S, preds); #round((sum(diag(tab))/sum(tab))*100,2)
#   
#   #define probabilities
#   mem_probs <- predict(mem_mod, newdata = target_dat, type = "probs")
#   S_mem <- c()
#   for (i in 1:nrow(mem_probs)) {
#     S_mem <- c(S_mem, sample(1:K, N, replace=T, prob=mem_probs[i,]))
#   }
#   
#   #assign study
#   new_mem <- target_dat %>%
#     slice(rep(1:n(), each=N)) %>%
#     mutate(S = S_mem,
#            tau_hat = NA)
#   new_feat <- new_mem %>%
#     dplyr::select(c(S, all_of(covars))) %>%
#     fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
#   
#   #predict CATE
#   target_pred <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
#   for (i in 1:nrow(target_pred)) {
#     new_mem$tau_hat[i] <- rnorm(1, mean=target_pred$predictions[i],
#                                 sd = sqrt(target_pred$variance.estimates[i]))
#   }
#   
#   #create confidence intervals
#   cis <- new_mem %>%
#     group_by(sex, smstat, weight, age, age2, madrs, tau) %>%
#     summarise(mean = mean(tau_hat),
#               sd = sd(tau_hat)) %>%
#     mutate(lower = mean - 2*sd,
#            upper = mean + 2*sd)
#   
#   return(cis)
# }
# 
# 
# #### OPTION 3: WITHIN-FOREST ####
# 
# impute_default <- function(K, target_dat, tau_forest, covars) {
#   
#   #default method: https://grf-labs.github.io/grf/REFERENCE.html#missing-values
#   #assign study
#   #we don't need to replicate because we will get the same prediction each time
#   new_default <- target_dat %>%
#     mutate(S = factor(NA, levels=1:K))
#   new_feat <- new_default %>%
#     dplyr::select(c(S, all_of(covars))) %>%
#     fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T, ignore_na = T)
#   
#   #predict CATE
#   cate_default <- predict(tau_forest, newdata = new_feat, estimate.variance = T)
#   new_default$mean <- cate_default$predictions 
#   new_default$sd <- sqrt(cate_default$variance.estimates)
#   
#   #create confidence intervals
#   cis <- new_default %>%
#     mutate(lower = mean - 2*sd,
#            upper = mean + 2*sd)
#   
#   return(cis)
# }

