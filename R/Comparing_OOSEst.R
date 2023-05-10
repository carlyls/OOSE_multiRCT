### Trying Methods for OOS Estimation in Random Effects Meta-Analysis and Causal Forests ###

library(tidyverse)
library(lme4)
library(rsample)
library(multcomp)
library(MASS)
library(grf)
library(nnet)

source("R/MDD_Generation_OOSEst.R")
source("R/MA_OOSEst.R")
source("R/Bootstrap_OOSEst.R")


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
compare_oos <- function(K=10, n_mean=500, n_sd=0, n_target=100, covars_fix="age", covars_rand="age",
                        lin=T, eps_study_m=0.05, eps_study_tau=0.05, eps_study_inter=0.05, 
                        distribution="same", target_dist="same") {
  
  
  ## Simulate training and target (OOS) data
  sim_dat <- gen_mdd(K, n_mean, n_sd, n_target, covars_fix, covars_rand, lin,
                     eps_study_m, eps_study_tau, eps_study_inter, 
                     distribution, target_dist)
  train_dat <- sim_dat[["train_dat"]]
  target_dat <- sim_dat[["target_dat"]]
  
  
  ## Mixed effects model: Correct
  #change for scenario with age^2
  if ("age2" %in% covars_fix) { 
    main_eff <- "Y ~ madrs + sex + age2 + W + "
  } else { 
    main_eff <- "Y ~ madrs + sex + age + W + "
  }
  
  formula <- as.formula(paste0(main_eff, 
                               paste("W", covars_fix, sep=":", collapse=" + "),
                               " + (W + ",
                               paste("W", covars_rand, sep=":", collapse=" + "),
                               " | S)"))
  mod <- lmer(formula, data=train_dat,
              control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000)))
  sum <- summary(mod)
  
  #manual PI
  res <- matrix_var(mod)
  manual_train <- manual_pi_train(train_dat, res, K, covars_fix, covars_rand)
  manual_target <- manual_pi_target(target_dat, res, K, covars_fix, covars_rand)
  manual_res <- assess_interval(manual_train, manual_target)
  
  
  ## Mixed effects model: Incorrect
  formula_wrong <- as.formula(paste0(main_eff, 
                                     paste("W", covars_fix, sep=":", collapse=" + "),
                                     " + (W | S)"))
  mod_wrong <- lmer(formula_wrong, data=train_dat,
                    control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000)))
  sum_wrong <- summary(mod_wrong)
  
  #manual PI
  res_wrong <- matrix_var(mod_wrong)
  manual_train_wrong <- manual_pi_train(train_dat, res_wrong, K, covars_fix, c())
  manual_target_wrong <- manual_pi_target(target_dat, res_wrong, K, covars_fix, c())
  manual_res_wrong <- assess_interval(manual_train_wrong, manual_target_wrong)
  
  
  ## Causal Forest
  if ("age2" %in% covars_fix) { 
    covars <- c("sex", "smstat", "weight", "age2", "madrs")
  } else { 
    covars <- c("sex", "smstat", "weight", "age", "madrs")
  }
  
  feat <- dplyr::select(train_dat, c(S, all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  
  tau_forest <- grf::causal_forest(X=feat, Y=train_dat$Y, W=train_dat$W, 
                                   num.threads=3, honesty=T, num.trees=1000)
  tau_hat <- predict(tau_forest, estimate.variance=T)
  
  #causal forest CI - training
  cf_train <- cf_ci(train_dat, tau_hat)
  
  #causal forest CI - target
  #random
  rand_target <- impute_rand(1000, target_dat, tau_forest, covars)
  #study membership model
  mem_target <- impute_mem(1000, train_dat, target_dat, tau_forest, covars)
  #within-forest default
  default_target <- impute_default(K, target_dat, tau_forest, covars)
  
  #calculate mean and CIs for individuals and assess accuracy
  rand_res <- assess_interval(cf_train, rand_target)
  mem_res <- assess_interval(cf_train, mem_target)
  default_res <- assess_interval(cf_train, default_target)
  
  
  ## Save results
  #data frame of parameters
  params <- data.frame(K=K, n_mean=n_mean, n_sd=n_sd, n_target=n_target, 
                       covars_fix=paste(covars_fix, collapse = ", "), 
                       covars_rand=paste(covars_rand, collapse = ", "), lin=lin,
                       eps_study_m=eps_study_m, eps_study_tau=eps_study_tau, 
                       eps_study_inter=paste(eps_study_inter, collapse = ", "),
                       distribution=distribution, target_dist=target_dist)
  
  #data frame of results
  all_res <- cbind(manual_res, manual_res_wrong, rand_res, mem_res, default_res) %>%
    data.frame() %>%
    rownames_to_column("Metric") %>%
    cbind(params)
  
  return(all_res)
}





#### UNUSED ####

# #confidence interval
# glht_train <- glht_ci(train_dat, mod)
# glht_target <- glht_ci(target_dat, mod)
# glht_res <- assess_interval(glht_train, glht_target)

# #confidence interval - incorrectly specified
# glht_train_wrong <- glht_ci(train_dat, mod_wrong)
# glht_target_wrong <- glht_ci(target_dat, mod_wrong)
# glht_res_wrong <- assess_interval(glht_train_wrong, glht_target_wrong)

# #bootstrap PI
# boot_train <- boot_pi(train_dat, mod, covars_fix, covars_rand)
# boot_target <- boot_pi(target_dat, mod, covars_fix, covars_rand)
# boot_res <- assess_interval(boot_train, boot_target)
# 
# #bootstrap PI - incorrectly specified
# boot_train_wrong <- boot_pi(train_dat, mod_wrong, covars_fix, c())
# boot_target_wrong <- boot_pi(target_dat, mod_wrong, covars_fix, c())
# boot_res_wrong <- assess_interval(boot_train_wrong, boot_target_wrong)
