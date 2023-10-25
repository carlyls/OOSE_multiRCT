### Trying Methods for OOS Estimation in Random Effects Meta-Analysis and Causal Forests ###

# library(tidyverse)
# library(lme4)
# library(rsample)
# library(multcomp)
# library(MASS)
# library(grf)
# library(nnet)


#check results for all methods ####
assess_interval <- function(train_dat, target_dat) {
  
  #calculate mean absolute bias
  train_bias <- mean(abs(train_dat$mean - train_dat$tau))
  target_bias <- mean(abs(target_dat$mean - target_dat$tau))
  
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
  
  return(c(train_bias = train_bias, target_bias = target_bias,
           train_mse = train_mse, target_mse = target_mse,
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
  cf_train <- cf_pi_train(train_dat, tau_hat)
  
  #causal forest CI - target
  cf_target <- cf_pi_target(K, target_dat, tau_forest, covars)

  rm(tau_forest)
  
  #calculate mean and CIs for individuals and assess accuracy
  cf_res <- assess_interval(cf_train, cf_target)
  
  
  ## Adaptive Causal Forest
  tau_forest_a <- grf::causal_forest(X=feat, Y=train_dat$Y, W=train_dat$W, 
                                   num.threads=3, honesty=F, num.trees=1000)
  tau_hat_a <- predict(tau_forest_a, estimate.variance=T)
  
  #causal forest CI - training
  cf_train_a <- cf_pi_train(train_dat, tau_hat_a)
  
  #causal forest CI - target
  cf_target_a <- cf_pi_target(K, target_dat, tau_forest_a, covars)
  
  rm(tau_forest_a)
  
  #calculate mean and CIs for individuals and assess accuracy
  cf_a_res <- assess_interval(cf_train_a, cf_target_a)
  
  
  ## BART: S-learner
  #use covariates from above
  #update features to include W
  feat <- dplyr::select(train_dat, c(W, S, all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  
  #include counterfactual covariates (swap control and treatment)
  feat_cf <- feat %>%
    mutate(W = as.numeric(W == 0))
  y <- as.numeric(train_dat$Y)
  
  #run bart
  sbart <- dbarts::bart(x.train=as.matrix(feat), y.train=y, x.test=as.matrix(feat_cf), keeptrees=T)
  
  #S-BART credible interval - training
  sb_train <- sbart_ci(train_dat, sbart)
  
  #S-BART credible interval - target
  sb_target <- sbart_target(K, target_dat, sbart, covars)
  
  rm(sbart)
  
  #calculate mean and CIs for individuals and assess accuracy
  sb_res <- assess_interval(sb_train, sb_target)
  
  
  ## BART: T-learner
  #use covariates from above
  m1_setup <- tlearn_setup(train_dat, covars, w=1)
  m0_setup <- tlearn_setup(train_dat, covars, w=0)
  
  #run bart for m1
  tbart1 <- dbarts::bart(x.train = as.matrix(m1_setup[["feat"]]), y.train = m1_setup[["y"]],
                         x.test = as.matrix(m0_setup[["feat"]]), keeptrees = T)
  
  #run bart for m0
  tbart0 <- dbarts::bart(x.train = as.matrix(m0_setup[["feat"]]), y.train = m0_setup[["y"]],
                         x.test = as.matrix(m1_setup[["feat"]]), keeptrees = T)
  
  #T-BART credible interval - training
  tb_train <- tbart_ci(train_dat, tbart1, tbart0)
  
  #T-BART credible interval - target
  tb_target <- tbart_target(K, target_dat, tbart1, tbart0, covars)
  
  rm(tbart1)
  rm(tbart0)
  
  #calculate mean and CIs for individuals and assess accuracy
  tb_res <- assess_interval(tb_train, tb_target)
  
  
  ## Save results
  #data frame of parameters
  params <- data.frame(K=K, n_mean=n_mean, n_sd=n_sd, n_target=n_target, 
                       covars_fix=paste(covars_fix, collapse = ", "), 
                       covars_rand=paste(covars_rand, collapse = ", "), lin=lin,
                       eps_study_m=eps_study_m, eps_study_tau=eps_study_tau, 
                       eps_study_inter=paste(eps_study_inter, collapse = ", "),
                       distribution=distribution, target_dist=target_dist)
  
  #data frame of results
  all_res <- cbind(manual_res, manual_res_wrong, cf_res, cf_a_res, sb_res, tb_res) %>%
    data.frame() %>%
    rownames_to_column("Metric") %>%
    cbind(params)
  
  #individual prediction errors from each method
  ipe <- target_dat %>%
    mutate(manual_ipe = manual_target$tau - manual_target$mean,
           manual_wrong_ipe = manual_target_wrong$tau - manual_target_wrong$mean,
           cf_ipe = cf_target$tau - cf_target$mean,
           cf_a_ipe = cf_target_a$tau - cf_target_a$mean,
           sb_ipe = sb_target$tau - sb_target$mean,
           tb_ipe = tb_target$tau - tb_target$mean)
  
  
  return(list(all_res=all_res, ipe=ipe))
}

