# Testing variance calculation for S-learner

library(tidyverse)
library(dbarts)

source("R/BART_OOSEst.R")
source("R/MDD_Generation_OOSEst.R")
source("R/Comparing_OOSEst.R")


# set up parameters
K <- 10
n_mean <- 200
n_sd <- 0
n_target <- 100
honesty <- T
covars_fix <- "age"
covars_rand <- "age"
lin <- T
distribution <- "same"
target_dist <- "same"

eps_opt <- list(list(eps_study_m=0.05, eps_study_tau=0.05, eps_study_inter=0.05),
               list(eps_study_m=1, eps_study_tau=0.5, eps_study_inter=0.05))

settings <- expand.grid(moderators=c(1:2),
                        iteration = c(1:100))

#run iterations
res <- data.frame()
for (i in 1:nrow(settings)) {

  #Setup
  moderators <- settings$moderators[i]
  eps_study_m <- eps_opt[[moderators]]$eps_study_m
  eps_study_tau <- eps_opt[[moderators]]$eps_study_tau
  eps_study_inter <- eps_opt[[moderators]]$eps_study_inter
  iteration <- settings$iteration[i]
  seed <- i
  
  set.seed(seed)
  #generate data
  sim_dat <- gen_mdd(K, n_mean, n_sd, n_target, covars_fix, covars_rand, lin,
                     eps_study_m, eps_study_tau, eps_study_inter, 
                     distribution, target_dist)
  train_dat <- sim_dat[["train_dat"]]
  target_dat <- sim_dat[["target_dat"]]
  
  #set up variables
  if ("age2" %in% covars_fix) { 
    covars <- c("sex", "smstat", "weight", "age2", "madrs")
  } else { 
    covars <- c("sex", "smstat", "weight", "age", "madrs")
  }
  
  #update features to include W
  feat <- dplyr::select(train_dat, c(W, S, all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  
  #include counterfactual covariates (swap control and treatment)
  feat_cf <- feat %>%
    mutate(W = as.numeric(W == 0))
  y <- as.numeric(train_dat$Y)
  
  #run bart
  sbart <- dbarts::bart(x.train=as.matrix(feat), y.train=y, x.test=as.matrix(feat_cf), keeptrees=T)
  
  #S-BART credible interval - training - pairwise difference = F
  sb_train <- sbart_ci(train_dat, sbart)
  #S-BART credible interval - target - pairwise difference = F
  sb_target <- sbart_target(K, target_dat, sbart, covars)
  
  #S-BART credible interval - training - pairwise difference = T
  sb_train_p <- sbart_ci(train_dat, sbart, pairwise_diff=T)
  #S-BART credible interval - target - pairwise difference = T
  sb_target_p <- sbart_target(K, target_dat, sbart, covars, pairwise_diff=T)
  
  rm(sbart)
  
  #calculate mean and CIs for individuals and assess accuracy
  sb_res <- c(assess_interval(sb_train, sb_target), 
              avg_train_var = mean(sb_train$var), 
              avg_target_var = mean(sb_target$var))
  sb_res_p <- c(assess_interval(sb_train_p, sb_target_p), 
                avg_train_var = mean(sb_train_p$var), 
                avg_target_var = mean(sb_target_p$var))
  
  #report results
  iter_res <- cbind(sb_res, sb_res_p) %>%
    data.frame() %>%
    rownames_to_column("Metric") %>%
    mutate(moderators = moderators)
  
  res <- bind_rows(res, iter_res)

}
#saveRDS(res, "BART_vartest.RDS")

#check results
#View(res)
colnames(res) <- c("Metric", "S_NonPair", "S_Pair")

#mse - should always be equal
mse <- filter(res, grepl("mse", Metric)==T)
sum(round(mse$S_NonPair, 10) != round(mse$S_Pair, 10)) #0 is what we want to see

#coverage - want close to 95%
res %>%
  filter(grepl("coverage", Metric)==T) %>%
  pivot_longer(cols=c(S_NonPair, S_Pair), names_to = "Method", values_to = "Coverage") %>%
  ggplot(aes(x=Metric, color=Method, y=Coverage)) +
  geom_boxplot()

#variance
res %>%
  filter(grepl("var", Metric)==T) %>%
  pivot_longer(cols=c(S_NonPair, S_Pair), names_to = "Method", values_to = "Variance") %>%
  ggplot(aes(x=Method, y=Variance)) +
  geom_boxplot() +
  facet_wrap(~Metric, scales = "free")
