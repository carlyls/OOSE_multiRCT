
library(tidyverse)
library(lme4)
library(rsample)
library(multcomp)
library(MASS)
library(grf)
library(dbarts)
library(fastDummies)

source("R/MDD_Generation_OOSEst.R")
source("R/MA_OOSEst.R")
source("R/CF_OOSEst.R")
source("R/BART_OOSEst.R")
source("R/Comparing_OOSEst.R")

# set up parameters
K <- 10
n_mean <- 500
n_sd <- 0
target_sample <- readRDS("Data/target_sample.RDS")

mods <- list(list(covars_fix="age", covars_rand="age", lin=T,
                  eps_study_m=0.05, eps_study_tau=0.05, eps_study_inter=0.05),
             list(covars_fix="age", covars_rand="age", lin=T,
                  eps_study_m=1, eps_study_tau=0.05, eps_study_inter=0.05),
             list(covars_fix="age", covars_rand="age", lin=T,
                  eps_study_m=1, eps_study_tau=0.5, eps_study_inter=0.05),
             list(covars_fix="age", covars_rand="age", lin=T,
                  eps_study_m=1, eps_study_tau=1, eps_study_inter=0.5),
             list(covars_fix=c("age", "madrs"), covars_rand=c("age", "madrs"), 
                  lin=T, eps_study_m=0.05, eps_study_tau=0.05, eps_study_inter=c(0.05,0.05)),
             list(covars_fix=c("age", "madrs"), covars_rand=c("age", "madrs"), 
                  lin=T, eps_study_m=1, eps_study_tau=0.5, eps_study_inter=c(0.5,0.05)),
             list(covars_fix="age", covars_rand="age", lin=F,
                  eps_study_m=0.05, eps_study_tau=0.05, eps_study_inter=0.05),
             list(covars_fix="age", covars_rand="age", lin=F,
                  eps_study_m=1, eps_study_tau=0.05, eps_study_inter=0.05),
             list(covars_fix="age", covars_rand="age", lin=F,
                  eps_study_m=1, eps_study_tau=0.5, eps_study_inter=0.05),
             list(covars_fix="age", covars_rand="age", lin=F,
                  eps_study_m=1, eps_study_tau=1, eps_study_inter=0.5))

settings <- expand.grid(moderators = c(1:length(mods)),
                        distribution = c("same", "varying_madrs", "separate_age"),
                        iteration = c(1:100))

#set row of settings and define parameters
i=7

moderators <- settings$moderators[i]
covars_fix <- mods[[moderators]]$covars_fix
covars_rand <- mods[[moderators]]$covars_rand
eps_study_m <- mods[[moderators]]$eps_study_m
eps_study_tau <- mods[[moderators]]$eps_study_tau
eps_study_inter <- mods[[moderators]]$eps_study_inter
lin <- mods[[moderators]]$lin

distribution <- settings$distribution[i]
iteration <- settings$iteration[i]
seed <- i

#run main function
set.seed(seed)

## Simulate training and target (OOS) data
sim_dat <- gen_mdd(K, n_mean, n_sd, covars_fix, covars_rand, lin,
                   eps_study_m, eps_study_tau, eps_study_inter, 
                   distribution, target_sample)
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

#save target dataset individual results
target_res <- target_metrics(manual_target, "MA") 

##### LOOK THROUGH RESULTS
View(target_res)
ggplot(target_res, aes(x=age, y=coverage)) + geom_point()
ggplot(target_res, aes(x=age, y=abs_bias)) + geom_point()
ggplot(target_res, aes(x=age, y=length)) + geom_point()

#note: looks much more random with madrs

#true vs est
target_res %>% 
  rename(True=tau, Estimated=mean) %>%
  pivot_longer(cols=c(True,Estimated),names_to="Type", values_to="CATE") %>% 
  ggplot(aes(x=age, y=CATE)) + 
  geom_errorbar(aes(ymax=upper, ymin=lower), color="lightgrey") +
  geom_point(aes(color=Type))
