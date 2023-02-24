## Trying Methods for OOS Estimation in Random Effects Meta-Analysis ##

library(tidyverse)
library(rsample)
library(grf)
library(fastDummies)
library(nnet)
library(lme4)

source("MDD_Simulation_OOSEst.R")

#simulate data
sim_dat <- gen_mdd(K=6, n_mean=200, n_sd=0, eps_study_m=0.05, eps_study_tau=0.01, 
                   distribution="same", target_dist="same")
train_dat <- sim_dat[["train_dat"]]
target_dat <- sim_dat[["target_dat"]]

covars <- c("sex", "smstat", "weight", "age", "madrs")

#fit typical model (no heterogeneity)
mod_avg <- lmer(Y ~ W + sex + age + madrs + 
                  (1 + W | S), data=train_dat)
summary(mod_avg)

#prediction interval for treatment effect
mu <- fixef(mod_avg)["W"]
se <- summary(mod_avg)$coef[,2]["W"]
tau2 <- 1.004e-05 #how do I pull this out with code?

pred_lower <- round(mu + qt(.025, 4)*sqrt(tau2+se^2),2)
pred_upper <- round(mu - qt(.025, 4)*sqrt(tau2+se^2),2)
print(paste0("(",pred_lower,", ",pred_upper,")"))

#add in heterogeneity in fixed effects
mod_het <- lmer(Y ~ W + sex + age + madrs + W:age + W:madrs +
                  (1 + W | S), data=train_dat)
summary(mod_het)

#add in heterogeneity in fixed and random effects
mod_ranhet <- lmer(Y ~ W + sex + age + madrs + W:age +
                  (1 + W + W:age | S), data=train_dat)
summary(mod_ranhet)

