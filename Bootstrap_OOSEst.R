### Bootstrapping Confidence Intervals for Out of Sample CATE Estimation ###

library(tidyverse)
library(causalToolbox)
library(rsample)
library(rpart)
library(ranger)
library(glmnet)
library(grf)
library(fastDummies)
library(lme4)
library(nnet)

source("Comparing_methods_functions.R")
source("MDD_Simulation_OOSEst.R")


## Simulate dataset
sim_dat <- gen_mdd(K=6, n_mean=200, n_sd=0, 
                   scenario="1a", distribution="same")

covars <- c("sex", "smstat", "weight", "age", "madrs")
feat <- select(sim_dat, c(S,all_of(covars))) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)

feat_nostudy <- select(sim_dat, all_of(covars))
tr <- sim_dat$W
y <- sim_dat$Y
tau_true <- sim_dat$tau


## Fit models (for now, causal forest with pooling with trial indicator)
tau_forest <- causal_forest(X=feat, Y=y, W=tr, num.threads=3, honesty=F, num.trees=1000)
tau_hat <- predict(tau_forest)$predictions


## Define new sample grid
N <- 100
test_dat <- expand.grid(W = c(0, 1),
                        sex = c(0, 1),
                        smstat = c(0, 1),
                        weight = seq(45, 130, by=10),
                        age = seq(18, 75, by=10),
                        madrs = seq(20, 50, by=5))



#### OPTION 1: COMPLETELY RANDOM ####

#assign study
new_rand <- test_dat%>%
  slice(rep(1:n(), each=N)) %>%
  mutate(S = sample(1:K, nrow(new)*N, replace = T)) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)

#predict CATE
new_rand$tau_hat <- predict(tau_forest, newdata = new_rand)$predictions

#create confidence intervals
new_rand_sum <- new_rand %>%
  group_by(X1,X2,X3,X4,X5) %>%
  summarise(mean = mean(tau_hat),
            sd = sd(tau_hat)) %>%
  mutate(lower = mean + qt(.025, df=N-1)*sd/sqrt(N),
         upper = mean + qt(.975, df=N-1)*sd/sqrt(N))



#### OPTION 2: STUDY MEMBERSHIP MODEL ####

#should I include the outcome? Then will need to simulate Y in the new dataset
mem_mod <- multinom(S ~ X1 + X2 + X3 + X4 + X5, data=sim_dat)
summary(mem_mod)

#can check this model using a few things
#round(fitted(mem_mod), 2) #looks at probabilities in each class
#check accuracy
#preds <- predict(mem_mod, newdata = sim_dat, "class")
#tab <- table(sim_dat$S, preds)
#round((sum(diag(tab))/sum(tab))*100,2)

#define probabilities
mem_probs <- predict(mem_mod, newdata = new, type = "probs")
S_mem <- c()
for (i in 1:nrow(mem_probs)) {
  S_mem <- c(S_mem, sample(1:K, N, replace=T, prob=mem_probs[i,]))
}

#assign study
new_mem <- test_dat%>%
  slice(rep(1:n(), each=N)) %>%
  mutate(S = S_mem) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)

#predict CATE
new_mem$tau_hat <- predict(tau_forest, newdata = new_mem)$predictions

#create confidence intervals
new_mem_sum <- new_mem %>%
  group_by(X1,X2,X3,X4,X5) %>%
  summarise(mean = mean(tau_hat),
            sd = sd(tau_hat)) %>%
  mutate(lower = mean + qt(.025, df=N-1)*sd/sqrt(N),
         upper = mean + qt(.975, df=N-1)*sd/sqrt(N))



#### OPTION 3: WITHIN-FOREST ####

#default method: https://grf-labs.github.io/grf/REFERENCE.html#missing-values

#assign study
#we don't need to replicate because we will get the same prediction each time
new_default <- test_dat%>%
  mutate(S = factor(NA, levels=1:10)) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T, ignore_na = T)

#predict CATE
cate_default <- predict(tau_forest, newdata = new_default, estimate.variance = T)
new_default$mean <- cate_default$predictions 
new_default$sd <- sqrt(cate_default$variance.estimates)

#create confidence intervals
new_default_sum <- new_default %>%
  mutate(lower = mean + qt(.025, df=N-1)*sd/sqrt(N),
         upper = mean + qt(.975, df=N-1)*sd/sqrt(N))



#### OPTION 4: WITHIN-FOREST RANDOM SAMPLING ####



#### OPTION 5: GRADE OF MEMBERSHIP MODEL ####




##################################

## Comparing with tau_true
#MSE

#CI coverage