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
source("Simulation_MLOptions.R")


## Simulate dataset
K <- 10
sim_dat <- gen_data(K=K, n_mean=500, n_sd=0, study_mean=0, study_inter_mean=0,
                    study_sd=.5, study_inter_sd=0, scenario="1a")

covars <- grep("^X", names(sim_dat), value=TRUE)
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
new <- expand.grid(X1=seq(0, 3, by=1),
                   X2=seq(0, 3, by=1),
                   X3=seq(0, 3, by=1),
                   X4=seq(0, 3, by=1),
                   X5=seq(0, 3, by=1)) 



#### OPTION 1: COMPLETELY RANDOM ####

#assign study
new_rand <- new %>%
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
new_mem <- new %>%
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




##################################

## Comparing with tau_true
#MSE

#CI coverage