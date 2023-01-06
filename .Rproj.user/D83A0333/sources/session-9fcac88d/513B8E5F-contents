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


## Predict on out of sample group
#define new sample grid
N <- 1000
new <- expand.grid(X1=seq(0, 3, by=.5),
                   X2=seq(0, 3, by=.5),
                   X3=seq(0, 3, by=.5),
                   X4=seq(0, 3, by=.5),
                   X5=seq(0, 3, by=.5)) %>%
  slice(rep(1:n(), each=N))

#add study indicator
new$S <- sample(1:K, nrow(new), replace = T)

#make dummy variable
new_feat <- new %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)

#predict CATE
new_feat$tau_hat <- predict(tau_forest, newdata = new_feat)$predictions


## Create confidence intervals
new_sum <- new %>%
  group_by(X1,X2,X3,X4,X5) %>%
  summarise(mean = mean(tau_hat),
            sd = sd(tau_hat)) %>%
  mutate(lower = mean + qt(.025, df=N-1)*sd/sqrt(N),
         upper = mean + qt(.975, df=N-1)*sd/sqrt(N))

