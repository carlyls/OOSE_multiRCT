### Running Simulations for Methods to Estimate on OOS Individuals ###

library(tidyverse)
library(lme4)
library(rsample)
library(multcomp)
library(MASS)
#library(pimeta)
#library(merTools)

source("R/MDD_Generation_OOSEst.R")
source("R/MA_OOSEst.R")

# set up parameters
K <- 10
n_mean <- 200
n_sd <- 0
n_target <- 100

settings <- expand.grid(eps_study_m = c(0.05, 1),
                        eps_combo = c("0.05 0", "0.5 0.05", "0.5 0.5",
                                      "1 0.05", "1 0.5", "1 1"),
                        distribution = c("same", "varying_madrs", "halfdiff_madrsage", "separate_age"),
                        target_dist = c("same", "different"),
                        iteration = c(1:100)) %>%
  separate(eps_combo, into=c("eps_study_tau", "eps_study_age"), sep=" ") %>%
  mutate(eps_study_tau = as.numeric(eps_study_tau),
         eps_study_age = as.numeric(eps_study_age))

#sets the row of the settings that you will use
i=as.numeric(Sys.getenv('SGE_TASK_ID'))

eps_study_m <- settings$eps_study_m[i]
eps_study_tau <- settings$eps_study_tau[i]
eps_study_age <- settings$eps_study_age[i]
distribution <- settings$distribution[i]
target_dist <- settings$target_dist[i]
iteration <- settings$iteration[i]
seed <- i

#now code
set.seed(seed)
results <- compare_oos(K=K, n_mean=n_mean, n_sd=n_sd, n_target=n_target, eps_study_m=eps_study_m, 
                       eps_study_tau=eps_study_tau, eps_study_age=eps_study_age,
                       distribution=distribution, target_dist=target_dist)
save(results, file=paste(paste("results",seed,iteration,K,n_mean,n_sd,n_target,eps_study_m,eps_study_tau,
                               eps_study_age,distribution,target_dist,sep = "_"),
                         ".Rdata",sep=""))

