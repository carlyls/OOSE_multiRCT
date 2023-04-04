### Running Simulations for Methods to Estimate on OOS Individuals ###

library(tidyverse)
library(rsample)
library(grf)
library(fastDummies)
library(nnet)

source("R/MDD_Generation_OOSEst.R")
source("R/MA_OOSEst.R")

# set up data
N <- 100
K <- 10
n_mean <- 200
n_sd <- 0
eps_target <- 0

settings <- expand.grid(eps_combo = c("0.05 0.05 0.05", "0.05 1 0.05", "0.05 1 1",
                                           "1 1 1", "1 3 1"),
                        distribution = c("same", "varying_madrs", "halfdiff_madrsage", "separate_age"),
                        target_dist = c("same", "upweight", "different"),
                        iteration = c(1:100)) %>%
  separate(eps_combo, into=c("eps_study_m", "eps_study_tau", "eps_study_age"), sep=" ") %>%
  mutate(eps_study_m = as.numeric(eps_study_m),
         eps_study_tau = as.numeric(eps_study_tau),
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
results <- compare_oos(N=N, K=K, n_mean=n_mean, n_sd=n_sd, eps_study_m=eps_study_m, 
                       eps_study_tau=eps_study_tau, eps_study_age=eps_study_age,
                       distribution=distribution, target_dist=target_dist, eps_target=eps_target)
save(results, file=paste(paste("results",seed,N,K,n_mean,n_sd,eps_study_m,eps_study_tau,
                               eps_study_age,distribution,target_dist,eps_target,iteration,sep = "_"),
                         ".Rdata",sep=""))

