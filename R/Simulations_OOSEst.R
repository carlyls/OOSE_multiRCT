### Running Simulations for Methods to Estimate on OOS Individuals ###

library(tidyverse)
library(rsample)
library(grf)
library(fastDummies)
library(nnet)

source("R/Bootstrap_OOSEst.R")
source("R/MDD_Generation_OOSEst.R")

# set up data
N <- 100
K <- 6
n_sd <- 0
#covars <- c("sex", "smstat", "weight", "age", "madrs")

settings <- expand.grid(n_mean = c(200, 500),
                        eps_combo = c("0.05 0.01", "1 0.01", "1 0.05",
                                           "1 1", "2 1"),
                        distribution = c("same", "varying_madrs", "halfdiff_madrsage", "separate_age"),
                        target_dist = c("same", "upweight", "different"),
                        iteration = c(1:500)) %>%
  separate(eps_combo, into=c("eps_study_m", "eps_study_tau"), sep=" ") %>%
  mutate(eps_study_m = as.numeric(eps_study_m),
         eps_study_tau = as.numeric(eps_study_tau))

#sets the row of the settings that you will use
i=as.numeric(Sys.getenv('SGE_TASK_ID'))

n_mean <- settings$n_mean[i]
eps_study_m <- settings$eps_study_m[i]
eps_study_tau <- settings$eps_study_tau[i]
distribution <- settings$distribution[i]
target_dist <- settings$target_dist[i]
iteration <- settings$iteration[i]
seed <- i

#now code
set.seed(seed)
results <- compare_oos(N=N, K=K, n_mean=n_mean, n_sd=n_sd, eps_study_m=eps_study_m, 
                       eps_study_tau=eps_study_tau, distribution=distribution, target_dist=target_dist)
save(results, file=paste(paste("results",seed,N,K,n_mean,n_sd,eps_study_sd,
                               eps_study_tau,distribution,target_dist,sep = "_"),".Rdata",sep=""))

