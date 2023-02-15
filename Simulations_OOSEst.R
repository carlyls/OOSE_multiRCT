### Running Simulations for Methods to Estimate on OOS Individuals ###

library(tidyverse)
library(rsample)
library(grf)
library(fastDummies)
library(nnet)

source("Bootstrap_OOSEst.R")
source("MDD_Simulation_OOSEst.R")

# set up data
N <- 100
K <- 6
n_sd <- 0
#covars <- c("sex", "smstat", "weight", "age", "madrs")

settings <- expand.grid(n_mean = c(200, 500),
                        scenario_combo = c("simple 0.01", "simple 0.05", "simple 1",
                                           "simple 3", "linear 0"),
                        distribution = c("same", "varying_madrs", "halfdiff_madrsage", "separate_age"),
                        test_dist = c("same", "upweight", "different"),
                        iteration = c(1:500)) %>%
  separate(scenario_combo, into=c("scenario", "eps_study_sd"), sep=" ") %>%
  mutate(eps_study_sd = as.numeric(eps_study_sd))

#sets the row of the settings that you will use
i=as.numeric(Sys.getenv('SGE_TASK_ID'))

n_mean <- settings$n_mean[i]
eps_study_sd <- settings$eps_study_sd[i]
scenario <- settings$scenario[i]
distribution <- settings$distribution[i]
test_dist <- settings$test_dist[i]
iteration <- settings$iteration[i]
seed <- i

#now code
set.seed(seed)
results <- compare_oos(N=N, K=K, n_mean=n_mean, n_sd=n_sd, eps_study_sd=eps_study_sd, 
                       scenario=scenario, distribution=distribution, test_dist=test_dist)
save(results, file=paste(paste("results",seed,N,K,n_mean,n_sd,eps_study_sd,
                               scenario,distribution,test_dist,sep = "_"),".Rdata",sep=""))

