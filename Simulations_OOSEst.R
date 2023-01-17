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
test_dat <- expand.grid(W = c(0, 1),
                        sex = c(0, 1),
                        smstat = c(0, 1),
                        weight = seq(45, 130, by=10),
                        age = seq(18, 75, by=10),
                        madrs = seq(20, 50, by=5))

settings <- expand.grid(n_mean = c(200, 500),
                        scenario = c("linear", "nonlinear"),
                        distribution = c("same", "varying_madrs", 
                                         "halfdiff_madrsage", "separate_age"),
                        test_scenario = c("random"),
                        iteration = c(1:500))

#sets the row of the settings that you will use
i=as.numeric(Sys.getenv('SGE_TASK_ID'))

n_mean <- settings$n_mean[i]
scenario <- settings$scenario[i]
distribution <- settings$distribution[i]
test_scenario <- settings$test_scenario[i]
iteration <- settings$iteration[i]
seed <- i

#now code
set.seed(seed)
results <- compare_oos(N=N, K=K, n_mean=n_mean, n_sd=n_sd, scenario=scenario, 
                       distribution=distribution, test_dat=test_dat, test_scenario=test_scenario)
save(results, file=paste(paste("results",seed,N,K,n_mean,n_sd,scenario,distribution,test_scenario,
                               sep = "_"),".Rdata",sep=""))

