### Running Simulations for Methods to Estimate on OOS Individuals ###

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
n_target <- 100
honesty <- T

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
               list(covars_fix=c("age2", "sex"), covars_rand=c("age2"), 
                    lin=T, eps_study_m=1, eps_study_tau=0.5, eps_study_inter=0.05),
               list(covars_fix=c("age2", "sex"), covars_rand=c("age2"), 
                    lin=T, eps_study_m=1, eps_study_tau=0.5, eps_study_inter=0.5),
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
                        target_dist = c("same", "different"),
                        iteration = c(1:100))

#set row of settings and define parameters
i=as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

moderators <- settings$moderators[i]
covars_fix <- mods[[moderators]]$covars_fix
covars_rand <- mods[[moderators]]$covars_rand
eps_study_m <- mods[[moderators]]$eps_study_m
eps_study_tau <- mods[[moderators]]$eps_study_tau
eps_study_inter <- mods[[moderators]]$eps_study_inter
lin <- mods[[moderators]]$lin

distribution <- settings$distribution[i]
target_dist <- settings$target_dist[i]
iteration <- settings$iteration[i]
seed <- i

#run main function
set.seed(seed)
results <- compare_oos(K=K, n_mean=n_mean, n_sd=n_sd, n_target=n_target, covars_fix=covars_fix,
                       covars_rand=covars_rand, lin=lin, eps_study_m=eps_study_m, eps_study_tau=eps_study_tau, 
                       eps_study_inter=eps_study_inter, distribution=distribution, target_dist=target_dist)
save(results, file=paste(paste("results",seed,iteration,K,n_mean,n_sd,n_target,"modset",moderators,
                               lin,eps_study_m,eps_study_tau,mean(eps_study_inter),distribution,
                               target_dist,sep = "_"),
                         ".Rdata",sep=""))

