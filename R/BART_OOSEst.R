## Fitting BART to Simulated Data

library(tidyverse)
library(BayesTree)
library(dbarts)

source("R/MDD_Generation_OOSEst.R")
source("R/MA_OOSEst.R")
source("R/Bootstrap_OOSEst.R")
source("R/Comparing_OOSEst.R")

#generate dataset using defaults
K <- 5
sim_dat <- gen_mdd(K=5, n_mean=100)
train_dat <- sim_dat[["train_dat"]]
target_dat <- sim_dat[["target_dat"]]


#bart ###########

#define covariates (ignore age^2 for now)
covars <- c("sex","smstat","weight","age","madrs")  #how to treat S as categorical?
ncovars <- length(covars)

#s-learner #####
feat <- dplyr::select(train_dat, c(W, S, all_of(covars))) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)

feat_cf <- feat %>%
  mutate(W = as.numeric(W == 0)) #swap control and treatment for test data

y <- as.numeric(train_dat$Y)

#run bart
set.seed(2)
sbart <- dbarts::bart(x.train=as.matrix(feat), y.train=y, x.test=as.matrix(feat_cf), keeptrees=T)
#save.image()

#check convergence
#plot(sbart$sigma)

sbart_ci <- function(train_dat, sbart) {
  
  #get means and variance for each person
  w <- train_dat$W
  means <- apply(sbart$yhat.train, 2, mean)
  means_cf <- apply(sbart$yhat.test, 2, mean)
  vars <- apply(sbart$yhat.train, 2, var)
  vars_cf <- apply(sbart$yhat.test, 2, var)
  
  #estimate cate and interval
  mu1 <- w*means + (1-w)*means_cf
  mu0 <- (1-w)*means + w*means_cf
  means_cate <- mu1 - mu0
  vars_cate <- vars + vars_cf
  
  #add to dataframe
  cis <- train_dat %>%
    mutate(mean = means_cate,
           lower = means_cate - 1.96*sqrt(vars_cate),
           upper = means_cate + 1.96*sqrt(vars_cate))
  
  return(cis)
}

# #results
# #get means and variance for each person
# w <- train_dat$W
# means <- apply(sbart$yhat.train, 2, mean)
# means_cf <- apply(sbart$yhat.test, 2, mean)
# vars <- apply(sbart$yhat.train, 2, var)
# vars_cf <- apply(sbart$yhat.test, 2, var)
# #cate
# mu1 <- w*means + (1-w)*means_cf
# mu0 <- (1-w)*means + w*means_cf
# means_cate <- mu1 - mu0
# vars_cate <- vars + vars_cf
# 
# #add to dataframe
# res <- train_dat %>%
#   mutate(mean = means_cate,
#          lower = means_cate - 1.96*sqrt(vars_cate),
#          upper = means_cate + 1.96*sqrt(vars_cate))

## now need to do the same but get intervals for the target data
#use predict function on target data

sbart_target <- function(K, target_dat, sbart)
#set up one row per study for all rows of target data
new_dat <- target_dat %>%
  slice(rep(1:n(), each=K)) %>%
  mutate(S = rep(1:K, nrow(target_dat)))
new_feat <- new_dat %>%
  dplyr::select(c(W, S, all_of(covars))) %>%
  fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
#define counterfactual
new_feat_cf <- new_feat %>%
  mutate(W = as.numeric(W == 0))

#predict on target data
#get means and variance for each person
w_target <- target_dat$W
target_pred <- predict(sbart, new_feat)
target_pred_cf <- predict(sbart, new_feat_cf)

#calculate differences across all posterior draws
#get means and variance for each person
means_cate_target <- lower_cate_target <- upper_cate_target <- c()
for (i in 1:nrow(target_dat)) {
  inds <- (K*(i-1)+1):(K*i)
  pred <- c(target_pred[,inds])
  pred_cf <- c(target_pred_cf[,inds])
  
  mu1_target <- w_target[i]*mean(pred) + (1-w_target[i])*mean(pred_cf)
  mu0_target <- (1-w_target[i])*mean(pred) + w_target[i]*mean(pred_cf)
  mean_cate_target <- mu1_target - mu0_target
  var_cate_target <- var(pred) + var(pred_cf)
  
  means_cate_target <- c(means_cate_target, mean_cate_target)
  lower_cate_target <- c(lower_cate_target, 
                         mean_cate_target - 1.96*sqrt(var_cate_target))
  upper_cate_target <- c(upper_cate_target, 
                         mean_cate_target + 1.96*sqrt(var_cate_target))
}

target_res <- target_dat %>%
  mutate(mean = means_cate_target,
         lower = lower_cate_target,
         upper = upper_cate_target)
