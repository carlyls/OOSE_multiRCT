## Fitting BART to Simulated Data

# library(tidyverse)
# library(BayesTree)
# library(dbarts)

#generate dataset using defaults
# K <- 5
# sim_dat <- gen_mdd(K=5, n_mean=100)
# train_dat <- sim_dat[["train_dat"]]
# target_dat <- sim_dat[["target_dat"]]


#bart ###########

#define covariates (ignore age^2 for now)
# covars <- c("sex","smstat","weight","age","madrs")  #how to treat S as categorical?
# ncovars <- length(covars)

#t-learner #####
#define two model setups
tlearn_setup <- function(train_dat, covars, w) {
  
  mod_dat <- train_dat %>% filter(W == w)
  feat <- mod_dat %>%
    dplyr::select(c(S, all_of(covars))) %>%
    fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
  y <- as.numeric(mod_dat$Y)
  
  return(list(mod_dat=mod_dat, feat=feat, y=y))
}

# #m1
# m1_setup <- tlearn_setup(train_dat, covars, w=1)
# tbart1 <- dbarts::bart(x.train = as.matrix(m1_setup[["feat"]]), y.train = m1_setup[["y"]],
#                        x.test = as.matrix(m0_setup[["feat"]]), keeptrees = T)
# 
# #m0
# m0_setup <- tlearn_setup(train_dat, covars, w=0)
# tbart0 <- dbarts::bart(x.train = as.matrix(m0_setup[["feat"]]), y.train = m0_setup[["y"]],
#                        x.test = as.matrix(m1_setup[["feat"]]), keeptrees = T)

#training data
tbart_ci <- function(train_dat, tbart1, tbart0) {
  
  #treated people: Y1 - m0(X1), so use training outcome from tbart1 and testing outcome from tbart0
  #control people: m1(X0) - Y0, so use testing outcome from tbart1 and training outcome from tbart0
  cate1 <- tbart1$yhat.train.mean - tbart0$yhat.test.mean
  cate0 <- tbart1$yhat.test.mean - tbart0$yhat.train.mean
  var1 <- apply(tbart1$yhat.train, 2, var) + apply(tbart0$yhat.test, 2, var)
  var0 <- apply(tbart1$yhat.test, 2, var) + apply(tbart0$yhat.train, 2, var)
  
  #add to dataframe
  cis <- filter(train_dat, W == 1) %>%
    bind_rows(filter(train_dat, W == 0)) %>%
    mutate(mean = c(cate1, cate0),
           lower = c(cate1 - 1.96*sqrt(var1), cate0 - 1.96*sqrt(var0)),
           upper = c(cate1 + 1.96*sqrt(var1), cate0 + 1.96*sqrt(var0)))
  
  #reorder to match target data
  cis_ord <- train_dat %>%
    left_join(cis, by = c(names(train_dat)))
  
  return(cis_ord)
}

#target data
tbart_target <- function(K, target_dat, tbart1, tbart0, covars) {
  
  #set up one row per study for all rows of target data
  new_dat <- target_dat %>%
    slice(rep(1:n(), each=K)) %>%
    mutate(S = rep(1:K, nrow(target_dat)))
  new_feat1 <- tlearn_setup(new_dat, covars, w=1)[["feat"]]
  new_feat0 <- tlearn_setup(new_dat, covars, w=0)[["feat"]]
  
  #predict on target data
  target_pred1 <- predict(tbart1, new_feat1) #treated individual Y1
  target_pred1_cf <- predict(tbart0, new_feat1) #treated individual Y0
  target_pred0_cf <- predict(tbart1, new_feat0) #control individual Y1
  target_pred0 <- predict(tbart0, new_feat0) #control individual Y0
  y1 <- cbind(target_pred1, target_pred0_cf)
  y0 <- cbind(target_pred1_cf, target_pred0)
  
  #calculate differences across all posterior draws
  #get means and variance for each person
  means_cate_target <- lower_cate_target <- upper_cate_target <- c()
  for (i in 1:nrow(target_dat)) {
    inds <- (K*(i-1)+1):(K*i)
    pred1 <- c(y1[,inds])
    pred0 <- c(y0[,inds])
    
    mean_cate_target <- mean(pred1) - mean(pred0)
    var_cate_target <- var(pred1) + var(pred0)
    
    means_cate_target <- c(means_cate_target, mean_cate_target)
    lower_cate_target <- c(lower_cate_target, 
                           mean_cate_target - 1.96*sqrt(var_cate_target))
    upper_cate_target <- c(upper_cate_target, 
                           mean_cate_target + 1.96*sqrt(var_cate_target))
  }
  
  #add to target data
  cis <- filter(target_dat, W == 1) %>%
    bind_rows(filter(target_dat, W == 0)) %>%
    mutate(mean = means_cate_target,
           lower = lower_cate_target,
           upper = upper_cate_target)
  
  #reorder to match target data
  cis_ord <- target_dat %>%
    left_join(cis, by = c(names(target_dat)))
  
  return(cis_ord)
}



#s-learner #####
# feat <- dplyr::select(train_dat, c(W, S, all_of(covars))) %>%
#   fastDummies::dummy_cols(select_columns="S", remove_selected_columns=T)
# 
# feat_cf <- feat %>%
#   mutate(W = as.numeric(W == 0)) #swap control and treatment for test data
# 
# y <- as.numeric(train_dat$Y)

#run bart
# set.seed(2)
# sbart <- dbarts::bart(x.train=as.matrix(feat), y.train=y, x.test=as.matrix(feat_cf), keeptrees=T)
#save.image()

#check convergence
#plot(sbart$sigma)

## S-learner for training data
sbart_ci <- function(train_dat, sbart, pairwise_diff=F) {
  
  if (pairwise_diff == F) {
    
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
             var = vars_cate,
             lower = means_cate - 1.96*sqrt(vars_cate),
             upper = means_cate + 1.96*sqrt(vars_cate))
    
  } else {
    
    #get pairwise differences
    w <- train_dat$W
    w_fac <- ifelse(w == 1, 1, -1)
    pairwise_diffs <- sweep(sbart$yhat.train - sbart$yhat.test, #some of these are the wrong order
                            2, w_fac, "*") #multiply by -1 if the person was in control group
    
    #estimate mean and variance
    means_cate <- apply(pairwise_diffs, 2, mean)
    lower_cate <- apply(pairwise_diffs, 2, quantile, probs=.025)
    upper_cate <- apply(pairwise_diffs, 2, quantile, probs=.975)
    
    #add to dataframe
    cis <- train_dat %>%
      mutate(mean = means_cate,
             lower = lower_cate,
             upper = upper_cate)
  }
  
  return(cis)
}

## S-learner for target data
#between and within study variance
sbart_target <- function(K, target_dat, sbart, covars) {
  
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
  target_pred <- predict(sbart, new_feat)
  target_pred_cf <- predict(sbart, new_feat_cf)
  
  #get mean and variance across posterior for each person-study combination
  new_dat <- new_dat %>%
    mutate(mean_obs = apply(target_pred, 2, mean),
           mean_cf = apply(target_pred_cf, 2, mean),
           var_obs = apply(target_pred, 2, var),
           var_cf = apply(target_pred_cf, 2, var),
           w_fac_new = ifelse(W == 0, -1, 1),
           mean_cate = w_fac_new*(mean_obs - mean_cf),
           var_cate = var_obs + var_cf)
  
  #average across studies
  cis <- new_dat %>%
    group_by(W, sex, smstat, weight, age, age2, madrs, Y, tau) %>%
    summarise(mean = mean(mean_cate),
              var_within = mean(var_cate),
              var_btwn = var(mean_cate),
              n_K = n()) %>%
    mutate(var_tot = var_within + var_btwn,
           sd = sqrt(var_tot),
           lower = mean - qt(.975, K-2)*sd,
           upper = mean + qt(.975, K-2)*sd)
  
  #reorder to match target data
  cis_ord <- target_dat %>%
    left_join(cis, by = c("W","sex","smstat","weight","age",
                          "age2","madrs","Y","tau"))
  
  return(cis_ord)
}


#original approach: do not decompose variance
sbart_rand <- function(K, target_dat, sbart, covars, pairwise_diff=F) {
  
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
  target_pred <- predict(sbart, new_feat)
  target_pred_cf <- predict(sbart, new_feat_cf)
  
  if (pairwise_diff == F) {
    
    w_target <- target_dat$W
    #calculate differences across all posterior draws
    #get means and variance for each person
    means_cate_target <- vars_cate_target <- lower_cate_target <- upper_cate_target <- c()
    for (i in 1:nrow(target_dat)) {
      inds <- (K*(i-1)+1):(K*i)
      pred <- c(target_pred[,inds])
      pred_cf <- c(target_pred_cf[,inds])
      
      mu1_target <- w_target[i]*mean(pred) + (1-w_target[i])*mean(pred_cf)
      mu0_target <- (1-w_target[i])*mean(pred) + w_target[i]*mean(pred_cf)
      mean_cate_target <- mu1_target - mu0_target
      var_cate_target <- var(pred) + var(pred_cf)
      
      means_cate_target <- c(means_cate_target, mean_cate_target)
      vars_cate_target <- c(vars_cate_target, var_cate_target)
      lower_cate_target <- c(lower_cate_target, 
                             mean_cate_target - qt(.975, K-2)*sqrt(var_cate_target))
      upper_cate_target <- c(upper_cate_target, 
                             mean_cate_target + qt(.975, K-2)*sqrt(var_cate_target))
    }
    
  } else {
    
    w_target_fac <- ifelse(new_dat$W == 1, 1, -1)
    pairwise_diffs <- sweep(target_pred - target_pred_cf, #some of these are the wrong order
                            2, w_target_fac, "*") #multiply by -1 if the person was in control group
    
    means_cate_target <- lower_cate_target <- upper_cate_target <- c()
    for (i in 1:nrow(target_dat)) {
      inds <- (K*(i-1)+1):(K*i)
      pairwise_diffs_inds <- c(pairwise_diffs[,inds])
      
      means_cate_target <- c(means_cate_target, mean(pairwise_diffs_inds))
      lower_cate_target <- c(lower_cate_target, 
                             quantile(pairwise_diffs_inds, probs=.025))
      upper_cate_target <- c(upper_cate_target, 
                             quantile(pairwise_diffs_inds, probs=.975))
    }
    
  }

  #add to target data
  cis <- target_dat %>%
    mutate(mean = means_cate_target,
           lower = lower_cate_target,
           upper = upper_cate_target)
  
  return(cis)
}

 