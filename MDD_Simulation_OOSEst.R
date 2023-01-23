### Simulation Approach for OOS Estimation ###
### Based on MDD Data ###

library(tidyverse)
library(rsample)
library(truncnorm)

#interior function
add_agemadrs <- function(dat, n, k, distribution) {
  
  #add age and madrs
  if (distribution == "same") {
    dat <- dat %>%
      mutate(age = rtruncnorm(n=n, a=18, b=75, mean=45, sd=10),
             madrs = rtruncnorm(n=n, a=26, b=60, mean=31, sd=4.1))
    
  } else if (distribution == "varying_madrs") {
    dat <- dat %>%
      mutate(age = rtruncnorm(n=n, a=18, b=75, mean=45, sd=10),
             madrs = rtruncnorm(n=n, a=26-k*1.5, b=60-k*1.5, mean=31-k*1.5, sd=4.1))
    
  } else if (distribution == "halfdiff_madrsage") {
    if (k%%2 == 0 ) {
      dat <- dat %>%
        mutate(age = rtruncnorm(n=n, a=30, b=75, mean=50, sd=10),
               madrs = rtruncnorm(n=n, a=30, b=60, mean=40, sd=4.1))
    } else {
      dat <- dat %>%
        mutate(age = rtruncnorm(n=n, a=18, b=75, mean=45, sd=10),
               madrs = rtruncnorm(n=n, a=26, b=60, mean=31, sd=4.1))
    }
    
  } else if (distribution == "separate_age") {
    ages <- seq(18,75,by=(75-18)/10)
    dat <- dat %>%
      mutate(age = runif(n=n, min=ages[k], max=ages[k+1]),
             madrs = rtruncnorm(n=n, a=26, b=60, mean=31, sd=4.1))
  }
  
  return(dat)
}

#main function
gen_mdd <- function (K=6, n_mean=200, n_sd=0, scenario="linear", 
                      distribution="same", test_dist="same") {
  
  #training data
  train_dat <- data.frame()
  n_study <- floor(rnorm(K, mean=n_mean, sd=n_sd))
  
  for (k in 1:K) {
    n <- n_study[k]
    
    #sample covariates
    dat <- data.frame(
      #age = rtruncnorm(n=n, a=18, b=75, mean=45, sd=10),
      sex = rbinom(n=n, size=1, prob=.65),
      smstat = rbinom(n=n, size=1, prob=.3),
      weight = rtruncnorm(n=n, a=45, b=140, mean=80, sd=15),
      #madrs = rtruncnorm(n=n, a=26, b=60, mean=31, sd=4.1),
      W = rbinom(n=n, size=1, prob=.5),
      S = rep(k, n),
      id = seq(1, n),
      eps = rnorm(n, mean=0, sd=.05),
      eps_study = rnorm(1, mean=0, sd=.05)
    )
    
    dat <- add_agemadrs(dat, n, k, distribution) %>%
      mutate(weight = round(weight, 2),
             age = round(age, 2),
             madrs = round(madrs, 0))
    
    train_dat <- bind_rows(train_dat, dat)
  }
  
  #testing data
  if (test_dist == "same") {
    test_dat <- train_dat[sample(nrow(train_dat), 100),] %>%
      mutate(S = NA,
             eps_study = rnorm(100, mean=0, sd=.05))
  } else if (test_dist == "upweight") {
    train_weight <- train_dat %>%
      mutate(weight = ifelse(S %in% c(3, 5), 3, 1))
    test_dat <- train_weight[sample(nrow(train_weight), 100, prob=train_weight$weight),] %>%
      mutate(S = NA,
             eps_study = rnorm(100, mean=0, sd=.05)) %>%
      select(-weight)
  } else if (test_dist == "different") {
    test_dat <- data.frame(
      sex = rbinom(n=100, size=1, prob=.65),
      smstat = rbinom(n=100, size=1, prob=.3),
      weight = rtruncnorm(n=100, a=45, b=140, mean=80, sd=15),
      W = rbinom(n=100, size=1, prob=.5),
      S = NA,
      id = seq(1,100),
      eps = rnorm(n=100, mean=0, sd=.05),
      eps_study = rnorm(n=100, mean=0, sd=.05),
      age = rtruncnorm(n=100, a=18, b=40, mean=25, sd=10),
      madrs = rtruncnorm(n=100, a=16, b=40, mean=25, sd=4.1)
    )
  }
  
  #tau and Y
  if (scenario == "simple") {
    train_dat <- train_dat %>% 
      mutate(m = -0.02*age - 0.87*madrs - 0.15*sex,
             tau = -8.5 + 0.07*age + 0.20*madrs + eps_study)
    test_dat <- test_dat %>% 
      mutate(m = -0.02*age - 0.87*madrs - 0.15*sex,
             tau = -8.5 + 0.07*age + 0.20*madrs + eps_study)
  }
  
  if (scenario == "linear") {
    study_main <- rnorm(K, mean=-10, sd=-0.5)
    study_inter <- rnorm(K, mean=0.25, sd=0.05)
    study_tau <- rnorm(K, mean=3, sd=0.01)
    train_dat <- train_dat %>% 
      mutate(study_main = study_main[S],
             study_inter = study_inter[S],
             study_tau = study_tau[S],
             m = 10.7 + study_main - 0.02*age - 0.87*madrs -
               0.15*sex + study_inter*madrs,
             tau = -8.5 + 0.07*age + 0.20*madrs + study_tau) %>%
      select(-study_main, -study_inter, -study_tau)
    
    test_main <- rnorm(1, mean=-10, sd=-0.5)
    test_inter <- rnorm(1, mean=0.25, sd=0.05)
    test_tau <- rnorm(1, mean=3, sd=0.01)
    test_dat <- test_dat %>%
      mutate(test_main = test_main + rnorm(nrow(test_dat), mean=0, sd=0.01),
             test_inter = test_inter + rnorm(nrow(test_dat), mean=0, sd=0.01),
             test_tau = test_tau + rnorm(nrow(test_dat), mean=0, sd=0.01),
             m = 10.7 + test_main - 0.02*age - 0.87*madrs -
               0.15*sex + test_inter*madrs,
             tau = -8.5 + 0.07*age + 0.20*madrs + test_tau) %>%
      select(-test_main, -test_inter, -test_tau)
  }
  
  train_dat <- train_dat %>%
    mutate(Y = m + W*tau + eps,
           S = factor(S)) %>%
    select(S, id, W, sex, smstat, weight, age, madrs, Y, tau)
  
  test_dat <- test_dat %>%
    mutate(Y = m + W*tau + eps) %>%
    select(S, id, W, sex, smstat, weight, age, madrs, Y, tau)
  return(list(train_dat=train_dat, test_dat=test_dat))
}












## other tau options
if (scenario == "linear") {
  study_main <- runif(K, min=-14, max=-7)
  study_inter <- runif(K, min=0.1, max=0.5)
  study_tau <- runif(K, min=2.5, max=3.5)
  train_dat <- train_dat %>% 
    mutate(study_main = study_main[S],
           study_inter = study_inter[S],
           study_tau = study_tau[S],
           m = 10.7 + study_main - 0.02*age - 0.87*madrs -
             0.15*sex + study_inter*madrs,
           tau = -8.5 + 0.07*age + 0.20*madrs + study_tau) %>%
    select(-study_main, -study_inter, -study_tau)
  
  test_main <- runif(1, min=-14, max=-7)
  test_inter <- runif(1, min=0.1, max=0.5)
  test_tau <- runif(1, min=2.5, max=3.5)
  test_dat <- test_dat %>%
    mutate(test_main = test_main + rnorm(nrow(test_dat), mean=0, sd=0.1),
           test_inter = test_inter + rnorm(nrow(test_dat), mean=0, sd=0.1),
           test_tau = test_tau + rnorm(nrow(test_dat), mean=0, sd=0.1),
           m = 10.7 + test_main - 0.02*age - 0.87*madrs -
             0.15*sex + test_inter*madrs,
           tau = -8.5 + 0.07*age + 0.20*madrs + test_tau) %>%
    select(-test_main, -test_inter, -test_tau)
}

if (scenario == "nonlinear") {
  study_tau <- runif(K, min=2.5, max=3.5)
  train_dat <- train_dat %>% 
    mutate(study_tau = study_tau[S],
           m = 0,
           tau = (study_tau/(1+exp(-1/12*age)))*(study_tau/(1+exp(-12*madrs)))) %>%
    select(-study_tau)
  
  test_tau <- runif(1, min=2.5, max=3.5)
  test_dat <- test_dat %>%
    mutate(test_tau = test_tau + rnorm(nrow(test_dat), mean=0, sd=0.1),
           m = 0,
           tau = (test_tau/(1+exp(-1/12*age)))*(test_tau/(1+exp(-12*madrs)))) %>%
    select(-test_tau)
}

