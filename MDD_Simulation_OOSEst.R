### Simulation Approach for OOS Estimation ###
### Based on MDD Data ###

library(tidyverse)
library(rsample)
library(truncnorm)

#interior function
add_study <- function(dat, K, n, distribution) {
  if (distribution == "same") {
    dat <- dat %>%
      mutate(S = sample(rep(1:K, n), replace=F))
    
    
    ## DEFINE CATEGORICAL MODEL FOR S and assign based on probabilities for the below options
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
  
  eps_m <- rnorm(K, mean=0, sd=eps_study_sd)
  eps_tau <- rnorm(K, mean=0, sd=eps_study_sd)
  
  dat <- dat %>%
    arrange(S) %>%
    mutate(W = rbinom(n=n*K, size=1, prob=.5), ### FIX THIS SO THAT IT'S p=0.5 WITHIN EACH STUDY
           eps_m = eps_m[S],
           eps_tau = eps_tau[S])
  
  return(dat)
}

#main function
gen_mdd <- function (K=6, n=200, eps_study_m=0.05, eps_study_tau=0.01, 
                      distribution="same", test_dist="same") {
  
  #training data
  train_dat <- data.frame()
  
  #simulate joint distribution
  mu <- c(age=44.8971, sex=0.6784, smstat=0.3043, weight=79.0253, madrs_base=31.4088)
  Sigma <- data.frame(age=c(165.6471, 0.2448, -0.5180, 1.6408, -0.9666),
                      sex=c(0.2448, 0.2183, -0.0218, -1.9030, 0.1380),
                      smstat=c(-0.5180, -0.0218, 0.2118, -0.1429, 0.1155),
                      weight=c(1.6408, -1.9030, -0.1428, 452.6100, -7.6864),
                      madrs=c(-0.9666, 0.1380, 0.1155, -7.6864, 17.5343),
                      row.names=c("age","sex","smstat","weight","madrs"))
  dat <- MASS::mvrnorm(n=n*K, mu=mu, Sigma=Sigma) %>%
    as.data.frame() %>%
    mutate(sex = ifelse(sex > 1-0.6784, 1, 0),
           smstat = ifelse(smstat > 1-0.3043, 1, 0),
           eps = rnorm(n*K, mean=0, sd=.05))
  
  train_dat <- add_study(dat, K, n, distribution) %>%
    mutate()
  
  #testing data
  if (test_dist == "same") {
    test_dat <- train_dat[sample(nrow(train_dat), 100),] %>%
      select(-S, -id)
  } else if (test_dist == "upweight") {
    train_weight <- train_dat %>%
      mutate(study_weight = ifelse(S %in% c(3, 5), 3, 1))
    test_dat <- train_weight[sample(nrow(train_weight), 100, prob=train_weight$study_weight),] %>%
      select(-study_weight, -S, -id)
  } else if (test_dist == "different") {
    ## MAKE THIS MVRNORM
    test_dat <- data.frame(
      sex = rbinom(n=100, size=1, prob=.65),
      smstat = rbinom(n=100, size=1, prob=.3),
      weight = rtruncnorm(n=100, a=45, b=140, mean=80, sd=15),
      W = rbinom(n=100, size=1, prob=.5),
      eps = rnorm(n=100, mean=0, sd=.05),
      age = rtruncnorm(n=100, a=18, b=40, mean=25, sd=10), #younger
      madrs = rtruncnorm(n=100, a=16, b=40, mean=25, sd=4.1) #less severe
    )
  }
  
  #m and tau
  train_dat <- train_dat %>% 
    mutate(m = -0.02*age - 0.7*madrs - 0.15*sex + eps_m,
           tau = -8.5 + 0.07*age + 0.20*madrs + eps_tau)
  test_dat <- test_dat %>% 
    mutate(m = -0.02*age - 0.7*madrs - 0.15*sex,
           tau = -8.5 + 0.07*age + 0.20*madrs)
  
  #outcome Y
  train_dat <- train_dat %>%
    mutate(Y = round(m + W*tau + eps, 0),
           S = factor(S)) %>%
    select(S, W, sex, smstat, weight, age, madrs, Y, tau)
  
  test_dat <- test_dat %>%
    mutate(Y = round(m + W*tau + eps, 0)) %>%
    select(W, sex, smstat, weight, age, madrs, Y, tau)
  
  return(list(train_dat=train_dat, test_dat=test_dat))
}






##other option for simulating covariates
# add_agemadrs <- function(dat, n, k, distribution) {
# 
#   #add age and madrs
#   if (distribution == "same") {
#     dat <- dat %>%
#       mutate(age = rtruncnorm(n=n, a=18, b=75, mean=45, sd=10),
#              madrs = rtruncnorm(n=n, a=26, b=60, mean=31, sd=4.1))
# 
#   } else if (distribution == "varying_madrs") {
#     dat <- dat %>%
#       mutate(age = rtruncnorm(n=n, a=18, b=75, mean=45, sd=10),
#              madrs = rtruncnorm(n=n, a=26-k*1.5, b=60-k*1.5, mean=31-k*1.5, sd=4.1))
# 
#   } else if (distribution == "halfdiff_madrsage") {
#     if (k%%2 == 0 ) {
#       dat <- dat %>%
#         mutate(age = rtruncnorm(n=n, a=30, b=75, mean=50, sd=10),
#                madrs = rtruncnorm(n=n, a=30, b=60, mean=40, sd=4.1))
#     } else {
#       dat <- dat %>%
#         mutate(age = rtruncnorm(n=n, a=18, b=75, mean=45, sd=10),
#                madrs = rtruncnorm(n=n, a=26, b=60, mean=31, sd=4.1))
#     }
# 
#   } else if (distribution == "separate_age") {
#     ages <- seq(18,75,by=(75-18)/10)
#     dat <- dat %>%
#       mutate(age = runif(n=n, min=ages[k], max=ages[k+1]),
#              madrs = rtruncnorm(n=n, a=26, b=60, mean=31, sd=4.1))
#   }
# 
#   return(dat)
# }


## other tau options
# if (scenario == "linear") {
#   study_main <- rnorm(K, mean=-10, sd=1)
#   study_inter <- rnorm(K, mean=0.25, sd=0.05)
#   study_tau <- rnorm(K, mean=3, sd=0.05)
#   train_dat <- train_dat %>% 
#     mutate(study_main = study_main[S],
#            study_inter = study_inter[S],
#            study_tau = study_tau[S],
#            m = 10.7 + study_main - 0.02*age - 0.87*madrs -
#              0.15*sex + study_inter*madrs,
#            tau = -10.5 + 0.07*age + 0.20*madrs + study_tau) %>%
#     select(-study_main, -study_inter, -study_tau)
#   
#   test_main <- rnorm(1, mean=-10, sd=1)
#   test_inter <- rnorm(1, mean=0.25, sd=0.05)
#   test_tau <- rnorm(1, mean=3, sd=0.05)
#   test_dat <- test_dat %>%
#     mutate(m = 10.7 + test_main - 0.02*age - 0.87*madrs -
#              0.15*sex + test_inter*madrs,
#            tau = -10.5 + 0.07*age + 0.20*madrs + test_tau)
# }

# if (scenario == "linear") {
#   study_main <- runif(K, min=-14, max=-7)
#   study_inter <- runif(K, min=0.1, max=0.5)
#   study_tau <- runif(K, min=2.5, max=3.5)
#   train_dat <- train_dat %>% 
#     mutate(study_main = study_main[S],
#            study_inter = study_inter[S],
#            study_tau = study_tau[S],
#            m = 10.7 + study_main - 0.02*age - 0.87*madrs -
#              0.15*sex + study_inter*madrs,
#            tau = -8.5 + 0.07*age + 0.20*madrs + study_tau) %>%
#     select(-study_main, -study_inter, -study_tau)
#   
#   test_main <- runif(1, min=-14, max=-7)
#   test_inter <- runif(1, min=0.1, max=0.5)
#   test_tau <- runif(1, min=2.5, max=3.5)
#   test_dat <- test_dat %>%
#     mutate(test_main = test_main + rnorm(nrow(test_dat), mean=0, sd=0.1),
#            test_inter = test_inter + rnorm(nrow(test_dat), mean=0, sd=0.1),
#            test_tau = test_tau + rnorm(nrow(test_dat), mean=0, sd=0.1),
#            m = 10.7 + test_main - 0.02*age - 0.87*madrs -
#              0.15*sex + test_inter*madrs,
#            tau = -8.5 + 0.07*age + 0.20*madrs + test_tau) %>%
#     select(-test_main, -test_inter, -test_tau)
# }

# if (scenario == "nonlinear") {
#   study_tau <- runif(K, min=2.5, max=3.5)
#   train_dat <- train_dat %>% 
#     mutate(study_tau = study_tau[S],
#            m = 0,
#            tau = (study_tau/(1+exp(-1/12*age)))*(study_tau/(1+exp(-12*madrs)))) %>%
#     select(-study_tau)
#   
#   test_tau <- runif(1, min=2.5, max=3.5)
#   test_dat <- test_dat %>%
#     mutate(test_tau = test_tau + rnorm(nrow(test_dat), mean=0, sd=0.1),
#            m = 0,
#            tau = (test_tau/(1+exp(-1/12*age)))*(test_tau/(1+exp(-12*madrs)))) %>%
#     select(-test_tau)
# }
# 
