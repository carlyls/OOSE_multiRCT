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
gen_data <- function (K=6, n_mean=200, n_sd=0, scenario="1a", 
                      distribution="same") {
  
  all_dat <- data.frame()
  n_study <- floor(rnorm(K, mean=n_mean, sd=n_sd))
  
  for (k in 1:K) {
    n <- n_study[k]
    
    #sample covariates
    dat <- data.frame(
      #age = rtruncnorm(n=n, a=18, b=75, mean=45, sd=10),
      sex = rbinom(n=n, size=1, prob=.65),
      smstat = rbinom(n=n, size=1, prob=.3),
      weight = rlnorm(n=n, meanlog=4.37, sdlog=.25),
      #madrs = rtruncnorm(n=n, a=26, b=60, mean=31, sd=4.1),
      W = rbinom(n=n, size=1, prob=.5),
      S = rep(k, n),
      id = seq(1, n),
      eps = rnorm(n, mean=0, sd=.05)
    )
    
    dat <- add_agemadrs(dat, n, k, distribution)
    
    all_dat <- bind_rows(all_dat, dat)
  }

  #tau and Y
  if (scenario == "1a") {
    study_main <- runif(K, min=-14, max=-7)
    study_inter <- runif(K, min=0.1, max=0.5)
    study_tau <- runif(K, min=2.5, max=3.5)
    all_dat <- all_dat %>% 
      mutate(study_main = study_main[S],
             study_inter = study_inter[S],
             study_tau = study_tau[S],
             m = 10.7 - study_main - 0.02*age - 0.87*madrs -
               0.15*sex + study_inter*madrs,
             tau = -8.5 + 0.07*age + 0.20*madrs + study_tau) %>%
      select(-study_main, -study_inter, -study_tau)
  }
  
  if (scenario == "1b") {
    all_dat$m <- 0
    study_tau <- runif(K, min=2.5, max=3.5)
    all_dat$tau <- (2/(1+exp(-12*(all_dat$X1-1/2))))*(2/(1+exp(-12*(all_dat$X2-1/2)))) +
      study_main[all_dat$S] + study_inter[all_dat$S]*all_dat$X1
    all_dat <- all_dat %>% 
      mutate(study_tau = study_tau[S],
             m = 0,
             tau = (study_tau/(1+exp(-1/12*age)))*(study_tau/(1+exp(-12*madrs)))) %>%
      select(-study_tau)
  }
  
  all_dat <- all_dat %>%
    mutate(Y = m + W*tau + eps) %>%
    select(-c(eps, m)) %>%
    mutate(S = factor(S)) %>%
    relocate(S, id, W, sex, smstat, weight, age, madrs, Y, tau)
  return(all_dat)
}


