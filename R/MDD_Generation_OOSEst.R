### Simulation Approach for OOS Estimation: Based on MDD Data ###


# library(tidyverse)
# library(rsample)
# library(lme4)
# library(multcomp)
# library(MASS)
# library(grf)
# library(nnet)


#generate target sample - JUST ONCE
#use standardized covariates, same distribution as the main training data (distribution = "same")
# mu <- c(age=0, sex=0.6784, smstat=0.3043, weight=0, madrs=0)
# Sigma <- data.frame(age=c(1, 0.0190, -0.0402, 0.0060, -0.0179),
#                     sex=c(0.0190, 0.2183, -0.0218, -0.0894, 0.0330),
#                     smstat=c(-0.0402, -0.0218, 0.2118, -0.0067, 0.0276),
#                     weight=c(0.0060, -0.0894, -0.0067, 1, -0.0863),
#                     madrs=c(-0.0179, 0.0330, 0.0276, -0.0863, 1),
#                     row.names=c("age","sex","smstat","weight","madrs"))
# set.seed(1)
# target_sample <- MASS::mvrnorm(n=100, mu=mu, Sigma=Sigma) %>%
#   as.data.frame() %>%
#   mutate(W = rbinom(n=100, size=1, prob=.5),
#          sex = ifelse(sex > 1-0.6784, 1, 0),
#          smstat = ifelse(smstat > 1-0.3043, 1, 0),
#          age2 = age^2,
#          eps = rnorm(n=100, mean=0, sd=.05),
#          id = seq(1:100)) %>%
#   dplyr::select(id, W, sex, smstat, age, age2, weight, madrs, eps)
# saveRDS(target_sample, "Data/target_sample.RDS")


#interior function
sample_dist <- function(k, n, Sigma, eps_study_m, eps_study_tau, 
                        eps_study_inter, covars_rand, distribution) {
  
  #define mu based on distribution input
  if (distribution == "same") {
    mu <- c(age=0, sex=0.6784, smstat=0.3043, weight=0, madrs=0)

    } else if (distribution == "varying_madrs") {
      mu <- c(age=0, sex=0.6784, smstat=0.3043, weight=0, madrs=0-k*0.2)
    
    } else if (distribution == "halfdiff_madrsage") {
      if (k%%2 == 0 ) {
        mu <- c(age=1, sex=0.6784, smstat=0.3043, weight=0, madrs=1)
      } else {
        mu <- c(age=0, sex=0.6784, smstat=0.3043, weight=0, madrs=0)
      }
    
    } else if (distribution == "separate_age") {
      if (k == 1) {
        mu <- c(age=1, sex=0.6784, smstat=0.3043, weight=0, madrs=0)
      } else {
        mu <- c(age=-1, sex=0.6784, smstat=0.3043, weight=0, madrs=0)
      }
  
    }
  
  #define random slopes for moderators
  eps_inter <- matrix(nrow=n, ncol=length(covars_rand))
  colnames(eps_inter) <- paste0("eps_",covars_rand)
  for (i in 1:length(covars_rand)) {
    eps_inter[,i] <- rnorm(n=1, mean=0, sd=eps_study_inter[i])
  }
  
  #simulate data
  dat <- MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma) %>%
    as.data.frame() %>%
    mutate(S = k,
           W = rbinom(n=n, size=1, prob=.5),
           sex = ifelse(sex > 1-0.6784, 1, 0),
           smstat = ifelse(smstat > 1-0.3043, 1, 0),
           age2 = age^2,
           eps = rnorm(n=n, mean=0, sd=.05),
           eps_m = rnorm(n=1, mean=0, sd=eps_study_m),
           eps_tau = rnorm(n=1, mean=0, sd=eps_study_tau)) %>%
    bind_cols(eps_inter)
  
  return(dat)
}


#main function
gen_mdd <- function (K=10, n_mean=500, n_sd=0, covars_fix="age", covars_rand="age",
                     lin=T, eps_study_m=0.05, eps_study_tau=0.05, eps_study_inter=0.05,
                     distribution="same", target_sample) {
  
  #training data
  train_dat <- data.frame()
  n_study <- floor(rnorm(K, mean=n_mean, sd=n_sd))
  
  #define covariance matrix
  Sigma <- data.frame(age=c(1, 0.0190, -0.0402, 0.0060, -0.0179),
                      sex=c(0.0190, 0.2183, -0.0218, -0.0894, 0.0330),
                      smstat=c(-0.0402, -0.0218, 0.2118, -0.0067, 0.0276),
                      weight=c(0.0060, -0.0894, -0.0067, 1, -0.0863),
                      madrs=c(-0.0179, 0.0330, 0.0276, -0.0863, 1),
                      row.names=c("age","sex","smstat","weight","madrs"))
  
  for (k in 1:K) {
    n <- n_study[k]
    
    #sample
    dat <- sample_dist(k, n, Sigma, eps_study_m, eps_study_tau, 
                       eps_study_inter, covars_rand, distribution)
    train_dat <- bind_rows(train_dat, dat)
    
  }
  
  #add m and tau
  if (length(covars_fix) == 1 & length(covars_rand) == 1) { #age only moderators
    
    if (lin == T) { #linear cate
      train_dat <- train_dat %>% 
        mutate(m = (-17.40 + eps_m) - 0.13*age - 2.05*madrs - 0.11*sex,
               tau = (2.505 + eps_tau) + (0.82 + eps_age)*age)
      target_dat <- target_sample %>% 
        mutate(m = (-17.40) - 0.13*age - 2.05*madrs - 0.11*sex,
               tau = (2.505) + (0.82)*age)
      
    } else { #nonlinear cate
      train_dat <- train_dat %>% 
        mutate(m = (-17.52 + eps_m) - 0.08*age,
               tau = (2.2 + eps_tau)*exp((.35 + eps_age)*age))
      target_dat <- target_sample %>% 
        mutate(m = (-17.52) - 0.08*age,
               tau = (2.2)*exp((.35)*age))
      
    }
  } else if (length(covars_fix) == 2 & length(covars_rand) == 2) { #age and madrs moderators
    
    train_dat <- train_dat %>% 
      mutate(m = (-17.29 + eps_m) - 0.20*age - 2.63*madrs - 0.14*sex,
             tau = (2.506 + eps_tau) + (0.91 + eps_age)*age + (0.84 + eps_madrs)*madrs)
    target_dat <- target_sample %>% 
      mutate(m = (-17.29) - 0.20*age - 2.63*madrs - 0.14*sex,
             tau = (2.506) + (0.91)*age + (0.84)*madrs)
      
  } else if (length(covars_fix) == 2 & length(covars_rand) == 1) { #age^2 and sex (fixed) moderators
    
    train_dat <- train_dat %>% 
      mutate(m = (-17.16 + eps_m) + 0.16*age2 - 2.07*madrs - 0.43*sex,
             tau = (2.14 + eps_tau) + (-0.16 + eps_age2)*age2 + (0.49)*sex)
    target_dat <- target_sample %>% 
      mutate(m = (-17.16) + 0.16*age2 - 2.07*madrs - 0.43*sex,
             tau = (2.14) + (-0.16)*age2 + (0.49)*sex)
    
  } 
  
  #outcome Y
  train_dat <- train_dat %>%
    mutate(Y = m + W*tau + eps,
           S = factor(S)) %>%
    dplyr::select(S, W, sex, smstat, age, age2, weight, madrs, Y, tau)
  
  target_dat <- target_dat %>%
    mutate(Y = m + W*tau + eps) %>%
    dplyr::select(id, W, sex, smstat, age, age2, weight, madrs, Y, tau)
  
  return(list(train_dat=train_dat, target_dat=target_dat))
}
