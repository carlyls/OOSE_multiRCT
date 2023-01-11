### Simulation Approach for OOS Estimation ###
### Based on MDD Data ###

library(tidyverse)
library(rsample)

## Data Generation Function
## Tau based on Tan et al., 2021; Wager and Athey, 2018; Kunzel et al., 2019; MDD Vivli Data

gen_data <- function (K, n_mean, n_sd, study_mean, study_inter_mean,
                      study_sd, study_inter_sd, scenario, distribution,
                      sd=sqrt(0.5)) {
  
  all_dat <- data.frame()
  
  n_study <- floor(rnorm(K, mean=n_mean, sd=n_sd))
  
  if (scenario %in% c("1a","1b")) {
    study_main <- rnorm(K, mean=study_mean, sd=study_sd)
    study_inter <- rnorm(K, mean=study_inter_mean, sd=study_inter_sd)
  }
  
  for (k in 1:K) {
    
    n <- n_study[k]
    
    #sample covariates
    dat <- data.frame(matrix(rnorm(n*ncovar), nrow=n, ncol=ncovar))
    colnames(dat) <- paste0("X", seq(1,ncovar))
    
    #treatment
    dat$W <- rbinom(n, size=1, prob=0.5)
    
    #study and id
    dat$S <- rep(k, n)
    dat$id <- seq(1, n)
    
    #noise
    dat$eps <- rnorm(n, mean=0, sd=sd)
    
    all_dat <- bind_rows(all_dat, dat)
  }
  
  #tau and Y
  if (scenario == "1a") {
    all_dat$m <- all_dat$X1/2 + all_dat$X2 + all_dat$X3 + all_dat$X4 +
      study_main[all_dat$S] + study_inter[all_dat$S]*all_dat$X1
    all_dat$tau <- all_dat$X1*(all_dat$X1>0) + study_main[all_dat$S] + 
      study_inter[all_dat$S]*all_dat$X1
  }
  
  if (scenario == "1b") {
    all_dat$m <- 0
    all_dat$tau <- (2/(1+exp(-12*(all_dat$X1-1/2))))*(2/(1+exp(-12*(all_dat$X2-1/2)))) +
      study_main[all_dat$S] + study_inter[all_dat$S]*all_dat$X1
  }
  
  if (scenario == "2") {
    all_dat$m <- all_dat$X1/2 + all_dat$X2 + all_dat$X3 + all_dat$X4
    all_dat$tau <- ifelse(all_dat$S %in% c(1:4), (2/(1+exp(-12*(all_dat$X1-1/2))))*(2/(1+exp(-12*(all_dat$X2-1/2)))),
                          ifelse(all_dat$S %in% c(5:8), all_dat$X1*(all_dat$X1>0), 0))
  }
  
  all_dat$Y <- all_dat$m + (2*all_dat$W-1)/2*all_dat$tau + all_dat$eps
  
  all_dat <- all_dat %>%
    select(-eps,-m) %>%
    mutate(S = factor(S)) %>%
    relocate(S, id, W, X1, X2, X3, X4, X5, Y, tau)
  return(all_dat)
}


