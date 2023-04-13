## Try out 0 matrix idea

sim_dat <- gen_mdd(K, n_mean, n_sd, n_target, eps_study_m, eps_study_tau, 
                   eps_study_age, distribution, target_dist)
train_dat <- sim_dat[["train_dat"]]
target_dat <- sim_dat[["target_dat"]]

mod <- lmer(Y ~  madrs + sex + W*age +
              (W + W:age | S), data=train_dat,
            control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000)))
sum <- summary(mod)

# train_dat <- train_dat %>%
#   mutate(intercept = 1,
#          W_age = W*age)
# target_dat <- target_dat %>%
#   mutate(intercept = 1,
#          W_age = W*age)
# 
# ## Fit mixed effects models
# #correct
# mod <- lmer(Y ~  0 + intercept + madrs + sex + W + age + W_age +
#               (0 + intercept + W + W_age | S), data=train_dat,
#             control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000)))
# sum <- summary(mod)
# 
# #prediction interval by hand ####
# beta <- fixef(mod) #beta-hat
# u <- ranef(mod)$S %>% t() %>% matrix()
# 
# fc <- vcov(mod) #covariance matrix of fixed effects
# rc <- Matrix::bdiag(VarCorr(mod))
# 
# X <- train_dat %>%
#   mutate(W = 1,
#          W_age = age,
#          intercept = 0,
#          madrs = 0,
#          sex = 0,
#          age = 0)
# 
# p <- predict(mod, X)


covars_fix <- c("age")
covars_rand <- c("age")
formula <- as.formula(paste0("Y ~ madrs + sex + ", 
                             paste("W", covars_fix, sep="*", collapse=" + "),
                             " + (W + ",
                             paste("W", covars_rand, sep=":", collapse=" + "),
                             " | S)"))
mod <- lmer(formula, data=train_dat,
            control=lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000)))
sum <- summary(mod)


#get variance components of model
matrix_var <- function(mod) {
  
  #calculate variances
  fc <- vcov(mod) #covariance matrix of fixed effects
  rc <- Matrix::bdiag(VarCorr(mod)) #variance and covariances of random effects
  
  #fixed effects
  beta <- fixef(mod)[grep("W", names(fixef(mod)))] %>% matrix() #beta-hat
  var_beta <- fc[grep("W", rownames(fc)),
                 grep("W", colnames(fc))] #Var(beta-hat)
  
  #random effects
  u <- ranef(mod)$S[grep("W", colnames(ranef(mod)$S))] %>% t() %>% matrix()
  var_rand <- rc[grep("W", rownames(rc)),
                 grep("W", colnames(rc))]
  
  return(list(beta=beta, var_beta=var_beta, u=u, var_rand=var_rand))
}

#add pis to dataset
manual_pi <- function(df, mod, K, covars_fix, covars_rand) {
  
  #get variances
  res <- matrix_var(mod)
  beta <- res$beta
  var_beta <- res$var_beta
  u <- res$u
  var_rand <- res$var_rand
  
  #calculate theta-hats
  X <- df %>%
    dplyr::select(W, all_of(covars_fix)) %>%
    mutate(W = 1) %>%
    as.matrix()
  Z <- df %>%
    dplyr::select(W, all_of(covars_rand)) %>%
    mutate(W = 1) %>%
    as.matrix()
  
  if ("S" %in% names(df)) { #training data structure is different
    Zlist <- list()
    var_rand_train <- rep(list(var_rand), K)
    var_rand <- do.call("bdiag", var_rand_train) %>% as.matrix()
    
    for (s in 1:K) {
      Z_S <- Z[df$S==s,]
      Zlist[[s]] <- Z_S
    }
    Z <- do.call("bdiag", Zlist) %>% as.matrix()
    
    mean_theta <- X %*% beta + Z %*% u %>% c()
    
  } else { #don't know u's for target data
    
    mean_theta <- X %*% beta %>% c() #theta-hat
  }
  
  vcov_theta <- X %*% var_beta %*% t(X) + Z %*% var_rand %*% t(Z)
  var_theta <- diag(vcov_theta) #Var(theta-hat)
  
  #prediction interval
  pred_lower <- mean_theta + qt(.025, K-2)*sqrt(var_theta)
  pred_upper <- mean_theta - qt(.025, K-2)*sqrt(var_theta)
  
  df <- df %>%
    mutate(lower = pred_lower,
           mean = mean_theta,
           upper = pred_upper)
  
  return(df)
}
