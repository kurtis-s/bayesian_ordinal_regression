rm(list=ls())

library(arm)

set.seed(83299)

source("helper_funcs.R")

J <- 3 # Number of categories
P <- 3 # Number of covariates, beta
true_delta <- c(.1, 1)
true_beta <- c(1.5, -1, 0)

n_obs <- 5000
true_coefs <- list(beta=true_beta, delta=true_delta)
Xobs <- matrix(runif(n_obs*P, min=-1, max=1), ncol=P)
Yobs <- gen_obs(true_coefs)
dat <- data.frame(Y=Yobs, X=Xobs)

ordinal_data <- list(dat=dat, beta=true_beta, delta=true_delta)
saveRDS(ordinal_data, file="ordinal_dat.rds")