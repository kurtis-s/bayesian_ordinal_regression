rm(list=ls())

library(arm)

set.seed(83299)

source("helper_funcs.R")

J <- 3 # Number of categories
P <- 3 # Number of covariates, beta
true_delta <- c(.1, 1)
true_beta <- c(1.5, -1, 0)

n_obs <- 5000
X <- matrix(runif(n_obs*P, min=-1, max=1), ncol=P)
eta <- get_eta(true_beta, X)
# column 1 is delta1 - eta, column 2 is delta2 - eta, etc.
linear_part <- get_linear_part(true_delta, eta)
true_gamma <- get_gamma_from_linear_part(linear_part)
true_pi <- get_pi_from_gamma(true_gamma)

Y <- vector(length=n_obs)
for(i in 1:n_obs) {
    category <- rmultinom(1, 1, prob=true_pi[i,])
    Y[i] <- which(category==1)
}

dat <- data.frame(Y=Y, X=X)

ordinal_data <- list(dat=dat, beta=true_beta, delta=true_delta,
                     gamma=true_gamma, pi=true_pi)
saveRDS(ordinal_data, file="ordinal_dat.rds")