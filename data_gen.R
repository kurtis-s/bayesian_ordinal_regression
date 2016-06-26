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
eta <- get_eta(X, true_beta)
# column 1 is delta1 - eta, column 2 is delta2 - eta, etc.
linear_part <- sapply(true_delta, function(delta) delta - eta)
true_gamma <- cbind(invlogit(linear_part), 1)
true_pi <- t(apply(true_gamma, 1,
                 function(gamma_vec) c(gamma_vec[1], diff(gamma_vec))))

Y <- vector(length=n_obs)
for(i in 1:n_obs) {
    category <- rmultinom(1, 1, prob=true_pi[i,])
    Y[i] <- which(category==1)
}

dat <- data.frame(Y=Y, X=X)
clm(factor(Y) ~ X, data=dat)
