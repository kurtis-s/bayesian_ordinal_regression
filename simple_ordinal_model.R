rm(list=ls())

library(ordinal)
library(coda)

set.seed(8392921)
source("helper_funcs.R")

ordinal_dat <- readRDS(file="ordinal_dat.rds")
dat <- ordinal_dat$dat

# Frequentist estimates ---------------------------------------------------
freq_model <- clm(factor(Y) ~ ., data=dat)
summary(freq_model)

# Bayesian estimates ------------------------------------------------------
Pdim <- 3
Jdim <- nlevels(factor(dat$Y))
Xobs <- as.matrix(dat[,1 + 1:Pdim])
Yobs <- model.matrix( ~ factor(Y) - 1, data=dat)

log_likelihood <- function(pi_probs) {
    #' @param pi_probs N x J matrix of category probabilities
    #'
    #' @return log-likelihood of the observations up to a constant

    # TODO: Need to add restriction for the ordering of delta
    return(sum(Yobs * log(pi_probs)))
}

get_coefficient_proposal <- function(coefficient) {
    #' Get a new proposed coefficient value for the MCMC
    return(coefficient + rnorm(1, 0, sd=.1))
}

metropolis_chooser <- function(param_cur, param_prop, log_lik_cur, log_lik_prop) {
    #' Return current or proposed value for Metropolis-Hastings step
    #'
    #' @param param_cur current parameter value
    #' @param param_prop proposed parameter value
    acceptance_ratio <- exp(log_lik_prop - log_lik_cur)
    u <- runif(1)
    if(acceptance_ratio > u) {
        ret <- param_prop
    }
    else {
        ret <- param_cur
    }

    return(ret)
}

get_coef_proposal <- function(coefs, idx, coef_name) {
    #' Get a new proposal for the specified vector (for use in M-H step)
    #'
    #' @param coefs list of the different coefficient vectors
    #' @param coef_name string indicating which coefficients to update.  The
    #'   string should be a key in the coefs vector
    coefs[[coef_name]][idx] <- get_coefficient_proposal(coefs[[coef_name]][idx])

    return(coefs)
}

sample_coef <- function(coefs, coef_name) {
    #' Sample a new coefficient vector
    #'
    #' @param coefs list of the different coefficient vectors
    #' @param coef_name string indicating which coefficients to update.  The
    #'   string should be a key in the coefs vector (e.g. "delta", or "beta")
    #'
    #' @return list of updated coefficient vectors
    Kdim <- length(coefs[[coef_name]])
    for(k in 1:Kdim) {
        pi_cur <- get_pi_from_coefs(coefs, Xobs)
        log_lik_cur <- log_likelihood(pi_cur)

        coefs_prop <- get_coef_proposal(coefs, k, coef_name)
        pi_prop <- get_pi_from_coefs(coefs_prop, Xobs)
        log_lik_prop_p <- log_likelihood(pi_prop)

        new_coefs <- metropolis_chooser(coefs, coefs_prop, log_lik_cur, log_lik_prop_p)
        coefs <- new_coefs
    }

    return(coefs)
}

sampler <- function(nsamps, nburn, coefs_init) {
    beta_samps <- matrix(nrow=nsamps, ncol=Pdim)
    delta_samps <- matrix(nrow=nsamps, ncol=Jdim - 1)
    coefs <- coefs_init
    for(b in -(nburn+1):nsamps) {
        coefs <- sample_coef(coefs, "beta")
        coefs <- sample_coef(coefs, "delta")
        beta_samps[b,] <- coefs$beta
        delta_samps[b,] <- coefs$delta
    }

    return(list(beta=beta_samps, delta=delta_samps))
}
nsamps <- 1000
nburn <- 1000
samps1 <- sampler(nsamps, nburn, list(beta=ordinal_dat$beta + .5, delta=ordinal_dat$delta))
samps2 <- sampler(nsamps, nburn, list(beta=rep(0, 3), delta=ordinal_dat$delta))
chain1 <- mcmc(cbind(samps1$beta, samps1$delta))
chain2 <- mcmc(cbind(samps2$beta, samps2$delta))
mcmc_chains <- mcmc.list(chain1, chain2)
gelman.diag(mcmc_chains)
