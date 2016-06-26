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
Xobs <- as.matrix(dat[,1 + 1:Pdim])
Yobs <- model.matrix( ~ factor(Y) - 1, data=dat)

log_likelihood <- function(pi_probs) {
    #' @param pi_probs N x J matrix of category probabilities
    #'
    #' @return log-likelihood of the observations up to a constant

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

sample_beta <- function(beta, delta) {
    for(p in 1:Pdim) {
        pi_cur <- get_pi_from_delta_beta(delta, beta, Xobs)
        log_lik_cur <- log_likelihood(pi_cur)

        beta_p_proposal <- get_coefficient_proposal(beta[p])
        beta_proposal <- beta
        beta_proposal[p] <- beta_p_proposal

        pi_prop <- get_pi_from_delta_beta(delta, beta_proposal, Xobs)
        log_lik_prop_p <- log_likelihood(pi_prop)

        new_beta_p <- metropolis_chooser(beta[p], beta_p_proposal, log_lik_cur, log_lik_prop_p)
        beta[p] <- new_beta_p
    }

    return(beta)
}

sampler <- function(nsamps, nburn, beta_init) {
    beta_samps <- matrix(nrow=nsamps, ncol=Pdim)
    beta <- ordinal_dat$beta + .5
    for(b in -(nburn+1):nsamps) {
        beta <- sample_beta(beta, ordinal_dat$delta)
        beta_samps[b,] <- beta
    }

    return(list(beta=beta_samps))
}
nsamps <- 1000
nburn <- 1000
samps1 <- sampler(nsamps, nburn, original_dat$beta + .5)
samps2 <- sampler(nsamps, nburn, rep(0, 3))
chain1 <- mcmc(samps1$beta)
chain2 <- mcmc(samps2$beta)
mcmc_chains <- mcmc.list(chain1, chain2)
gelman.diag(mcmc_chains)