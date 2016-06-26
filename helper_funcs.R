get_eta <- function(beta, X) {
    return(X %*% beta)
}

get_pi_from_gamma <- function(gamma) {
    #' @param gamma N x J matrix of cumulative probabilities
    #'
    #' @return pi N x J matrix of category probabilities
    ret <- t(rbind(gamma[,1], diff(t(gamma))))

    return(ret)
}

get_pi_from_gamma2 <- function(gamma) {
    #' @param gamma N x J matrix of cumulative probabilities
    #'
    #' @return pi N x J matrix of category probabilities
    diff(t(gamma))
}

get_pi_from_coefs <- function(coefs, X) {
    beta <- coefs$beta
    delta <- coefs$delta

    eta <- get_eta(beta, X)
    linear_part <- get_linear_part(delta, eta)
    gamma <- get_gamma_from_linear_part(linear_part)
    pi <- get_pi_from_gamma(gamma)

    return(pi)
}

get_gamma_from_linear_part <- function(linear_part) {
    gamma <- cbind(invlogit(linear_part), 1)

    return(gamma)
}

get_linear_part <- function(delta, eta) {
    #' @param delta vector of length J-1
    #' @param eta matrix X %*% beta

    # column 1 is delta1 - eta, column 2 is delta2 - eta, etc.
    linear_part <- sapply(delta, function(deltaj) deltaj - eta)

    return(linear_part)
}