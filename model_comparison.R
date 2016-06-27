rm(list=ls())

get_gg_criterion_G_and_P <- function(model) {
    Yobs <- model.matrix( ~ factor(Y), dat=model$dat)
    predicted_category_probs <- sapply(1:3,
        function(category) colMeans(model$post_preds==category))
    predicted_category_vars <- sapply(1:3,
        function(category) apply(model$post_preds==category, 2, var))

    G <- sum((Yobs - predicted_category_probs)^2) # Goodness of fit
    P <- sum(predicted_category_vars) # Penalty term

    return(list(G=G, P=P))
}

gg_criterion <- function(k, model) {
    g_and_p <- get_gg_criterion_G_and_P(model)
    G <- g_and_p$G
    P <- g_and_p$P
    D <- (k/(k+1))*G + P

    return(D)
}

simple_ordinal_model <- readRDS(file="simple_ordinal_samps.rds")
sapply(c(.01, .1, 1, 10, 100, 1000), gg_criterion, simple_ordinal_model)
