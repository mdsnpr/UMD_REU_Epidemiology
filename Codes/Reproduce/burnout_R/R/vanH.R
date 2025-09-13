##' Persistence probability via van Herwaarden (1997) approximation
##'
##' @seealso \code{\link{burnout_prob_vanH}}, \code{\link{P1_prob_MS}},
##'     \code{\link{P1_prob}}
##'
##' @param persist return persistence probability rather than burnout probability
##' @inheritParams P1_prob_other
##' @inheritParams stats::integrate
##'
##' @importFrom stats integrate
##'
##' @references
##' \insertRef{Ballard2016}{burnout}
##' 
##' \insertRef{vanH1997}{burnout}
##'
##' @export
##'
P1_prob_vanH <- function(R0, epsilon, k=1, N=10^6, subdivisions=1000L,
                         persist = FALSE, tiny=0,
                         ... ) {

    P1 <- P1_prob_other(R0 = R0, epsilon = epsilon,
                        burnout_prob_fun = burnout_prob_vanH, k = k, N = N,
                        subdivisions = subdivisions, tiny = tiny, ... )
    return(P1)
}

##' Burnout probability based on van Herwaarden (1997)
##'
##' @param debug print debugging info?
##' @details
##' See equation (8) of Ballard et al (2016)
##'
##' @seealso \code{\link{P1_prob_vanH}}
##'
##' @inheritParams P1_prob_vanH
##' @param tiny amount which to avoid integration limits (where
##'     integrand blow up)
##' 
##' @importFrom stats integrate
##'
##' @references
##' \insertRef{Ballard2016}{burnout}
##' 
##' \insertRef{vanH1997}{burnout}
##' 
##' @export
##'
burnout_prob_vanH <- function( R0, epsilon, N=10^6,
                              subdivisions=1000L,
                              tiny=0,
                              debug = FALSE,
                              persist = FALSE,
                              ... ) {

    dfun <- function(x) {
        if (debug) cat(x, get(x), "\n")
    }
    ## choose units such that gamma+mu = 1, i.e., mean time infected
    ## period is 1:
    beta <- R0
    mu <- epsilon
    gamma <- 1-epsilon

    bog <- beta/gamma
    x1A <- (-1/bog)*W0(-bog*exp(-bog))

    dfun("x1A")
    integrand <- function(s) {
        (x1A/(1-x1A)) * gamma*(s - s*log(s) - 1) /
            ((beta*s^2)*(1-s+(1/bog)*log(s))) +
            1/(s-x1A)
    }
    messy.integral <-
        try(stats::integrate(f=integrand, lower=x1A+tiny, upper=1-tiny,
                             subdivisions=subdivisions, ...)$value)

    ## if (debug) curve(integrand(x), from = x1A+tiny, to = 1-tiny)
    
    dfun("messy.integral")
    
    ## FIX: this destroys the automatic vectorization:
    if ("try-error" %in% class(messy.integral)) return(NA)

    C3 <- -log(-beta*x1A / (beta*x1A - gamma)) - messy.integral

    dfun("C3")
    
    K <- (1/mu)*exp(((1/mu)*(beta*x1A + (beta-gamma)*log(1-x1A))) + C3)

    logK <- -log(mu) + (((1/mu)*(beta*x1A + (beta-gamma)*log(1-x1A))) + C3)

    dfun("K")
    

    x1 <- (beta-gamma-mu)/mu
    lnum <- logK + log(N) + 2*log(mu) + (x1*log(beta/mu) - beta/mu)
    ldenom <- log(gamma+mu)  + lgamma(x1)
    log_p <- -1*exp(lnum-ldenom)
    p0 <- if (!persist) exp(log_p) else -expm1(log_p)

    ## Gamma <- base::gamma
    ## p0 <- exp(-K*N*mu^2*(beta/mu)^((beta-gamma-mu)/mu) * exp(-beta/mu) /
    ##           ((gamma+mu)*Gamma((beta-gamma-mu)/mu))
    ##           )

    return(p0)
}
