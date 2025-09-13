##' Persistence probability via Meerson and Sasorov (2009) approximation
##'
##' @seealso \code{\link{burnout_prob_MS}}, \code{\link{P1_prob_MS}},
##'     \code{\link{P1_prob}}
##' 
##' @inheritParams P1_prob_other
##'
##' @importFrom stats integrate
##'
##' @references
##' \insertRef{Meerson2009}{burnout}
##' 
##' @export
##'
P1_prob_MS <- function( R0, epsilon, k=1, N=10^6, subdivisions=1000L, tiny=0, ... ) {

    P1 <- P1_prob_other(R0 = R0, epsilon = epsilon,
                        burnout_prob_fun = burnout_prob_MS, k = k, N = N,
                        subdivisions = subdivisions, tiny=tiny, ... )
    return(P1)
}

##' Burnout probability based on Meerson and Sasorov (2009)
##'
##' @details
##' See equation (9) of Ballard et al (2016)
##'
##' The derivation assumes \eqn{N S_0 \gg 1}, where \eqn{S_0} is
##' defined in equation (9) of Ballard et al (2016).
##'
##' @seealso \code{\link{P1_prob_MS}}
##'
##' @inheritParams P1_prob_MS
##' @param debug print debugging info?
##' 
##' @importFrom stats integrate
##'
##' @references
##' \insertRef{Meerson2009}{burnout}
##' 
##' @export
##'
burnout_prob_MS <- function( R0, epsilon, N=10^6, subdivisions=1000L, debug = FALSE, tiny=0, ... ) {

    dfun <- function(x) {
        if (debug) cat(x, get(x), "\n")
    }

    ## choose units such that gamma+mu = 1, i.e., mean time infected is 1:
    beta <- R0
    mu <- epsilon
    gamma <- 1-epsilon

    K <- beta/mu
    delta <- 1 - 1/R0
    xm <- (-1/R0) * W0(-R0*exp(-R0)) - 1
    dfun(xm)
    
    integrand <- function(s) {
        s*(s+delta) / ((1+s)^2 * (s - (1-delta)*log(1+s))) -
            xm / ((1+xm)*(s-xm))
    }

    i_lwr <- integrand(0)
    dfun(i_lwr)

    i_upr <- integrand(xm)
    dfun(i_upr)
    
    Q1 <- try(stats::integrate(f=integrand, lower=0+tiny, upper=xm-tiny,
                                    subdivisions=subdivisions, ...)$value,
              silent = TRUE)
    ## FIX: this destroys the automatic vectorization:
    if (inherits(Q1, "try-error")) return(NA)

    ym <- (delta+xm)*xm/(1+xm)*(-xm/delta)^(K*delta) *
        exp(K*(xm+delta) - (1+1/xm)*Q1)

    logC <- log(ym) + log(delta) - log(2*pi) -log1p( -delta)
    C <- exp(logC)
    
    S0 <- C * sqrt( 2*pi / (K*delta) )

    p0 <- exp(-N * S0)

    return(p0)
}
