##' Peak prevalence
##'
##' Approximation of peak prevalence in the SIR model with vital dynamics.
##'
##' The peak prevalence is
##'
##' \deqn{
##'   y_{\rm max} \approx 1-
##'     \frac{1}{{\cal R}_{0}}\big(1 +\ln{{\cal R}_{0}}\big)
##'     - \varepsilon
##'   \Big(1-\frac{1}{{\cal R}_0}\Big)
##'   \frac{(\ln{{\cal R}_0})\big/\big(1-\frac{1}{{\cal R}_0}\big) - 1}
##'   {(\ln{{\cal R}_0})\big/\big(1-\frac{1}{{\cal R}_0}\big) + {\cal R}_0} \,.
##' }
##'
##' The expression is exact for \eqn{\varepsilon = 0}.
##'
##' @seealso \code{\link{peak_prev_nvd}}
##'
##' @inheritParams fizzle_prob
##' @param epsilon mean infectious period (\eqn{1/(\gamma+\mu)}) as a
##'     proportion of mean lifetime (\eqn{1/\mu}; \eqn{\varepsilon =
##'     0} corresponds to infinite lifetime, i.e., no mortality)
##' @param eta mean infectious period (\eqn{1/(\gamma+\mu)}) as a
##'     proportion of mean duration of immunity (\eqn{1/\delta};
##'     \eqn{\eta = 0} corresponds to permanent immunity
##' @param alpha antigenic evolution enhancement factor for effective
##'     immunity waning; rate of immune decay is \eqn{\delta(1+\alpha
##'     Y)}.
##' @param xi initial susceptible proportion \eqn{x_{\rm i}}
##' @param yi initial infective proportion \eqn{y_{\rm i}}
##'
##' @return real number between 0 and 1
##' @export
##'
##' @examples
##' peak_prev(R0=2)
##' peak_prev(R0=2, epsilon=c(0, 0.001, 0.01, 0.1, 0.9))
##' peak_prev(R0=c(2,4), epsilon=0.1)
##' peak_prev(R0=2, epsilon=0.001, eta=c(0.01, 0.1))
##' peak_prev(R0=c(2,4), epsilon=c(0,0.1)) # FIX: fails when both are vectors
##'
peak_prev <- function( R0, epsilon=0, eta.i=0, eta.a=0, xi=1, yi=0 ) {
    pc <- 1 - 1/R0 # p_crit
    lnR <- log(R0)
    ##ymax <- 1 - (1/R0)*(1 + lnR) - # orig version assuming (xi,yi)=(1,0)
    ymax <- peak_prev_nvd(R0=R0, xi=xi, yi=yi) -
        ## FIX: the corrections assume (xi,yi)=(1,0):
        (epsilon+eta.i+eta.a)*pc*(lnR/pc - 1)/(lnR/pc + R0) +
        ((eta.i+eta.a)/R0)*(lnR - pc)
    return(ymax)
}

##' Peak prevalence
##'
##' Exact peak prevalence in the SIR model with \emph{no vital dynamics} (nvd).
##'
##' \deqn{
##'   y_{\rm max} = 1 -
##'     \frac{1}{{\cal R}_{0}}\big(1 +\ln{{\cal R}_{0}}\big) \,.
##' }
##'
##' @seealso \code{\link{peak_prev}}, \code{\link{x_in}}
##'
##' @param ... (ignored) for compatibility 
##' with the more accurate \code{\link{peak_prev}} function
##' that takes additional arguments
##'
##' @inheritParams peak_prev
##'
##' @return real number between 0 and 1
##' @export
##'
peak_prev_nvd <- function( R0, xi=1, yi=0, ... ) {
    ##ymax <- 1 - (1/R0)*(1 + log(R0)) # orig version assuming (xi,yi)=(1,0)
    ymax <- yi+xi - (1/R0)*(1 + log(R0*xi))
    return(ymax)
}
