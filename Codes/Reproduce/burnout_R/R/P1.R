##' Persistence probability
##'
##' Probability of surviving the trough after the first epidemic.
##'
##' If \eqn{p_k} is the probability of \emph{not}
##' \code{\link[=fizzle_prob]{fizzling}} if \eqn{k} individuals are
##' initially infected, then the probability of burning out
##' conditional on not fizzling is
##'
##' \deqn{ p_k q^{Ny^\star} \,,}
##'
##' where \eqn{q} is Kendall's \eqn{q}, the expected total population
##' size is \eqn{N}, and the equilibrium prevalence is
##'
##' \deqn{y^\star = \varepsilon\Big(1 - \frac{1}{{\cal R}_0}\Big) \,.}
##'
##' Hence, the probability of \emph{either} fizzling before a first
##' major outbreak \emph{or} burning out after a first major outbreak
##' is
##'
##' \deqn{
##'      (1-p_k) + p_k q^{Ny^\star} =
##'      1 - p_k\big( 1 - q^{Ny^\star} \big) \,.
##' }
##'
##' Consequently, the probability of \emph{persisting} beyond the
##' first epidemic is
##'
##' \deqn{
##'     {{\cal P}_1}({\cal R}_0,\varepsilon,k,N) =
##'       p_k\big( 1 - q^{Ny^\star} \big) \,.
##' }
##'
##' @seealso \code{\link{fizzle_prob}}, \code{\link{x_in}},
##'     \code{\link{llig}}, \code{\link{q_approx}}, \code{\link{q_exact}}
##'
##' @inheritParams peak_prev
##' @inheritParams q_approx
##' @inheritParams fizzle_prob
##' @param N population size
##' @param q_fun function to compute Kendall's \eqn{q}
##' @param ystar prevalence at boundary layer (equilibrium prevalence
##'     by default)
##'
##' @return real number between 0 and 1
##' @importFrom Rdpack reprompt
##'
##' @references
##' \insertRef{Kendall1948b}{burnout}
##'
##' @export
##'
##' @examples
##' P1_prob(R0=2, epsilon=0.01)
##' P1_prob(R0=2, epsilon=0)
##' P1_prob(R0=c(2,4), epsilon=0.01)
##' P1_prob(R0=2, epsilon=c(0.01, 0.1))
##' curve(P1_prob(x,epsilon=0.001), from=1.01, to=2, las=1, n=1001, ylim=c(0,1))
##' curve(P1_prob(x,epsilon=0.01), from=1.01, to=2, las=1, add=TRUE, col="magenta", n=1001)
##' curve(P1_prob(x,epsilon=0.1), from=1.01, to=2, las=1, add=TRUE, col="cyan", n=1001)
##' xvals <- exp(seq(log(1.01), log(64), length=1001))
##' P1vals <- P1_prob(xvals,epsilon=0.01)
##' ## results are complex:
##' max(Im(P1vals))
##' P1vals <- Re(P1vals)
##' plot(xvals, P1vals, las=1, col="blue", type="l", lwd=2, log="x")
##' cat("max(P1) = ", max(P1vals), "\n")
##'
P1_prob <- function(R0, epsilon,
              eta.i=0, eta.a=0, k=1, N=10^6,
              xin = x_in(R0=R0, epsilon=epsilon, eta.i=eta.i, eta.a=eta.a),
              q_fun = q_approx,
              ystar = eqm_prev(R0=R0, epsilon=epsilon, eta.i=eta.i, eta.a=eta.a)
              ) {
    fizz <- fizzle_prob(R0, k)
    ## pk = probability of not fizzling:
    notfizz <- 1 - fizz
    ## burnout probability (conditional on not fizzling):
    burn <- burnout_prob(R0=R0, epsilon=epsilon,
                         eta.i=eta.i, eta.a=eta.a,
                         N=N, xin = xin,
                         q_fun=q_fun, ystar=ystar)
    ## probability of not fizzling and then burning out:
    notfizz.and.burn <- notfizz * burn
    ## probability of either fizzling or burning out:
    fizz.or.burn <- fizz + notfizz.and.burn
    ## persist after neither fizzling nor burning out:
    P1 <- 1 - fizz.or.burn # P1 <- pk * (1 - q^(N*ystar))
    return(P1)
}

##' Burnout probability based on Kendall's \eqn{q}
##'
##' \deqn{
##'     q^{Ny^\star}
##' }
##'
##' @seealso \code{\link{fizzle_prob}}, \code{\link{x_in}},
##'     \code{\link{llig}}, \code{\link{q_approx}}, \code{\link{q_exact}}
##'
##' @inheritParams P1_prob
##'
##' @return real number between 0 and 1
##' @importFrom Rdpack reprompt
##'
##' @references
##' \insertRef{Kendall1948b}{burnout}
##'
##' @export
burnout_prob <- function(R0, epsilon, eta.i=0, eta.a=0, N=10^6,
                         xin = x_in(R0=R0, epsilon=epsilon, eta.i=eta.i, eta.a=eta.a),
                         q_fun = q_approx,
                         ystar = eqm_prev(R0=R0, epsilon=epsilon, eta.i=eta.i, eta.a=eta.a)
                         ) {
    ## Kendall's q:
    q <- q_fun(R0, epsilon, eta.i, eta.a, xin)
    ## burnout probability (conditional on not fizzling):
    return( q^(N*ystar) )
}

##' Burnout probability using higher order corner/boundary
##' approximation of \eqn{x_{\rm in}}
##'
##' This is just a convenient wrapper of \code{\link{burnout_prob}} to
##' use \code{\link{x_in_hocb}} rather than \code{\link{x_in}}.
##'
##' @seealso \code{\link{burnout_prob}}
##'
##' @inheritParams burnout_prob
##'
##' @return real number between 0 and 1
##' @importFrom Rdpack reprompt
##'
##' @references
##' \insertRef{Kendall1948b}{burnout}
##'
##' @export
burnout_prob_hocb <- function(R0, epsilon, N=10^6,
                              xin = x_in_hocb(R0,epsilon),
                              q_fun = q_approx) {
    return( burnout_prob(R0, epsilon, N=N,
                         xin = xin,
                         q_fun = q_fun) )
}

##' Persistence probability via other approximations
##'
##' This computes \eqn{{\cal P}_1} using the burnout probability
##' approximations of van Herwaarden (1997) or Meerson and Sasorov
##' (2009).
##'
##' @seealso \code{\link{burnout_prob_MS}},
##'     \code{\link{burnout_prob_vanH}}, \code{\link{P1_prob}}
##'
##' @inheritParams P1_prob
##' @inheritParams stats::integrate
##' @param burnout_prob_fun function with which to calculate burnout
##'     probability (either \code{\link{burnout_prob_vanH}} or
##'     \code{\link{burnout_prob_MS}})
##' @param tiny avoid integration limits by this amount
##' @param ... additional arguments pass to
##'     \code{\link[stats]{integrate}}
##'
##' @importFrom stats integrate
##'
##' @export
##'
P1_prob_other <- function(R0, epsilon, burnout_prob_fun, k=1, N=10^6,
                          subdivisions=1000L, tiny=0, ... ) {

    ## Ballard et al (2016) use the notation p0:
    p0 <- burnout_prob_fun(R0 = R0, epsilon = epsilon, N = N,
                            subdivisions = subdivisions, tiny=tiny, ...)

    ## FIX: This is identical to the code in P1_prob() except that the
    ##      probability of burning out conditional on not fizzling
    ##      (p0) is calculated above via van H or MS's formulae.
    fizz <- fizzle_prob(R0=R0, k=k)
    ## pk = probability of not fizzling:
    notfizz <- 1 - fizz
    ## probability of not fizzling and then burning out:
    notfizz.and.burn <- notfizz * p0
    ## probability of either fizzling or burning out:
    fizz.or.burn <- fizz + notfizz.and.burn
    ## persist after neither fizzling nor burning out:
    P1 <- 1 - fizz.or.burn
    return(P1)
}
