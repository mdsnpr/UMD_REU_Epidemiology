##' Persistence probability after \eqn{m}'th epidemic wave
##'
##' Probability of surviving the trough after the \eqn{m}'th epidemic.
##'
##' \deqn{
##'     {{\cal P}_m}({\cal R}_0,\varepsilon,k,N) =
##'       p_k\prod_{j=1}^m\big( 1 - q(x_{{\rm in},j})^{Ny^\star} \big) \,.
##' }
##'
##' @seealso \code{\link{P1_prob}}, \code{\link{fizzle_prob}},
##'     \code{\link{x_in}}, \code{\link{llig}},
##'     \code{\link{q_approx}}, \code{\link{q_exact}}
##'
##' @inheritParams P1_prob
##' @param m epidemic wave number
##'
##' @return real number between 0 and 1
##' @importFrom Rdpack reprompt
##'
##' @note This function is recursive and ultimate calls
##'     \code{\link{P1_prob}}.
##'
##' @references
##' \insertRef{Kendall1948b}{burnout}
##' 
##' @export
##'
##' @examples
##' \dontrun{
##' Pm_prob(R0=2, epsilon=0.01, m=2)
##' Pm_prob(R0=2, epsilon=0, m=2)
##' Pm_prob(R0=c(2,4), epsilon=0.01, m=2)
##' Pm_prob(R0=2, epsilon=c(0.01, 0.1), m=2)
##' curve(Pm_prob(x,epsilon=0.001, m=2), from=1.01, to=2, las=1, n=1001, ylim=c(0,1))
##' curve(Pm_prob(x,epsilon=0.01, m=2), from=1.01, to=2, las=1, add=TRUE, col="magenta", n=1001)
##' curve(Pm_prob(x,epsilon=0.1, m=2), from=1.01, to=2, las=1, add=TRUE, col="cyan", n=1001)
##' xvals <- exp(seq(log(1.01), log(64), length=1001))
##' Pmvals <- Pm_prob(xvals,epsilon=0.01, m=2)
##' ## results are complex:
##' max(Im(Pmvals))
##' Pmvals <- Re(Pmvals)
##' plot(xvals, Pmvals, las=1, col="blue", type="l", lwd=2, log="x")
##' cat("max(Pm) = ", max(Pmvals), "\n")
##' }
##'
Pm_prob <- function(R0, epsilon, m, k=1, N=10^6,
                    xin = x_in(R0,epsilon),
                    q_fun = q_approx,
                    ystar = epsilon * (1 - 1/R0)
                    ) {
    if (m == 1) return(P1_prob(R0, epsilon = epsilon,
                               k = k, N = N, xin = xin,
                               q_fun = q_fun, ystar = ystar))
    ## we know m > 1:
    xim <- x_i(R0, j=m)
    xinm <- x_in(R0, epsilon, xi=xim)
    ## burnout probability (conditional on not fizzling):
    burn <- burnout_prob(R0=R0, epsilon=epsilon, N=N, xin = xinm,
                         q_fun=q_fun, ystar=ystar)
    message("m = ", m, ", xim = ", xim, ", xinm = ", xinm,
            ", burn = ", burn)
    return((1-burn) *
           Pm_prob(R0, epsilon, m-1, k, N, xin, q_fun, ystar))
}
