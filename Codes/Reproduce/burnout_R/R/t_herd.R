
##' Duration of herd immunity
##'
##' The time \eqn{t_{\rm H}} (after the first epidemic) from when the
##' proportion susceptible drops below the herd immunity threshold
##' (defined via \eqn{X(t_{\rm H}) = 1/{\cal R}_0}) until the
##' proportion susceptible next rises above the threshold.
##'
##' \deqn{
##' 	t_{\rm H} =  \frac{1}{\varepsilon}
##'     \ln\Big( \frac{1-x_{\rm in}}{1 - \frac{1}{{\cal R}_0}} \Big)
##' }
##'
##' @inheritParams P1_prob
##' @inheritParams fizzle_prob
##'
##' @examples
##' \dontrun{
##' op <- par(mfrow = c(1,2))
##' plot_t_herd(col = c("darkred", "darkgreen", "darkblue"), main="")
##' epsilonseq <- seq(0,0.03,length=1001)
##' plot(epsilonseq, t_herd(R0=1.5,epsilonseq), type="l", lwd=2, bty="L",
##'      las=1, col="darkred", ylim=c(0,60),
##'      xlab = expression(epsilon), ylab = expression(t[H]))
##' lines(epsilonseq, t_herd(R0=3,epsilonseq), lwd=2, col="darkgreen")
##' lines(epsilonseq, t_herd(R0=4.5,epsilonseq), lwd=2, col="darkblue")
##' par(mfrow = op)
##' title(main=latex2exp::TeX("Duration of herd immunity $t_{H}$"))
##' }
##' @export
##' 
t_herd <- function(R0, epsilon, xin = x_in(R0,epsilon)) {
    ##return((1/epsilon) * log((1-xin)/(1-1/R0)))
    return((1/epsilon) * (log(1-xin) - log(1-1/R0)))
}

##' Plot herd immunity duration
##'
##' @inheritParams t_herd
##' @param x_in_fun function to use to estimate \eqn{x_{\rm in}}
##' @param col,lwd,ylim,lty,main,xlab,ylab,... see \code{\link{graphical parameters}}
##' @param show.legend if \code{TRUE} then show legend
##' @param add if \code{TRUE} then add to existing plot
##'
##' @seealso \code{\link{t_herd}}
##'
##' @importFrom graphics title
##'
##' @details It is best to use \code{\link{x_in_exact}} when
##'     calculating \code{t_herd}.
##'
##' @export
##'
##' @examples
##' ## plot tH with x_in as dotted under tH with x_in_exact as solid:
##' plot_t_herd(x_in_fun = x_in, lty="dotted")
##' plot_t_herd(add=TRUE)
##'
plot_t_herd <- function(R0 = exp(seq(log(1.001),log(5),length=1001)),
                        epsilon = c(0.01, 0.02, 0.03)
                      , x_in_fun = x_in_exact
                      , col = 1:length(epsilon)
                        ##col = c("darkred", "darkgreen", "darkblue")
                      , lwd=2
                        ##, log="x"
                      , lty="solid"
                      , ylim=c(0,60)
                      , xlab=expression(R[0])
                      , ylab=expression(t[H])
                      , main=latex2exp::TeX("Duration of herd immunity $t_{H}$")
                      , show.legend = TRUE
                      , add=FALSE
                      , ...
                        ) {
    tH_fun <- function(R0,epsilon) {
        t_herd(R0,epsilon,xin=x_in_fun(R0,epsilon))
    }
    if (!add) {
        plot(R0, tH_fun(R0,epsilon=epsilon[1]), type="n", lwd=lwd, bty="L",
             las=1, col=col[1], ylim=ylim,
             xlab = xlab, ylab = ylab, ...)
        title(main = main)
    }
    for (iepsilon in seq_along(epsilon)) {
        lines(R0, tH_fun(R0,epsilon=epsilon[iepsilon]), lwd=2, col=col[iepsilon], lty=lty)
    }
    if (show.legend) legend("topright", bty="n",
                            title=expression(epsilon),
                            legend = epsilon, col = col, lwd=lwd)
}
