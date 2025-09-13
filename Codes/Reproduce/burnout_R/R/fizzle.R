##' Fizzle probability
##' 
##' Probability that a disease will go extinct before causing a
##' substantial outbreak.
##'
##' The fizzle probability is
##' 
##' \deqn{1 - p_k =
##'   \left\{
##'   \begin{array}{rl}
##'      1\,, & 0 \le {\cal R}_0 \le 1, \\
##'      \left(\frac{1}{{\cal R}_0}\right)^{k} \,, & 1\le{\cal R}_0.
##'   \end{array}
##'   \right.
##' }
##'
##' @param R0 basic reproduction number (\eqn{{\cal R}_0})
##' @param k initial number of infected individuals (\eqn{k})
##'
##' @return real number between 0 and 1
##' @export
##'
##' @seealso \code{\link{fizzle_prob}}
##'
##' @examples
##' fizzle_prob(R0=c(1/2,2))
##' fizzle_prob(R0=2, k=2)
##' fizzle_prob(R0=c(2,4,8))
##' fizzle_prob(R0=2, k=c(1,2))  # FIX: fails with vector k
##' 
fizzle_prob <- function(R0, k=1) {
    ifelse( R0 < 1, 1, (1/R0)^k )
}

##' Not fizzle probability
##' 
##' Probability that a disease will \emph{not} go extinct before
##' causing a substantial outbreak.
##'
##' The probability of not fizzling is
##' 
##' \deqn{p_k =
##'   \left\{
##'   \begin{array}{rl}
##'      0\,, & 0 \le {\cal R}_0 \le 1, \\
##'      1 - \left(\frac{1}{{\cal R}_0}\right)^{k} \,, & 1\le{\cal R}_0.
##'   \end{array}
##'   \right.
##' }
##'
##' @seealso \code{\link{fizzle_prob}}
##'
##' @inheritParams fizzle_prob
##'
##' @return real number between 0 and 1
##' @export
##'
##' @examples
##' not_fizzle_prob(R0=c(1/2,2))
##' not_fizzle_prob(R0=2, k=2)
##' not_fizzle_prob(R0=c(2,4,8))
##' not_fizzle_prob(R0=2, k=c(1,2))  # FIX: fails with vector k
##' 
not_fizzle_prob <- function(R0, k=1) {
    1 - fizzle_prob(R0, k)
}

##' Fizzle time
##'
##' Time \eqn{t_\delta} for which probability of fizzle after
##' \eqn{t_\delta} is \eqn{< \delta}
##'
##' \deqn{
##' 	t_{\delta} =  \frac{1}{{\cal R}_{0}-1} \ln\left(\frac{(1-\delta)^{-\frac{1}{k}}-\frac{1}{{\cal R}_{0}}}{(1-\delta)^{-\frac{1}{k}}-1}\right).
##' }
##'
##' @inheritParams P1_prob
##' @inheritParams fizzle_prob
##' @param delta probability threshold
##' @seealso plot_fizzle_time
##' @aliases t_delta
##'
##' @export
##' 
##' @examples
##' fizzle_time(R0 = 2, k = 1, delta = 0.01)
##' plot_fizzle_time()
##' 
fizzle_time <- function(R0, k, delta) {
    tmp <- (1-delta)^(-(1/k))
    fizz.time <- 1/(R0-1) * log((tmp - 1/R0) / (tmp - 1))
    return(fizz.time)
}

##' Plot fizzle time
##'
##' @inheritParams fizzle_time
##' @inheritParams base::plot
##' @param col,lwd,log,... see \code{\link{graphical parameters}}
##'
##' @seealso \code{\link{fizzle_time}}
##'
##' @importFrom graphics title
##' 
##' @export
##'
##' @examples
##' plot_fizzle_time(delta = 0.01)
##'
## FIX: this would be better using facet_wrap if there is more than one delta value
##
plot_fizzle_time <- function(R0 = exp(seq(log(1.001),log(2),length=1001)),
                             k = 1:3, delta = 0.1
                           , col = 1:length(k)
                             ##col = c("darkred", "darkgreen", "darkblue")
                           , lwd=2
                           , log="x"
                           , ...
                             ) {
    plot(R0, fizzle_time(R0,k=k[1],delta=delta), type="l", lwd=lwd, bty="L",
         las=1, log=log, col=col[1], ylim=c(0,30),
         xlab = expression(R[0]), ylab = expression(t[delta]), ...)
    title(main = latex2exp::TeX(sprintf(
          "Time $t_\\delta$ for which probability of fizzle after $t_\\delta$ is $< \\delta = %g$", delta)))
    for (ik in 2:length(k)) {
        lines(R0, fizzle_time(R0,k=k[ik],delta=delta), lwd=2, col=col[ik])
    }
    legend("topright", bty="n", title=expression(I[0]),
           legend = k, col = col, lwd=lwd)
}
