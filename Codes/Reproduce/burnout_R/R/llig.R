##' log lower incomplete gamma function
##'
##' Compute the logarithm of the lower incomplete gamma function via
##' \code{\link[stats]{pgamma}}.
##'
##' The \emph{lower incomplete gamma function} ([NIST
##' 8.2.1](https://dlmf.nist.gov/8.2#E1)) is
##'
##' \deqn{\texttt{lig}(a,z) = \int_0^z t^{a-1}e^{-t}\,{\rm d}t \,,
##' \qquad \Re{a}>0.}
##'
##' We are interested only in real arguments, and we write
##' \eqn{\texttt{lig}(a,x)} rather than the usual \eqn{\gamma(a,x)}
##' for this function to avoid confusion with the standard notation
##' for the recovery rate in the SIR model.
##'
##' Since \eqn{\Gamma(a) = \texttt{lig}(a,\infty)}, the
##' \emph{normalized} area under the curve \eqn{t^{a-1}e^{-t}} ([NIST
##' 8.2.4](https://dlmf.nist.gov/8.2#E4)) is
##' 
##' \deqn{P(a,x) = \frac{\texttt{lig}(a,x)}{\Gamma(a)} .}
##' 
##' This is the cumulative distribution of the Gamma distribution
##' (\code{\link[stats]{pgamma}}), which we can calculate and express
##' on the log scale via
##' 
##' \deqn{\fcolorbox{red}{blue}{\code{pgamma(shape=a, q=x, lower.tail = TRUE, log.p = TRUE)}}}
##' 
##' Since \eqn{\log(\Gamma(a))} can be calculated directly on the
##' log scale via \code{\link[base]{lgamma}}, we can do the entire
##' calculation on the log scale via
##'
##' \deqn{\texttt{llig}(a,x) = \log\big(P(a,x)\big) + \log\big(\Gamma(a)\big)}
##' 
##' This approach allows us to compute \code{llig} for very large
##' values of \code{a} and \code{x} without overflow.
##' 
##' If you need to convert to the linear scale and have values of
##' \eqn{\texttt{lig}(a,x) > 10^{308}} you're going to be in trouble
##' (although there are packages such as
##' \code{\link[Brobdingnag]{Brobdingnag}} that handle these kinds of
##' large numbers.  For our purposes with the \code{burnout} package,
##' we can use logspace addition and subtraction and avoid ever
##' working with excessively large numbers on the linear scale.
##'
##' There are other implementations of the lower incomplete gamma
##' function, which fail for large arguments (e.g.,
##' \code{\link[expint]{gammainc}(a,x)},
##' \code{\link[gsl]{gamma_inc}(a,x)}).
##'
##' @note Champredon _et al_ (2018) \doi{10.1137/18M1186411}
##'     used the notation \eqn{{\cal G}} for \code{lig}.
##'
##' @param a shape parameter (non-negative real number)
##' @param x quantile parameter (non-negative real number)
##'
##' @importFrom stats pgamma
##'
##' @return real number
##' @export
##'
##' @seealso \code{\link[stats]{pgamma}}, \code{\link[base]{lgamma}}
##'
##' @aliases log_lig
##' 
##' @examples
##' llig( a = 1, x = 1 )
##' llig( a = 10^200, x = 10^200 )
##' xseq <- seq(0,10,length=100)
##' aval <- 10
##' plot(xseq, llig(a=aval, x=xseq), type="o", pch=21, bg="darkred",
##'      cex=0.5, las=1)
##' abline(h = lgamma(aval), col="cyan")
##' 
llig <- function(a, x) {
    pg <- try(pgamma(shape = a, q = x, lower.tail = TRUE, log.p = TRUE),
              silent = TRUE) # we'll spew error message below
    if ("try-error" %in% class(pg)) {
        ## provide helpful info rather than just crashing:
        message("\nllig: shape = a = ", paste(a,collapse=", "),
                ",\n      q = x = ", paste(x,collapse=", "))
        if (any(Im(a) != 0)) message("llig: a has non-zero imaginary part")
        if (any(Im(x) != 0)) message("llig: x has non-zero imaginary part")
        if (all(Im(x) == 0) && any(x < 0)) message("llig: x is negative")
        ## now go ahead and crash with default error message:
        ## pgamma(shape = a, q = x, lower.tail = TRUE, log.p = TRUE)
        ## actually, instead return the result using the real parts
        ## of a and x as input:
        pg <- pgamma(shape = Re(a), q = Re(x), lower.tail = TRUE,
                     log.p = TRUE)
        return(pg + lgamma(Re(a)))
    }
    return(pg + lgamma(a))
}
