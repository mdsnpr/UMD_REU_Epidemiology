##' Compare functions of \eqn{{\cal R}_0} and \eqn{\varepsilon}
##'
##' @details The default \code{R0} value of \code{NULL} yields a
##'     sensible sprinkling of \code{nR0} \eqn{{\cal R}_0} values.
##'     Typical uses are to compare exact vs approximate
##'     calculations of Kendall's \eqn{q}, or different approximations
##'     of \eqn{x_{\rm in}}.
##'
##' @inheritParams peak_prev
##' @inheritParams P1_prob
##' @param fun1,fun2 functions of \code{R0} and \code{epsilon} to be
##'     compared
##' @param Rmin,Rmax minimum and maximum values of \eqn{{\cal R}_0}
##' @param nR0,nepsilon number of values of \code{R0} and
##'     \code{epsilon}
##' @param show.progress logical: if \code{TRUE} then spew each
##'     element of the data frame to \code{stdout} as it is computed
##' @param ... additional parameters passed to \code{fun1} and
##'     \code{fun2}
##'
##' @seealso \code{\link{q_exact}}, \code{\link{q_approx}},
##'     \code{\link{x_in}}
##'
##' @export
##'
##' @return data frame with \eqn{\code{nR0*nepsilon}} rows
##'
##' @examples
##' (cf <- compare_funs(nR0=101))
## too slow:
## (cf <- compare_funs(x_in_exact, x_in, epsilon=0.01, Rmax=2))
##
compare_funs <- function(fun1 = q_exact, fun2 = q_approx,
                         R0 = NULL, Rmin = 1.001, Rmax = 64, nR0 = 101,
                         epsilon = 10^(-(4:1)), nepsilon = length(epsilon),
                         N = NULL,
                         ...,
                         show.progress=FALSE) {

    if (is.null(R0)) {
        if (Rmax > 8) {
            ## sprinkle R0 values where they are needed to get a smooth curve
            n1 <- round(0.2*nR0)
            n2 <- round(0.1*nR0)
            R0 <- c(seq(Rmin,1.04,length=n1),
                    seq(1.04,1.1,length=n2),
                    seq(1.1,2,length=n2),
                    exp(seq(log(2),log(8),length=n2)),
                    seq(8,Rmax,length=(nR0-3*n2-n1)))
        } else {
            R0 <- exp(seq(log(Rmin), log(Rmax), length=nR0))
        }
    } else {
        nR0 <- length(R0)
    }

    raw <- expand.grid(R0, epsilon)
    names(raw) <- c("R0", "epsilon")
    ## FIX: The following would be much cleaner if dplyr::mutate were
    ##      applied to the raw data frame, but the loop allows me to
    ##      print the table as far as it can be be computed, which is
    ##      helpful.
    nn <- nrow(raw)
    f1 <- f2 <- rep(NA,nn)
    f1name <- deparse(substitute(fun1))
    f2name <- deparse(substitute(fun2))
    if (show.progress) cat(sprintf("%s\t%s\t%s\t%s\t%s\n",
                                   "i", "R0", "epsilon", f1name, f2name))
    for (i in 1:nn) {
        r <- raw[i,"R0"]
        e <- raw[i,"epsilon"]
        if (show.progress) cat(sprintf("%d\t%g\t%g", i, r, e))
        if (is.null(N)) {
            f1i <- try(fun1(R0=r, epsilon=e, ...))
        } else {
            f1i <- try(fun1(R0=r, epsilon=e, N=N, ...))
        }
        if ("try-error" %in% class(f1i)) f1i <- NA
        if (show.progress) cat(sprintf("\t%g", f1i))
        if (is.null(N)) {
            f2i <- try(fun2(R0=r, epsilon=e, ...))
        } else {
            f2i <- try(fun2(R0=r, epsilon=e, N=N, ...))
        }
        if ("try-error" %in% class(f2i)) f2i <- NA
        if (show.progress) cat(sprintf("\t%g\n", f2i))
        f1[i] <- f1i
        f2[i] <- f2i
    }
    ## FIX: not sure why this is failing...
    ## nIme <- sum(Im(qexact) != 0)
    ## nIma <- sum(Im(qapprox) != 0)
    ## if (nIme != 0) warning(nIme, " complex values of qexact")
    ## if (nIma != 0) warning(nIma, " complex values of qapprox")
    f1 <- Re(f1) # FIX: avoid this hack
    f2 <- Re(f2) # FIX: avoid this hack
    dd <- cbind(raw, f1, f2)
    class(dd) <- c("compare_funs", "data.frame")
    ## attr(dd,"nIme") <- nIme
    ## attr(dd,"nIma") <- nIma
    attr(dd,"N") <- N
    ## save names of functions that were compared:
    attr(dd,"f1name") <- deparse(substitute(fun1))
    attr(dd,"f2name") <- deparse(substitute(fun2))
    return(dd)
}
