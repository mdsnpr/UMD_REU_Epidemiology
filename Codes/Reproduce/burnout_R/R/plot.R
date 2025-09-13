##' Plot persistence probability
##'
##' @inheritParams P1_prob
##' @param Rmin,Rmax minimum and maximum value of \eqn{{\cal R}_0}
##' @param show.N,show.epsilon if \code{TRUE}, display parameter value in legend
##' @param P1_fun function that computes \eqn{{\cal P}_1}
##' @param add if \code{TRUE}, add to existing plot
##' @param col,lwd,log standard \link{graphical parameters}
##' @param ... additional parameters passed to \code{\link{plot}}
##'
##' @importFrom graphics axis
##' @importFrom graphics legend
##' @importFrom graphics lines
##' @importFrom graphics par
##' @importFrom latex2exp TeX
##'
##' @seealso P1_prob
##' 
##' @export
##'
##' @examples
##' op <- par(mfrow = c(2,2))
##' plot_P1(epsilon=0.01, N=10^4)
##' plot_P1(epsilon=0.01, N=10^5)
##' plot_P1(epsilon=0.01, N=10^6)
##' plot_P1(epsilon=0.01, N=10^7)
##' par(mfrow = op)
##' 
plot_P1 <- function(epsilon = 0.01, N = 10^6,
                    ##npts = 1001,
                    Rmin = 1.001,
                    Rmax = 64,
                    P1_fun = P1_prob, add=FALSE,
                    col="black", lwd=2, log="x",
                    show.N = TRUE,
                    show.epsilon = FALSE,
                    ... ) {

    ## naive R0 values:
    ##xvals <- exp(seq(log(Rmin), log(Rmax), length=npts))

    ## sprinkle R0 values where they are needed to get a smooth curve
    Rseq <- c(seq(Rmin,1.04,length=200),
              seq(1.04,1.1,length=100),
              seq(1.1,2,length=100),
              seq(2,Rmax,length=1000))
    
    P1vals <- P1_fun(Rseq, epsilon=epsilon, N=N)
    if (any(Im(P1vals) != 0)) {
        warning("results are complex")
        print(summary(Im(P1vals)))
        message("max(Im(P1)) = ", max(Im(P1vals)), "\n")
        P1vals <- Re(P1vals)
        cat("max(Re(P1)) = ", max(P1vals), "\n")
    }
    if (add) {
        lines(Rseq, P1vals, las=1, col=col, type="l", lwd=lwd,
              ...)
    } else {
        plot(Rseq, P1vals, las=1, col=col, type="l", lwd=lwd, log=log,
             bty="L", xaxs="i", xaxt="n",
             xlim=c(2^(-1/8),Rmax), ylim=c(10^(-4.5),1),
             xlab = expression(R[0]),
             ylab = "persistence probability",
             ...)
        xticks <- 2^(0:6)
        axis(side=1, at=xticks, labels=xticks)
        legend.text <- c()
        N.string <- sprintf("$N = 10^{%g}$", log10(N))
        epsilon.string <- sprintf("$\\varepsilon = %g$", epsilon)
        if (show.N)
            legend.text <- c(legend.text, latex2exp::TeX(N.string))
        if (show.epsilon)
            legend.text <- c(legend.text, latex2exp::TeX(epsilon.string))
        legend("topleft", legend = legend.text, bty="n")
    }
}

##' Plot comparison of two functions.
##'
##' \code{\link{plot}} method for \code{\link{compare_funs}} objects.
##' Typical use: exact vs approximate \eqn{q}.
##'
##' @inheritParams plot_P1
##' @param x \code{\link{compare_funs}} object, i.e., a
##'     \code{\link{data.frame}} with particular structure
##' @param show.plots logical: if \code{FALSE} then return plots in a
##'     list without displaying them
##'
##' @seealso \code{\link{compare_funs}}, \code{\link{q_exact}},
##'     \code{\link{q_approx}}, \code{\link{x_in}}
##'
##' @import dplyr
##' @import ggplot2
##
## The following avoids check() from complaining "no visible binding
## for global variable ‘type’".
## https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887/2
## https://www.r-bloggers.com/2019/08/no-visible-binding-for-global-variable/
##' @importFrom rlang .data
##' 
##' @export
plot.compare_funs <- function(x, show.plots=TRUE, ...) {

    ## make title
    N <- attr(x,"N")
    title.text <- if (is.null(N)) "" else sprintf("$N = 10^{%d}$",
                                                  log10(N))

    ## scatter plot coloured by epsilon:
    scatter.plot <- (x %>% ggplot()
        + geom_abline(slope=1, colour="magenta")
        + geom_point(aes(x=.data$f1, y=.data$f2, colour=.data$epsilon),
                     size = 0.5, alpha=0.5)
        + ggtitle(latex2exp::TeX(paste("scatter plot:", title.text)))
    )

    title.text <- paste0(title.text, ", facet by $\\epsilon$")
    f1name <- attr(x,"f1name")
    f2name <- attr(x,"f2name")
    line.plot <- (x %>% ggplot()
        + geom_point(aes(x=.data$R0, y=.data$f2), colour = "black")
        + geom_point(aes(x=.data$R0, y=.data$f1), colour = "red", size=0.25)
        + facet_wrap(~.data$epsilon, scales="free_y")
        + scale_x_continuous(trans='log2')
        + ggtitle(latex2exp::TeX(paste("line plot:", title.text, " (red = ", f1name, ", black = ", f2name, ")")))
    )
    ##FIX: allowing for 3rd data column... not working as expected...
    ##if ("f0" %in% names(x)) line.plot  <- line.plot + geom_point(aes(x=.data$R0, y=.data$f0), colour = "blue", size=0.1)

    relative.plot <- (x %>% ggplot()
        + geom_point(aes(x=.data$R0, y=(.data$f2-.data$f1)/.data$f1),
                         colour="black")
        + facet_wrap(~.data$epsilon, scales="free_y")
        + scale_x_continuous(trans='log2')
        + ggtitle(latex2exp::TeX(paste("relative plot:", title.text)))
    )

    if (show.plots) {
        print(scatter.plot)
        print(line.plot)
        print(relative.plot)
    }
    out <- list(scatter = scatter.plot, line = line.plot,
                relative = relative.plot)
    return(invisible(out))
}
