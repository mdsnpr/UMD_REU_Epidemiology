##' Equilibrium prevalence \eqn{y^\star}
##'
##' @details If there is no effective waning from antigenic evolution,
##'     i.e., \eqn{\alpha=0}, then
##' \deqn{
##' 	y^\star = \frac{\varepsilon+\eta}{1+\eta}
##'       \big( 1 - \frac{1}{{\mathcal R}_0} \big)
##' }
##'
##' @inheritParams peak_prev
##' @export
##'
eqm_prev <- function(R0, epsilon, eta.i=0, eta.a=0) {
    if (eta.a == 0) {
        yeqm <- (epsilon+eta.i)/(1+eta.i)*(1 - 1/R0)
    } else {
        pc <- 1 - 1/R0
        ##tmp <- 1 + eta.i - 2*eta.a*pc ## SIGN ERROR
        tmp <- 1 + eta.i + 2*eta.a*pc
        yeqm <-tmp/(2*eta.a) - (1/(2*eta.a)) * sqrt(
            tmp^2 - 4*eta.a*pc*(epsilon + eta.i + eta.a*pc) )
    }
    return(yeqm)
}

##' Minimum \eqn{{\mathcal R}_0} for which \eqn{x_{\rm in}}
##' approximation is valid
##'
##' @details
##' \deqn{
##' 	{\mathcal R}_0 >
##'     -\frac{W_{-1}\left(-(1-\varepsilon)
##'      e^{-(1-\varepsilon)}\right)}{1-\varepsilon}
##' }
##'
##' @inheritParams x_in
##' @export
##'
R0_min_for_x_in <- function(epsilon) {
    mome <- -(1-epsilon) # mome = minus-one-minus-epsilon
    R0min <- Wm1(mome * exp(mome))/mome
    return(R0min)
}

##' Susceptibles at boundary layer (boundary/boundary layer
##' approximation) (scalar version)
##'
##' Fraction of hosts that are susceptible when the trajectory enters
##' the boundary layer at the end of the first epidemic.
##'
##' \deqn{
##' 	x_{\rm in} = -\frac{1}{{\cal R}_{0}} W_{0}\left(-{\cal R}_{0} e^{-{\cal R}_{0}\left(1-y^{\star}\right)}\right) + \varepsilon\, e^{{\cal R}_{0} y^{\star}}\big(E_{1}({\cal R}_{0} y^{\star}) - E_{1}({\cal R}_{0} y_{\rm max}) \big) \,.
##' }
##'
##' Here, \eqn{W_{0}} denotes the principal branch of the Lambert
##' \eqn{W} function
##' (\insertCite{Olver2010;textual}{burnout}, \insertCite{Corless1996;textual}{burnout},
##' [NIST 4.13](https://dlmf.nist.gov/4.13)),
##' \eqn{E_{1}(x) =
##' \int_{x}^{\infty} \frac{e^{-t}}{t}\,{\rm d}t} is the exponential
##' integral function (\insertCite{Olver2010;textual}{burnout}, [NIST
##' 6.2](https://dlmf.nist.gov/6.2)),
##' \eqn{y^\star=\varepsilon(1-\frac{1}{{\cal R}_0})}
##' is the equilibrium prevalence,
##' and \eqn{y_{\rm max}} is the \link[=peak_prev]{peak prevalence}.
##'
##' @note The initial conditions \eqn{(x_{\rm i},y_{\rm i}} enter the
##'     boundary layer component only via \eqn{y_{\rm max}}.
##'     \eqn{x_{\rm i}} also appears in the KM term involving
##'     \eqn{W_0}.
##'
##' @seealso \code{\link{peak_prev}}, \code{\link[gsl]{lambert_W0}},
##'     \code{\link[expint]{expint_E1}}, , \code{\link{x_in_cb}},
##'     \code{\link{x_in_exact}}
##'
##' @inheritParams peak_prev
##' @param peakprev_fun function of \eqn{{\cal R}_0} and
##'     \eqn{\varepsilon} to use to compute peak prevalence
##' @param ... additional arguments are ignored
##'
##' @return real number between 0 and 1
##' @importFrom expint expint_E1
##' @importFrom Rdpack reprompt
##'
##' @references
##' \insertRef{Corless1996}{burnout}
##'
##' \insertRef{Olver2010}{burnout}
##'
##' @export
##'
##' @examples
##' x_in(R0=2, epsilon=0.01)
##' x_in(R0=2, epsilon=0)
##' x_in(R0=c(2,4), epsilon=0.01)
##' x_in(R0=2, epsilon=c(0.01, 0.1))
##' curve(x_in(x,epsilon=0.001), from=1.01, to=5, las=1, n=1001)
##' curve(x_in(x,epsilon=0.01), from=1.01, to=5, las=1, add=TRUE, col="magenta", n=1001)
##' curve(x_in(x,epsilon=0.1), from=1.01, to=5, las=1, add=TRUE, col="cyan", n=1001)
##'
x_in_scalar <- function(R0, epsilon, eta.i=0, eta.a=0, xi = 1, yi=0,
                        ## FIX: using nvd version for now because version with
                        ##      with vital dynamics is currently wrong unless at DFE:
                        ##peakprev_fun = peak_prev,
                        peakprev_fun = peak_prev_nvd,
                        ...) {
    R0min <- R0_min_for_x_in(epsilon+eta.i+eta.a) # FIX: I set eta=eta.i+eta.a: OK?
    ## R0_min_for_x_in(1) yields NaN, hence:
    if (!is.finite(R0min) || R0 <= R0min) return(NA)
    yeqm <- eqm_prev(R0=R0, epsilon=epsilon, eta.i=eta.i, eta.a=eta.a) # equilibrium prevalence
    ##ymax <- peakprev_fun(R0=R0, epsilon=epsilon, eta.i=eta.i, eta.a=eta.a, xi=xi, yi=yi) # peak prevalence
    ## x_in expression in Fadeout-SIRS.tex assume ymax is the nvd version:
    ymax <- peakprev_fun(R0=R0, epsilon=0, eta.i=0, eta.a=0, xi=xi, yi=yi) # peak prevalence
    E1 <- expint::expint_E1
    xin <- ifelse (yeqm < ymax
                   ## see eq:xin and eq:X^0y.eps.etas in SIRiaS.tex:
                 , ## first term is the outer solution, i.e., KM
                   ## phase plane solution, inverted via Lambert W:
                   -(1/R0)*W0(-R0*xi*exp(-R0*(xi+yi-yeqm))) +
                   (epsilon+eta.i+eta.a)*exp(R0*yeqm)*(E1(R0*yeqm) - E1(R0*ymax)) -
                   (eta.i/R0)*(1 - exp(-R0*(ymax-yeqm))) +
                   (eta.a/R0)*((1 - exp(-R0*(ymax-yeqm)))/R0 +
                      exp(-R0*(ymax-yeqm))*(2-ymax) - (2-yeqm))
                 , ## epsilon correction is garbage so ignore it (first order
                   ## correction isn't good enough in this part of parameter
                   ## space where yeqm >= ymax):
                   -(1/R0)*W0(-R0*xi*exp(R0*(yeqm-xi)))
                   )
    return(xin)
}
##
## from original code for Xb_approx in Xb.R:
##  Xb <- -(1/R)*lambertW(-R*exp(R*(Yb-1)), b=0) +
##               eps*exp(R*Yb)*(expint_E1(R*Yb) - expint_E1(R*Ymax))

##' Susceptibles at boundary layer (vector version)
##'
##' @rdname x_in_scalar
##' @export
##'
x_in <- Vectorize(x_in_scalar, vectorize.args = c("R0", "epsilon", "eta.i", "eta.a", "xi"))

## not sure why checks are picking this up ?
utils::globalVariables("peak_prev")
utils::globalVariables("peak_prev_nvd")


##' Crude "Kermack-McKendrick" version of \eqn{x_{\rm in}}
##'
##' FIX: this is the same as \code{\link{x_in_crude}}, but coded
##' differently.
##'
##' This is just the case where \eqn{\varepsilon=0},
##' \deqn{x_{\rm in,KM}({\mathcal R}_0) = x_{\rm in}({\mathcal R}_0, 0)}
##'
##' @inheritParams x_in
##' @seealso \code{\link{x_in}}, \code{\link{x_in_crude}}
##'
##' @export
x_in_KM <- function(R0, ...) {
    xin <- x_in(R0, epsilon=0, ...)
    return(xin)
}

##' Susceptibles at boundary layer (crude approximation)
##'
##' @details
##' Argument \code{peakprev_fun} is ignored but is present so that the argument
##' list is the same as for \code{\link{x_in}}.
##'
##' @inheritParams x_in
##' @export
##'
x_in_crude <- function(R0, epsilon, peakprev_fun = NULL, ...) {
    yeqm <- eqm_prev(R0, epsilon, eta.i=0, eta.a=0) # equilibrium prevalence
    xin <- -(1/R0)*W0(-R0*exp(R0*(yeqm-1)))
    return(xin)
}

##' Susceptibles at boundary layer (corner/boundary layer approximation)
##'
##' @details
##' \deqn{
##' 	x_{\rm in} = 1 + \left(1-\frac{1}{{\mathcal R}_{0}}\right) W_{-1}\left(-\frac{1-x_{\rm f}}{1-\frac{1}{{\mathcal R}_{0}}} e^{-\frac{1-x_{\rm f}}{1-\frac{1}{{\mathcal R}_{0}}}}\left(\frac{y_{\rm max}}{y^*}\right)^{\frac{\varepsilon}{{\mathcal R}_{0}}\frac{1}{1-\frac{1}{{\mathcal R}_{0}}}}\right)
##' }
##'
##' @inheritParams x_in
##' @seealso \code{\link{x_in}}, \code{\link{x_in_cb}}, \code{\link{x_in_exact}}
##' @export
##'
x_in_cb <- function(R0, epsilon, peakprev_fun = peak_prev, ...) {
    yeqm <- eqm_prev(R0, epsilon, eta.i=0, eta.a=0) # equilibrium prevalence
    ymax <- peakprev_fun(R0=R0, epsilon=epsilon) # peak prevalence
    pc <- 1 - 1/R0 # p_crit
    Z <- final_size(R0)
    xf <- 1-Z
    xin <- ifelse (yeqm < ymax
                 , 1 + pc * Wm1(-(Z/pc)*exp(-(Z/pc)) * (ymax/yeqm)^((epsilon/R0)/pc))
                 , # else epsilon correction is garbage so ignore it
                   1 + pc * Wm1(-(Z/pc)*exp(-(Z/pc)))
                   )
    return(xin)
}

##' \eqn{{\tilde Y}_1(x_{\rm f})} (scalar version)
##'
##' @details
##' \deqn{
##' \tilde{Y}_{1}(x_{\rm f}) =  \int_{x_{\rm f}}^{1} \left[\left(\frac{1}{{\mathcal R}_{0}t} -1\right) \frac{1-t}{{\mathcal R}_{0}\,t\,Y_{0}(t)} + \frac{1-x_{\rm f}}{{\mathcal R}_{0}x_{\rm f}}\frac{1}{t-x_{\rm f}}\right] dt
##' }
##' where
##' \deqn{
##' 	Y_{0}(x) = 1-x+\frac{1}{{\mathcal R}_{0}}\ln{x}
##' }
##'
##' @note \code{subdivisions} and \code{...} are arguments to
##'     \code{\link[stats]{integrate}}.  The function \code{f} to
##'     which \code{...} is passed is the integrand.
##'
##' @inheritParams peak_prev
##' @inheritParams stats::integrate
##'
##' @param xf final size (\eqn{0<x_{\rm f}<1})
##' @param tiny amount by which to avoid integration limits
##'
##' @importFrom stats integrate
##' @seealso \code{\link{Ytilde_1}}
##'
##' @export
Ytilde_1_scalar <- function(xf, R0, tiny=1e-12, subdivisions=1000L, ...) {
    Y0 <- function(x) {1 - x + (1/R0)*log(x)}
    integrand <- function(t) {
        ## using tiny adjustment instead of cases:
        ##out <- ifelse(
        ##    t==1,
        ##    (1/(R0*t) - 1)/(R0*t) + (1-xf)/(R0*xf)/(t-xf),
        ##ifelse(
        ##    t==xf,
        ##    (1/(R0*t) - 1)*(1-t)/(R0*t*Y0(t)) + 1/(R0*xf),
        ##    (1/(R0*t) - 1)*(1-t)/(R0*t*Y0(t)) + (1-xf)/(R0*xf)/(t-xf)
        ##)
        ##)
        ##return(out)
        return( (1/(R0*t) - 1)*(1-t)/(R0*t*Y0(t)) - (1-xf)/(R0*xf)/(t-xf) )
    }
    the.integral <-
        try(stats::integrate(f=integrand, lower=xf+tiny, upper=1-tiny
                             , subdivisions=subdivisions
                             , ...)$value, silent = TRUE)
    if (inherits(the.integral, "try-error")) {
        warning("Ytilde_1: try error with xf = ", xf, "; returning NA")
        return(NA)
    }
    return(-the.integral)
}

##' \eqn{{\tilde Y}_1(x_{\rm f})} (vector version)
##'
##' @inheritParams Ytilde_1_scalar
##' @seealso \code{\link{Ytilde_1_scalar}}
##'
##' @export
##' @examples
##' R0 <- 2
##' xf <- 1 - final_size(R0)
##' Ytilde_1(xf,R0)
##'
Ytilde_1  <- Vectorize(Ytilde_1_scalar)

##' Standard final size
##'
##' for a newly invading infectious disease (so \eqn{I_0\to0}, \eqn{S_0\to1}).
##' \deqn{
##'   Z = 1 + \frac{1}{{\mathcal R}_{0}} W_{0}\left(-{\mathcal R}_{0} e^{-{\mathcal R}_{0}}\right)
##' }
##' @inheritParams x_in
##' @export
##' @examples
##' final_size(R0 = 2)
final_size <- function(R0) {1 + (1/R0)*W0(-R0*exp(-R0))}

## complement of final size (to avoid numeric problems), 1-Z
comp_final_size <- function(R0) {-(1/R0)*W0(-R0*exp(-R0))}

##' Susceptibles at boundary layer (higher order corner/boundary layer approximation)
##'
##' @details
##' \deqn{
##' 	x_{\rm in} = 1+\Big(1-\frac{1}{{\mathcal R}_{0}}\Big)  W_{-1}\left[-\frac{1-x_{\rm f}}{1-\frac{1}{{\mathcal R}_{0}}}\left(\frac{\left(\frac{1}{x_{\rm f}}-1\right)\left(\frac{1}{{\mathcal R}_{0}}-x_{\rm f}\right)}{y^*}\right)^{\frac{\varepsilon}{{\mathcal R}_{0}}\frac{1}{1-\frac{1}{{\mathcal R}_{0}}}} \exp\left({-\frac{1-x_{\rm f}}{1-\frac{1}{{\mathcal R}_{0}}}-\varepsilon\frac{1}{\left(\frac{1}{x_{\rm f}}-1\right)\left(1-\frac{1}{{\mathcal R}_{0}}\right)}\tilde{Y}_{1}(x_{\rm f})}\right)\right]
##' }
##' where
##' \deqn{
##' \tilde{Y}_{1}(x_{\rm f}) =  -\int_{x_{\rm f}}^{1} \left[\Big(\frac{1}{{\mathcal R}_{0}\,t} -1\Big) \frac{1-t}{{\mathcal R}_{0}\,t\,Y_{0}(t)} - \frac{1-x_{\rm f}}{{\mathcal R}_{0}\,x_{\rm f}}\frac{1}{t-x_{\rm f}}\right] {\rm d}t \,.
##' }
##' \deqn{
##' Y_{0}(x) = 1-x+\frac{1}{{\mathcal R}_{0}}\ln{x}
##' }
##' and \eqn{x_{\rm f} = 1 - Z}, where \eqn{Z} is the standard \code{\link{final_size}}.
##' @inheritParams x_in
##' @param debug print debugging output?
##' @seealso \code{\link{x_in}}, \code{\link{x_in_cb}}, \code{\link{x_in_exact}}
##' @export
##'
x_in_hocb <- function(R0, epsilon, peakprev_fun = peak_prev, debug = FALSE, ...) {

    dfun <- function(x) {
        if (debug) cat(x, get(x), "\n")
    }

    yeqm <- eqm_prev(R0, epsilon, eta.i=0, eta.a=0) # equilibrium prevalence [aka y_in]
    dfun("yeqm")
    ##ymax <- peakprev_fun(R0, epsilon) # peak prevalence [not needed for hocb]
    pc <- 1 - 1/R0 # p_crit
    dfun("pc")
    Z <- final_size(R0)
    dfun("Z")
    xf <- comp_final_size(R0)
    dfun("xf")
    pc_arg0 <- ((1/R0-xf)/xf)
    dfun("pc_arg0")
    pc_arg1 <- ((Z/yeqm)*pc_arg0)
    dfun("pc_arg1")
    pc_pow <- ((epsilon/R0)/pc)
    dfun("pc_pow")
    pc_arg <- -(Z/pc)*pc_arg1^pc_pow*
        exp(-(Z/pc) - epsilon*(xf/(Z*pc))*Ytilde_1(xf,R0, tiny = 0))
    dfun("pc_arg")
    xin <- 1 + pc * Wm1(pc_arg)
    dfun("xin")
    return(xin)
}

##' Plot \eqn{x_{\rm in}}
##'
##' @inheritParams x_in
##' @inheritParams compare_funs
##'
##' @param xvar string indicating whether \code{R0} or \code{epsilon}
##'     is the desired independent variable
##' @param col,log,lwd,xlim,ylim,lty,main,... see \code{\link{graphical parameters}}
##' @param add if \code{TRUE} then add to existing plot
##' @param xin_fun function of \code{R0} and \code{epsilon} that
##'     returns susceptible proportion at entry to boundary layer
##'
##' @seealso \code{\link{x_in}}, \code{\link{x_in_exact}}
##'
##' @importFrom graphics title
##'
##' @export
##'
##' @examples
##' \dontrun{
##' op <- par(mfrow = c(1,2))
##' plot_x_in(epsilon = 0.01)
##' plot_x_in(epsilon = 0.1)
##' par(mfrow = op)
##' op <- par(mfrow = c(1,3))
##' plot_x_in(R0 = 2, epsilon = seq(0,1,length=101), xvar="epsilon", log="")
##' plot_x_in(R0 = 5, epsilon = seq(0,1,length=101), xvar="epsilon", log="")
##' plot_x_in(R0 = 20, epsilon = seq(0,1,length=101), xvar="epsilon", log="")
##' par(mfrow = op)
##' }
##'
plot_x_in <- function(Rmin=1.0001, Rmax=20,
                      ##R0 = exp(seq(log(Rmin),log(Rmax),length=1001)),
                      R0 = seq(Rmin,Rmax,length=1001),
                      epsilon = c(0.01, 0.02, 0.03),
                      xvar = "R0"
                    , xin_fun = x_in
                    , col = 1:length(epsilon)
                      ##col = c("darkred", "darkgreen", "darkblue")
                    , lwd=2
                    , lty="solid"
                    , log="x"
                    , main="Susceptible proportion at\nentry to boundary layer"
                    , xlim = if (xvar == "R0") c(1,Rmax) else c(0,1)
                    , ylim = if (xvar == "R0") c(0,1) else c(0,max(1/R0))
                    , add=FALSE
                    , ...
                      ) {
    if (xvar == "R0") {
        cat("plot_x_in: Rmin = ", Rmin, ", Rmax = ", Rmax, "\n")
        xx <- R0
    } else {
        cat("plot_x_in: plotting as a function of epsilon...\n")
        xx <- epsilon
        if (length(R0) != 1) stop("only one R0 value allowed if xvar is epsilon")
    }

    ## show naive approx (eqm susceptible proportion) as dashed line first:
    plot(xx,
         if (xvar == "R0") 1/R0 else rep(1/R0,length(xx)),
         type="n", log=log, lwd=lwd/2, bty="L", lty="dashed", las=1,
         xlim=xlim, ylim=ylim, xaxs="i", yaxs="i",
         xlab = if (xvar == "R0") expression(R[0]) else expression(epsilon),
         ylab = expression(x[i][n]), ...)
    title(main = main)
    lines(xx,
         if (xvar == "R0") 1/R0 else rep(1/R0,length(xx)),
         lwd=lwd/2, lty="dashed", ...)

    if (xvar == "R0") {
        ## curves with peak prevalence estimate from approximation of integral:
        for (iepsilon in 1:length(epsilon)) {
            ## estimate obtained by taking x_in to be x_f from final size formula:
            lines(R0, x_in_crude(R0,epsilon=epsilon[iepsilon]), lwd=lwd/2, col="pink")
            ## better estimate
            lines(R0, xin_fun(R0,epsilon=epsilon[iepsilon]), lwd=lwd, col=col[iepsilon], ...)
        }
        ## plot with dotted curves when using the crude peak prevalence
        ## ignoring vital dynamics:
        for (iepsilon in 1:length(epsilon)) {
            lines(R0, xin_fun(R0,epsilon=epsilon[iepsilon], peakprev_fun=peak_prev_nvd),
                  lwd=lwd,
                  col=col[iepsilon],
                  lty="dotted")
        }
        legend("topright", bty="n", title=expression(epsilon),
               legend = epsilon, col = col, lwd=lwd)
    } else {
        ## estimate obtained by taking x_in to be x_f from final size formula:
        lines(xx, x_in_crude(R0,xx), lwd=lwd/2, col="pink")
        ## better estimate
        lines(xx, xin_fun(R0,xx), lwd=lwd, col=col)
        lines(xx, xin_fun(R0,xx, peakprev_fun=peak_prev_nvd),
              lwd=lwd, col=col, lty="dotted")
        title(sub = latex2exp::TeX(sprintf("$R_0 = %g$", R0)))
    }
}

##' Plot all approximations to \eqn{x_{\rm in}}
##'
##' @inheritParams plot_x_in
##' @inheritParams compare_funs
##' @param col,lwd,xlim,ylim,log,... see \code{\link{graphical parameters}}
##'
##' @seealso \code{\link{plot_x_in}}, \code{\link{x_in}},
##'     \code{\link{x_in_exact}}
##'
##' @importFrom graphics title
##'
##' @export
##'
plot_all_x_in <- function(Rmin=1.0001, Rmax=20,
                      ##R0 = exp(seq(log(Rmin),log(Rmax),length=1001)),
                      R0 = seq(Rmin,Rmax,length=1001),
                      epsilon = c(0.01, 0.02, 0.03),
                      xvar = "R0"
                    , col = 1:length(epsilon)
                      ##col = c("darkred", "darkgreen", "darkblue")
                    , lwd=2
                    , log="x"
                    , xlim = if (xvar == "R0") c(1,Rmax) else c(0,1)
                    , ylim = if (xvar == "R0") c(0,1) else c(0,max(1/R0))
                    , ...
                      ) {
    omf <- par(mfrow = c(2,2))
    on.exit(par(mfrow = omf))

    ##plot_x_in(xin_fun = x_in_exact) # too slow
    plot_x_in(xin_fun = x_in_crude, ...)
    legend("top", bty="n", legend="x_in_crude")
    plot_x_in(xin_fun = x_in, ...)
    legend("top", bty="n", legend="x_in")
    plot_x_in(xin_fun = x_in_cb, ...)
    legend("top", bty="n", legend="x_in_cb")
    plot_x_in(xin_fun = x_in_hocb, ...)
    legend("top", bty="n", legend="x_in_hocb")
}

