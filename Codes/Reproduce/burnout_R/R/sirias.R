
##' Vector field for SIRiaS model
##'
##' @param tau time in natural units
##' @param vars state variable list (\eqn{S}, \eqn{\log{I}}, \eqn{R})
##' @param parms parameter list
##' @export
sirias_vf <- function(tau, vars, parms) {
    ## vars:
    S <- vars[1]
    I <- exp(vars[2])
    R <- vars[3]
    ## parms:
    beta <- parms[["beta"]]
    gamma <- parms[["gamma"]]
    mu <- parms[["mu"]]
    delta.i <- parms[["delta.i"]]
    delta.a <- parms[["delta.a"]]
    ## gradient:
    dS <- mu * (1-S) - beta*S*I + (delta.i + delta.a*R)*R
    dlogI <- beta*S - (gamma + mu)
    dR <- gamma*I - mu*R - (delta.i + delta.a*R)*R
    return(list(c(dS=dS, dlogI=dlogI, dR=dR)))
}

##' SIRiaS trajectory computation via deSolve
##'
##' @export
##' @inheritParams x_in_exact
##' 
##' @importFrom deSolve lsodar
##' @importFrom utils tail
sirias_traj <- function(R0, epsilon,
                        eta.i=0, eta.a=0,
                        I0=1e-6,
                        t_int = NULL,
                        tmin = 20,
                        maxit = 10,
                        nt = 1000) {

    ## reasonable guess at an approximate max time for a specified R0
    ## (given default I0)
    if (is.null(t_int)) t_int <- max(100/(R0-1), tmin)

    ## precompute useful factors, as well as hard-coding vars[i]
    ## lookup below rather than using with(as.list(...)), for speed
    herdImm <- if (R0 > 1) (1 - 1/R0) else 0

    ## epsilon == mu/(gamma + mu)
    ## -> epsilon*(gamma+mu) == mu
    ## -> gamma*epsilon = mu*(1-epsilon)
    ## -> mu = gamma*epsilon/(1-epsilon)

    ## Ihat <- epsilon * herdImm ## WRONG IF IMMUNITY DECAYS
    Ihat <- eqm_prev(R0, epsilon, eta.i, eta.a)

    gamma <- 1
    mu <- gamma*epsilon/(1-epsilon)
    beta <- R0*(gamma+mu)
    delta.i <- eta.i*(gamma+mu)
    delta.a <- eta.a*(gamma+mu)

    parms <- list(R0 = R0, epsilon = epsilon,
                 beta = beta, gamma = gamma, mu = mu,
                 delta.i = delta.i,
                 delta.a = delta.a)

    print(as.data.frame(parms))

    dd <- deSolve::lsodar(
                       y = c(S=1-I0, logI = log(I0), R=0),
                       ## ?? what should length.out be?
                       times = seq(0, t_int, length.out = nt),
                       func = sirias_vf,
                       parms = parms
                          )

    class(dd) <- c("sirias_traj", class(dd))

    attr(dd, "Ihat") <- Ihat
    attr(dd, "parms") <- parms
    attr(dd, "args") <- list(
        R0 = R0, epsilon = epsilon,
        eta.i = eta.i, eta.a = eta.a,
        I0 = I0,
        t_int = t_int,
        tmin = tmin,
        maxit = maxit,
        nt = nt
    )

    return(dd)
}

#' @export
plot.sirias_traj <- function(dd, ...
                           , las = 1
                           , type = "l"
                           , lwd = 3
                           , log = "y"
                           , col.xin.exact = "red"
                           , col.ystar = "red"
                           , col.xin = "blue"
                           , x.legend = "bottomleft"
                             ) {

    S <- dd[,"S"]
    logI <- dd[,"logI"]
    I <- exp(logI)
    Ihat <- attr(dd, "Ihat")

    plot(S, I,
         log = log, las = las, type = type, lwd = lwd, ...)

    abline(h = Ihat, col = "red")
    
    with(attr(dd,"args"), {
        xin.exact <- 
            x_in_exact(R0=R0, epsilon=epsilon, eta.i=eta.i, eta.a=eta.a,
                       I0=I0, t_int=t_int, tmin=tmin, maxit=maxit, nt=nt)
        abline(v = xin.exact, col = col.xin.exact, lwd = lwd/2)
        xin <- x_in(R0=R0, epsilon=epsilon, eta.i=eta.i, eta.a=eta.a)
        abline(v = xin, col = col.xin, lwd = lwd*1.5, lty = "dotted")
        parms.string <- paste0("$R_0 = ", R0, ", \\epsilon = ", epsilon,
                           ", \\eta_i = ", eta.i, ", \\eta_a = ", eta.a,
                           ", y_* = ", Ihat, "$")
        title(main = latex2exp::TeX(parms.string))
        relerr <- (xin - xin.exact)/xin.exact
        legend.title <- "$\\frac{\\Delta x_{{i}{n}}}{x_{{i}{n}}}$"
        legend(x.legend, bty = "n", title = latex2exp::TeX(legend.title),
               legend = signif(relerr,5))
    })

}
