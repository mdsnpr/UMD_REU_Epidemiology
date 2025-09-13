
##' Susceptibles at boundary layer (exact) (scalar version)
##'
##' @export
##' @inheritParams x_in
##' @param I0 initial prevalence (proportion)
##' @param nt number of time steps
##' @param tmin minimum integration time
##' @param maxit max number of tries to double integration time
##' @param t_int integration time (default is max(tmin, 100/(R0-1)))
##' @param return_traj (logical) return trajectory?
##' 
##' @importFrom deSolve lsodar
##' @importFrom utils tail
x_in_exact_scalar <- function(R0, epsilon,
                              eta.i=0, eta.a=0,
                              I0=1e-6,
                              t_int = NULL,
                              tmin = 20,
                              maxit = 10,
                              nt = 100,
                              return_traj = FALSE) {

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

    ## Ihat <- epsilon * herdImm ## WRONG IF IMMUNITY DECAYS!!
    Ihat <- eqm_prev(R0, epsilon, eta.i, eta.a)

    gamma <- 1
    mu <- gamma*epsilon/(1-epsilon)
    beta <- R0*(gamma+mu)
    delta.i <- eta.i*(gamma+mu)
    delta.a <- eta.a*(gamma+mu)

    sirgrad <- function(tau, vars, parms) {
        S <- vars[1]
        I <- exp(vars[2])
        R <- vars[3]
        dS <- mu * (1-S) - beta*S*I + (delta.i + delta.a*R)*R
        dlogI <- beta*S - (gamma + mu)
        dR <- gamma*I - mu*R - (delta.i + delta.a*R)*R
        return(list(c(dS=dS, dlogI=dlogI, dR=dR)))
    }

    eventfun <- function(t, x, p) { ..nroot <<- ..nroot + 1; return(x) }
    ## could also condition eventfun on (I = I^* AND deriv is negative)
    rfunc <- function(tau, vars, parms) {
        ## criterion for I == Ihat
        Ieqm <- exp(vars[2]) - Ihat
        ## eventfun is always called once initially (why?), so we want
        ## to make sure this criterion is non-zero until after the first root of Ieqm
        Ieqm2 <- if (..nroot <= 1) -1 else Ieqm
        return(c(Ieqm, Ieqm2))
    }

    ok <- FALSE

    it <- 0
    
    while (!ok && it < maxit) {
        ..nroot <- 0
        dd <- deSolve::lsodar(y = c(S=1-I0, logI = log(I0), R=0),
                              ## ?? what should length.out be?
                              times = seq(0, t_int, length.out = nt),
                              func = sirgrad,
                              parms = list(R0 = R0, epsilon = epsilon,
                                           delta.i = delta.i,
                                           delta.a = delta.a),
                              rootfunc = rfunc,
                              events = list(func = eventfun, root = TRUE,
                                            ## terminate on criterion 2 (-> I==Ihat for second time)
                                            terminalroot = 2))
        ok <- max(dd[, "time"] < t_int)
        it <- it + 1
    }
    
    if (return_traj) return(dd)
    if (max(dd[, "time"]) == t_int && it == maxit) stop(sprintf("reached maximum time (%1.1f): increase maxit?",
                                                               t_int))
    tail(dd[,"S"], 1)
}

##'
##' Susceptibles at boundary layer (exact) (vector version)
##'
##' @rdname x_in_exact_scalar
##' @export
##' 
x_in_exact <- Vectorize(x_in_exact_scalar)


