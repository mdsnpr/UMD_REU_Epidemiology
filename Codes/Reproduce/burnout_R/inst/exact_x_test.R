library(burnout)
library(deSolve)
library(microbenchmark)

do_benchmark <- FALSE  ## a little slow (~1 minute)

## utility plotting function
pfun <- function(d, parms) {
    par(las = 1, bty = "l")
    plot(logI ~ time, data = as.data.frame(d), type ="l", log = "x")
    Ihat <- with(as.list(parms), epsilon*(1-1/R0))
    abline(h = log(Ihat), col = 2)
}

## checking required tmax scaling as R0 approaches 1
## (from sirr, hardcoded there)
R0vec <- c(1.5, 1.1, 1.05, 1.01, 1.005)
tmax <- c(200, 300, 1000, 2000, 20000)

## ... looks like tmax = 100/(R0-1) is a good approximation
plot(R0vec-1, tmax, log = "xy")
curve(100/x, add= TRUE)


## testing, basic picture
Sval <- x_in_exact_scalar(R0=5 , epsilon = 0.1)
dd <- x_in_exact_scalar(R0=5 , epsilon = 0.1,
                         nt = 1000,
                         return_traj = TRUE)
pfun(dd, parms = list(R0=5, epsilon = 0.1))

## compare (time and values) vs sirr

if (!require("sirr")) stop()

data("HtypeTable", package = "sirr")
data("HtypeReferences", package = "sirr")
sirr_x_in_exact_scalar <- function(R0, epsilon, I0=1e-6, ...) {
    m <- sirr::create_SIRmodel(R0=R0, eps=epsilon)
    ## FIX: tmax=100 works for a large range of R0, but is still too
    ## short as R0 --> 1.  I need to have a stopcrit that ends when
    ## the 2nd crossing occurs.
    tmax <- 100
    if (R0 < 1.5) tmax <- 200
    if (R0 < 1.1) tmax <- 300
    if (R0 < 1.05) tmax <- 1000
    if (R0 < 1.01) tmax <- 2000
    if (R0 < 1.005) tmax <- 20000
    mts <- sirr::compute_SIRts(m, saveRoots=TRUE, inits=c(S=1-I0), tmax=tmax)
    dd <- attr(mts, "rootPoints")
    Ihat.pts <- dd[dd$type == "Ieqm", ]
    xin <- Ihat.pts[2,"S"] # 2nd crossing
    ##Yb <- Ihat.pts[1,"I"]
    ##out <- c(Xb=Xb, Yb=Yb)
    return(xin)
}

## this takes a minute or so ...
if (do_benchmark) {
    mm <- microbenchmark(
        sirr = burnout:::x_in_exact_scalar(R0 = 5, epsilon = 0.1),
        new = burnout:::x_in_exact_scalar2(R0=5 , epsilon = 0.1))
    ## approx 200 ms (sir) vs 2 ms (new)
    print(mm)
}

b1 <- sirr_x_in_exact_scalar(R0 = 5, epsilon = 0.1)
b2 <- burnout:::x_in_exact_scalar(R0 = 5, epsilon = 0.1)
stopifnot(all.equal(b1, b2, tolerance = 5e-6))

## recapitulate SIR integration, draw some pictures
epsilon <- 0.1
R0 <- 5

herdImm <- if (R0 > 1) (1 - 1/R0) else 0
Ihat <- epsilon * herdImm

I0 <- 1e-6
Ihat <- epsilon * herdImm
gamma <- 1
mu <- gamma*epsilon/(1-epsilon)
beta <- R0*(gamma+mu)


sirgrad <- function(tau, vars, parms) {
    S <- vars[1]
    I <- exp(vars[2])
    dS <- mu * (1-S) - beta*S*I
    dlogI <- beta*S - (gamma + mu)
    return(list(c(dS=dS, dlogI=dlogI)))
}

dd2 <- as.data.frame(ode(y = c(S = 1 - I0, logI = log(I0)), 
                         times = seq(0, 20, length.out = 501), func = sirgrad))

par(mfrow=c(1,2), bty = "l", las = 1)
plot(logI ~ time, data = dd2, type = "l",
     ylim = c(-3.2,-0.5),
     xlim = c(0,10))
abline(v = tail(dd[,"time"], 1), lty = 2)
abline(h = log(Ihat), lty = 2)
logIhat_c <- tail(dd[, "logI"], 1)
abline(h = logIhat_c, lty = 3, col = 2)
stopifnot(all.equal(log(Ihat), logIhat_c))
plot(S ~ time, data = dd2, type = "l",
     xlim = c(3, 8),
     ylim = c(0.1, 0.3))
abline(v = tail(dd[,"time"], 1), lty = 2)
abline(h = b1, col = 4, lty = 2)
abline(h = b2, col = 5, lty = 2)
legend("topleft",
       col = 4:5, lty = 2,
       legend = c("sirr", "new"))
