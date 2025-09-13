##' burnout: A package for computing infectious disease persistence probabilities
##'
##' The burnout package uses analytical approximations to compute the
##' probabilities of fizzle and burnout in the stochastic SIR epidemic
##' model.
##' 
##' @section Key functions:
##'
##' The primary functions are \code{\link{P1_prob}},
##'     \code{\link{fizzle_prob}}, \code{\link{fizzle_time}}, and
##'     \code{\link{t_herd}}.  There are also functions to compare
##'     approximations of \eqn{{\cal P}_1} and component
##'     calculations.
##'
##' Most of these functions depend on the (real branches of the)
##' Lambert \eqn{W} function (\eqn{W_0} or \eqn{W_{-1}}).  Lambert
##' \eqn{W} is implemented in several R packages.  By default we use
##' the implementation in \pkg{\link{gsl}}.  You can change to a
##' different implementation via \code{options(burnout.lambertW =
##' "emdbook")}.
##' 
##' @name burnout
"_PACKAGE"

