
`burnout` is an R package for calculating the probability of epidemic
burnout (or persistence) after a major outbreak.  It was used to create the figures in

   Parsons TL, Bolker BM, Dushoff J, Earn DJD (2024),
   "The probability of epidemic burnout in the stochastic SIR model with vital dynamics",
   PNAS, Volume: 121, Issue: 5, DOI: [10.1073/pnas.2313708120](https://doi.org/10.1073/pnas.2313708120)

The `burnout` package requires R ≥ 4.2.0 in order for the documentation to be displayed as desired.  As mentioned [here](https://cran.r-project.org/doc/manuals/r-devel/NEWS.html), 4.2.0 automatically supports [KaTeX](https://katex.org/docs/support_table.html) and [MathJax](https://www.mathjax.org/), which I use liberally in [roxygen](https://roxygen2.r-lib.org/) documentation.

----

_The remainder of this file currently contains notes for the
authors. We will clean this up when time permits._

- LaTeX issues: `\\mathscr` (`mathrsfs` package?, 
   Debugging: `R CMD Rdconv man/burnout.Rd`
   
```
R CMD Rd2pdf burnout/man/burnout.Rd
```

# generating graphs from Parsons _et al_ (2024)

Generate the analytical curves in Figure 3 via
```
example(plot_P1)
```

Generate Figures B1 and B2 via
```
example(plot_x_in)
```
Note that the green/blue curves (i.e., crude/better approximations of
peak prevalence `ymax`) in Todd's plots correspond to dotted/solid
curves in this plot.

# comparing approximations

The function `compare_funs` creates a data frame with values of two
functions of `R0` and `epsilon` at a grid of points.  It returns an
object of class `compare_funs`, for which there is a `plot` method.
So, for example, the following compares the exact and approximate
expressions for Kendalls's q:
```
cf <- compare_funs()
plot(cf)
```

To compare the crude and better approximations of `xin`, use
```
cxinc <- compare_funs(x_in, x_in_crude, Rmax=20)
plot(cxinc)
```
In the above plot (with `Rmax = 20`) the bottom two panels (for
`epsilon = 0.01` and `0.1`) correspond to Figure B1 in the paper.

To compare the approximate and exact `xin`, use
```
cxine <- compare_funs(x_in, x_in_exact)
plot(cxine)
```
or, to look a little more closely use, for example,
```
cxine <- compare_funs(x_in, x_in_exact, Rmax=8, epsilon=0.01)
plot(cxine)
```

To compare van Herwaarden's (1997) approximation to ours:
```
ch <- compare_funs(P1_prob, P1_prob_vanH)
plot(ch)
```
_Note, however, that the above works only because of a `try` catch inside `vanH_prob`.  We actually need to resolve this._

A better test is to compare the burnout probabilities (conditional on
not fizzling), as opposed to the persistence probabiities, because
`P1` is derived from the burnout probability identically in all cases.
The formulae that are different are the components from burning out
conditional on not fizzling:
```
cb <- compare_funs(burnout_prob, burnout_prob_vanH)
plot(cb)
```

To compare van Meerson and Sasorov's (2009) approximation to ours:
```
cms <- compare_funs(P1_prob, P1_prob_MS)
plot(cp)
```
