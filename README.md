
# poissonexact

Exact conditional inference for Poisson regression: exact conditional
maximum likelihood estimates, exact p-values, and exact confidence
intervals for categorical effects in Poisson regression models with an
optional offset and stratification. Nuisance parameters (stratum
intercepts and other model effects not under test) are eliminated by
conditioning on their sufficient statistics rather than relying on
asymptotic approximations. This matters when case counts are small or
sparse and standard Wald/likelihood-ratio inference is unreliable (rare
outcomes, small strata, vaccine/drug-safety signal detection, and
similar settings).

Three algorithms are provided for different problem shapes:

- `ep_closed_form()`: a fast closed-form method for a single two-level
  effect with strata as the only nuisance structure, following Guess,
  Lydick, Small and Miller (1987) and Guess and Thomas (1990).
- `ep_multi_shift()`: a multivariate shift algorithm for arbitrary
  categorical covariates and strata, generalizing Hirji, Mehta and
  Patel (1987) to bounded Poisson counts.
- `ep_network()`: a network algorithm offering more aggressive
  mid-computation pruning for larger problems, generalizing Mehta, Patel
  and Senchaudhuri (1992) to bounded Poisson counts.

`recompute_exact_inference()` lets a caller ask a *different* question
(a different null, significance level, or CI type/direction) about an
already-fitted result’s `OUTDIST` table without rebuilding it.

The conditioning approach mirrors that used by SAS’s `PROC GENMOD`
`EXACT` statement; users familiar with that procedure’s options and
output will find the terminology (`EXACTMAX`, `OFFSET`, `OUTDIST`)
intentionally similar.

## Installation

This package is internal to GSK and is not on CRAN. Install from the
internal GitHub mirror:

``` r
# install.packages("pak")
pak::pkg_install("gsk-tech/poissonexact")
```

## A minimal example

A two-arm, single-stratum comparison (the simplest case
`ep_closed_form()` is built for):

``` r
library(poissonexact)

df <- data.frame(
  group = c("treatment", "control"),
  cases = c(2, 8),
  person_time = c(100, 100)
)

fit <- ep_closed_form(
  formula = cases ~ group,
  data = df,
  offset = "person_time",
  offset_log = FALSE,
  ref_levels = list(group = "control"),
  verbose = FALSE
)
fit$exactparmest
```

See `vignette("Theory-Implementation")` for the full statistical
derivation and `vignette("Example")` for a realistic, stratified worked
example (including a Monte Carlo check of empirical coverage against the
usual asymptotic Wald approach).

## Status

This package is under active internal development (`Version: 0.1.0`).
See `NEWS.md` for what’s currently implemented and
`CHANGES-updated-package.md` / `REMAINING-POINTS.md` (in this staging
folder) for the current review status against the tech-team’s code
review.
