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

- [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md):
  a fast closed-form method for a single two-level effect with strata as
  the only nuisance structure, following Guess, Lydick, Small and
  Miller (1987) and Guess and Thomas (1990).
- [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md):
  a multivariate shift algorithm for arbitrary categorical covariates
  and strata, generalizing Hirji, Mehta and Patel (1987) to bounded
  Poisson counts.
- [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md):
  a network algorithm offering more aggressive mid-computation pruning
  for larger problems, generalizing Mehta, Patel and Senchaudhuri (1992)
  to bounded Poisson counts.

[`recompute_exact_inference()`](https://gsk-biostatistics.github.io/poissonexact/reference/recompute_exact_inference.md)
lets a caller ask a *different* question (a different null, significance
level, or CI type/direction) about an already-fitted result’s `OUTDIST`
table without rebuilding it.

The conditioning approach mirrors that used by SAS’s `PROC GENMOD`
`EXACT` statement; users familiar with that procedure’s options and
output will find the terminology (`EXACTMAX`, `OFFSET`, `OUTDIST`)
intentionally similar.

## Installation

Install the released version of {poissonexact} from CRAN:

``` r
install.packages("poissonexact")
```

Or install the development version from GitHub with:

``` r
# install.packages("pak")
pak::pak("GSK-Biostatistics/poissonexact")
```

## A minimal example

A two-arm, single-stratum comparison (the simplest case
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
is built for):

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

See
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
for the full statistical derivation and
[`vignette("Example")`](https://gsk-biostatistics.github.io/poissonexact/articles/Example.md)
for a realistic, stratified worked example (including a Monte Carlo
check of empirical coverage against the usual asymptotic Wald approach).

## Status

This package is under active internal development (`Version: 0.1.0`).
See `NEWS.md` for what’s currently implemented and
`CHANGES-updated-package.md` / `REMAINING-POINTS.md` (in this staging
folder) for the current review status against the tech-team’s code
review.
