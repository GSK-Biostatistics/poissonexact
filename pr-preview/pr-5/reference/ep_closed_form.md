# Exact conditional Poisson regression (fast, closed-form algorithm)

Computes the exact conditional maximum likelihood estimate, exact
confidence interval, and exact p-value for the coefficient of a single
two-level (treatment vs. control) factor in a Poisson regression with
log offset and, optionally, one or more stratification variables acting
purely as nuisance intercepts. This is the special case handled by
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
when the model contains exactly one exact effect with two levels and no
other categorical covariates; in that special case, the exact
distribution can be obtained in closed form, stratum by stratum,
following Guess, Lydick, Small and Miller (1987) and Guess and Thomas
(1990), instead of by the general per-observation network algorithm. See
the accompanying theory note for why this is faster in this scenario,
and use
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
whenever the model contains additional categorical covariates to adjust
for, or more than one exact effect.

## Usage

``` r
ep_closed_form(
  formula,
  data,
  offset = NULL,
  offset_log = TRUE,
  strata = NULL,
  ref_levels = list(),
  alpha = 0.05,
  ci_type = c("two-sided", "one-sided"),
  verbose = TRUE
)
```

## Arguments

- formula:

  A model formula of the form `response ~ group`, where `group` is the
  two-level factor of interest and is the formula's only right-hand-side
  term. Strata and the offset are supplied separately (as in
  [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)),
  not as terms in this formula.

- data:

  A data frame. One row per observation; rows belonging to the same
  stratum and the same level of `group` are summed internally, so `data`
  may be individual-level or already aggregated.

- offset:

  Name of the person-time / exposure variable in `data`, or `NULL`
  (default). If `NULL`, every observation is treated as contributing one
  unit of exposure (offset = 0 on the log scale), matching
  [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md).
  Note that, unlike the general algorithms, supplying `offset` here
  never imposes a per-observation cap on any observation's count (see
  Details).

- offset_log:

  If `TRUE` (default, matching
  [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)),
  `offset` already holds log(exposure). If `FALSE`, `offset` holds raw
  exposure and is logged internally.

- strata:

  Optional character vector naming one or more stratification variables
  in `data`. If `NULL` (default), all observations form a single
  stratum, i.e. an unstratified two-sample exact Poisson test.

- ref_levels:

  Optional named list giving the reference (control) level of the
  formula's group variable, e.g. `list(trt = "0")`; the other level is
  treated as the treatment group. Defaults to R's usual first level.

- alpha:

  Type I error rate for the confidence interval (default 0.05).

- ci_type:

  `"two-sided"` or `"one-sided"`.

- verbose:

  If `TRUE` (default), prints a short note identifying this as the
  closed-form two-group/strata-only algorithm.

## Value

An object of class `c("ep_closed_form", "ep_fit")`: a list with elements
`exactparmest` (a one-row data frame with `Parameter`, `Estimate` (the
coefficient beta = log rate ratio, not the rate ratio), `StdError`,
`Alpha`, `LowerCL`, `UpperCL` (beta scale), `PValue`, and `SuffStat`),
`OUTDIST` (the exact distribution of the sufficient statistic), and
`call`.

## Details

Unlike
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
and
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md),
this function has no `exactmax` argument and never imposes a
per-observation cap beyond what conditioning on the stratum total
already implies (SAS's own "Rule 2"; see the theory vignette, Sections
5.2 and 10). In particular, supplying `offset` here does *not* trigger
the offset-implied per-observation cap ("Rule 3": \\S_i = \max(\lceil
e^{o_i} \rceil, 1)\\) that the general algorithms apply automatically
whenever `offset` is supplied without `exactmax`. When that cap would
bind, small exposures relative to the observed counts, typically in
unstratified or lightly-stratified two-group data, `ep_closed_form()`
and the general algorithms are, by design, answering subtly different
questions, and can legitimately disagree on the estimate, confidence
interval, and p-value even given identical call arguments. This is
expected, documented behaviour, not a bug in either algorithm; see the
theory vignette for a worked numerical example and
`tests/testthat/test-cross-algorithm-agreement.R` for a regression test
that locks in the expected divergence.

## References

Guess HA, Lydick EG, Small RD, Miller LP, Exact binomial confidence
intervals for the relative risk in follow-up studies with sparsely
stratified incidence density data, American Journal of Epidemiology,
Volume 125, Issue 2, February 1987, Pages 340-347,
https://doi.org/10.1093/oxfordjournals.aje.a114535

Guess HA, Thomas JE. A rapidly converging algorithm for exact binomial
confidence intervals about the relative risk in follow-up studies with
stratified incidence-density data. Epidemiology. 1990 Jan;1(1):75-7.
doi: 10.1097/00001648-199001000-00016. PMID: 2081244.

## Examples

``` r
df <- data.frame(
  cases = c(3, 32, 1, 19, 2, 24, 4, 40),
  fu = c(100, 100, 120, 120, 150, 150, 170, 170),
  trt = c(1, 0, 1, 0, 1, 0, 1, 0),
  strata = c(1, 1, 2, 2, 3, 3, 4, 4)
)
ep_closed_form(
  formula = cases ~ trt,
  data = df,
  offset = "fu",
  offset_log = FALSE,
  strata = "strata"
)
#> 
#> Closed-form two-group/strata-only exact Poisson algorithm
#>                 (Guess & Thomas, 1990): 4 stratum(strata), 10 treatment cases,
#>                 115 control cases.
#> 
#> Exact conditional Poisson, two groups + strata (closed-form Guess & Thomas algorithm)
#> -------------------------------------------------------------------------------------
#>  Parameter  Estimate  StdError Alpha   LowerCL   UpperCL      PValue SuffStat
#>       trt1 -2.442347 0.3296902  0.05 -3.203615 -1.797111 9.11799e-24       10
```
