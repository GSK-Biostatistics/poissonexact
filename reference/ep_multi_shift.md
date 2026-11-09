# Exact conditional Poisson regression (multivariate-shift algorithm)

Exact conditional inference for one or more categorical effects in a
Poisson regression with a log link, an optional offset, and optional
stratification. Nuisance parameters (stratum intercepts and every model
effect not listed in `exact_vars`) are eliminated by conditioning on
their sufficient statistics, and the resulting exact conditional
distribution of each requested effect's sufficient statistic is obtained
by the multivariate-shift algorithm: a per-stratum dynamic program
combined with explicit cross-stratum convolution, generalizing the
multivariate shift algorithm of Hirji, Mehta, and Patel (1987) for exact
logistic regression to bounded Poisson counts, following the
`EXACTMAX`/`OFFSET` rules documented for SAS `PROC GENMOD`'s `EXACT`
statement.

## Usage

``` r
ep_multi_shift(
  formula,
  data,
  exact_vars,
  offset = NULL,
  offset_log = TRUE,
  exactmax = NULL,
  strata = NULL,
  ref_levels = list(),
  alpha = 0.05,
  ci_type = c("two-sided", "one-sided"),
  chunk_size = 4000L,
  condition_on = c("model", "exact"),
  verbose = TRUE
)
```

## Arguments

- formula:

  A model formula.

- data:

  A data frame.

- exact_vars:

  Character vector of model *terms* to compute exact results for. Each
  entry must name a term that appears in `formula`, exactly as SAS's
  `EXACT` statement requires. An entry can be a main effect (`"region"`)
  or an interaction of any order, named by joining every variable in it
  with `*` or `:` (`"region:gender"` and `"region*gender"` both work,
  and the variables may be listed in either order). Interaction columns
  are matched via the model's own
  [`terms()`](https://rdrr.io/r/stats/terms.html)/`assign` structure,
  not by name-guessing, so a term is never confused with a different
  term that happens to share a text prefix (e.g. requesting `"region"`
  in a model that also contains `region:gender` selects only region's
  own main-effect columns). Requesting a term does *not* automatically
  request any other term built from the same variable(s): a main effect
  and its interactions are entirely independent entries, so testing
  `region` while leaving `region:gender` untested keeps the interaction
  as an ordinary nuisance covariate (conditioned on under
  `condition_on = "model"`, never tested); to test the interaction
  itself, add it as its own entry, e.g.
  `exact_vars = c("region", "region:gender")`.

- offset:

  Name of an offset variable in `data`, or `NULL` (default). If `NULL`,
  every observation is treated as contributing one unit of exposure
  (offset = 0 on the log scale).

- offset_log:

  If `TRUE` (default), `offset` already holds log(exposure); if `FALSE`,
  it holds raw exposure.

- exactmax:

  `NULL` (not specified), `TRUE` (specified without a variable), or a
  character variable name (SAS `EXACTMAX=`).

- strata:

  Optional character vector of stratification variable names in `data`.

- ref_levels:

  Optional named list of reference levels for categorical variables.

- alpha:

  Significance level for confidence intervals.

- ci_type:

  `"two-sided"` (default) or `"one-sided"`, matching
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)'s
  argument of the same name. For `"one-sided"`, only the confidence
  limit on the informative side is computed for each parameter (the
  other is left at its boundary, `-Inf` or `Inf`); which side is
  informative is determined per parameter from the two one-sided tail
  probabilities at the null (see `two_sided_from_tails()`'s
  `direction`), exactly as in
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md).

- chunk_size:

  Maximum number of rows processed at once while combining two strata's
  coefficient tables. Combining strata requires pairing every row of one
  stratum's table with every row of the other's; `chunk_size` bounds how
  many such pairs are materialized simultaneously, trading memory for
  speed. The default is generous for typical problems; lower it if
  combining strata runs out of memory, or raise it if memory is
  plentiful and speed matters more.

- condition_on:

  `"model"` (default; condition on every other model effect, matching
  SAS) or `"exact"` (condition only on the other requested
  `exact_vars`).

- verbose:

  If `TRUE` (default), prints progress and which SAS rule determined the
  per-observation caps.

## Value

An object of class `c("ep_multi_shift", "ep_fit")`: a list with elements
`exactparmest` (a one-row data frame with `Parameter`, `Estimate` (the
coefficient beta = log rate ratio, not the rate ratio), `StdError`,
`Alpha`, `LowerCL`, `UpperCL` (beta scale), `PValue`, and `SuffStat`),
`OUTDIST` (the exact distribution of the sufficient statistic), and
`call`.

## Details

`formula` may contain interaction terms (`a*b`, `a:b`, or higher-order
interactions) freely, in two distinct roles. An interaction can be a
plain adjustment covariate, left out of `exact_vars` entirely, in which
case it is conditioned on like any other nuisance term (under the
default `condition_on = "model"`) but never itself tested. Or an
interaction can be the effect of interest, including an interaction of
the effect of interest with another covariate, by listing it as its own
entry in `exact_vars`; see that argument's documentation for the exact
naming rules. The two roles are independent: listing a main effect in
`exact_vars` never implicitly tests its interactions, and listing an
interaction never implicitly tests the main effects it is built from,
each term to be tested must be named explicitly.

## References

Hirji KF, Mehta CR, Patel NR. Computing Distributions for Exact Logistic
Regression. JASA, 82(400):1110-1117, 1987.

Mehta CR, Patel NR. Exact logistic regression: theory and examples.
Statistics in Medicine, 14(19):2143-2160, 1995.

## Examples

``` r
df <- data.frame(
  group = c(1L, 1L, 1L, 0L, 0L, 0L),
  age = c(50L, 60L, 70L, 50L, 60L, 70L),
  subjects = c(3492, 2141, 1711, 3525, 2166, 1724),
  cases = c(3, 2, 1, 10, 7, 5),
  futime = c(50, 60, 40, 51, 62, 42)
)

df$lnfutime <- log(df$futime)

ep_multi_shift(
  formula = cases ~ group,
  data = df,
  offset = "lnfutime",
  offset_log = TRUE,
  strata = "age",
  exact_vars = "group",
  ref_levels = list(group="0", age="50"),
  alpha = 0.05,
  verbose = TRUE,
  ci_type = "two-sided"
)
#> 
#> SAS rule (3): OFFSET specified and EXACTMAX not specified =>
#>                 S_i = max(ceil(exp(o_i)),1).
#> 
#> === EXACT for parameter column: group1 ===
#>   stratum 50: n=2, M=13
#>   stratum 60: n=2, M=9
#>   stratum 70: n=2, M=6
#>   after combining 2/3 strata: joint rows = 16
#>   after combining 3/3 strata: joint rows = 29
#> 
#> Exact conditional Poisson (categorical + STRATA; multivariate-shift algorithm)
#> ------------------------------------------------------------------------------
#>  Parameter Estimate  StdError Alpha   LowerCL   UpperCL      PValue SuffStat
#>     group1 -1.26913 0.4605711  0.05 -2.372642 -0.335739 0.004763324        6
```
