# Exact conditional Poisson regression (network algorithm)

A structurally different, more aggressively pruned alternative to
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md),
generalizing the network algorithm of Mehta, Patel, and Senchaudhuri
(1992), SAS's own `EXACTOPTIONS METHOD=NETWORK`, to Poisson counts. This
algorithm avoids constructing separate per-stratum tables and
cross-joining them. Instead, it processes all observations across every
stratum in a single, continuous forward pass. To optimize performance,
it instantly discards any partial state the moment it becomes
mathematically impossible to hit the observed nuisance total (either
because it cannot reach the target or will inevitably overshoot it).
These mid-stratum pruning boundaries are calculated just once per
stratum using a closed-form backward induction, bypassing the need to
wait for stratum boundaries or explicit cross-joins. Every argument, cap
rule, interaction-naming convention, and output format is identical to
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md),
this function reuses the same data-preparation, cap-rule,
stratum-ordering, `compute_stratum_max_add()`, CI, p-value, and CMLE
helpers; only the table-construction strategy differs. See
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)'s
documentation for the full description of every argument. Unlike
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md),
this function has no `chunk_size` argument: it never materializes a
cross-join between two strata's tables for `chunk_size` to bound.

## Usage

``` r
ep_network(
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

- condition_on:

  `"model"` (default; condition on every other model effect, matching
  SAS) or `"exact"` (condition only on the other requested
  `exact_vars`).

- verbose:

  If `TRUE` (default), prints progress and which SAS rule determined the
  per-observation caps.

## Value

An object of class `c("ep_network", "ep_multi_shift", "ep_fit")`, in the
same format as
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md).

## References

Mehta CR, Patel NR, Senchaudhuri P. Exact Stratified Linear Rank Tests
for Ordered Categorical and Binary Data. Journal of Computational and
Graphical Statistics, 1(1):21-40, 1992.

## Examples

``` r
df <- data.frame(
  group = c( 1L,  1L,  1L,  0L,  0L,  0L),
  age      = c(50L, 60L, 70L, 50L, 60L, 70L),
  subjects = c(3492, 2141, 1711, 3525, 2166, 1724),
  cases    = c(   3,    2,    1,   10,   7,   5),
  futime   = c(50, 60, 40, 51, 62, 42)
)

df$lnfutime <- log(df$futime)

ep_network(
  formula    = cases ~ group,
  data       = df,
  offset     = "lnfutime",
  offset_log = TRUE,
  strata     = "age",
  exact_vars = "group",
  ref_levels = list(group="0", age="50"),
  alpha      = 0.05,
  verbose    = TRUE,
  ci_type = "two-sided"
)
#> 
#> SAS rule (3): OFFSET specified and EXACTMAX not specified =>
#>                 S_i = max(ceil(exp(o_i)),1).
#> 
#> === EXACT for parameter column: group1 (network algorithm) ===
#>   stratum 70: n=2, M=6
#>   stratum 60: n=2, M=9
#>   stratum 50: n=2, M=13
#>   after finishing stratum 1/3 (70): live states = 7
#>   after finishing stratum 2/3 (60): live states = 16
#>   after finishing stratum 3/3 (50): live states = 29
#> 
#> Exact conditional Poisson (categorical + STRATA; network algorithm)
#> -------------------------------------------------------------------
#>  Parameter Estimate  StdError Alpha   LowerCL   UpperCL      PValue SuffStat
#>     group1 -1.26913 0.4605711  0.05 -2.372642 -0.335739 0.004763324        6
```
