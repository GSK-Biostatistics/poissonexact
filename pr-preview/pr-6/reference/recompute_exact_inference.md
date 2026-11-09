# Recompute a p-value and/or confidence interval from an already-fitted exact distribution

Answers a *different* question about a parameter already fit by
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md),
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md),
or
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md);
a different null value `beta0`, a different significance level `alpha`,
or a different `ci_type`/`direction`, directly from that fit's `OUTDIST`
table, without rebuilding the joint distribution. This is possible, and
exact (not an approximation), because every algorithm's reference
distribution is the same one-dimensional exponential-tilt family,
\\\Pr\_\beta(T = t) \propto C(t)\\ e^{\beta t}\\, which sums to one at
every value of \\\beta\\: re-tilting an already-built `C(t)` table at a
new `beta0`, or re-deriving a confidence limit at a new `alpha`, only
ever needs that one table, not a new enumeration over the model's
stratification and nuisance structure. Building the joint table in the
first place is the expensive part of every algorithm in this package;
this function is the cheap part, made reusable.

## Usage

``` r
recompute_exact_inference(
  fit,
  parm,
  beta0 = 0,
  alpha = 0.05,
  ci_type = c("two-sided", "one-sided"),
  direction = c("upper", "lower")
)
```

## Arguments

- fit:

  A fitted object returned by
  [`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md),
  [`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md),
  or
  [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
  (any object of class `"ep_fit"`), i.e. anything with an `OUTDIST` list
  and an `exactparmest` data frame in the format those three functions
  return.

- parm:

  The name of one parameter already present in
  `fit$OUTDIST`/`fit$exactparmest$Parameter` (e.g. `"armvaccinated"`);
  matches the column name `exactparmest` reports for that term.

- beta0:

  The null value of `beta` (the log rate ratio) to test against; `0`
  (the default) reproduces the same null the original fit itself tested.
  Supply a non-zero value to test a different margin, e.g. `log(0.5)`
  for a 50% rate-ratio margin.

- alpha:

  Significance level for the confidence interval; need not match the
  `alpha` the original fit used.

- ci_type:

  `"two-sided"` (default) or `"one-sided"`; need not match the original
  fit's own `ci_type`.

- direction:

  For `ci_type = "one-sided"`: `"upper"` reports only the upper
  confidence limit (lower limit fixed at `-Inf`); `"lower"` reports only
  the lower limit (upper limit fixed at `Inf`). Ignored for
  `ci_type = "two-sided"`.

## Value

A list with:

- `beta0`: the null value tested, echoed back.

- `p_lower`: the lower tail probability at `beta0`, \\\Pr\_{\beta_0}(T
  \le t\_{obs})\\.

- `p_upper`: the upper tail probability at `beta0`, \\\Pr\_{\beta_0}(T
  \ge t\_{obs})\\.

- `pvalue_two_sided`: twice the smaller one-sided tail probability,
  capped at 1 (the package's two-sided exact p-value convention,
  matching SAS's `PROC GENMOD` `EXACT` statement).

- `pvalue_one_sided`: `min(p_lower, p_upper)`.

- `alpha`: the significance level used, echoed back.

- `LowerCL`, `UpperCL`: the recomputed confidence limits, on the same
  scale (`beta`, the log rate ratio) as `fit$exactparmest$Estimate`.

## See also

[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md),
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md),
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md),
and
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)'s
discussion of reusing an already-built exact distribution.

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

fit <- ep_multi_shift(
  formula = cases ~ group,
  data = df,
  offset = "lnfutime",
  offset_log = TRUE,
  strata = "age",
  exact_vars = "group",
  ref_levels = list(group = "0", age = "50"),
  verbose = FALSE
)

# Same question the fit itself already answered (beta0 = 0, alpha = 0.05):
# pvalue_two_sided/LowerCL/UpperCL should match fit$exactparmest exactly.
recompute_exact_inference(fit, "group1")
#> $beta0
#> [1] 0
#> 
#> $p_lower
#> [1] 0.002381662
#> 
#> $p_upper
#> [1] 0.9993993
#> 
#> $pvalue_two_sided
#> [1] 0.004763324
#> 
#> $pvalue_one_sided
#> [1] 0.002381662
#> 
#> $alpha
#> [1] 0.05
#> 
#> $LowerCL
#> [1] -2.372642
#> 
#> $UpperCL
#> [1] -0.335739
#> 

# A different confidence level, reusing the same OUTDIST table:
recompute_exact_inference(fit, "group1", alpha = 0.10)
#> $beta0
#> [1] 0
#> 
#> $p_lower
#> [1] 0.002381662
#> 
#> $p_upper
#> [1] 0.9993993
#> 
#> $pvalue_two_sided
#> [1] 0.004763324
#> 
#> $pvalue_one_sided
#> [1] 0.002381662
#> 
#> $alpha
#> [1] 0.1
#> 
#> $LowerCL
#> [1] -2.193102
#> 
#> $UpperCL
#> [1] -0.4607849
#> 
```
