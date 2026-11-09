# A Low-Incidence Vaccine Efficacy Example

``` r

library(poissonexact)
```

## Why this vignette exists

The companion vignette,
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md),
derives *why* exact conditional inference works and *how* each of this
package’s three algorithms computes it. This vignette is the “so what”:
a single, realistic scenario where the difference between exact
inference and the usual asymptotic (Wald) approach is not academic.

Vaccine efficacy trials are close to the worst case for large-sample
approximations, for exactly the reasons
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 1 lists: the outcome (a confirmed case) is rare, the trial is
often broken into several strata (age band, region, risk group) each
with only a handful of events, and a genuinely effective vaccine pushes
the observed rate ratio toward the edge of what a handful of events can
distinguish from zero. That is precisely the regime this package
targets.

This vignette does three things:

1.  Simulates a small, stratified, low-incidence vaccine efficacy trial.
2.  Analyses it with
    [`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
    and, for comparison, with
    [`glm()`](https://rdrr.io/r/stats/glm.html)’s ordinary asymptotic
    (Wald) inference, fit to the *same* model, using the equivalence
    between stratifying and adjusting for a fully saturated factor
    established in
    [`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
    Section 4.1-4.2.
3.  Shows how an already-computed exact distribution can be reused to
    answer a *different* question (a different null hypothesis, or a
    different confidence level) without repeating the combinatorial work
    that built it in the first place
    ([`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
    Sections 6.5 and 8.4).

A small Monte Carlo study at the end repeats the simulation many times
to check, empirically, how often each method’s nominal 95% interval
actually contains the true value, rather than relying on a single lucky
(or unlucky) draw to make the point.

## The scenario

Consider a trial with two arms (vaccine, placebo), each followed for one
year, stratified into three age bands with different baseline attack
rates: older participants are assumed to be at higher baseline risk of
the outcome. The vaccine is assumed to reduce the incidence rate by the
same proportion in every stratum: a true vaccine efficacy of 70% (a true
rate ratio of 0.30), exactly the “constant across strata” assumption
that every algorithm in this package makes for the effect under test
([`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 2.4).

This trial is deliberately sized small, a few hundred participants per
stratum rather than a few thousand, because the point of this vignette
is sharpest exactly where events are scarce. As the next section shows,
a trial this size can (and, with the seed fixed below, does) come back
with *zero* observed vaccine-arm cases in every single stratum by chance
alone, even though the vaccine is genuinely 70% effective. That outcome
is not a contrived pathology; it is a realistic, if unlucky, draw from a
small low-incidence trial, and it is exactly the situation where the
choice of analysis method matters most.

``` r

age_grp <- factor(
  c("18-40", "41-64", "65+"),
  levels = c("18-40", "41-64", "65+")
)

true_ve <- 0.70 # true vaccine efficacy, constant across strata
true_rr <- 1 - true_ve # true rate ratio, vaccine vs. placebo
true_beta <- log(true_rr) # true log rate ratio (the target of inference)

trial_design <- data.frame(
  age_grp = rep(age_grp, each = 2),
  arm = factor(
    rep(c("placebo", "vaccine"), times = 3),
    levels = c("placebo", "vaccine")
  ),
  n = c(300, 300, 250, 250, 150, 150),
  person_years = c(300, 300, 250, 250, 150, 150),
  placebo_rate = rep(c(0.0060, 0.0100, 0.0180), each = 2)
)
trial_design$true_rate <- with(
  trial_design,
  ifelse(
    arm == "vaccine",
    placebo_rate * true_rr,
    placebo_rate
  )
)
trial_design$expected_cases <- with(trial_design, person_years * true_rate)
trial_design
#>   age_grp     arm   n person_years placebo_rate true_rate expected_cases
#> 1   18-40 placebo 300          300        0.006    0.0060           1.80
#> 2   18-40 vaccine 300          300        0.006    0.0018           0.54
#> 3   41-64 placebo 250          250        0.010    0.0100           2.50
#> 4   41-64 vaccine 250          250        0.010    0.0030           0.75
#> 5     65+ placebo 150          150        0.018    0.0180           2.70
#> 6     65+ vaccine 150          150        0.018    0.0054           0.81
```

Even before simulating anything, notice how small the *expected* case
counts already are, especially in the vaccine arm: under 1 expected case
in every single stratum (0.54, 0.75, and 0.81), and under 3 even in the
larger placebo arm. This is the sparse-data regime
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 1 motivates specifically so that the vaccine arm’s expected
counts sit close enough to zero that an all-zero outcome is a real,
non-negligible possibility rather than a one-in-a-million curiosity.

## Simulating the trial

Actual case counts are simulated from these expected values with a fixed
seed, so this document reproduces identically every time it is rendered.

``` r

set.seed(2026)
trial_design$cases <- rpois(nrow(trial_design), trial_design$expected_cases)
trial_design[, c(
  "age_grp",
  "arm",
  "n",
  "person_years",
  "expected_cases",
  "cases"
)]
#>   age_grp     arm   n person_years expected_cases cases
#> 1   18-40 placebo 300          300           1.80     2
#> 2   18-40 vaccine 300          300           0.54     0
#> 3   41-64 placebo 250          250           2.50     1
#> 4   41-64 vaccine 250          250           0.75     0
#> 5     65+ placebo 150          150           2.70     3
#> 6     65+ vaccine 150          150           0.81     0
```

With expected counts this low, it would not be at all unusual for one or
more cells above to come out as 0, 1, or 2 cases by chance alone. With
the seed fixed above, that is exactly what happened here: every
vaccine-arm cell came back with **zero** observed cases, while the
placebo arm shows a handful of events in every stratum. A vaccine-arm
total of zero is not a data error and not a sign the simulation is
broken; with expected vaccine-arm counts of 0.54, 0.75, and 0.81
(summing to a total expected count of only 2.1), an all-zero outcome
across all three strata has roughly a 1-in-8 chance of happening on any
given trial (e^{-2.1} \approx 0.12), purely from Poisson variation, even
though the vaccine is genuinely working. This is precisely where
asymptotic Wald inference is at its least trustworthy, and precisely
where exact conditional inference has no such weakness, because it never
relies on a large-sample approximation in the first place. The next two
sections show what each method actually does with this data.

## Recommended: aggregate repeated covariate patterns before fitting

Every algorithm in this package
([`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md),
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md),
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md))
processes its input one row at a time when building the exact
conditional distribution. If your raw data has multiple rows that share
the exact same covariate pattern (every column that enters the model,
plus `strata`, plus (if used) the `EXACTMAX` variable), summing those
rows into one before fitting is a free performance win, **not an
approximation**: it never changes the exact result.

The reason this is exact, not just a convenient shortcut, is the same
sufficiency argument the algorithms themselves are built on
([`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)):
independent Poisson counts that share the same rate parameter sum to a
Poisson with the combined exposure (the standard Poisson convolution
identity). Concretely, if two rows are identical on every column the
model conditions on and differ only in their own `cases`/exposure,
replacing them with one row (`cases` summed, exposure summed, or, if
`offset_log = TRUE`, the log-offset combined as
`log(sum(exp(log_offset)))`) gives the identical sufficient statistics
the dynamic program would have derived from the two original rows. This
matters in practice because the dynamic program’s cost grows with the
number of rows processed, not the number of underlying subjects;
collapsing, say, patient-level records down to one row per distinct
covariate pattern can turn a slow fit into an instant one with no change
to any p-value or confidence interval.

One caveat: this equivalence is only *exactly* guaranteed under two of
the package’s three SAS-style per-observation cap rules:
`EXACTMAX = var` (as long as every merged row shares the same `EXACTMAX`
value) and the no-offset/no-`EXACTMAX` default (whose cap is always the
whole stratum total regardless of how many rows it’s split across).
Under the third rule (`offset` supplied without `EXACTMAX`, where each
row’s own cap is derived from its own exposure), pre-aggregating changes
*which* computational truncation bound gets enumerated against, since
that per-row cap doesn’t combine linearly under merging. In practice
this essentially never changes the answer (the cap is a generous bound
the true distribution almost never approaches), but it isn’t on the same
provably-exact footing as the other two rules; if that distinction
matters for a particular analysis, verify by comparing the aggregated
and unaggregated fits directly, as below.

The example below is one row per `age_grp` x `arm` cell already (the
smallest this data naturally gets), so there is nothing left to
aggregate in `trial_design` itself. To make the point concretely, this
recreates the same trial as if it had instead arrived as separate
per-site records: three recruiting sites contributing to every cell,
each with its own slice of that cell’s person-years and (simulated)
cases, and confirms that aggregating back down to one row per cell
before fitting gives an *identical* result to fitting the per-site data
directly:

``` r

set.seed(20260819)
n_sites <- 3

# Same 6 cells as trial_design, but split into 3 site-level rows each,
# with person_years randomly divided across sites (summing back to the
# original cell total) and cases simulated per site.
trial_design_by_site <- do.call(
  rbind, lapply(seq_len(nrow(trial_design)), function(i) {
  cell <- trial_design[i, ]
  py_split <- as.numeric(
    rmultinom(1, cell$person_years, prob = rep(1, n_sites)))
  data.frame(
    age_grp = cell$age_grp,
    arm = cell$arm,
    site = paste0("site", seq_len(n_sites)),
    person_years = py_split,
    cases = rpois(n_sites, py_split * cell$placebo_rate *
      ifelse(cell$arm == "vaccine", true_rr, 1))
  )
}))
trial_design_by_site <-
  trial_design_by_site[trial_design_by_site$person_years > 0, ]
nrow(trial_design_by_site)
#> [1] 18

# Aggregate back down to one row per (age_grp, arm) pattern: sum cases and
# person_years within each group (no log() needed here, since this
# vignette's offset_log = FALSE throughout; person_years is already the
# raw exposure).
trial_design_agg <- aggregate(
  cbind(cases, person_years) ~ age_grp + arm,
  data = trial_design_by_site,
  FUN = sum
)

fit_by_site <- ep_closed_form(
  formula = cases ~ arm,
  data = trial_design_by_site,
  offset = "person_years",
  offset_log = FALSE,
  strata = "age_grp",
  ref_levels = list(arm = "placebo"),
  verbose = FALSE
)

fit_agg <- ep_closed_form(
  formula = cases ~ arm,
  data = trial_design_agg,
  offset = "person_years",
  offset_log = FALSE,
  strata = "age_grp",
  ref_levels = list(arm = "placebo"),
  verbose = FALSE
)

# Identical estimate, CI, and p-value either way; aggregating first only
# changed how many rows the dynamic program had to process, not the answer.
all.equal(
  fit_by_site$exactparmest[, c("Estimate", "LowerCL", "UpperCL", "PValue")],
  fit_agg$exactparmest[, c("Estimate", "LowerCL", "UpperCL", "PValue")]
)
#> [1] TRUE
```

This particular demonstration uses
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md),
which has no per-observation cap of any kind (see
[`?ep_closed_form`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)’s
Details), so aggregation is exact here unconditionally, independent of
the Rule 3 caveat above.
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
and
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
follow the same convolution argument for the same reason, but do go
through the `EXACTMAX`/offset cap rules described above; the caveat
there only bites when `offset` is supplied without `EXACTMAX`, and even
then only in the rare case where that heuristic cap is actually close to
binding.

The rest of this vignette continues to use `trial_design` as originally
simulated (already one row per cell); the recommendation is simply to do
the equivalent of the
[`aggregate()`](https://rdrr.io/r/stats/aggregate.html) step above on
your own raw data before calling any of this package’s fitting
functions, whenever your data has more rows than distinct covariate
patterns.

## Analysing it the exact way

[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
is the right tool here: one two-level effect of interest (`arm`), with
age strata as the only other structure in the model
([`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 6.1).

``` r

fit_exact <- ep_closed_form(
  formula = cases ~ arm,
  data = trial_design,
  offset = "person_years",
  offset_log = FALSE,
  strata = "age_grp",
  ref_levels = list(arm = "placebo"),
  verbose = FALSE
)
fit_exact
#> 
#> Exact conditional Poisson, two groups + strata (closed-form Guess & Thomas algorithm)
#> -------------------------------------------------------------------------------------
#>   Parameter Estimate StdError Alpha LowerCL    UpperCL  PValue SuffStat
#>  armvaccine     -Inf       NA  0.05    -Inf -0.1633296 0.03125        0
```

The observed vaccine-arm total (`SuffStat`) sits at 0, the very edge of
what this statistic can attain, since it cannot go any lower.
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Sections 8.3-8.4 describe exactly this situation: when the observed
value is at the boundary of its attainable range, the corresponding
confidence limit is set to \mp\infty directly rather than estimated, and
here that shows up as `Estimate = -Inf` and `LowerCL = -Inf`. This is
not a failure or a placeholder; it is the exact, honest answer to “how
small could the true rate ratio be and still be consistent with zero
observed events?”: arbitrarily small, so the interval says so explicitly
instead of printing a number that looks precise but isn’t. The *upper*
limit and the p-value are still perfectly well defined and informative,
as the next chunk shows.

`Estimate` is on the \beta = \log(\text{rate ratio}) scale used
throughout
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md).
Converting to the rate-ratio and vaccine-efficacy scales:

``` r

rr_exact <- exp(fit_exact$exactparmest$Estimate)
rr_ci_exact <- exp(c(
  fit_exact$exactparmest$LowerCL,
  fit_exact$exactparmest$UpperCL
))

ve_exact <- 1 - rr_exact
ve_ci_exact <- rev(1 - rr_ci_exact)

cat(sprintf(
  paste("Exact:  RR = %.3f (95%% CI %.3f, %.3f)",
  "   VE = %.1f%% (95%% CI %.1f%%, %.1f%%)   p = %.4f\n"),
  rr_exact,
  rr_ci_exact[1],
  rr_ci_exact[2],
  100 * ve_exact,
  100 * ve_ci_exact[1],
  100 * ve_ci_exact[2],
  fit_exact$exactparmest$PValue
))
#> Exact:  RR = 0.000 (95% CI 0.000, 0.849)    VE = 100.0% (95% CI 15.1%, 100.0%)   p = 0.0312
```

## Analysing the *same* model with `glm()`

For a fair comparison, [`glm()`](https://rdrr.io/r/stats/glm.html) needs
to fit the identical model. Age strata, here, are a single categorical
variable with nothing else in the model, so
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 4.1-4.2 guarantees that adjusting for the full set of `age_grp`
dummy columns is an exact reparametrisation of stratifying by `age_grp`:
the same model, just written the ordinary
[`glm()`](https://rdrr.io/r/stats/glm.html) way instead of the `strata=`
way.

``` r

fit_glm <- glm(
  cases ~ arm + age_grp + offset(log(person_years)),
  family = poisson(link = "log"),
  data = trial_design
)
summary(fit_glm)$coefficients
#>                 Estimate   Std. Error       z value     Pr(>|z|)
#> (Intercept)   -5.0106353 7.071068e-01 -7.0861083893 1.379358e-12
#> armvaccine   -23.9108390 3.854339e+04 -0.0006203617 9.995050e-01
#> age_grp41-64  -0.5108256 1.224745e+00 -0.4170873753 6.766145e-01
#> age_grp65+     1.0986123 9.128709e-01  1.2034694650 2.287947e-01
```

Look closely at the `armvaccine` row: the standard error is *four orders
of magnitude larger* than the estimate itself (\mathrm{SE} \approx
38{,}543 against an estimate of \approx -23.9). That imbalance is the
signature of quasi-complete separation: with zero observed events in
every vaccine-arm cell, the log-linear model wants to push the vaccine
coefficient to -\infty, and the optimizer stops only because it ran out
of room to keep moving in a fixed number of iterations, not because it
found a well-determined answer. This is exactly the failure mode this
section set out to demonstrate, and the next chunk checks for the
diagnostic a user would normally rely on to catch it:

``` r

cat("fit_glm$converged:", fit_glm$converged, "\n")
#> fit_glm$converged: TRUE

withCallingHandlers(
  glm(
    cases ~ arm + age_grp + offset(log(person_years)),
    family = poisson(link = "log"),
    data = trial_design
  ),
  warning = function(w) cat("Warning raised:", conditionMessage(w), "\n")
)
#> 
#> Call:  glm(formula = cases ~ arm + age_grp + offset(log(person_years)), 
#>     family = poisson(link = "log"), data = trial_design)
#> 
#> Coefficients:
#>  (Intercept)    armvaccine  age_grp41-64    age_grp65+  
#>      -5.0106      -23.9108       -0.5108        1.0986  
#> 
#> Degrees of Freedom: 5 Total (i.e. Null);  2 Residual
#> Null Deviance:       10.87 
#> Residual Deviance: 4.953e-10     AIC: 15.61
cat("(if no \"Warning raised\" line appears above, glm() raised none)\n")
#> (if no "Warning raised" line appears above, glm() raised none)
```

`fit_glm$converged` reports `TRUE`, and
[`glm()`](https://rdrr.io/r/stats/glm.html) raises **no warning at all**
for this fit. There is no `glm.fit: fitted rates numerically 0 occurred`
message, no `algorithm did not converge` message, nothing: by every
signal [`glm()`](https://rdrr.io/r/stats/glm.html) itself exposes, this
fit looks completely unremarkable. A user who checks for warnings before
trusting a model, a common and generally reasonable habit, would see
nothing here and proceed with unwarranted confidence. The next two
chunks show what that unwarranted confidence would actually produce.

The Wald confidence interval and p-value, the usual large-sample
approximation
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 1 warns about, are computed directly from the estimate and its
asymptotic standard error, exactly the textbook construction (this uses
the plain Wald interval, not a profile-likelihood interval, to isolate
the specific approximation being compared against):

``` r

co_glm <- summary(fit_glm)$coefficients["armvaccine", ]
z975 <- qnorm(0.975)

beta_glm <- co_glm["Estimate"]
wald_ci_beta <- beta_glm + c(-1, 1) * z975 * co_glm["Std. Error"]
pval_glm <- co_glm["Pr(>|z|)"]

rr_glm <- exp(beta_glm)
rr_ci_glm <- exp(wald_ci_beta)
ve_glm <- 1 - rr_glm
ve_ci_glm <- rev(1 - rr_ci_glm)

cat(sprintf(
  paste("Wald:   RR = %.3f (95%% CI %.3f, %.3f)",
  "   VE = %.1f%% (95%% CI %.1f%%, %.1f%%)   p = %.4f\n"),
  rr_glm,
  rr_ci_glm[1],
  rr_ci_glm[2],
  100 * ve_glm,
  100 * ve_ci_glm[1],
  100 * ve_ci_glm[2],
  pval_glm
))
#> Wald:   RR = 0.000 (95% CI 0.000, Inf)    VE = 100.0% (95% CI -Inf%, 100.0%)   p = 0.9995
```

## Side by side

``` r

comparison <- data.frame(
  Method = c("Exact (ep_closed_form)", "Wald (glm)"),
  RR = c(rr_exact, rr_glm),
  RR_lo = c(rr_ci_exact[1], rr_ci_glm[1]),
  RR_hi = c(rr_ci_exact[2], rr_ci_glm[2]),
  `VE %` = 100 * c(ve_exact, ve_glm),
  `VE % lo` = 100 * c(ve_ci_exact[1], ve_ci_glm[1]),
  `VE % hi` = 100 * c(ve_ci_exact[2], ve_ci_glm[2]),
  PValue = c(fit_exact$exactparmest$PValue, pval_glm),
  check.names = FALSE
)
knitr::kable(comparison, digits = 4)
```

|          | Method                 |  RR | RR_lo |  RR_hi | VE % | VE % lo | VE % hi | PValue |
|:---------|:-----------------------|----:|------:|-------:|-----:|--------:|--------:|-------:|
|          | Exact (ep_closed_form) |   0 |     0 | 0.8493 |  100 | 15.0689 |     100 | 0.0312 |
| Estimate | Wald (glm)             |   0 |     0 |    Inf |  100 |    -Inf |     100 | 0.9995 |

The two point estimates land in the same place: both report essentially
100% efficacy, since both are estimating the same quantity from the same
all-zero vaccine-arm data. If a reader only glanced at the point
estimates, this would look like two methods agreeing closely. The two
rows tell completely different stories the moment you look past the
point estimate:

- The exact interval, `(0.000, 0.849)` on the rate-ratio scale,
  `(15.1%, 100.0%)` on the efficacy scale, is wide on one side (as it
  should be, given zero events) but still genuinely informative: it
  rules out rate ratios above 0.849, i.e., efficacy below 15%, at the
  95% level, and its p-value of `0.0312` correctly flags this as a
  significant result.
- The Wald interval is `(0.000, Inf)` on the rate-ratio scale and
  `(-Inf%, 100.0%)` on the efficacy scale: literally the entire range of
  values either quantity could ever take. It rules out nothing. Its
  p-value, `0.9995`, says the opposite of the exact method’s conclusion:
  no evidence of any effect at all. And, as the previous section showed,
  [`glm()`](https://rdrr.io/r/stats/glm.html) reached this completely
  uninformative, wrong-direction answer while reporting
  `converged = TRUE` and raising no warning whatsoever.

This is the central point of this example: a confidence interval
spanning every value a parameter could possibly take is every bit as
much a sign of a broken fit as an error message would be, but nothing
about [`glm()`](https://rdrr.io/r/stats/glm.html)’s own output says so.
Checking `summary(fit_glm)` for a plausible-looking table of
coefficients, or checking
[`warnings()`](https://rdrr.io/r/base/warnings.html) for a diagnostic
message, both miss it. The width of the interval itself, and the sheer
size of the standard error relative to the estimate, are the only
signals available, and they are signals a user has to know to look for.
The exact interval never has this problem, not because it happens to
give a narrower answer here, but because its width is determined
entirely by what the data can and cannot rule out, with no large-sample
approximation to silently break down in the first place. Rerunning the
simulation chunk above with a different seed will sometimes produce a
nonzero (if still small) vaccine-arm count, in which case both methods
will report finite, comparable intervals; this is exactly the comparison
the Monte Carlo study at the end of this vignette runs systematically,
rather than leaving it to one draw of the dice.

## Cross-checking with the general algorithms

[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 4.1-4.2 is not just a theoretical curiosity here; it is directly
checkable. Adjusting for `age_grp` as an ordinary covariate (instead of
stratifying by it) should give *exactly* the same exact inference,
because age is the only nuisance structure in the model and the
adjustment is fully saturated:

``` r

fit_multishift <- ep_multi_shift(
  formula = cases ~ arm + age_grp,
  data = trial_design,
  offset = "person_years",
  offset_log = FALSE,
  exact_vars = "arm",
  ref_levels = list(arm = "placebo"),
  verbose = FALSE
)
fit_multishift$exactparmest
#>    Parameter Estimate StdError Alpha LowerCL    UpperCL  PValue SuffStat
#> 1 armvaccine     -Inf       NA  0.05    -Inf -0.1633296 0.03125        0
```

``` r

fit_network <- ep_network(
  formula = cases ~ arm + age_grp,
  data = trial_design,
  offset = "person_years",
  offset_log = FALSE,
  exact_vars = "arm",
  ref_levels = list(arm = "placebo"),
  verbose = FALSE
)
fit_network$exactparmest
#>    Parameter Estimate StdError Alpha LowerCL    UpperCL  PValue SuffStat
#> 1 armvaccine     -Inf       NA  0.05    -Inf -0.1633296 0.03125        0
```

``` r

all.equal(
  fit_exact$exactparmest[, c("Estimate", "LowerCL", "UpperCL", "PValue")],
  fit_multishift$exactparmest[, c("Estimate", "LowerCL", "UpperCL", "PValue")],
  tolerance = 1e-6
)
#> [1] TRUE
all.equal(
  fit_multishift$exactparmest[, c("Estimate", "LowerCL", "UpperCL", "PValue")],
  fit_network$exactparmest[, c("Estimate", "LowerCL", "UpperCL", "PValue")],
  tolerance = 1e-6
)
#> [1] TRUE
```

All three algorithms (the closed-form special case and both fully
general algorithms) agree, exactly as
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Sections 4 and 10.3 predict. This is a useful habit generally, not just
for this example: running an unfamiliar model through more than one
algorithm and comparing `exactparmest` is cheap insurance before relying
on a result for anything important.

## Reusing an already-computed exact distribution

Here is the part of this vignette that answers a different question:
once
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
has run, is every subsequent question about the same effect stuck
re-running the whole enumeration?

[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 6.5 derives the fact that makes the answer no: because the
exponential-tilt family (3.3)/(6.3) sums to one at *every* value of
\beta, a tail probability, and therefore a p-value against any null, or
a confidence limit at any \alpha, can be recomputed directly from the
*already-built* coefficient table, without rebuilding it.
`fit_exact$OUTDIST` *is* that table:

``` r

parm <- "armvaccine"
dist1d <- fit_exact$OUTDIST[[parm]]
head(dist1d)
#>   SuffStat     prob LogWeight
#> 1        0 0.015625  0.000000
#> 2        1 0.093750  1.791759
#> 3        2 0.234375  2.708050
#> 4        3 0.312500  2.995732
#> 5        4 0.234375  2.708050
#> 6        5 0.093750  1.791759
t_obs <- fit_exact$exactparmest$SuffStat[
  fit_exact$exactparmest$Parameter == parm
]
t_obs
#> [1] 0
```

`SuffStat` is the attainable values of the sufficient statistic T (the
vaccine-arm case count, summed over strata); `LogWeight` is \log C(t)
(Section 3.4); `prob` is the null (\beta = 0) probability of each value,
already normalized. Turning this table into a p-value or a confidence
interval at a *different* null value or significance level is exactly
what
[`recompute_exact_inference()`](https://gsk-biostatistics.github.io/poissonexact/reference/recompute_exact_inference.md)
does: it is built from the same tilt-probability/root-finding machinery
every fitting function already uses internally, exported as one small
convenience wrapper, (see
[`?recompute_exact_inference`](https://gsk-biostatistics.github.io/poissonexact/reference/recompute_exact_inference.md)).

### A different null hypothesis: is efficacy above 50%, not just above 0%?

The original fit tested the usual null of no effect (\mathrm{RR}=1,
\beta=0). A regulator might instead want to know whether this trial’s
data are consistent with a *pre-specified superiority margin*, for
example, whether the observed efficacy is significantly above 50%,
i.e. H_0: \mathrm{RR}=0.50 against H_1: \mathrm{RR}\<0.50 (equivalently,
true VE exceeds 50%):

``` r

beta0_margin <- log(1 - 0.50) # beta corresponding to a 50% efficacy margin

margin_test <- recompute_exact_inference(fit_exact, parm, beta0 = beta0_margin)

cat(sprintf(
  "One-sided p-value for H0: VE = 50%% vs H1: VE > 50%%:  p = %.4f\n",
  margin_test$p_lower
))
#> One-sided p-value for H0: VE = 50% vs H1: VE > 50%:  p = 0.0878
```

Fewer vaccine-arm cases than the margin implies gives a small `p_lower`
(evidence that true efficacy exceeds 50%); this is read directly off the
lower tail of the *same* `dist1d` table used for the original test,
evaluated at a different \beta_0: no new dynamic program, no new
convolution, no new network pass.

### A different confidence level: 90% instead of 95%

Likewise, a 90% interval instead of the original 95% one is a different
\alpha applied to the same table:

``` r

ci_90 <- recompute_exact_inference(fit_exact, parm, alpha = 0.10)

rr_ci_90 <- exp(c(ci_90$LowerCL, ci_90$UpperCL))
ve_ci_90 <- rev(1 - rr_ci_90)

cat(sprintf(
  "90%% CI:  RR (%.3f, %.3f)   VE (%.1f%%, %.1f%%)\n",
  rr_ci_90[1],
  rr_ci_90[2],
  100 * ve_ci_90[1],
  100 * ve_ci_90[2]
))
#> 90% CI:  RR (0.000, 0.648)   VE (35.2%, 100.0%)
```

Both of these (a different null, and a different confidence level) are
recomputed from `fit_exact$OUTDIST` alone. For this two-stratum-effect
example the combinatorial work behind the original fit was already
cheap, so the saving is mostly conceptual. But for a general model with
several categorical covariates, run through
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
or
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md),
building the joint coefficient table is exactly the expensive part
([`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 10); this is the part that answering a second question about the
same effect never has to repeat.

## A Monte Carlo check: does the nominal coverage actually hold?

Rather than rely on a single simulated trial to make the point, this
section repeats the simulation many times and refits both ways each
time. Every replicate is kept, including every all-zero-vaccine-arm
replicate (roughly a 1-in-8 chance each time, e^{-2.1} \approx 0.122),
because the point of this section is to measure what actually happens
across the full range of outcomes this design can produce, boundary
cases included, rather than to report a tidy number that only describes
the well-behaved fits.

Two things make a plain confidence-interval-width comparison a poor fit
for that goal here. First, a replicate at the boundary gives the exact
method an `-Inf` lower limit, which forces any width computed from it to
be `Inf`; averaging widths across replicates that include even one such
case is not informative. Second, and more importantly,
[`glm()`](https://rdrr.io/r/stats/glm.html)’s own `$converged` flag and
its warnings are both supposed to be the tools for catching a bad fit,
and the single worked example above already showed both can stay silent.
So instead of interval width, this simulation tracks the diagnostics a
user would actually have available: whether each method’s interval
excludes the null in the correct direction (statistical significance for
a genuine vaccine effect), whether
[`glm()`](https://rdrr.io/r/stats/glm.html) itself reports having
converged, and whether its standard error is a plausible number at all.

``` r

n_sims <- 1000
set.seed(20260731)

# On the beta (log rate ratio) scale, a well-behaved fit in this design never
# produces a standard error above about 1.5; a fit that has hit the all-zero
# vaccine-arm boundary produces one in the tens of thousands (see the single
# worked example earlier in this vignette). A threshold of 10 sits with a
# wide safety margin between the two and cleanly separates them.
se_large_cutoff <- 10

one_rep <- function() {
  d <- trial_design
  d$cases <- rpois(nrow(d), d$expected_cases)

  fit_i <- ep_closed_form(
    cases ~ arm,
    data = d,
    offset = "person_years",
    offset_log = FALSE,
    strata = "age_grp",
    ref_levels = list(arm = "placebo"),
    verbose = FALSE
  )
  exact_row <- fit_i$exactparmest
  exact_covers <- (exact_row$LowerCL <= true_beta) &&
    (true_beta <= exact_row$UpperCL)
  # "Significant" here means the entire CI for beta sits below 0, i.e., the
  # rate ratio's CI sits below 1, i.e., the VE% CI's lower limit sits above 0%.
  exact_sig <- exact_row$UpperCL < 0

  ws <- character(0)
  glm_i <- withCallingHandlers(
    glm(
      cases ~ arm + age_grp + offset(log(person_years)),
      family = poisson,
      data = d
    ),
    warning = function(w) {
      ws <<- c(ws, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  co_i <- summary(glm_i)$coefficients["armvaccine", ]
  se_i <- co_i["Std. Error"]
  lo_i <- co_i["Estimate"] - z975 * se_i
  hi_i <- co_i["Estimate"] + z975 * se_i
  wald_covers <- (lo_i <= true_beta) && (true_beta <= hi_i)
  wald_sig <- hi_i < 0

  data.frame(
    exact_covers = exact_covers,
    exact_sig = exact_sig,
    wald_covers = wald_covers,
    wald_sig = wald_sig,
    glm_not_converged = !isTRUE(glm_i$converged),
    se_large = se_i > se_large_cutoff,
    glm_flagged = length(ws) > 0
  )
}

mc_results <- do.call(rbind, lapply(seq_len(n_sims), function(i) one_rep()))

mc_summary <- data.frame(
  Method = c("Exact (closed-form)", "Wald (glm)"),
  `Nominal level` = "95%",
  `Empirical coverage` = c(
    mean(mc_results$exact_covers),
    mean(mc_results$wald_covers)
  ),
  `% VE significant` = c(
    100 * mean(mc_results$exact_sig),
    100 * mean(mc_results$wald_sig)
  ),
  `% glm() did not converge` = c(NA, 100 * mean(mc_results$glm_not_converged)),
  `% SE very large` = c(NA, 100 * mean(mc_results$se_large)),
  `% glm() flagged a warning` = c(NA, 100 * mean(mc_results$glm_flagged)),
  check.names = FALSE
)
knitr::kable(mc_summary, digits = 3)
```

| Method | Nominal level | Empirical coverage | % VE significant | % glm() did not converge | % SE very large | % glm() flagged a warning |
|:---|:---|---:|---:|---:|---:|---:|
| Exact (closed-form) | 95% | 0.983 | 25.5 | NA | NA | NA |
| Wald (glm) | 95% | 0.968 | 17.0 | 0 | 12.5 | 3.3 |

Every replicate is included in every column above, boundary and
non-boundary alike, so this table is a direct empirical answer to “what
actually happens across 1000 trials of this size,” not a description of
only the trials that happened to behave nicely.

**`% glm() did not converge` is 0.0%.**
[`glm()`](https://rdrr.io/r/stats/glm.html)’s own convergence flag never
once reports a problem in this experiment, including on every replicate
where the vaccine arm came back with zero events. Checking `$converged`
is not a safety net here: it is uniformly, uninformatively `TRUE`.

**`% SE very large` is 12.5%,** and this is not a coincidence: it
matches the all-zero-vaccine-arm boundary rate almost exactly, because
that boundary is precisely what drives the standard error to an
implausible size in this design. Unlike `$converged`, a standard error
of the estimate’s own magnitude is at least, in principle, something a
careful user could notice; the trouble, as the single worked example
demonstrated, is that nothing in
[`glm()`](https://rdrr.io/r/stats/glm.html)’s default output forces that
noticing to happen.

**`% glm() flagged a warning` is only 3.3%,** far short of the 12.5% of
replicates that actually have a blown-up standard error. Most of the
replicates driving that 12.5% figure produce no diagnostic message
whatsoever; the warning some replicates do trigger
(`glm.fit: fitted rates numerically 0 occurred`) depends on the fitted
rate underflowing to exact floating-point zero, which happens on some
but not all of the boundary replicates, essentially by numerical
accident.

**`% VE significant`** is where the practical cost of this shows up most
directly: the exact method flags a significant vaccine effect on 25.5%
of replicates, against only 17.0% for the Wald method, even though both
are analysing the identical simulated data every time. The gap is
concentrated exactly among the boundary replicates: with zero
vaccine-arm events, the exact method’s `UpperCL` is still a small,
finite, sub-zero number (as the single worked example showed, `UpperCL`
\approx -0.16, correctly significant), while the Wald method’s own
`hi_i` is a huge positive number on those same replicates, so it never
registers as significant at all.
[`glm()`](https://rdrr.io/r/stats/glm.html) does not merely produce an
uninformative interval on a boundary replicate; it actively loses the
power to detect a real effect exactly where the data are sparsest, which
is exactly where a trial’s analysis plan needs that power the most.

## A more general model: adjusting for a covariate on top of stratification

Every fit so far has used `age_grp` purely as `strata`, the one special
case
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
handles, where strata are the *only* nuisance structure in the model
([`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 6.1). Real studies often need more than that: an ordinary
categorical covariate to adjust for, in addition to stratification, that
is not itself under test. Suppose this trial had also recorded each
participant’s baseline risk status (`risk_group`: `"standard"` or
`"high"`), a real confounder (high-risk participants have a higher
baseline incidence rate regardless of arm or age) that the analysis
should adjust for without testing it. That is exactly the model
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
cannot fit (it only accepts a single two-level effect, with strata as
the sole other structure); this is what
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
and
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
are for.

``` r

set.seed(20260819)

design_adj <- expand.grid(
  age_grp = age_grp,
  arm = c("placebo", "vaccine"),
  risk_group = factor(c("standard", "high"), levels = c("standard", "high")),
  KEEP.OUT.ATTRS = FALSE
)
design_adj$arm <- factor(design_adj$arm, levels = c("placebo", "vaccine"))

placebo_rate_by_age <- c(
  "18-40" = 0.0060,
  "41-64" = 0.0100,
  "65+"   = 0.0180
)
person_years_by_age <- c(
  "18-40" = 2000,
  "41-64" = 1500,
  "65+"   = 500
)
risk_mult  <- c(standard = 1.0, high = 1.6)
risk_share <- c(standard = 0.65, high = 0.35)

design_adj$placebo_rate <-
  placebo_rate_by_age[as.character(design_adj$age_grp)]
design_adj$person_years <-
  person_years_by_age[as.character(design_adj$age_grp)] *
  risk_share[as.character(design_adj$risk_group)]
design_adj$true_rate <- with(
  design_adj,
  placebo_rate * risk_mult[as.character(risk_group)] *
    ifelse(arm == "vaccine", true_rr, 1)
)
design_adj$expected_cases <- design_adj$person_years * design_adj$true_rate
design_adj$cases <- rpois(nrow(design_adj), design_adj$expected_cases)

design_adj[, c(
  "age_grp",
  "arm",
  "risk_group",
  "person_years",
  "expected_cases",
  "cases"
)]
#>    age_grp     arm risk_group person_years expected_cases cases
#> 1    18-40 placebo   standard         1300          7.800     9
#> 2    41-64 placebo   standard          975          9.750     6
#> 3      65+ placebo   standard          325          5.850     7
#> 4    18-40 vaccine   standard         1300          2.340     3
#> 5    41-64 vaccine   standard          975          2.925     1
#> 6      65+ vaccine   standard          325          1.755     0
#> 7    18-40 placebo       high          700          6.720     7
#> 8    41-64 placebo       high          525          8.400     4
#> 9      65+ placebo       high          175          5.040     3
#> 10   18-40 vaccine       high          700          2.016     0
#> 11   41-64 vaccine       high          525          2.520     3
#> 12     65+ vaccine       high          175          1.512     1
```

The vaccine effect (`arm`) is still constant across every stratum and
every risk group, by construction (true VE = 70%, exactly as before);
only the *baseline* rate differs by `risk_group`, which is exactly what
“adjusting for a confounder” means: something that shifts everyone’s
baseline risk without itself being the effect of interest, and without
modifying how strongly `arm` acts. `exact_vars = "arm"` tests only the
vaccine effect; `risk_group` enters `formula` like any other covariate
and is conditioned on (by default, `condition_on = "model"`) rather than
tested:

``` r

fit_adj <- ep_multi_shift(
  formula = cases ~ arm + risk_group,
  data = design_adj,
  offset = "person_years",
  offset_log = FALSE,
  strata = "age_grp",
  exact_vars = "arm",
  ref_levels = list(arm = "placebo", risk_group = "standard"),
  verbose = FALSE
)
fit_adj$exactparmest
#>    Parameter  Estimate StdError Alpha   LowerCL    UpperCL      PValue SuffStat
#> 1 armvaccine -1.504077 0.390868  0.05 -2.416557 -0.7211634 2.54486e-05        8
```

``` r

rr_adj <- exp(fit_adj$exactparmest$Estimate)
rr_ci_adj <- exp(c(fit_adj$exactparmest$LowerCL, fit_adj$exactparmest$UpperCL))
ve_adj <- 1 - rr_adj
ve_ci_adj <- rev(1 - rr_ci_adj)

cat(sprintf(
  paste("Exact (adjusted for risk_group): RR = %.3f (95%% CI %.3f, %.3f)",
   "  VE = %.1f%% (95%% CI %.1f%%, %.1f%%)   p = %.4f\n"),
  rr_adj,
  rr_ci_adj[1],
  rr_ci_adj[2],
  100 * ve_adj,
  100 * ve_ci_adj[1],
  100 * ve_ci_adj[2],
  fit_adj$exactparmest$PValue
))
#> Exact (adjusted for risk_group): RR = 0.222 (95% CI 0.089, 0.486)   VE = 77.8% (95% CI 51.4%, 91.1%)   p = 0.0000
```

As in every other example in this vignette,
[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
must agree exactly with
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)
on this same model: both algorithms answer the identical question, just
by different computational routes
([`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Sections 4 and 10.3):

``` r

fit_adj_net <- ep_network(
  formula = cases ~ arm + risk_group,
  data = design_adj,
  offset = "person_years",
  offset_log = FALSE,
  strata = "age_grp",
  exact_vars = "arm",
  ref_levels = list(arm = "placebo", risk_group = "standard"),
  verbose = FALSE
)
all.equal(
  fit_adj$exactparmest[, c("Estimate", "LowerCL", "UpperCL", "PValue")],
  fit_adj_net$exactparmest[, c("Estimate", "LowerCL", "UpperCL", "PValue")]
)
#> [1] TRUE
```

For comparison, the asymptotic
[`glm()`](https://rdrr.io/r/stats/glm.html) fit of the identical model
(`arm`, `age_grp`, and `risk_group` all as ordinary covariates, exactly
the same saturation argument
[`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
Section 4.1-4.2 makes for strata alone):

``` r

glm_adj <- glm(
  cases ~ arm + age_grp + risk_group + offset(log(person_years)),
  family = poisson(link = "log"),
  data = design_adj
)
summary(glm_adj)$coefficients
#>                   Estimate Std. Error      z value     Pr(>|z|)
#> (Intercept)    -4.95244436  0.2709545 -18.27777443 1.243997e-74
#> armvaccine     -1.50407740  0.3908673  -3.84805136 1.190611e-04
#> age_grp41-64   -0.01769958  0.3522216  -0.05025125 9.599222e-01
#> age_grp65+      0.83975065  0.3788675   2.21647606 2.665892e-02
#> risk_grouphigh  0.25131443  0.3066219   0.81962309 4.124310e-01
```

Finally, note what happens if
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
is asked to fit this same model: it refuses, rather than silently
ignoring `risk_group` and returning an unadjusted (and therefore biased)
estimate:

``` r

ep_closed_form(
  formula = cases ~ arm + risk_group,
  data = design_adj,
  offset = "person_years",
  offset_log = FALSE,
  strata = "age_grp",
  ref_levels = list(arm = "placebo", risk_group = "standard"),
  verbose = FALSE
)
#> Error:
#> ! `formula` must contain exactly one right-hand-side term (the
#>                 effect to test), with strata and the offset supplied separately
#>                 via the `strata`/`offset` arguments, not as additional terms in
#>                 `formula`; found term(s): arm, risk_group. Use `ep_multi_shift()` or
#>                 `ep_network()` for models with additional categorical
#>                 covariates.
```

That error message is the package telling you which function to use: any
time there’s more structure to condition on than strata alone,
[`ep_closed_form()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_closed_form.md)
is the wrong tool, and
[`ep_multi_shift()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_multi_shift.md)/[`ep_network()`](https://gsk-biostatistics.github.io/poissonexact/reference/ep_network.md)
are what this vignette’s earlier sections on reusing an already-built
distribution and recomputing inference at a different null/alpha apply
to just as well.

## References

- See
  [`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
  for the full derivation of every method used here, and for the
  complete reference list (Guess, Lydick, Small and Miller 1987; Guess
  and Thomas 1990; Hirji, Mehta and Patel 1987; Mehta, Patel and
  Senchaudhuri 1992; and others).
- ICH E9(R1), *Addendum on Estimands and Sensitivity Analysis in
  Clinical Trials*, for the convention (used in
  [`vignette("Theory-Implementation")`](https://gsk-biostatistics.github.io/poissonexact/articles/Theory-Implementation.md)
  Section 4.7 and followed here) of carrying a trial’s
  randomization-stratification factors into the primary analysis as
  strata.
