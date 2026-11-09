# poissonexact (development version)

## Exported functions

* `ep_closed_form()`: fast closed-form exact inference for a
  single two-level effect with strata as the only nuisance structure
  (Guess, Lydick, Small and Miller 1987; Guess and Thomas 1990).
* `ep_multi_shift()`: multivariate shift algorithm for
  arbitrary categorical covariates and strata (generalizing Hirji, Mehta
  and Patel 1987 to bounded Poisson counts).
* `ep_network()`: network algorithm with more aggressive
  mid-computation pruning for larger problems (generalizing Mehta, Patel
  and Senchaudhuri 1992 to bounded Poisson counts).
* `recompute_exact_inference()`: recomputes a p-value and/or confidence
  interval from an already-fitted object's `OUTDIST` table at a different
  null value, significance level, `ci_type`, or `direction`, without
  rebuilding the joint distribution.
