# poissonexact: Exact Conditional Inference for Poisson Regression

Compute exact conditional maximum likelihood estimates, exact p-values,
and exact confidence intervals for categorical effects in Poisson
regression models with an optional offset and stratification. Nuisance
parameters (stratum intercepts and other model effects not under test)
are eliminated by conditioning on their sufficient statistics rather
than relying on asymptotic approximations, which matters when case
counts are small or sparse and standard Wald/likelihood-ratio inference
is unreliable. Three algorithms are provided for different problem
shapes: a fast closed-form method for a single two-level effect with
strata as the only nuisance structure, following Guess, Lydick, Small
and Miller (1987)
[doi:10.1093/oxfordjournals.aje.a114535](https://doi.org/10.1093/oxfordjournals.aje.a114535)
and Guess and Thomas (1990)
[doi:10.1097/00001648-199001000-00016](https://doi.org/10.1097/00001648-199001000-00016)
; a multivariate shift algorithm for arbitrary categorical covariates
and strata, generalizing Hirji, Mehta and Patel (1987)
[doi:10.2307/2289388](https://doi.org/10.2307/2289388) to bounded
Poisson counts; and a network algorithm offering more aggressive
mid-computation pruning for larger problems, generalizing Mehta, Patel
and Senchaudhuri (1992)
[doi:10.1080/10618600.1992.10474574](https://doi.org/10.1080/10618600.1992.10474574)
to bounded Poisson counts.

## See also

Useful links:

- <https://github.com/GSK-Biostatistics/poissonexact>

- <https://gsk-biostatistics.github.io/poissonexact/>

- Report bugs at
  <https://github.com/GSK-Biostatistics/poissonexact/issues>

## Author

**Maintainer**: Ali Charkhi <ali.x.charkhi@gsk.com>

Authors:

- Ali Charkhi <ali.x.charkhi@gsk.com>

Other contributors:

- Dragoș Moldovan-Grünfeld <dragos.v.moldovan-grunfeld@gsk.com>
  \[contributor\]

- GlaxoSmithKline Research & Development Limited \[copyright holder,
  funder\]
