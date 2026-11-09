# ==============================================================================
# Shared fixtures for the test suite ----
# ==============================================================================

#' Reorder every element of a list-of-vectors "table" by one or more of
#' its own elements, ascending
#'
#' The package's internal joint-table representation is now a plain named
#' list of equal-length vectors (see `.tbl_nrow()`/`.tbl_subset()` in
#' `ep-utils.R`).
#'
#' @param tbl A named list of equal-length vectors.
#' @param cols Character vector of element names to sort by, in priority
#'   order.
#' @return `tbl`, with every element reordered the same way.
order_tbl <- function(tbl, cols) {
    ord <- do.call(order, tbl[cols])
    lapply(tbl, function(col) col[ord])
}
# ------------------------------------------------------------------------
# test1 ----
# the small, hand-checkable 3-stratum dataset used throughout the
# theory vignette's own worked examples and in the user's SAS-validated
# tests.R (age = 50/60/70; group = treatment/control).
# ------------------------------------------------------------------------
test1 <- data.frame(
    group = c(1L, 1L, 1L, 0L, 0L, 0L),
    age = c(50L, 60L, 70L, 50L, 60L, 70L),
    subjects = c(3492L, 2141L, 1711L, 3525L, 2166L, 1724L),
    cases = c(3L, 2L, 1L, 10L, 7L, 5L),
    futime = c(50, 60, 40, 51, 62, 42),
    all = c(1L, 1L, 1L, 1L, 1L, 1L)
)
test1$lnfutime <- log(test1$futime)

# ------------------------------------------------------------------------
# test2 ----
# a larger synthetic dataset with a full 2 x 4 x 2 x 4 factorial
# design (group x age x gender x region), Poisson counts simulated from
# known coefficients with a fixed seed, reproducible, and rich enough to
# exercise multiple categorical covariates, interactions, and every
# EXACTMAX/OFFSET cap rule. This mirrors the second scenario in the user's
# own tests.R exactly (same seed, same coefficients), so results here
# should match whatever SAS output was obtained from that same data.
# ------------------------------------------------------------------------
.make_test2 <- function(seed = 420) {
    group <- c(0, 1)
    age <- c(50, 60, 70, 80)
    gender <- c(0, 1)
    region <- c(1, 2, 3, 4)

    result <- expand.grid(
        group = group,
        age = age,
        gender = gender,
        region = region
    )
    result$age <- as.factor(result$age)
    result$gender <- as.factor(result$gender)
    result$region <- as.factor(result$region)

    set.seed(seed)
    result$futime <- round(runif(nrow(result), min = 5, max = 8), 1)

    beta_0 <- -1.5
    beta_group <- -0.4
    beta_age60 <- 0.1
    beta_age70 <- 0.2
    beta_age80 <- 0.3
    beta_gender1 <- -0.1
    beta_region2 <- 0.1
    beta_region3 <- -0.2
    beta_region4 <- 0.2

    coefs <- c(
        "(Intercept)" = beta_0,
        group = beta_group,
        age60 = beta_age60,
        age70 = beta_age70,
        age80 = beta_age80,
        gender1 = beta_gender1,
        region2 = beta_region2,
        region3 = beta_region3,
        region4 = beta_region4
    )

    x_mm <- model.matrix(~ group + age + gender + region, data = result)
    log_rate <- x_mm %*% coefs
    expected_cases <- exp(log_rate) * result$futime / 2
    result$cases <- rpois(nrow(result), lambda = expected_cases)
    result$lnfutime <- log(result$futime)

    result$group <- relevel(factor(result$group), ref = "0")
    result$region <- relevel(factor(result$region), ref = "1")
    result$age <- relevel(factor(result$age), ref = "50")
    result$gender <- relevel(factor(result$gender), ref = "0")

    result
}
test2 <- .make_test2()

# ------------------------------------------------------------------------
# test3 ----
# a small, fully deterministic multi-covariate dataset used ONLY
# for cross-algorithm *agreement* checks (multi_shift vs. network) that
# need more than one adjustment covariate and/or an interaction term, but
# don't need to match any specific SAS-validated number. test2 must stay
# exactly as constructed above because test-sas-validated-scenarios.R's
# fit21-fit25 are pinned to it (real SAS output was obtained from that
# exact data); this fixture exists so the *other*, purely structural
# agreement tests don't have to pay test2's combinatorial cost.
#
# The runtime of the exact algorithms is driven far more by the total
# case count (the conditioning statistic, which sets the size of the
# combinatorial state space) than by row count or factor cardinality
# alone. So instead of shrinking test2's design (2x4x2x4) proportionally,
# test3 keeps region and age at 3 levels each (group and gender stay at
# 2), a 2x3x2x3 factorial (group x age x gender x region, 36 rows),
# and fixes a small, hand-picked case vector (total = 28 cases, values up
# to 3, not just 0/1) rather than simulating from rpois, so the total
# conditioning mass is small and known by construction, not just
# "probably small for this seed". Margins are deliberately kept
# non-degenerate (every level of every factor has at least some mass, and
# both interaction dummy columns have a real signal) so the multi-
# covariate and interaction tests remain meaningful:
#   group total   = 26 (group1) vs 2 (group0)
#   region total  =  7 (region1), 9 (region2), 12 (region3)
#   gender total  = 16 (gender0) vs 12 (gender1)
#   age total     = 14 (age50), 9 (age60), 5 (age70)
#   group1:region2 interaction total = 9
#   group1:region3 interaction total = 11
#
# futime is a constant 5: with offset = "lnfutime" and no
# exactmax, .setup_exact_model() falls back to SAS rule 3
# (si_from_offset(), ep-common.R), which caps every row's
# count at S_i = ceiling(futime). That cap must stay comfortably above
# the largest single-row case count used below (currently 3) or
# .joint_dist_stratum_dp_vec() errors with "Observed y_i exceeds S_i for
# at least one row.", if the case vector is edited to include a larger
# value, futime must be raised to at least that value too.
# ------------------------------------------------------------------------
# fmt: skip
test3 <- data.frame(
    group = c(
        0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L,
        0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 1L
    ),
    age = c(
        50L, 50L, 60L, 60L, 70L, 70L, 50L, 50L, 60L, 60L, 70L, 70L,
        50L, 50L, 60L, 60L, 70L, 70L, 50L, 50L, 60L, 60L, 70L, 70L,
        50L, 50L, 60L, 60L, 70L, 70L, 50L, 50L, 60L, 60L, 70L, 70L
    ),
    gender = c(
        0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L,
        1L, 1L, 1L, 1L, 1L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 1L
    ),
    region = c(
        1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L,
        2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L
    ),
    cases = c(
        0L, 1L, 0L, 2L, 0L, 1L, 1L, 1L, 0L, 1L, 0L, 0L, 0L, 2L, 0L, 1L, 0L, 2L,
        0L, 2L, 0L, 1L, 0L, 1L, 1L, 3L, 0L, 2L, 0L, 1L, 0L, 3L, 0L, 2L, 0L, 0L
    )
)
test3$futime <- 4
test3$lnfutime <- log(test3$futime)
test3$group <- relevel(factor(test3$group), ref = "0")
test3$region <- relevel(factor(test3$region), ref = "1")
test3$gender <- relevel(factor(test3$gender), ref = "0")
test3$age <- relevel(factor(test3$age), ref = "50")

# ------------------------------------------------------------------------
# toy2 ----
# a fully hand-derivable, minimal 2-stratum, single-observation-per-
# arm dataset used as an *analytic* correctness anchor throughout the
# suite (not just an algorithm-agreement check). With one case per stratum
# split symmetrically (stratum 1's case falls in the treatment arm,
# stratum 2's case falls in the control arm) and equal, unit exposure
# everywhere, the per-stratum polynomial is (1 + theta) for each stratum,
# so the whole-sample coefficient sequence is the coefficients of
# (1+theta)^2 = 1 + 2*theta + theta^2, i.e. c_0 = 1, c_1 = 2, c_2 = 1,
# exactly the shape of a Binomial(2, 0.5) distribution. At beta = 0 this
# gives Pr(K=0) = Pr(K=2) = 1/4, Pr(K=1) = 1/2, with the observed
# K = 1 landing exactly on the (symmetric) mean, so by hand:
#   - the conditional MLE (beta_hat) is exactly 0 (mean_t(0) == t_obs == 1)
#   - the two-sided exact p-value is exactly 1
#     (p_lower = p_upper = Pr(K<=1) = Pr(K>=1) = 3/4;
#      two-sided = min(1, 2*min(3/4,3/4)) = 1)
#   - the two-sided exact CI is symmetric about 0 (lo == -hi), because the
#     coefficient sequence is itself symmetric (c_0 == c_2) under the
#     treatment/control relabeling theta -> 1/theta.
# See test-ep-closed-form.R for the full derivation reused
# across multiple test files.
# ------------------------------------------------------------------------
toy2 <- data.frame(
    group = c(1L, 0L, 1L, 0L),
    strat = c(1L, 1L, 2L, 2L),
    cases = c(1L, 0L, 0L, 1L)
)

# ------------------------------------------------------------------------
# Small helper ----
# compare two ep_fit objects' numeric results, matched by Parameter name
# (so row order need not match), used by the cross-algorithm-agreement tests.
# ------------------------------------------------------------------------
expect_exactparmest_equal <- function(a, b, tolerance = 1e-6) {
    cols <- c(
        "Estimate",
        "StdError",
        "LowerCL",
        "UpperCL",
        "PValue",
        "SuffStat"
    )
    ea <- a[order(a$Parameter), c("Parameter", cols)]
    eb <- b[order(b$Parameter), c("Parameter", cols)]
    testthat::expect_identical(ea$Parameter, eb$Parameter)
    for (col in cols) {
        testthat::expect_equal(ea[[col]], eb[[col]], tolerance = tolerance)
    }
}
