# ==============================================================================
# Tests for ep-closed-form.R; the fast, closed-form
# (Guess-Lydick-Small-Miller 1987 / Guess-Thomas 1990) algorithm for a
# single two-level exact effect with strata as the only nuisance
# structure. This is one of the three exported main functions, so it is
# tested extensively: the internal coefficient/cumulative-distribution
# building blocks in isolation, then the exported
# ep_closed_form() end-to-end, anchored on two independently
# verifiable numeric examples:
#
#   * `toy2` (see helper-data.R): a 2-stratum, symmetric-coefficient
#     (c_k = 1, 2, 1) example, hand-derived to give Estimate = 0,
#     PValue = 1, and a CI symmetric about 0.
#
#   * `binom5`: a single-stratum example with unequal treatment/control
#     exposure (n = 2, m = 1), which reduces exactly to a
#     Binomial(5, p = n/(n+m) = 2/3) distribution at beta = 0; so its
#     tail probabilities can be cross-checked directly against R's own
#     `stats::pbinom()`, entirely independently of this package's code,
#     and its exact CMLE works out to the clean closed form
#     beta_hat = log(3/4) (derived in the comments below).
# ==============================================================================

# ------------------------------------------------------------------------
# compute_ck_natlog ----
# ------------------------------------------------------------------------

test_that("compute_ck_natlog reduces to binomial coefficients (one stratum)", {
    # t = x + y = 5, n = m = 1 => log_multiplier = lchoose(5, r) exactly,
    # i.e. C(k) = choose(5, k).
    ck <- compute_ck_natlog(x = 3, y = 2, n = 1, m = 1, max_k = 5)
    expect_identical(ck$k, as.double(0:5))
    expect_equal(exp(ck$log_ck), choose(5, 0:5), tolerance = 1e-10)
})

test_that("compute_ck_natlog matches toy2 hand-derived (1, 2, 1) sequence", {
    # Two strata, each contributing polynomial (1 + theta); their product is
    # '1 + 2*theta + theta^2'.
    ck <- compute_ck_natlog(
        x = c(1, 0),
        y = c(0, 1),
        n = c(1, 1),
        m = c(1, 1),
        max_k = 2
    )
    expect_identical(ck$k, as.double(0:2))
    expect_equal(exp(ck$log_ck), c(1, 2, 1), tolerance = 1e-10)
})

test_that("compute_ck_natlog incorporates unequal exposure via (n/m)^r", {
    # Single stratum, t = 5, n = 2, m = 1: C(k) = choose(5, k) * 2^k.
    ck <- compute_ck_natlog(x = 3, y = 2, n = 2, m = 1, max_k = 5)
    expect_equal(exp(ck$log_ck), choose(5, 0:5) * 2^(0:5), tolerance = 1e-8)
})

test_that("compute_ck_natlog truncated at max_k agrees with full table", {
    full <- compute_ck_natlog(x = 3, y = 2, n = 2, m = 1, max_k = 5)
    trunc <- compute_ck_natlog(x = 3, y = 2, n = 2, m = 1)
    expect_identical(trunc$k, as.double(0:3))
    expect_equal(trunc$log_ck, full$log_ck[1:4], tolerance = 1e-10)
})

test_that("compute_ck_natlog convolves more than two strata correctly", {
    # Three equal-exposure strata each contributing (1 + theta):
    # overall contribution (1 + theta)^3 = 1 + 3*theta + 3*theta^2 + theta^3.
    ck <- compute_ck_natlog(
        x = c(1, 1, 1),
        y = c(0, 0, 0),
        n = c(1, 1, 1),
        m = c(1, 1, 1),
        max_k = 3
    )
    expect_equal(exp(ck$log_ck), c(1, 3, 3, 1), tolerance = 1e-8)
})

# ------------------------------------------------------------------------
# cumulative_dist_natlog ----
# ------------------------------------------------------------------------

test_that("cumulative_dist_natlog at beta=0 matches pbinom, p=0.5", {
    x <- 3
    y <- 2
    n <- 1
    m <- 1
    ck <- compute_ck_natlog(x, y, n, m, max_k = 5)
    cum <- cumulative_dist_natlog(x, y, n, m, ck, beta = 0, equal = TRUE)
    expect_equal(cum, stats::pbinom(5, 5, 0.5), tolerance = 1e-10)
})

test_that("cumulative_dist_natlog reproduces pbinom() for unequal exposure", {
    x <- 3
    y <- 2
    n <- 2
    m <- 1
    ck <- compute_ck_natlog(x, y, n, m, max_k = 5) # full table
    p <- n / (n + m) # = 2/3, the Binomial "success" probability at beta = 0

    cum_le3 <- cumulative_dist_natlog(x, y, n, m, ck, beta = 0, equal = TRUE)
    expect_equal(cum_le3, stats::pbinom(5, 5, p), tolerance = 1e-10)

    # Restrict to the truncated table (k = 0..3, matching k_obs = sum(x)):
    # cumulative_dist_natlog's `equal = TRUE`/`FALSE` cases correspond to
    # P(K <= k_obs) / P(K < k_obs) using only that truncated table, which
    # should equal stats::pbinom() evaluated at the *truncation point* 3,
    # thanks to the analytic normalizing factor `gg`.
    ck_trunc <- compute_ck_natlog(x, y, n, m) # default max_k = sum(x) = 3
    expect_equal(
        cumulative_dist_natlog(x, y, n, m, ck_trunc, beta = 0, equal = TRUE),
        stats::pbinom(3, 5, p),
        tolerance = 1e-10
    )
    expect_equal(
        cumulative_dist_natlog(x, y, n, m, ck_trunc, beta = 0, equal = FALSE),
        stats::pbinom(2, 5, p),
        tolerance = 1e-10
    )
})

test_that("cumulative_dist_natlog at nonzero beta matches pbinom()", {
    x <- 3
    y <- 2
    n <- 2
    m <- 1
    beta <- log(2)
    theta <- exp(beta)
    p_beta <- theta * n / (theta * n + m) # implied Binom prob at this beta
    ck_trunc <- compute_ck_natlog(x, y, n, m)
    expect_equal(
        cumulative_dist_natlog(x, y, n, m, ck_trunc, beta = beta, equal = TRUE),
        stats::pbinom(3, 5, p_beta),
        tolerance = 1e-8
    )
})

# ------------------------------------------------------------------------
# ep_closed_form(): the toy2 anchor (symmetric, strata as nuisance only) ----
# ------------------------------------------------------------------------

test_that("ep_closed_form reproduces every hand-derived toy2 quantity", {
    fit <- ep_closed_form(
        cases ~ group,
        data = toy2,
        strata = "strat",
        alpha = 0.10,
        verbose = FALSE
    )
    expect_s3_class(fit, "ep_closed_form")
    expect_s3_class(fit, "ep_fit")

    row <- fit$exactparmest
    expect_identical(row$Parameter, "group1")
    expect_identical(row$SuffStat, 1)
    expect_equal(
        as.list(row[c("Estimate", "StdError", "PValue")]),
        list(Estimate = 0, StdError = sqrt(2), PValue = 1),
        tolerance = 1e-6
    )
    expect_equal(row$LowerCL, -row$UpperCL, tolerance = 1e-8)

    od <- fit$OUTDIST[["group1"]]
    expect_identical(as.integer(od$SuffStat), 0:2)
    expect_equal(od$prob, c(0.25, 0.5, 0.25), tolerance = 1e-10)
})

test_that("ep_closed_form's default offset (NULL) matches unit offset", {
    toy2_off <- toy2
    toy2_off$one <- 1
    fit_default <- ep_closed_form(
        cases ~ group,
        data = toy2,
        strata = "strat",
        verbose = FALSE
    )
    fit_explicit <- ep_closed_form(
        cases ~ group,
        data = toy2_off,
        strata = "strat",
        offset = "one",
        offset_log = FALSE,
        verbose = FALSE
    )
    expect_equal(
        fit_default$exactparmest,
        fit_explicit$exactparmest,
        tolerance = 1e-10
    )
})

test_that("ep_closed_form's ref_levels flips which level is the effect", {
    fit <- ep_closed_form(
        cases ~ group,
        data = toy2,
        strata = "strat",
        ref_levels = list(group = "1"),
        verbose = FALSE
    )
    expect_identical(fit$exactparmest$Parameter, "group0")
    # toy2's coefficient sequence is palindromic, so swapping which level is
    # "treatment" leaves the (symmetric) point estimate and p-value unchanged.
    expect_equal(fit$exactparmest$Estimate, 0, tolerance = 1e-8)
    expect_equal(fit$exactparmest$PValue, 1, tolerance = 1e-8)
})

test_that("ep_closed_form's one-sided p-value is the un-doubled tail", {
    fit_two <- ep_closed_form(
        cases ~ group,
        data = toy2,
        strata = "strat",
        ci_type = "two-sided",
        verbose = FALSE
    )
    fit_one <- ep_closed_form(
        cases ~ group,
        data = toy2,
        strata = "strat",
        ci_type = "one-sided",
        verbose = FALSE
    )
    expect_equal(fit_one$exactparmest$PValue, 0.75, tolerance = 1e-8)
    expect_equal(fit_two$exactparmest$PValue, 1, tolerance = 1e-8)
    expect_identical(fit_one$exactparmest$LowerCL, -Inf)
})

# ------------------------------------------------------------------------
# ep_closed_form(): the binom5 anchor ----
# (single stratum, unequal exposure, cross-checked against stats::pbinom()
# and a clean closed-form CMLE)
# ------------------------------------------------------------------------
#
# Data: one treatment row (group = 1, cases = 3, fu = 2) and one control
# row (group = 0, cases = 2, fu = 1). No strata. This is exactly the
# (x=3, y=2, n=2, m=1) example above: t = 5, p = n/(n+m) = 2/3 at beta = 0.
#
# Derivation of the exact CMLE: tilting by beta multiplies each
# coefficient C(t,r) = choose(5,r) * 2^r by exp(beta*r), so the combined
# natural parameter is (log(2) + beta); solving E[R] = k_obs = 3 under
# Binomial(5, q) forces q = 3/5, i.e. odds q/(1-q) = 3/2 = exp(log(2) +
# beta_hat), giving beta_hat = log(3/2) - log(2) = log(3/4) exactly. The
# corresponding SE is 1/sqrt(Var) = 1/sqrt(5 * (3/5) * (2/5)) = 1/sqrt(1.2).

binom5 <- data.frame(
    group = c(1, 0),
    cases = c(3, 2),
    fu = c(2, 1)
)

test_that("ep_closed_form's CMLE matches log(3/4), binom5 anchor", {
    fit <- ep_closed_form(
        cases ~ group,
        data = binom5,
        offset = "fu",
        offset_log = FALSE,
        verbose = FALSE
    )
    row <- fit$exactparmest
    expect_identical(row$SuffStat, 3)
    expect_equal(row$Estimate, log(3 / 4), tolerance = 1e-6)
    expect_equal(row$StdError, 1 / sqrt(1.2), tolerance = 1e-6)
})

test_that("ep_closed_form's p-value matches stats::pbinom(), binom5 anchor", {
    fit <- ep_closed_form(
        cases ~ group,
        data = binom5,
        offset = "fu",
        offset_log = FALSE,
        verbose = FALSE
    )
    p <- 2 / 3
    p_lower <- stats::pbinom(3, 5, p)
    p_upper <- 1 - stats::pbinom(2, 5, p)
    expect_equal(
        fit$exactparmest$PValue,
        min(1, 2 * min(p_lower, p_upper)),
        tolerance = 1e-8
    )
})

test_that("ep_closed_form's OUTDIST matches analytic choose(5,k)*2^k shape", {
    fit <- ep_closed_form(
        cases ~ group,
        data = binom5,
        offset = "fu",
        offset_log = FALSE,
        verbose = FALSE
    )
    od <- fit$OUTDIST[["group1"]]
    expect_identical(od$SuffStat, as.double(0:5))
    unnorm <- choose(5, 0:5) * 2^(0:5)
    expect_equal(od$prob, unnorm / sum(unnorm), tolerance = 1e-8)
})

# ------------------------------------------------------------------------
# ep_closed_form(): error paths ----
# ------------------------------------------------------------------------

test_that("ep_closed_form errors when formula contains more than one term", {
    expect_error(
        ep_closed_form(
            cases ~ group + strat,
            data = toy2,
            verbose = FALSE
        ),
        "exactly one right-hand-side term"
    )
})

test_that("ep_closed_form errors when the term resolves to >2 levels", {
    expect_error(
        ep_closed_form(
            cases ~ age,
            data = test1,
            verbose = FALSE
        ),
        "must have exactly two levels"
    )
})

test_that("ep_closed_form errors when a stratum has zero exposure, one arm", {
    bad <- data.frame(
        group = c(1, 0, 1, 0),
        strat = c(1, 1, 2, 2),
        cases = c(2, 1, 0, 3),
        fu = c(5, 5, 5, 0) # stratum 2's control arm has zero exposure
    )
    expect_error(
        ep_closed_form(
            cases ~ group,
            data = bad,
            offset = "fu",
            offset_log = FALSE,
            strata = "strat",
            verbose = FALSE
        ),
        "exposure must be strictly positive"
    )
})

test_that("ep_closed_form's zero-exposure check fires, no rows in one arm", {
    bad <- data.frame(
        group = c(1, 0, 1), # stratum 2 has a treatment row but no control row
        strat = c(1, 1, 2),
        cases = c(2, 1, 3),
        fu = c(5, 5, 5) # every individual exposure is strictly positive
    )
    expect_error(
        ep_closed_form(
            cases ~ group,
            data = bad,
            offset = "fu",
            offset_log = FALSE,
            strata = "strat",
            verbose = FALSE
        ),
        "Every stratum must contribute positive exposure"
    )
})

test_that("ep_closed_form reports '(none)' when formula has no terms", {
    d <- data.frame(
        cases = c(1L, 0L, 2L, 1L),
        group = c(1L, 0L, 1L, 0L)
    )
    expect_error(
        ep_closed_form(
            cases ~ 1,
            data = d,
            verbose = FALSE
        ),
        "\\(none\\)"
    )
})

test_that("ep_closed_form's verbose=TRUE prints the algorithm ID note", {
    expect_message(
        ep_closed_form(
            cases ~ group,
            data = toy2,
            strata = "strat",
            verbose = TRUE
        ),
        "Closed-form two-group"
    )
})

test_that("ep_closed_form's print method shows the closed-form label", {
    fit <- ep_closed_form(
        cases ~ group,
        data = toy2,
        strata = "strat",
        verbose = FALSE
    )
    expect_snapshot(print(fit))
})
