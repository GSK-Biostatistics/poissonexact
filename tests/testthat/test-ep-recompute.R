# ==============================================================================
# Tests for ep-recompute.R, the exported
# recompute_exact_inference() convenience wrapper.
#
# recompute_exact_inference() is built entirely from three internal
# pieces already tested directly elsewhere (tilt_probs(),
# two_sided_from_tails(), exact_ci_from_dist1d(), see
# test-ep-utils.R and test-ep-common.R), so these
# tests focus on:
#
#   * that it reproduces a fit's own exactparmest/OUTDIST exactly when
#     called with the same beta0/alpha/ci_type/direction the fit itself
#     used (the "same question" case);
#   * that changing alpha/beta0/ci_type actually recomputes something
#     different, and that the recomputed value agrees with a *fresh* fit
#     built with that different setting (proving the table is genuinely
#     reusable, not just echoed back);
#   * that it works unmodified against all three fitting functions'
#     output (ep_multi_shift(), ep_network(),
#     ep_closed_form()), since its documented contract is "any
#     object of class ep_fit";
#   * its own input validation.
#
# Anchored on `toy2` (see helper-data.R): a fully hand-derived 2-stratum
# example whose null (beta = 0) distribution is exactly (c_0, c_1, c_2) =
# (1, 2, 1), i.e. Pr(T=0) = Pr(T=2) = 0.25, Pr(T=1) = 0.5, with the
# observed T = 1 landing exactly on the symmetric mean. By hand:
#   probability of lower = Pr(T <= 1) = 0.25 + 0.5 = 0.75
#   probability of upper = Pr(T >= 1) = 0.5 + 0.25 = 0.75
#   two sided pvalue = min(1, 2 * min(0.75, 0.75)) = 1
#   one sided pvalue = 0.75
# (matching the PValue = 1 anchor already used throughout the suite).
# ==============================================================================

fit_toy2_shift <- ep_multi_shift(
    cases ~ group,
    data = toy2,
    strata = "strat",
    exact_vars = "group",
    verbose = FALSE
)

# ------------------------------------------------------------------------
# Reproduces the same question the original fit already answered ----
# ------------------------------------------------------------------------

test_that("recompute_exact_inference matches toy2 tails at beta0=0", {
    out <- recompute_exact_inference(fit_toy2_shift, "group1")

    expect_equal(
        out[c(
            "beta0",
            "alpha",
            "p_lower",
            "p_upper",
            "pvalue_two_sided",
            "pvalue_one_sided"
        )],
        list(
            beta0 = 0,
            alpha = 0.05,
            p_lower = 0.75,
            p_upper = 0.75,
            pvalue_two_sided = 1,
            pvalue_one_sided = 0.75
        ),
        tolerance = 1e-8
    )
})

test_that("recompute_exact_inference's defaults reproduce fit$exactparmest", {
    out <- recompute_exact_inference(fit_toy2_shift, "group1")
    row <- fit_toy2_shift$exactparmest

    expect_equal(
        out[c("pvalue_two_sided", "LowerCL", "UpperCL")],
        list(
            pvalue_two_sided = row$PValue,
            LowerCL = row$LowerCL,
            UpperCL = row$UpperCL
        ),
        tolerance = 1e-10
    )
    # toy2's coefficient sequence is symmetric, so the two-sided CI is
    # symmetric about 0.
    expect_equal(out$LowerCL, -out$UpperCL, tolerance = 1e-8)
})

# ------------------------------------------------------------------------
# Recomputing at a different alpha ----
# reuses the table but genuinely recomputes the interval, cross-checked
# against a *fresh* fit built at that alpha, not just against itself.
# ------------------------------------------------------------------------

test_that("recompute_exact_inference at a different alpha matches fresh fit", {
    out_90 <- recompute_exact_inference(fit_toy2_shift, "group1", alpha = 0.10)

    fit_90 <- ep_multi_shift(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        alpha = 0.10,
        verbose = FALSE
    )

    expect_equal(
        out_90[c("LowerCL", "UpperCL")],
        list(
            LowerCL = fit_90$exactparmest$LowerCL,
            UpperCL = fit_90$exactparmest$UpperCL
        ),
        tolerance = 1e-8
    )
    # The interval actually changed (wider alpha -> narrower interval), so
    # this isn't trivially passing because both sides are unchanged.
    expect_false(
        isTRUE(
            all.equal(
                out_90$UpperCL,
                fit_toy2_shift$exactparmest$UpperCL,
                tolerance = 1e-8
            )
        )
    )
})

test_that("recompute_exact_inference's one-sided CI matches a fresh fit", {
    # toy2's p_lower and p_upper are exactly tied (both 0.75), so
    # two_sided_from_tails() picks direction = "upper" (its documented
    # tie-breaking rule: direction = "upper" iff p_lower <= p_upper),
    # matching what a fit built with ci_type = "one-sided" would use
    # internally for this same data.
    out <- recompute_exact_inference(
        fit_toy2_shift,
        "group1",
        ci_type = "one-sided",
        direction = "upper"
    )

    fit_one_sided <- ep_multi_shift(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        ci_type = "one-sided",
        verbose = FALSE
    )

    expect_identical(out$LowerCL, -Inf)
    expect_identical(fit_one_sided$exactparmest$LowerCL, -Inf)
    expect_equal(
        out$UpperCL,
        fit_one_sided$exactparmest$UpperCL,
        tolerance = 1e-8
    )
})

# ------------------------------------------------------------------------
# Recomputing at a different beta0 (a different null / margin) ----
# is checked directly against the shared tilt_probs() primitive, not just
# against recompute_exact_inference() itself.
# ------------------------------------------------------------------------

test_that("recompute_exact_inference's beta0!=0 tails match tilt_probs()", {
    beta0 <- log(0.5)
    dist1d <- fit_toy2_shift$OUTDIST[["group1"]]
    t_obs <- fit_toy2_shift$exactparmest$SuffStat[
        fit_toy2_shift$exactparmest$Parameter == "group1"
    ]

    p0 <- tilt_probs(beta0, dist1d$SuffStat, dist1d$LogWeight)
    expected_p_lower <- sum(p0[dist1d$SuffStat <= t_obs])
    expected_p_upper <- sum(p0[dist1d$SuffStat >= t_obs])

    out <- recompute_exact_inference(fit_toy2_shift, "group1", beta0 = beta0)

    expect_equal(
        out[c("beta0", "p_lower", "p_upper")],
        list(
            beta0 = beta0,
            p_lower = expected_p_lower,
            p_upper = expected_p_upper
        ),
        tolerance = 1e-10
    )
})

test_that("recompute_exact_inference at beta0=0 explicit matches default", {
    out_default <- recompute_exact_inference(fit_toy2_shift, "group1")
    out_explicit <- recompute_exact_inference(
        fit_toy2_shift,
        "group1",
        beta0 = 0
    )
    expect_identical(out_default, out_explicit)
})

# ------------------------------------------------------------------------
# Works against every fitting function's output ----
# its documented contract is "any object of class ep_fit", not
# just ep_multi_shift()'s.
# ------------------------------------------------------------------------

test_that("recompute_exact_inference works against ep_network()'s output", {
    fit_network <- ep_network(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        verbose = FALSE
    )
    out <- recompute_exact_inference(fit_network, "group1")
    row <- fit_network$exactparmest

    expect_equal(
        out[c("pvalue_two_sided", "LowerCL", "UpperCL")],
        list(
            pvalue_two_sided = row$PValue,
            LowerCL = row$LowerCL,
            UpperCL = row$UpperCL
        ),
        tolerance = 1e-10
    )
})

test_that("recompute_exact_inference works against ep_closed_form()'s output", {
    fit_closed <- ep_closed_form(
        cases ~ group,
        data = toy2,
        strata = "strat",
        verbose = FALSE
    )
    out <- recompute_exact_inference(fit_closed, "group1")
    row <- fit_closed$exactparmest

    expect_equal(
        out[c("pvalue_two_sided", "LowerCL", "UpperCL")],
        list(
            pvalue_two_sided = row$PValue,
            LowerCL = row$LowerCL,
            UpperCL = row$UpperCL
        ),
        tolerance = 1e-10
    )
})

# ------------------------------------------------------------------------
# Input validation ----
# ------------------------------------------------------------------------

test_that("recompute_exact_inference errors on a fit that isn't an ep_fit", {
    expect_error(
        recompute_exact_inference(
            list(OUTDIST = list(), exactparmest = data.frame()),
            "group1"
        ),
        "ep_multi_shift"
    )
})

test_that("recompute_exact_inference errors on a non-scalar or empty parm", {
    expect_error(
        recompute_exact_inference(fit_toy2_shift, c("group1", "group1")),
        "single non-empty parameter name"
    )
    expect_error(
        recompute_exact_inference(fit_toy2_shift, ""),
        "single non-empty parameter name"
    )
})

test_that("recompute_exact_inference errors on parm absent from OUTDIST", {
    expect_error(
        recompute_exact_inference(fit_toy2_shift, "not_a_parameter"),
        "not_a_parameter.*not a parameter in fit\\$OUTDIST"
    )
    expect_error(
        recompute_exact_inference(fit_toy2_shift, "not_a_parameter"),
        "group1" # the "Available: ..." part of the message
    )
})

test_that("recompute_exact_inference errors on non-finite/non-scalar beta0", {
    expect_error(
        recompute_exact_inference(fit_toy2_shift, "group1", beta0 = NA_real_),
        "beta0.*finite"
    )
    expect_error(
        recompute_exact_inference(fit_toy2_shift, "group1", beta0 = Inf),
        "beta0.*finite"
    )
    expect_error(
        recompute_exact_inference(fit_toy2_shift, "group1", beta0 = c(0, 1)),
        "beta0.*finite"
    )
})

test_that("recompute_exact_inference errors on alpha outside (0, 1)", {
    expect_error(
        recompute_exact_inference(fit_toy2_shift, "group1", alpha = 0),
        "alpha.*between 0 and 1"
    )
    expect_error(
        recompute_exact_inference(fit_toy2_shift, "group1", alpha = 1),
        "alpha.*between 0 and 1"
    )
    expect_error(
        recompute_exact_inference(fit_toy2_shift, "group1", alpha = -0.05),
        "alpha.*between 0 and 1"
    )
})

test_that("recompute_exact_inference errors on invalid ci_type/direction", {
    expect_error(
        recompute_exact_inference(
            fit_toy2_shift,
            "group1",
            ci_type = "not-a-type"
        ),
        "`ci_type` must be one of"
    )
    expect_error(
        recompute_exact_inference(
            fit_toy2_shift,
            "group1",
            direction = "sideways"
        ),
        "`direction` must be one of"
    )
})

test_that("recompute_exact_inference errors when exactparmest has >1 row", {
    # A hand-built, deliberately malformed fit object; exercises the
    # defensive "expected exactly one exactparmest row" check, which no
    # well-formed fit from the package's own fitting functions can trigger
    # (each parameter appears exactly once in exactparmest by construction).
    fake_fit <- structure(
        list(
            OUTDIST = fit_toy2_shift$OUTDIST,
            exactparmest = rbind(
                fit_toy2_shift$exactparmest,
                fit_toy2_shift$exactparmest
            )
        ),
        class = c("ep_multi_shift", "ep_fit")
    )
    expect_error(
        recompute_exact_inference(fake_fit, "group1"),
        "Expected exactly one exactparmest row"
    )
})

# ------------------------------------------------------------------------
# .validate_recompute_args ----
# ------------------------------------------------------------------------

test_that(".validate_recompute_args passes silently on valid input", {
    expect_silent(
        .validate_recompute_args(
            fit_toy2_shift,
            "group1",
            beta0 = 0,
            alpha = 0.05
        )
    )
})

test_that(".recompute_lookup returns the parameter's dist1d and t_obs", {
    out <- .recompute_lookup(fit_toy2_shift, "group1")

    expect_type(out, "list")
    expect_named(out, c("dist1d", "t_obs"))
    expect_identical(out$dist1d, fit_toy2_shift$OUTDIST[["group1"]])
    expect_equal(
        out$t_obs,
        fit_toy2_shift$exactparmest$SuffStat[
            fit_toy2_shift$exactparmest$Parameter == "group1"
        ]
    )
})

test_that(".validate_recompute_args errors on a fit that isn't an ep_fit", {
    expect_error(
        .validate_recompute_args(list(), "group1", beta0 = 0, alpha = 0.05),
        "ep_multi_shift"
    )
})

test_that(".validate_recompute_args errors on a non-scalar or empty parm", {
    expect_error(
        .validate_recompute_args(
            fit_toy2_shift,
            c("group1", "group1"),
            beta0 = 0,
            alpha = 0.05
        ),
        "single non-empty parameter name"
    )
    expect_error(
        .validate_recompute_args(fit_toy2_shift, "", beta0 = 0, alpha = 0.05),
        "single non-empty parameter name"
    )
})

test_that(".validate_recompute_args errors on a parm absent from OUTDIST", {
    expect_error(
        .validate_recompute_args(
            fit_toy2_shift,
            "not_a_parameter",
            beta0 = 0,
            alpha = 0.05
        ),
        "not_a_parameter.*not a parameter in fit\\$OUTDIST"
    )
})

test_that(".validate_recompute_args errors on a non-finite or non-scalar beta0", {
    expect_error(
        .validate_recompute_args(
            fit_toy2_shift,
            "group1",
            beta0 = NA_real_,
            alpha = 0.05
        ),
        "beta0.*finite"
    )
    expect_error(
        .validate_recompute_args(
            fit_toy2_shift,
            "group1",
            beta0 = Inf,
            alpha = 0.05
        ),
        "beta0.*finite"
    )
    expect_error(
        .validate_recompute_args(
            fit_toy2_shift,
            "group1",
            beta0 = c(0, 1),
            alpha = 0.05
        ),
        "beta0.*finite"
    )
})

test_that(".validate_recompute_args delegates alpha checking to .check_alpha_arg", {
    # Every one of these is the exact same message .check_alpha_arg() itself
    # raises (see test-ep-validate.R); confirms the delegation,
    # not a re-implementation, including the NA case that recompute_exact_
    # inference()'s own pre-refactor inline check did not handle cleanly.
    expect_error(
        .validate_recompute_args(
            fit_toy2_shift,
            "group1",
            beta0 = 0,
            alpha = 0
        ),
        "alpha.*between 0 and 1"
    )
    expect_error(
        .validate_recompute_args(
            fit_toy2_shift,
            "group1",
            beta0 = 0,
            alpha = 1
        ),
        "alpha.*between 0 and 1"
    )
    expect_error(
        .validate_recompute_args(
            fit_toy2_shift,
            "group1",
            beta0 = 0,
            alpha = NA_real_
        ),
        "alpha.*between 0 and 1"
    )
})

test_that(".validate_recompute_args errors when exactparmest has >1 row for parm", {
    fake_fit <- structure(
        list(
            OUTDIST = fit_toy2_shift$OUTDIST,
            exactparmest = rbind(
                fit_toy2_shift$exactparmest,
                fit_toy2_shift$exactparmest
            )
        ),
        class = c("ep_multi_shift", "ep_fit")
    )
    expect_error(
        .validate_recompute_args(fake_fit, "group1", beta0 = 0, alpha = 0.05),
        "Expected exactly one exactparmest row"
    )
})
