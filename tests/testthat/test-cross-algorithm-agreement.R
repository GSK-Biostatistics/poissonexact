# ==============================================================================
# Cross-algorithm agreement tests.
#
# The package's central correctness property is that all three algorithms
# solve the *same* exact conditional inference problem by different
# computational routes, so on any well-posed problem within their shared
# scope, they must agree exactly (up to numerical root-finding
# tolerance): ep_multi_shift() and ep_network()
# always (both handle the fully general categorical-covariates-plus-
# strata model), and ep_closed_form() too whenever a problem
# also happens to fall inside its own narrower special case (exactly one
# two-level exact effect, with strata as the only other structure, no
# other covariates in the mean model).
# ==============================================================================

# ------------------------------------------------------------------------
# test1: single two-level effect + strata only; ----
# all three algorithms are in scope and must all agree.
# ------------------------------------------------------------------------

test_that("all three algorithms agree on test1 with strata and a log-offset", {
    args_common <- list(
        formula = cases ~ group,
        data = test1,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = "age",
        exact_vars = "group",
        verbose = FALSE
    )
    fit_closed <- do.call(
        ep_closed_form,
        args_common[setdiff(names(args_common), "exact_vars")]
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)

    expect_exactparmest_equal(
        fit_closed$exactparmest,
        fit_shift$exactparmest,
        tolerance = 1e-6
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

test_that("all three algorithms agree on test1, no offset, unstratified", {
    args_common <- list(
        formula = cases ~ group,
        data = test1,
        exact_vars = "group",
        verbose = FALSE
    )
    fit_closed <- do.call(
        ep_closed_form,
        args_common[setdiff(names(args_common), "exact_vars")]
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)

    expect_exactparmest_equal(
        fit_closed$exactparmest,
        fit_shift$exactparmest,
        tolerance = 1e-6
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

test_that("multi_shift and network agree on test1 with EXACTMAX=subjects", {
    # ep_closed_form() has no exactmax argument (it always uses
    # SAS rule 2, the stratum total), so only the two general algorithms,
    # which both honour EXACTMAX=var identically, are compared here.
    args_common <- list(
        formula = cases ~ group,
        data = test1,
        strata = "age",
        exact_vars = "group",
        exactmax = "subjects",
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)

    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

test_that("all three algorithms agree on test1: ref_levels, one-sided CI", {
    args_common <- list(
        formula = cases ~ group,
        data = test1,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = "age",
        exact_vars = "group",
        ref_levels = list(group = "1"),
        ci_type = "one-sided",
        verbose = FALSE
    )
    fit_closed <- do.call(
        ep_closed_form,
        args_common[setdiff(names(args_common), "exact_vars")]
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)

    expect_exactparmest_equal(
        fit_closed$exactparmest,
        fit_shift$exactparmest,
        tolerance = 1e-6
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

# ------------------------------------------------------------------------
# test3: same single-effect-plus-strata special case ----
# reproduced by folding every other covariate (age, gender, region) into
# `strata` instead of the model formula, still within ep_closed_form()'s
# scope, so all three algorithms remain comparable. Uses the small,
# deterministic test3 fixture (see helper-data.R) instead of test2, since
# this is a structural agreement check, not a SAS-validated scenario,
# it runs in milliseconds instead of needing the slow-test gate.
# ------------------------------------------------------------------------

test_that("all three algorithms agree on test3, age/gender/region as strata", {
    args_common <- list(
        formula = cases ~ group,
        data = test3,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = c("age", "gender", "region"),
        exact_vars = "group",
        verbose = FALSE
    )
    fit_closed <- do.call(
        ep_closed_form,
        args_common[setdiff(names(args_common), "exact_vars")]
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)

    expect_exactparmest_equal(
        fit_closed$exactparmest,
        fit_shift$exactparmest,
        tolerance = 1e-6
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

# ------------------------------------------------------------------------
# test3: genuine multi-covariate models ----
# (region/gender/age as actual adjustment covariates, not strata),
# outside ep_closed_form()'s scope, so only the two general algorithms are
# compared, across several exact_vars/strata/condition_on combinations.
# Same rationale as above: these check algorithm *agreement*, not a
# specific SAS-matched value, so the small test3 fixture is sufficient
# and much faster than test2.
# ------------------------------------------------------------------------

test_that("multi_shift/network agree on test3, region covariate, age strata", {
    args_common <- list(
        formula = cases ~ group + region,
        data = test3,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = "age",
        exact_vars = "group",
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

test_that("multi_shift and network agree on test3 with multiple exact_vars", {
    args_common <- list(
        formula = cases ~ group + region + gender,
        data = test3,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = "age",
        exact_vars = c("group", "region"),
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

test_that("multi_shift/network agree on test3, condition_on = 'exact'", {
    args_common <- list(
        formula = cases ~ group + region + gender,
        data = test3,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = "age",
        exact_vars = c("group", "region"),
        condition_on = "exact",
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

test_that("multi_shift/network agree on test3, interaction exact effect", {
    args_common <- list(
        formula = cases ~ group * region,
        data = test3,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = c("age", "gender"),
        exact_vars = c("group", "region", "group:region"),
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})

test_that("multi_shift's result does not depend on chunk_size (test3)", {
    args_common <- list(
        formula = cases ~ group + region,
        data = test3,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = "age",
        exact_vars = "group",
        verbose = FALSE
    )
    fit_default <- do.call(ep_multi_shift, args_common)
    fit_chunk1 <- do.call(
        ep_multi_shift,
        c(args_common, list(chunk_size = 1L))
    )
    expect_exactparmest_equal(
        fit_default$exactparmest,
        fit_chunk1$exactparmest,
        tolerance = 1e-8
    )
})

# ------------------------------------------------------------------------
# EXPECTED divergence ----
# ep_closed_form() has no exactmax argument
# and never applies any per-observation cap (always SAS Rule 2, the
# stratum total; see theory.Rmd Sections 5.2, 6.4, 10). The two general
# algorithms, by contrast, apply Rule 3's offset-implied cap automatically
# whenever `offset` is supplied without `exactmax`. Whenever that cap
# actually binds, the closed-form algorithm and the two general algorithms
# are, by design, constraining the reference set differently, and can
# legitimately disagree on the estimate, CI, and p-value. This is not a
# bug; these two tests lock in the expected behaviour on both sides of that
# divide, using the exact fixture that first surfaced the difference.
# ------------------------------------------------------------------------

test_that("closedform legitimately disagrees with multi_shift/network
when the offset-implied cap (Rule 3) binds", {
    binom5 <- data.frame(group = c(1, 0), cases = c(3, 2), fu = c(4, 3))
    args_common <- list(
        formula = cases ~ group,
        data = binom5,
        offset = "fu",
        offset_log = FALSE,
        exact_vars = "group",
        verbose = FALSE
    )
    fit_closed <- do.call(
        ep_closed_form,
        args_common[setdiff(names(args_common), "exact_vars")]
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)

    # Mechanism, directly on the support: the stratum total is
    # M = 3 + 2 = 5. ep_closed_form() places no cap beyond that,
    # so its support for the group-1 count is the full 0:5. The general
    # algorithms cap each observation at ceiling(exp(log(fu))) (Rule 3).
    # For fu = 4 that is unambiguously 4. For fu = 3 it is *also* 4, not
    # the mathematically "clean" 3: in IEEE double precision,
    # exp(log(3)) == 3.0000000000000004, whose ceiling is 4, a
    # floating-point round-trip artifact of si_from_offset(), not a
    # decision the package makes about the data. So the true caps are
    # S_1 = S_0 = 4, which force the group-1 count y1 to satisfy y1 <= 4
    # and 5 - y1 <= 4 (i.e. y1 >= 1) simultaneously, restricting its
    # support to {1, 2, 3, 4} (see theory.Rmd, Section 10.2, for the full
    # derivation and the exact OUTDIST numbers this reproduces).
    expect_identical(as.integer(fit_closed$OUTDIST[["group1"]]$SuffStat), 0:5)
    expect_identical(fit_shift$OUTDIST[["group1"]]$SuffStat, 1:4)
    expect_identical(fit_net$OUTDIST[["group1"]]$SuffStat, 1:4)

    # The two general algorithms honour Rule 3 identically, so they still
    # agree with each other ...
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-8
    )

    # ... but the closed-form algorithm's unrestricted estimate is
    # genuinely different from theirs, not just numerically noisy.
    testthat::expect_gt(
        abs(fit_closed$exactparmest$Estimate - fit_shift$exactparmest$Estimate),
        0.01
    )
})

test_that("closedform agrees with others once Rule 3 cap is non-binding", {
    binom5_wide <- data.frame(group = c(1, 0), cases = c(3, 2), fu = c(5, 5))
    args_common <- list(
        formula = cases ~ group,
        data = binom5_wide,
        offset = "fu",
        offset_log = FALSE,
        exact_vars = "group",
        verbose = FALSE
    )
    fit_closed <- do.call(
        ep_closed_form,
        args_common[setdiff(names(args_common), "exact_vars")]
    )
    fit_shift <- do.call(ep_multi_shift, args_common)
    fit_net <- do.call(ep_network, args_common)

    # Both caps (S_1 = S_0 = 5) now exceed the stratum total (5), so Rule 3
    # excludes nothing beyond what Rule 2 already implies: every algorithm
    # builds the same, full 0:5 support, and all three agree exactly.
    expect_identical(as.integer(fit_closed$OUTDIST[["group1"]]$SuffStat), 0:5)
    expect_identical(fit_shift$OUTDIST[["group1"]]$SuffStat, 0:5)

    expect_exactparmest_equal(
        fit_closed$exactparmest,
        fit_shift$exactparmest,
        tolerance = 1e-6
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
})
