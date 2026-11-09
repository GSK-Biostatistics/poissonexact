# ==============================================================================
# SAS PROC GENMOD (EXACT statement)-validated scenarios  ----
# ==============================================================================

skip_if_no_sas <- function(sas_expected) {
    testthat::skip_if(
        all(vapply(sas_expected, function(v) all(is.na(v)), logical(1))),
        "SAS reference values not yet supplied for this scenario"
    )
}

expect_matches_sas <- function(row, sas_expected, tolerance = 0.0001) {
    skip_if_no_sas(sas_expected)
    for (nm in names(sas_expected)) {
        if (!all(is.na(sas_expected[[nm]]))) {
            # nolint start: object_usage_linter
            expect_equal(row[[nm]], sas_expected[[nm]], tolerance = tolerance)
            # nolint end: object_usage_linter
        }
    }
}

# A completely empty SAS reference block; fill in any subset of these
# fields (Estimate, StdError, PValue, LowerCL, UpperCL, SuffStat) once
# real PROC GENMOD EXACT output is available for a scenario.
sas_placeholder <- function(
    estimate = NA,
    stderror = NA,
    pvalue = NA,
    lowercl = NA,
    uppercl = NA
) {
    list(
        Estimate = estimate,
        StdError = stderror,
        PValue = pvalue,
        LowerCL = lowercl,
        UpperCL = uppercl
    )
}

# ------------------------------------------------------------------------
# test1-based scenarios (fit1, fit1_1, fit12, fit13, fit1_3, fit14)  ----
#
# Every one of these tests exact_vars = "group" on test1, so, regardless
# of strata/offset/formula differences between scenarios, the observed
# sufficient statistic for "group1" is always the same self-computable
# quantity: the total cases among the group == 1 rows.
# ------------------------------------------------------------------------

t1_group1_suffstat <- sum(test1$cases[test1$group == 1L]) # = 3 + 2 + 1 = 6

test_that("fit1: 3 algorithms agree on test1 with strata + offset", {
    args <- list(
        formula = cases ~ group,
        data = test1,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = "age",
        exact_vars = "group",
        ref_levels = list(group = "0", age = "50"),
        alpha = 0.05,
        ci_type = "two-sided",
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args)
    fit_net <- do.call(ep_network, args)
    fit_closed <- do.call(
        ep_closed_form,
        args[setdiff(names(args), "exact_vars")]
    )

    expect_identical(fit_shift$exactparmest$Parameter, "group1")
    expect_equal(fit_shift$exactparmest$SuffStat, t1_group1_suffstat)
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
    expect_equal(fit_closed$exactparmest$SuffStat, t1_group1_suffstat)
    expect_exactparmest_equal(
        fit_closed$exactparmest,
        fit_shift$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -1.2691,
        stderror = 0.4606,
        pvalue = 0.0048,
        lowercl = -2.3726,
        uppercl = -0.3357
    )
    fit_shift$exactparmest[1, -1] <- round(
        fit_shift$exactparmest[1, -1],
        digits = 4
    )
    expect_matches_sas(
        fit_shift$exactparmest[1, ],
        sas_result,
        tolerance = 1e-4
    )
})

test_that("fit12: all 3 algorithms agree on test1 with offset but no strata", {
    args <- list(
        formula = cases ~ group,
        data = test1,
        offset = "lnfutime",
        offset_log = TRUE,
        exact_vars = "group",
        ref_levels = list(group = "0", age = "50"),
        alpha = 0.05,
        ci_type = "two-sided",
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args)
    fit_net <- do.call(ep_network, args)
    fit_closed <- do.call(
        ep_closed_form,
        args[setdiff(names(args), "exact_vars")]
    )

    for (fit in list(fit_shift, fit_net, fit_closed)) {
        expect_identical(fit$exactparmest$Parameter, "group1")
        expect_equal(fit$exactparmest$SuffStat, t1_group1_suffstat)
    }
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_closed$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -1.2665,
        stderror = 0.4606,
        pvalue = 0.0049,
        lowercl = -2.3700,
        uppercl = -0.3331
    )
    fit_closed$exactparmest[1, -1] <- round(
        fit_closed$exactparmest[1, -1],
        digits = 4
    )
    expect_matches_sas(fit_closed$exactparmest[1, ], sas_result)
})

test_that("fit13: ep_multi_shift/ep_network agree on test1, no offset/strata", {
    args <- list(
        formula = cases ~ group,
        data = test1,
        exact_vars = "group",
        ref_levels = list(group = "0", age = "50"),
        alpha = 0.05,
        ci_type = "two-sided",
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args)
    fit_net <- do.call(ep_network, args)
    fit_closed <- do.call(
        ep_closed_form,
        args[setdiff(names(args), "exact_vars")]
    )

    expect_equal(fit_shift$exactparmest$SuffStat, t1_group1_suffstat)
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )

    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_closed$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -1.2993,
        stderror = 0.4606,
        pvalue = 0.0037,
        lowercl = -2.4028,
        uppercl = -0.3659
    )
    fit_shift$exactparmest[1, -1] <- round(
        fit_shift$exactparmest[1, -1],
        digits = 4
    )
    expect_matches_sas(fit_shift$exactparmest[1, ], sas_result)
})

test_that("fit14: ep_multi_shift/ep_network agree with age as cov, no strata", {
    args <- list(
        formula = cases ~ group + age,
        data = test1,
        exact_vars = "group",
        ref_levels = list(group = "0", age = "50"),
        alpha = 0.05,
        ci_type = "two-sided",
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, args)
    fit_net <- do.call(ep_network, args)

    expect_equal(fit_shift$exactparmest$SuffStat, t1_group1_suffstat)
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -1.2993,
        stderror = 0.4606,
        pvalue = 0.0037,
        lowercl = -2.4028,
        uppercl = -0.3659
    )
    fit_shift$exactparmest[1, -1] <- round(
        fit_shift$exactparmest[1, -1],
        digits = 4
    )
    expect_matches_sas(fit_shift$exactparmest[1, ], sas_result)
})

# ------------------------------------------------------------------------
# test2-based scenarios (fit21-fit25) ----
#
# Every parameter's SuffStat is likewise self-computable directly from
# test2 (defined with a fixed seed in helper-data.R), independent of
# which algorithm computed it: it is simply the total cases among the
# rows at that dummy column's non-reference level.
# ------------------------------------------------------------------------

t2_suffstat_for <- function(var, level) {
    sum(test2$cases[as.character(test2[[var]]) == level]) # nolint: object_usage_linter
}

test_that("fit21: 3 algos agree on test2 with age/gender/region as strata", {
    common <- list(
        formula = cases ~ group,
        data = test2,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = c("age", "gender", "region"),
        exact_vars = "group",
        ref_levels = list(
            group = "0",
            gender = "0",
            region = "1",
            age = "50"
        ),
        alpha = 0.05,
        verbose = FALSE
    )
    fit_shift <- do.call(
        ep_multi_shift,
        c(common, list(condition_on = "model"))
    )
    fit_net <- do.call(
        ep_network,
        c(common, list(condition_on = "model"))
    )
    fit_closed <- do.call(
        ep_closed_form,
        common[setdiff(names(common), "exact_vars")]
    )

    expected_suffstat <- t2_suffstat_for("group", "1")
    for (fit in list(fit_shift, fit_net, fit_closed)) {
        expect_identical(fit$exactparmest$Parameter, "group1")
        expect_equal(fit$exactparmest$SuffStat, expected_suffstat)
    }
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_closed$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -0.3959,
        stderror = 0.3354,
        pvalue = 0.3062,
        lowercl = -1.1259,
        uppercl = 0.3071
    )
    fit_shift$exactparmest[1, -1] <- round(
        fit_shift$exactparmest[1, -1],
        digits = 4
    )
    expect_matches_sas(fit_shift$exactparmest[1, ], sas_result)
})

test_that("fit22: ep_multi_shift/ep_network agree on test2 (gender as cov)", {
    common <- list(
        formula = cases ~ group + gender,
        data = test2,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = c("age", "region"),
        exact_vars = "group",
        ref_levels = list(
            group = "0",
            gender = "0",
            region = "1",
            age = "50"
        ),
        alpha = 0.05,
        ci_type = "two-sided",
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, common)
    fit_net <- do.call(ep_network, common)

    expect_equal(
        fit_shift$exactparmest$SuffStat,
        t2_suffstat_for("group", "1")
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -0.4002,
        stderror = 0.3352,
        pvalue = 0.2995,
        lowercl = -1.1297,
        uppercl = 0.3022
    )
    fit_shift$exactparmest[1, -1] <- round(
        fit_shift$exactparmest[1, -1],
        digits = 4
    )
    expect_matches_sas(fit_shift$exactparmest[1, ], sas_result)
})

test_that("fit23: ep_multi_shift/ep_network agree testing group AND gender", {
    common <- list(
        formula = cases ~ group + gender,
        data = test2,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = c("age", "region"),
        exact_vars = c("group", "gender"),
        ref_levels = list(
            group = "0",
            gender = "0",
            region = "1",
            age = "50"
        ),
        alpha = 0.05,
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, common)
    fit_net <- do.call(ep_network, common)

    for (fit in list(fit_shift, fit_net)) {
        expect_identical(fit$exactparmest$Parameter, c("group1", "gender1"))
        expect_equal(
            fit$exactparmest$SuffStat,
            c(t2_suffstat_for("group", "1"), t2_suffstat_for("gender", "1"))
        )
    }
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -0.4002,
        stderror = 0.3352,
        pvalue = 0.2995,
        lowercl = -1.1297,
        uppercl = 0.3022
    )
    fit_shift$exactparmest[1, -1] <- round(
        fit_shift$exactparmest[1, -1],
        digits = 4
    )

    expect_matches_sas(
        fit_shift$exactparmest[fit_shift$exactparmest$Parameter == "group1", ],
        sas_result
    )

    sas_result <- sas_placeholder(
        estimate = -0.2466,
        stderror = 0.3331,
        pvalue = 0.5648,
        lowercl = -0.9671,
        uppercl = 0.4550
    )
    fit_shift$exactparmest[2, -1] <- round(
        fit_shift$exactparmest[2, -1],
        digits = 4
    )
    expect_matches_sas(
        fit_shift$exactparmest[fit_shift$exactparmest$Parameter == "gender1", ],
        sas_result
    )
})

test_that("fit24: same as fit23 but with no offset; 2 algs agree", {
    common <- list(
        formula = cases ~ group + gender,
        data = test2,
        strata = c("age", "region"),
        exact_vars = c("group", "gender"),
        ref_levels = list(group = "0", gender = "0", region = "1", age = "50"),
        alpha = 0.05,
        verbose = FALSE
    )
    fit_shift <- do.call(ep_multi_shift, common)
    fit_net <- do.call(ep_network, common)

    expect_equal(
        fit_shift$exactparmest$SuffStat,
        c(t2_suffstat_for("group", "1"), t2_suffstat_for("gender", "1"))
    )
    expect_exactparmest_equal(
        fit_shift$exactparmest,
        fit_net$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -0.3830,
        stderror = 0.3348,
        pvalue = 0.3240,
        lowercl = -1.1118,
        uppercl = 0.3187
    )
    fit_shift$exactparmest[1, -1] <- round(
        fit_shift$exactparmest[1, -1],
        digits = 4
    )
    expect_matches_sas(
        fit_shift$exactparmest[fit_shift$exactparmest$Parameter == "group1", ],
        sas_result
    )

    sas_result <- sas_placeholder(
        estimate = -0.2719,
        stderror = 0.3318,
        pvalue = 0.5114,
        lowercl = -0.9897,
        uppercl = 0.4268
    )
    fit_shift$exactparmest[2, -1] <- round(
        fit_shift$exactparmest[2, -1],
        digits = 4
    )
    expect_matches_sas(
        fit_shift$exactparmest[fit_shift$exactparmest$Parameter == "gender1", ],
        sas_result
    )
})

test_that("fit25: ep_multi_shift/ep_network agree testing group and region", {
    common <- list(
        formula = cases ~ group + region,
        data = test2,
        offset = "lnfutime",
        offset_log = TRUE,
        strata = c("age", "gender"),
        exact_vars = c("group", "region"),
        ref_levels = list(
            group = "0",
            gender = "0",
            region = "1",
            age = "50"
        ),
        alpha = 0.05,
        verbose = FALSE
    )
    fit_net <- do.call(ep_network, common)
    fit_shift <- do.call(ep_multi_shift, common)

    expected_parms <- c("group1", "region2", "region3", "region4")
    expected_suffstat <- c(
        t2_suffstat_for("group", "1"),
        t2_suffstat_for("region", "2"),
        t2_suffstat_for("region", "3"),
        t2_suffstat_for("region", "4")
    )
    for (fit in list(fit_net, fit_shift)) {
        expect_identical(fit$exactparmest$Parameter, expected_parms)
        expect_equal(fit$exactparmest$SuffStat, expected_suffstat)
    }
    expect_exactparmest_equal(
        fit_net$exactparmest,
        fit_shift$exactparmest,
        tolerance = 1e-6
    )

    sas_result <- sas_placeholder(
        estimate = -0.3989,
        stderror = 0.3354,
        pvalue = 0.3017,
        lowercl = -1.1288,
        uppercl = 0.3040
    )
    fit_shift$exactparmest[1, -1] <- round(
        fit_shift$exactparmest[1, -1],
        digits = 4
    )
    expect_matches_sas(
        fit_shift$exactparmest[fit_shift$exactparmest$Parameter == "group1", ],
        sas_result
    )

    sas_result <- sas_placeholder(
        estimate = 0.3262,
        stderror = 0.4213,
        pvalue = 0.5696,
        lowercl = -0.5791,
        uppercl = 1.2629
    )
    fit_shift$exactparmest[2, -1] <- round(
        fit_shift$exactparmest[2, -1],
        digits = 4
    )
    expect_matches_sas(
        fit_shift$exactparmest[fit_shift$exactparmest$Parameter == "region2", ],
        sas_result
    )

    sas_result <- sas_placeholder(
        estimate = -0.4974,
        stderror = 0.5169,
        pvalue = 0.4719,
        lowercl = -1.7068,
        uppercl = 0.6147
    )
    fit_shift$exactparmest[3, -1] <- round(
        fit_shift$exactparmest[3, -1],
        digits = 4
    )
    expect_matches_sas(
        fit_shift$exactparmest[fit_shift$exactparmest$Parameter == "region3", ],
        sas_result
    )

    sas_result <- sas_placeholder(
        estimate = -0.1798,
        stderror = 0.4750,
        pvalue = 0.8883,
        lowercl = -1.2512,
        uppercl = 0.8563
    )
    fit_shift$exactparmest[4, -1] <- round(
        fit_shift$exactparmest[4, -1],
        digits = 4
    )
    expect_matches_sas(
        fit_shift$exactparmest[fit_shift$exactparmest$Parameter == "region4", ],
        sas_result
    )
})
