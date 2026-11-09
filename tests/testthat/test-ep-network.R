# ==============================================================================
# Tests for ep-network.R; the network algorithm (Mehta,
# Patel, Senchaudhuri 1992), generalized to Poisson counts: closed-form
# LP/SP backward induction, the continuous forward pass with mid-stratum
# pruning, and the exported ep_network(). One of the three
# exported main functions, so tested extensively: the internal building
# blocks with fully hand-traced examples (including a 2-step forward-pass
# trace that hand-verifies the LP/SP pruning logic against a direct
# combinatorial argument), then the end-to-end function anchored on the
# same toy2/binom5/multi_df examples used in
# test-ep-multivariate-shift.R, explicitly cross-checked for
# agreement with both ep_closed_form() and
# ep_multi_shift() wherever their special cases overlap.
# ==============================================================================

# ------------------------------------------------------------------------
# .stratum_suffix_bounds ----
# ------------------------------------------------------------------------

test_that(".stratum_suffix_bounds matches the hand-derived example", {
    xc_s <- matrix(c(1L, 0L, 1L), ncol = 1)
    si_s <- c(2L, 3L, 4L)
    out <- .stratum_suffix_bounds(xc_s, si_s, p_cond = 1)
    expect_equal(out$suffix1[, 1], c(6, 4, 4, 0))
    expect_equal(out$suffix0[, 1], c(3, 3, 0, 0))
})

test_that(".stratum_suffix_bounds returns all-zero bounds when p_cond = 0", {
    out <- .stratum_suffix_bounds(NULL, si_s = c(2L, 3L), p_cond = 0)
    expect_identical(dim(out$suffix1), c(3L, 0L))
    expect_identical(dim(out$suffix0), c(3L, 0L))
})

test_that(".stratum_suffix_bounds handles multiple nuisance columns", {
    xc_s <- matrix(c(1L, 0L, 0L, 1L), ncol = 2) # col1 = (1,0), col2 = (0,1)
    si_s <- c(3L, 5L)
    out <- .stratum_suffix_bounds(xc_s, si_s, p_cond = 2)
    # col 1: active = (T,F) -> contrib1=(3,0),contrib0=(0,5)
    expect_equal(out$suffix1[, 1], c(3, 0, 0))
    expect_equal(out$suffix0[, 1], c(5, 5, 0))
    # col 2: active = (F,T) -> contrib1=(0,5),contrib0=(3,0)
    expect_equal(out$suffix1[, 2], c(5, 5, 0))
    expect_equal(out$suffix0[, 2], c(3, 0, 0))
})

# ------------------------------------------------------------------------
# .order_strata_by_size ----
# ------------------------------------------------------------------------

# ordering is by ascending observed total
test_that(".order_strata_by_size orders strata by total, stable on ties", {
    slevels0 <- c("A", "B", "C")
    sid <- c("A", "A", "B", "C", "C", "C")
    y <- c(1, 2, 5, 1, 1, 1) # totals: A = 3, B = 5, C = 3
    out <- .order_strata_by_size(slevels0, sid, y)
    # A and C tie at 3; stable order() keeps A (appears first in slevels0)
    # ahead of C.
    expect_identical(out, c("A", "C", "B"))
})

test_that(".order_strata_by_size leaves a single stratum unchanged", {
    expect_identical(.order_strata_by_size("A", c("A", "A"), c(1, 2)), "A")
})

# ------------------------------------------------------------------------
# .compute_stratum_min_add ----
# ------------------------------------------------------------------------

test_that(".compute_stratum_min_add matches the hand-derived example", {
    y_s <- c(3, 4)
    xc_s <- matrix(c(1, 0), ncol = 1)
    si_s <- c(2, 5)
    # idx0 = {2}, cap0 = si_s[2] = 5, out = max(0, 7 - 5) = 2.
    expect_identical(.compute_stratum_min_add(y_s, xc_s, si_s), 2L)
})

# capacity is the inactive-observation capacity
test_that(".compute_stratum_min_add returns 0 when capacity covers the total", {
    y_s <- c(3, 4)
    xc_s <- matrix(c(1, 0), ncol = 1)
    si_s <- c(2, 10)
    expect_identical(.compute_stratum_min_add(y_s, xc_s, si_s), 0L)
})

test_that(".compute_stratum_min_add handles multiple nuisance columns", {
    y_s <- c(3, 4)
    xc_s <- matrix(c(1, 0, 0, 1), ncol = 2)
    si_s <- c(2, 5)
    # col1: idx0={2}, cap0=5, out=max(0,7-5)=2.
    # col2: idx0={1}, cap0=2, out=max(0,7-2)=5.
    expect_identical(.compute_stratum_min_add(y_s, xc_s, si_s), c(2L, 5L))
})

# ------------------------------------------------------------------------
# .network_stratum_slice ----
# ------------------------------------------------------------------------

test_that(".network_stratum_slice extends slice with min_add_s, suffix bound", {
    model <- list(
        sid = c("A", "A"),
        y = c(3L, 4L),
        # col1 = interest, col2 = nuisance (1,0)
        x = matrix(c(1, 0, 1, 0), ncol = 2),
        log_off = c(0, 0),
        use_rule1 = FALSE,
        use_rule2 = TRUE,
        use_rule3 = FALSE,
        exactmax_vec = NULL
    )
    info <- .network_stratum_slice(
        model,
        "A",
        jcol = 1,
        cond_cols_idx = 2,
        p_cond = 1,
        verbose = FALSE
    )

    # Base .stratum_slice() fields still present.
    expect_identical(info$y_s, c(3L, 4L))
    expect_identical(info$m_s_obs, 7L)

    # min_add_s: xc_s = (1,0), si_s (rule 2) = rep(7,2); idx0={2}, cap0=7,
    # and out = max(0, 7-7) = 0.
    expect_identical(info$min_add_s, 0L)

    # suffix bounds for xc_s=(1,0), si_s=(7,7): active=(T,F); contrib1=(7,0),
    # contrib0=(0,7); suffix1=(7,0,0); suffix0=(7,7,0).
    expect_equal(info$suffix1[, 1], c(7, 0, 0))
    expect_equal(info$suffix0[, 1], c(7, 7, 0))
})

# skipped to avoid needless computation cost
test_that(".network_stratum_slice skips min_add_s/bounds when p_cond = 0", {
    model <- list(
        sid = c("A", "A"),
        y = c(1L, 1L),
        x = matrix(c(1, 0), ncol = 1),
        log_off = c(0, 0),
        use_rule1 = FALSE,
        use_rule2 = TRUE,
        use_rule3 = FALSE,
        exactmax_vec = NULL
    )
    info <- .network_stratum_slice(
        model,
        "A",
        jcol = 1,
        cond_cols_idx = integer(0),
        p_cond = 0,
        verbose = FALSE
    )
    expect_identical(info$min_add_s, integer(0))
    expect_identical(dim(info$suffix1), c(3L, 0L))
})

# ------------------------------------------------------------------------
# .network_expand_observation / .network_finish_stratum ----
#
# A 2-observation, 1-stratum, hand trace with a nuisance column, verified
# against the direct combinatorial argument: y = (1, 1) (ms_here = 2),
# xj_s = (0, 0) (t_j irrelevant here), xc_s = (1, 0) (observation 1 is
# "active" for the nuisance coordinate, observation 2 is not),
# si_s = (2, 2), tn_obs = 1, no future strata (rem_future = min_future =
# 0). Since observation 2 can never affect the nuisance coordinate, the
# only way to hit cond_1 = tn_obs = 1 exactly is y1 = 1 (and then, to also
# hit ms_here = 2, y2 = 1); every other allocation is pruned by the LP/SP
# bounds strictly before the end of the stratum, exactly the single
# feasible table (y1=1, y2=1), weight 1/(1!1!) = 1.
# ------------------------------------------------------------------------

# result matches a hand-derived trace
test_that(".network_expand_observation prunes to the feasible state at r = 1", {
    info <- list(
        xj_s = c(0L, 0L),
        off_s = c(0, 0),
        si_s = c(2L, 2L),
        xc_s = matrix(c(1L, 0L), ncol = 1),
        suffix1 = matrix(c(2, 0, 0), ncol = 1),
        suffix0 = matrix(c(2, 2, 0), ncol = 1)
    )
    dt0 <- list(m_j = 0L, t_j = 0L, cond_1 = 0L, lw = 0.0)

    dt1 <- .network_expand_observation(
        dt0,
        r = 1,
        info,
        p_cond = 1,
        cond_cols = "cond_1",
        tn_obs = 1L,
        ms_here = 2L,
        rem_future = 0L,
        min_future = 0L
    )
    expect_identical(.tbl_nrow(dt1), 1L)
    expect_identical(dt1$m_j, 1L)
    expect_identical(dt1$cond_1, 1L)
    expect_equal(dt1$lw, 0)
})

test_that(".network_expand_observation continues the trace to r = 2", {
    info <- list(
        xj_s = c(0L, 0L),
        off_s = c(0, 0),
        si_s = c(2L, 2L),
        xc_s = matrix(c(1L, 0L), ncol = 1),
        suffix1 = matrix(c(2, 0, 0), ncol = 1),
        suffix0 = matrix(c(2, 2, 0), ncol = 1)
    )
    dt_after_r1 <- list(
        m_j = 1L,
        t_j = 0L,
        cond_1 = 1L,
        lw = 0.0
    )

    dt2 <- .network_expand_observation(
        dt_after_r1,
        r = 2,
        info,
        p_cond = 1,
        cond_cols = "cond_1",
        tn_obs = 1L,
        ms_here = 2L,
        rem_future = 0L,
        min_future = 0L
    )
    expect_identical(.tbl_nrow(dt2), 1L)
    expect_identical(dt2$m_j, 2L)
    expect_identical(dt2$cond_1, 1L)
    expect_equal(dt2$lw, 0)

    finished <- .network_finish_stratum(dt2, ms_here = 2L, cond_cols = "cond_1")
    expect_identical(finished$t_j, 0L)
    expect_identical(finished$cond_1, 1L)
    expect_equal(finished$lw, 0)
})

test_that(".network_finish_stratum filters to total, resets m_j, merges keys", {
    dt <- list(
        m_j = c(2L, 2L, 2L),
        t_j = c(0L, 0L, 1L),
        lw = log(c(1, 2, 4))
    )
    out <- .network_finish_stratum(dt, ms_here = 2L, cond_cols = character(0))
    out <- order_tbl(out, "t_j")
    expect_identical(out$t_j, c(0L, 1L))
    expect_equal(exp(out$lw), c(3, 4), tolerance = 1e-10)
    expect_true(all(out$m_j == 0L))
})

test_that(".network_finish_stratum returns zero rows when nothing matches", {
    dt <- list(
        m_j = c(1L, 3L),
        t_j = c(0L, 0L),
        lw = log(c(1, 1))
    )
    out <- .network_finish_stratum(dt, ms_here = 2L, cond_cols = character(0))
    expect_identical(.tbl_nrow(out), 0L)
})

# ------------------------------------------------------------------------
# ep_network(): end-to-end, anchored on toy2 ----
# ------------------------------------------------------------------------

test_that("ep_network reproduces every hand-derived toy2 quantity", {
    fit <- ep_network(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        alpha = 0.10,
        verbose = FALSE
    )
    expect_s3_class(fit, "ep_network")
    expect_s3_class(fit, "ep_multi_shift")
    row <- fit$exactparmest
    expect_identical(row$Parameter, "group1")
    expect_identical(row$SuffStat, 1L)
    expect_equal(
        as.list(row[c("Estimate", "StdError", "PValue")]),
        list(Estimate = 0, StdError = sqrt(2), PValue = 1),
        tolerance = 1e-6
    )
    expect_equal(row$LowerCL, -row$UpperCL, tolerance = 1e-8)

    od <- fit$OUTDIST[["group1"]]
    expect_identical(od$SuffStat, 0:2)
    expect_equal(od$prob, c(0.25, 0.5, 0.25), tolerance = 1e-10)
})

# ------------------------------------------------------------------------
# ep_network(): the multi_df combinatorial anchor directly ----
# ------------------------------------------------------------------------

# anchor is combinatorial
test_that("ep_network matches the hand-derived anchor, nuisance covariate", {
    multi_df <- data.frame(
        group = c(1L, 0L, 1L, 0L),
        region = c(1L, 1L, 0L, 0L),
        cases = c(1L, 0L, 0L, 1L)
    )
    fit <- ep_network(
        cases ~ group + region,
        data = multi_df,
        exact_vars = "group",
        verbose = FALSE
    )
    row <- fit$exactparmest
    expect_identical(row$SuffStat, 1L)
    expect_equal(row$Estimate, 0, tolerance = 1e-8)
    expect_equal(row$PValue, 1, tolerance = 1e-8)

    od <- fit$OUTDIST[["group1"]]
    expect_identical(od$SuffStat, 0:2)
    expect_equal(od$prob, c(0.25, 0.5, 0.25), tolerance = 1e-10)
})

# ------------------------------------------------------------------------
# ep_network(): condition_on, EXACTMAX, verbose, and errors ----
# ------------------------------------------------------------------------

# results are multi-parameter
test_that("ep_network's condition_on = 'exact'/'model' both return results", {
    multi_df <- data.frame(
        group = c(1L, 0L, 1L, 0L),
        region = c(1L, 1L, 0L, 0L),
        cases = c(1L, 0L, 0L, 1L)
    )
    fit_exact <- ep_network(
        cases ~ group + region,
        data = multi_df,
        exact_vars = c("group", "region"),
        condition_on = "exact",
        verbose = FALSE
    )
    fit_model <- ep_network(
        cases ~ group + region,
        data = multi_df,
        exact_vars = c("group", "region"),
        condition_on = "model",
        verbose = FALSE
    )
    for (fit in list(fit_exact, fit_model)) {
        expect_identical(fit$exactparmest$Parameter, c("group1", "region1"))
        expect_true(
            all(
                fit$exactparmest$PValue >= 0 & fit$exactparmest$PValue <= 1
            )
        )
        expect_identical(fit$exactparmest$SuffStat, c(1L, 1L))
    }
})

# ------------------------------------------------------------------------
# ep_network(): an actual interaction term as the exact effect ----
# hand-derived
#
# Same derivation as test-ep-multivariate-shift.R's own
# "group:region" anchor: testing exact_vars = "group:region" on multi_df
# (formula cases ~ group * region) conditions on both main effects
# (p_cond = 2), and re-enumerating the 10 configurations filtered to the
# observed (group total = 1, region total = 1) leaves exactly two
# configurations, giving log_ck = log(c(1, 1)) for the interaction's
# conditional distribution and an observed t_obs = 1 sitting at the
# upper boundary of that 2-state support (CMLE = +Inf, upper CI = Inf by
# definition).
# ------------------------------------------------------------------------

test_that("ep_network matches the interaction anchor (group:region)", {
    multi_df <- data.frame(
        group = c(1L, 0L, 1L, 0L),
        region = c(1L, 1L, 0L, 0L),
        cases = c(1L, 0L, 0L, 1L)
    )
    fit <- ep_network(
        cases ~ group * region,
        data = multi_df,
        exact_vars = "group:region",
        verbose = FALSE
    )
    row <- fit$exactparmest
    expect_identical(row$Parameter, "group1:region1")
    expect_identical(row$SuffStat, 1L)
    expect_identical(row$Estimate, Inf)
    expect_true(is.na(row$StdError))
    expect_equal(row$PValue, 1, tolerance = 1e-8)
    expect_identical(row$UpperCL, Inf)
    expect_true(is.finite(row$LowerCL))

    od <- fit$OUTDIST[["group1:region1"]]
    expect_identical(od$SuffStat, 0:1)
    expect_equal(od$prob, c(0.5, 0.5), tolerance = 1e-10)
})

test_that("ep_network agrees with ep_multi_shift on the interaction anchor", {
    multi_df <- data.frame(
        group = c(1L, 0L, 1L, 0L),
        region = c(1L, 1L, 0L, 0L),
        cases = c(1L, 0L, 0L, 1L)
    )
    fit_net <- ep_network(
        cases ~ group * region,
        data = multi_df,
        exact_vars = "group:region",
        verbose = FALSE
    )
    fit_shift <- ep_multi_shift(
        cases ~ group * region,
        data = multi_df,
        exact_vars = "group:region",
        verbose = FALSE
    )
    expect_exactparmest_equal(
        fit_net$exactparmest,
        fit_shift$exactparmest,
        tolerance = 1e-8
    )
})

test_that("ep_network never implicitly tests interaction, main effect only", {
    multi_df <- data.frame(
        group = c(1L, 0L, 1L, 0L),
        region = c(1L, 1L, 0L, 0L),
        cases = c(1L, 0L, 0L, 1L)
    )
    fit <- ep_network(
        cases ~ group * region,
        data = multi_df,
        exact_vars = "group",
        verbose = FALSE
    )
    expect_identical(fit$exactparmest$Parameter, "group1")
})

test_that("ep_network honours EXACTMAX=var (rule 1)", {
    fit <- ep_network(
        cases ~ group,
        data = test1,
        strata = "age",
        exact_vars = "group",
        exactmax = "subjects",
        verbose = FALSE
    )
    expect_identical(fit$exactparmest$Parameter, "group1")
    expect_true(is.finite(fit$exactparmest$PValue))
})

test_that("ep_network's verbose = TRUE prints per-param/stratum progress", {
    expect_snapshot(
        ep_network(
            cases ~ group,
            data = toy2,
            strata = "strat",
            exact_vars = "group",
            verbose = TRUE
        )
    )
})

test_that("ep_network propagates argument-validation errors", {
    expect_error(
        ep_network(
            cases ~ group,
            data = toy2,
            strata = "strat",
            exact_vars = character(0),
            verbose = FALSE
        ),
        "non-empty character vector"
    )
})

test_that("ep_network has no chunk_size argument", {
    expect_false("chunk_size" %in% names(formals(ep_network)))
})

# unexpected/defensive case
test_that("ep_network errors if a mid-stratum step empties the table", {
    # Defensive-only, same reasoning as the analogous multi_shift guards:
    # the LP/SP suffix bounds used by .network_expand_observation() are
    # valid (not merely approximate) bounds, and the true observed
    # allocation always satisfies them, so this mid-stratum "all states
    # pruned" check can never fire on real, valid input.
    testthat::local_mocked_bindings(
        .network_expand_observation = function(...) {
            list(
                m_j = integer(0),
                t_j = integer(0),
                cond_1 = integer(0),
                lw = numeric(0)
            )
        }
    )
    expect_error(
        ep_network(
            cases ~ group,
            data = toy2,
            strata = "strat",
            exact_vars = "group",
            verbose = FALSE
        ),
        "All states pruned mid-stratum"
    )
})

# unexpected/defensive case
test_that("ep_network errors if finishing a stratum empties the table", {
    # Same reasoning as above, for the guard that fires immediately after
    # a stratum's forward pass completes (rather than mid-stratum). Forced
    # here by mocking .network_finish_stratum() to return a zero-row
    # table, so the (real, correctly-computed) per-observation steps
    # leading up to it are unaffected.
    testthat::local_mocked_bindings(
        .network_finish_stratum = function(...) {
            list(
                m_j = integer(0),
                t_j = integer(0),
                cond_1 = integer(0),
                lw = numeric(0)
            )
        }
    )
    expect_error(
        ep_network(
            cases ~ group,
            data = toy2,
            strata = "strat",
            exact_vars = "group",
            verbose = FALSE
        ),
        "All states pruned after finishing stratum"
    )
})


test_that("ep_network's print method identifies the network algorithm", {
    fit <- ep_network(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        verbose = FALSE
    )
    expect_snapshot(print(fit))
})

# rather than silently returning a bogus fit
test_that("ep_network errors when a count exceeds its SAS-rule-3 cap", {
    # Stratum 2's control-arm row has futime = 1, so under SAS rule 3
    # (OFFSET supplied, no EXACTMAX) its cap is S_i = max(ceil(1), 1) = 1,
    # but the observed count for that row is 3, which violates 0 <= y_i <=
    # S_i. ep_multi_shift() already caught this (inside
    # .joint_dist_stratum_dp_vec()); ep_network() previously did
    # not, and returned a fit computed from an infeasible reference set
    # instead of raising an error.
    bad <- data.frame(
        group = c(1, 0, 1, 0),
        strat = c(1, 1, 2, 2),
        cases = c(2, 1, 0, 3),
        fu = c(5, 5, 5, 1)
    )
    expect_error(
        ep_network(
            cases ~ group,
            data = bad,
            offset = "fu",
            offset_log = FALSE,
            strata = "strat",
            exact_vars = "group",
            verbose = FALSE
        ),
        "exceed the per-observation cap"
    )
})

# ------------------------------------------------------------------------
# .network_one_parameter() directly ----
# the per-parameter worker that ep_network() now delegates its loop body to.
# As with .multishift_one_parameter() (see test-ep-multivariate-shift.R),
# its logic is already exercised indirectly through every
# ep_network() test above, including the verbose-output and
# mock-based pruning-error tests. These tests call the helper directly
# against a model built by hand with .setup_exact_model() and a stratum
# order built by hand with .order_strata_by_size(), to pin down its own
# contract independently of how the wrapper invokes it.
# ------------------------------------------------------------------------

test_that(".network_one_parameter agrees with ep_network (p_cond = 0)", {
    model <- .setup_exact_model(
        cases ~ group,
        data = toy2,
        offset = NULL,
        offset_log = TRUE,
        exactmax = NULL,
        strata = "strat",
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    slevels_ord <- .order_strata_by_size(model$slevels0, model$sid, model$y)
    jcol <- model$exact_cols[1]

    one <- .network_one_parameter(
        model,
        jcol,
        slevels_ord,
        condition_on = "model",
        alpha = 0.10,
        ci_type = "two-sided",
        verbose = FALSE
    )

    expect_type(one, "list")
    expect_named(one, c("parm", "outdist_entry", "result_row"))
    expect_identical(one$parm, "group1")

    fit <- ep_network(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        alpha = 0.10,
        verbose = FALSE
    )
    expect_identical(one$outdist_entry, fit$OUTDIST[["group1"]])
    expect_exactparmest_equal(
        as.data.frame(one$result_row),
        fit$exactparmest[fit$exactparmest$Parameter == "group1", ]
    )

    # Same hand-derived anchor as toy2's other tests: at beta = 0 the
    # two-sided p-value is exactly 1 and the MLE is exactly 0.
    expect_identical(one$result_row$SuffStat, 1L)
    expect_equal(one$result_row$Estimate, 0, tolerance = 1e-8)
    expect_equal(one$result_row$PValue, 1, tolerance = 1e-8)
})

# conditions on a nuisance covariate
test_that(".network_one_parameter agrees with ep_network, p_cond > 0", {
    multi_df <- data.frame(
        group = c(1L, 0L, 1L, 0L),
        region = c(1L, 1L, 0L, 0L),
        cases = c(1L, 0L, 0L, 1L)
    )
    for (cond in c("exact", "model")) {
        model <- .setup_exact_model(
            cases ~ group + region,
            data = multi_df,
            offset = NULL,
            offset_log = TRUE,
            exactmax = NULL,
            strata = NULL,
            exact_vars = c("group", "region"),
            ref_levels = list(),
            verbose = FALSE
        )
        slevels_ord <- .order_strata_by_size(model$slevels0, model$sid, model$y)
        jcol <- model$exact_cols[which(model$cn[model$exact_cols] == "group1")]

        one <- .network_one_parameter(
            model,
            jcol,
            slevels_ord,
            condition_on = cond,
            alpha = 0.05,
            ci_type = "two-sided",
            verbose = FALSE
        )
        expect_identical(one$parm, "group1")

        fit <- ep_network(
            cases ~ group + region,
            data = multi_df,
            exact_vars = c("group", "region"),
            condition_on = cond,
            verbose = FALSE
        )
        expect_identical(one$outdist_entry, fit$OUTDIST[["group1"]])
        expect_exactparmest_equal(
            as.data.frame(one$result_row),
            fit$exactparmest[fit$exactparmest$Parameter == "group1", ]
        )
        expect_identical(one$result_row$SuffStat, 1L)
    }
})

test_that(".network_one_parameter prints its parameter header when verbose", {
    model <- .setup_exact_model(
        cases ~ group,
        data = toy2,
        offset = NULL,
        offset_log = TRUE,
        exactmax = NULL,
        strata = "strat",
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    slevels_ord <- .order_strata_by_size(model$slevels0, model$sid, model$y)
    jcol <- model$exact_cols[1]

    expect_snapshot(
        .network_one_parameter(
            model,
            jcol,
            slevels_ord,
            condition_on = "model",
            alpha = 0.05,
            ci_type = "two-sided",
            verbose = TRUE
        )
    )
})

test_that(".network_one_parameter itself raises the mid-stratum pruning error", {
    model <- .setup_exact_model(
        cases ~ group,
        data = toy2,
        offset = NULL,
        offset_log = TRUE,
        exactmax = NULL,
        strata = "strat",
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    slevels_ord <- .order_strata_by_size(model$slevels0, model$sid, model$y)
    jcol <- model$exact_cols[1]

    testthat::local_mocked_bindings(
        .network_expand_observation = function(...) {
            list(m_j = integer(0), t_j = integer(0), lw = numeric(0))
        }
    )
    expect_error(
        .network_one_parameter(
            model,
            jcol,
            slevels_ord,
            condition_on = "model",
            alpha = 0.05,
            ci_type = "two-sided",
            verbose = FALSE
        ),
        "All states pruned mid-stratum"
    )
})

test_that(".network_one_parameter itself raises the end-of-stratum pruning error", {
    model <- .setup_exact_model(
        cases ~ group,
        data = toy2,
        offset = NULL,
        offset_log = TRUE,
        exactmax = NULL,
        strata = "strat",
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    slevels_ord <- .order_strata_by_size(model$slevels0, model$sid, model$y)
    jcol <- model$exact_cols[1]

    testthat::local_mocked_bindings(
        .network_finish_stratum = function(...) {
            list(
                m_j = integer(0),
                t_j = integer(0),
                lw = numeric(0)
            )
        }
    )
    expect_error(
        .network_one_parameter(
            model,
            jcol,
            slevels_ord,
            condition_on = "model",
            alpha = 0.05,
            ci_type = "two-sided",
            verbose = FALSE
        ),
        "All states pruned after finishing stratum"
    )
})
