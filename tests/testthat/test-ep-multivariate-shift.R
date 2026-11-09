# ==============================================================================
# Tests for ep-multivariate-shift.R; the per-stratum dynamic
# program, cross-stratum convolution, and the exported
# ep_multi_shift(). One of the three exported main functions,
# so tested extensively: the internal DP/convolution building blocks in
# isolation (with fully hand-traced examples, several cross-checked
# against independent combinatorial arguments), then the end-to-end
# function anchored on three independently derived numeric examples:
#
#   * `toy2` (2 strata, no nuisance covariates), same anchor used in
#     test-ep-closed-form.R, giving Estimate = 0, PValue = 1.
#   * `binom5`-equivalent (1 stratum, no nuisance covariates, unequal
#     exposure), cross-checked against ep_closed_form()'s own
#     result for the identical data, since both algorithms must agree
#     exactly on this special case.
#   * `multi_df` (1 stratum, WITH a nuisance covariate to condition on);
#     hand-derived from first principles via the combinatorial argument
#     that, conditional on 2 total cases distributed across 4 cells with
#     equal (unit) exposure, the joint distribution is a simple
#     enumeration over ways to place 2 indistinguishable case-units in 4
#     cells (see the derivation in comments below), which happens to
#     collapse to the same (1, 2, 1) shape once filtered to the observed
#     nuisance total.
# ==============================================================================

# ------------------------------------------------------------------------
# .joint_dist_stratum_dp_vec ----
# ------------------------------------------------------------------------

# 2-observation example
test_that(".joint_dist_stratum_dp_vec matches hand-derived no-nuisance case", {
    # y = c(1, 0), xj = c(1, 0), xc = NULL, log_offset = c(0, 0), si = c(1, 1).
    # Only one combination of (y_1, y_2) summing to 1 is feasible per
    # observation's own cap (si = 1 each): either (1, 0) or (0, 1), each with
    # weight 1/(1!0!) = 1. (1,0) contributes t = 1*1 + 0*0 = 1; (0,1)
    # contributes t = 0*1 + 1*0 = 0. Hence two states, t_j in {0, 1}, lw = 0
    # (log(1)) each.
    out <- .joint_dist_stratum_dp_vec(
        y = c(1L, 0L),
        xj = c(1L, 0L),
        xc = NULL,
        log_offset = c(0, 0),
        si = c(1L, 1L)
    )
    out <- order_tbl(out, "t_j")
    expect_identical(out$t_j, c(0L, 1L))
    expect_equal(out$lw, c(0, 0), tolerance = 1e-10)
})

# shape is (1,1/2,1), caps are per-observation
test_that(".joint_dist_stratum_dp_vec matches a 3-state shape, larger caps", {
    # y = c(1, 1) (ms = 2), xj = c(1, 0), xc = NULL, si = c(2, 2): three ways
    # to split 2 units between 2 observations, (0,2), (1,1), (2,0), with
    # weights 1/(0!2!) = 1/2, 1/(1!1!) = 1, 1/(2!0!) = 1/2 respectively, and
    # t = y1*1 + y2*0 = y1, giving t = 0 (weight 1/2), t = 1 (weight 1),
    # t = 2 (weight 1/2), i.e. log_ck = log(c(0.5, 1, 0.5)).
    out <- .joint_dist_stratum_dp_vec(
        y = c(1L, 1L),
        xj = c(1L, 0L),
        xc = NULL,
        log_offset = c(0, 0),
        si = c(2L, 2L)
    )
    out <- order_tbl(out, "t_j")
    expect_identical(out$t_j, c(0L, 1L, 2L))
    expect_equal(exp(out$lw), c(0.5, 1, 0.5), tolerance = 1e-10)
})

test_that(".joint_dist_stratum_dp_vec tracks a nuisance column correctly", {
    # Same setup as above but with a nuisance column xc equal to xj itself
    # reversed (xc = c(0, 1)), so cond_1 tallies observation 2's contribution
    # only. Enumerate: (y1,y2) = (0,2): t=0,cond=2*1=2, w=1/2.
    # (1,1): t=1, cond=1*1=1, w=1. (2,0): t=2, cond=0, w=1/2.
    out <- .joint_dist_stratum_dp_vec(
        y = c(1L, 1L),
        xj = c(1L, 0L),
        xc = matrix(c(0L, 1L), ncol = 1),
        log_offset = c(0, 0),
        si = c(2L, 2L)
    )
    out <- order_tbl(out, "t_j")
    expect_identical(out$t_j, c(0L, 1L, 2L))
    expect_identical(out$cond_1, c(2L, 1L, 0L))
    expect_equal(exp(out$lw), c(0.5, 1, 0.5), tolerance = 1e-10)
})

test_that(".joint_dist_stratum_dp_vec's ms == 0 case returns all-zero state", {
    out0 <- .joint_dist_stratum_dp_vec(
        y = c(0L, 0L),
        xj = c(1L, 1L),
        xc = NULL,
        log_offset = c(0, 0),
        si = c(3L, 3L)
    )
    expect_identical(out0$t_j, 0L)
    expect_equal(out0$lw, 0)

    out1 <- .joint_dist_stratum_dp_vec(
        y = c(0L, 0L),
        xj = c(1L, 1L),
        xc = matrix(c(1L, 2L), ncol = 1),
        log_offset = c(0, 0),
        si = c(3L, 3L)
    )
    expect_identical(out1$t_j, 0L)
    expect_identical(out1$cond_1, 0L)
})

# also covers negative counts
test_that(".joint_dist_stratum_dp_vec errors on NA/negative/bad offsets", {
    expect_error(
        .joint_dist_stratum_dp_vec(
            c(1L, NA),
            c(1L, 1L),
            NULL,
            c(0, 0),
            c(1L, 1L)
        ),
        "NA in"
    )
    expect_error(
        .joint_dist_stratum_dp_vec(
            c(1L, -1L),
            c(1L, 1L),
            NULL,
            c(0, 0),
            c(1L, 1L)
        ),
        "Negative counts"
    )
    expect_error(
        .joint_dist_stratum_dp_vec(
            c(1L, 1L),
            c(1L, 1L),
            NULL,
            c(0, Inf),
            c(1L, 1L)
        ),
        "Non-finite log_offset"
    )
})

test_that(".joint_dist_stratum_dp_vec errors when a count exceeds its cap", {
    expect_error(
        .joint_dist_stratum_dp_vec(
            c(2L, 0L),
            c(1L, 1L),
            NULL,
            c(0, 0),
            c(1L, 1L)
        ),
        "exceeds S_i"
    )
})

# Note: `.joint_dist_stratum_dp_vec()`'s `ms > sum(si)` check
# ("Stratum total M_s exceeds sum(S_i)") is unreachable given
# well-formed input: since the immediately preceding check already
# guarantees `y_i <= si_i` for every row (otherwise it stops first with
# "exceeds S_i for at least one row"), `sum(y) <= sum(si)` follows
# automatically. It is a defensive guard, not exercised here for the same
# reason `.stratum_slice()`'s callers never need to test it either.

# Rule A pruning against tn_obs
test_that(".joint_dist_stratum_dp_vec pruning can exhaust states mid-stratum", {
    # A single row with xc = 1 and tn_obs = -1: even the yy = 0 candidate
    # (cond_1 = 0) fails 0 <= -1, so every state is pruned immediately.
    expect_error(
        .joint_dist_stratum_dp_vec(
            y = 1L,
            xj = 1L,
            xc = matrix(1L, ncol = 1),
            log_offset = 0,
            si = 1L,
            tn_obs = -1L
        ),
        "All states pruned mid-stratum"
    )
})

# ------------------------------------------------------------------------
# .joint_dist_stratum_dp_vec's helper functions, tested directly ----
# ------------------------------------------------------------------------

test_that(".validate_stratum_dp_inputs passes valid input, errors otherwise", {
    expect_silent(.validate_stratum_dp_inputs(c(1L, 0L), c(0, 0), c(1L, 1L)))
    expect_error(
        .validate_stratum_dp_inputs(c(1L, NA), c(0, 0), c(1L, 1L)),
        "NA in"
    )
    expect_error(
        .validate_stratum_dp_inputs(c(1L, -1L), c(0, 0), c(1L, 1L)),
        "Negative counts"
    )
    expect_error(
        .validate_stratum_dp_inputs(c(1L, 1L), c(0, Inf), c(1L, 1L)),
        "Non-finite log_offset"
    )
})

test_that(".zero_total_stratum_table matches the ms == 0 dp_vec shortcut", {
    expect_identical(
        .zero_total_stratum_table(NULL),
        list(t_j = 0L, lw = 0.0)
    )
    expect_identical(
        .zero_total_stratum_table(matrix(c(1L, 2L), ncol = 1)),
        list(t_j = 0L, cond_1 = 0L, lw = 0.0)
    )
})

# also caps si and normalizes xc
test_that(".prepare_stratum_dp_inputs coerces si/xc, checks feasibility", {
    out <- .prepare_stratum_dp_inputs(
        y = c(1L, 1L),
        xj = c(1.0, 0.0),
        xc = NULL,
        si = c(5L, 5L),
        ms = 2L
    )
    expect_identical(out$si, c(2L, 2L))
    expect_identical(out$xj, c(1L, 0L))
    expect_identical(out$p_cond, 0L)
    expect_identical(dim(out$xc), c(2L, 0L))

    out_xc <- .prepare_stratum_dp_inputs(
        y = c(1L, 1L),
        xj = c(1L, 0L),
        xc = matrix(c(0.4, 1.0), ncol = 1),
        si = c(2L, 2L),
        ms = 2L
    )
    expect_identical(out_xc$p_cond, 1L)
    expect_identical(as.integer(out_xc$xc[, 1]), c(0L, 1L)) # round()ed

    expect_error(
        .prepare_stratum_dp_inputs(c(2L, 0L), c(1L, 1L), NULL, c(1L, 1L), 2L),
        "exceeds S_i"
    )
})

test_that(".prepare_stratum_dp_inputs errors when ms exceeds sum(si)", {
    # Unreachable through .joint_dist_stratum_dp_vec()'s public entry point:
    # `ms` is always computed there as `sum(y)`, and the immediately
    # preceding check already guarantees every row satisfies y_i <= si_i,
    # so sum(y) <= sum(si) follows automatically.
    expect_error(
        .prepare_stratum_dp_inputs(
            y = c(1L, 1L), # per-row check passes: y_i <= si_i for both rows
            xj = c(1L, 0L),
            xc = NULL,
            si = c(1L, 1L),
            ms = 5L # deliberately inconsistent with y (sum(y) = 2, not 5)
        ),
        "Stratum total M_s exceeds sum\\(S_i\\)"
    )
})

# xc is a zero-column matrix, not NULL
test_that(".prepare_stratum_dp_inputs handles an explicit zero-column xc", {
    # Distinct from the `xc = NULL` case already covered above: here the
    # caller passes an actual matrix with zero columns, exercising the
    # `else` branch of `if (ncol(xc) > 0L) ... else ...` rather than the
    # `is.null(xc)` branch.
    out <- .prepare_stratum_dp_inputs(
        y = c(1L, 1L),
        xj = c(1L, 0L),
        xc = matrix(numeric(0), nrow = 2, ncol = 0),
        si = c(2L, 2L),
        ms = 2L
    )
    expect_identical(out$p_cond, 0L)
    expect_identical(dim(out$xc), c(2L, 0L))
    expect_identical(out$cond_cols, character(0))
})

test_that(".dp_expand_observation matches one hand-traced step of the DP", {
    # Starting frontier: a single state (m_j=0, t_j=0, lw=0). One
    # observation with si = 1 (yy in {0,1}), xj_i = 1, log_offset_i = 0, no
    # nuisance columns: expect two resulting states, (m_j=0,t_j=0) and
    # (m_j=1,t_j=1), both lw = 0 (log(1/0!) = log(1/1!) = 0).
    state0 <- list(m_j = 0L, t_j = 0L, lw = 0.0)
    out <- .dp_expand_observation(
        state = state0,
        yy_vec = 0:1,
        xj_i = 1L,
        log_offset_i = 0,
        xc_i = numeric(0),
        lgam = lgamma(0:1 + 1),
        cond_cols = character(0),
        p_cond = 0L,
        ms = 1L,
        tn_obs = NULL
    )
    out <- order_tbl(out, "t_j")
    expect_identical(out$t_j, c(0L, 1L))
    expect_equal(out$lw, c(0, 0), tolerance = 1e-10)
})

# pruning uses Rule A
test_that(".dp_expand_observation prunes vs ms, errors if states exhausted", {
    state0 <- list(m_j = 0L, t_j = 0L, lw = 0.0)
    # ms = 0 means even the yy = 0 candidate produces m_j = 0, which is the
    # only surviving state (not pruned).
    out <- .dp_expand_observation(
        state = state0,
        yy_vec = 0:0,
        xj_i = 1L,
        log_offset_i = 0,
        xc_i = numeric(0),
        lgam = lgamma(0:0 + 1),
        cond_cols = character(0),
        p_cond = 0L,
        ms = 0L,
        tn_obs = NULL
    )
    expect_identical(out$m_j, 0L)

    # tn_obs = -1 with a nuisance column present: every candidate's cond_1
    # is >= 0, so every state fails `<= tn_obs` and the function errors.
    expect_error(
        .dp_expand_observation(
            state = list(m_j = 0L, t_j = 0L, cond_1 = 0L, lw = 0.0),
            yy_vec = 0:1,
            xj_i = 1L,
            log_offset_i = 0,
            xc_i = 1L,
            lgam = lgamma(0:1 + 1),
            cond_cols = "cond_1",
            p_cond = 1L,
            ms = 1L,
            tn_obs = -1L
        ),
        "All states pruned mid-stratum"
    )
})

test_that(".finalize_stratum_dp keeps m_j == ms, errors if nothing survives", {
    state <- list(m_j = c(0L, 1L, 2L), t_j = c(0L, 1L, 2L), lw = c(0, 0, 0))
    out <- .finalize_stratum_dp(state, ms = 1L, cond_cols = character(0))
    expect_identical(out$t_j, 1L)
    expect_null(out$m_j) # m_j itself is dropped from the finalized table

    expect_error(
        .finalize_stratum_dp(
            list(m_j = c(0L, 2L), t_j = c(0L, 2L), lw = c(0, 0)),
            ms = 1L,
            cond_cols = character(0)
        ),
        "No feasible states at m == ms"
    )
})

# ------------------------------------------------------------------------
# .order_joint_tables_by_size ----
# ------------------------------------------------------------------------

test_that(".order_joint_tables_by_size orders lists by ascending row count", {
    joint_list <- list(
        list(t_j = 1, lw = 0), # nrow 1
        list(t_j = c(1, 2, 3), lw = c(0, 0, 0)), # nrow 3
        list(t_j = c(1, 2), lw = c(0, 0)) # nrow 2
    )
    max_add_list <- list("a", "b", "c")
    out <- .order_joint_tables_by_size(joint_list, max_add_list)
    expect_equal(vapply(out$joint_list, .tbl_nrow, integer(1)), c(1L, 2L, 3L))
    expect_identical(out$max_add_list, list("a", "c", "b"))
})

test_that(".order_joint_tables_by_size leaves a single-element list as-is", {
    jl <- list(list(t_j = 1, lw = 0))
    out <- .order_joint_tables_by_size(jl, list("x"))
    expect_identical(out$joint_list, jl)
    expect_identical(out$max_add_list, list("x"))
})

# ------------------------------------------------------------------------
# .stratum_slice ----
# ------------------------------------------------------------------------

test_that(".stratum_slice extracts and computes rule-2 caps for a fake model", {
    model <- list(
        sid = c("A", "A", "B"),
        y = c(2L, 3L, 1L),
        x = matrix(c(1, 0, 1, 0, 1, 0), ncol = 2),
        log_off = c(0, 0, 0),
        use_rule1 = FALSE,
        use_rule2 = TRUE,
        use_rule3 = FALSE,
        exactmax_vec = NULL
    )
    info <- .stratum_slice(
        model,
        "A",
        jcol = 1,
        cond_cols_idx = 2,
        p_cond = 1,
        verbose = FALSE
    )
    expect_identical(info$idx, c(1L, 2L))
    expect_identical(info$y_s, c(2L, 3L))
    expect_identical(info$xj_s, c(1L, 0L))
    expect_identical(info$si_s, c(5L, 5L)) # rule 2: stratum total = 5 for both
    expect_identical(as.numeric(info$xc_s[, 1]), c(0, 1))
    expect_identical(info$m_s_obs, 5L)
    expect_identical(info$n_s, 2L)
    expect_identical(info$max_add_s, 5L)
})

test_that(".stratum_slice applies rule 1 (EXACTMAX=var) caps when requested", {
    model <- list(
        sid = c("A", "A"),
        y = c(1L, 1L),
        x = matrix(c(1, 0), ncol = 1),
        log_off = c(0, 0),
        use_rule1 = TRUE,
        use_rule2 = FALSE,
        use_rule3 = FALSE,
        exactmax_vec = c(9L, 4L)
    )
    info <- .stratum_slice(
        model,
        "A",
        jcol = 1,
        cond_cols_idx = integer(0),
        p_cond = 0,
        verbose = FALSE
    )
    expect_identical(info$si_s, c(9L, 4L))
    expect_null(info$xc_s)
    expect_identical(info$max_add_s, integer(0))
})

test_that(".stratum_slice applies rule 3 offset-derived caps when requested", {
    model <- list(
        sid = c("A", "A"),
        y = c(1L, 1L),
        x = matrix(c(1, 0), ncol = 1),
        log_off = log(c(3, 0.2)),
        use_rule1 = FALSE,
        use_rule2 = FALSE,
        use_rule3 = TRUE,
        exactmax_vec = NULL
    )
    info <- .stratum_slice(
        model,
        "A",
        jcol = 1,
        cond_cols_idx = integer(0),
        p_cond = 0,
        verbose = FALSE
    )
    expect_identical(info$si_s, si_from_offset(log(c(3, 0.2))))
})

test_that(".stratum_slice's verbose = TRUE prints stratum size and total", {
    model <- list(
        sid = c("A", "A"),
        y = c(2L, 3L),
        x = matrix(c(1, 0), ncol = 1),
        log_off = c(0, 0),
        use_rule1 = FALSE,
        use_rule2 = TRUE,
        use_rule3 = FALSE,
        exactmax_vec = NULL
    )
    expect_message(
        .stratum_slice(
            model,
            "A",
            jcol = 1,
            cond_cols_idx = integer(0),
            p_cond = 0,
            verbose = TRUE
        ),
        "stratum A: n=2, M=5"
    )
})

# ------------------------------------------------------------------------
# .convolve_joint_pnuis_vec ----
# ------------------------------------------------------------------------

# toy2 is (1,2,1), no nuisance columns
test_that(".convolve_joint_pnuis_vec reproduces the toy2 convolution", {
    d1 <- list(t_j = c(0L, 1L), lw = log(c(1, 1)))
    d2 <- list(t_j = c(0L, 1L), lw = log(c(1, 1)))
    out <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 0,
        tn_obs = integer(0),
        remaining_max = integer(0)
    )
    out <- order_tbl(out, "t_j")
    expect_identical(out$t_j, c(0L, 1L, 2L))
    expect_equal(exp(out$lw), c(1, 2, 1), tolerance = 1e-10)
})

# checked against asymmetric tables
test_that(".convolve_joint_pnuis_vec matches manual poly multiplication", {
    # (1 + 2x + x^2) * (1 + 3x) = 1 + 5x + 7x^2 + 3x^3 (worked out by hand
    # in the file header / test comments).
    d1 <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    d2 <- list(t_j = 0:1, lw = log(c(1, 3)))
    out <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 0,
        tn_obs = integer(0),
        remaining_max = integer(0)
    )
    out <- order_tbl(out, "t_j")
    expect_identical(out$t_j, 0:3)
    expect_equal(exp(out$lw), c(1, 5, 7, 3), tolerance = 1e-8)
})

test_that(".convolve_joint_pnuis_vec is invariant to chunk_size", {
    d1 <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    d2 <- list(t_j = 0:1, lw = log(c(1, 3)))
    out_small_chunk <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 0,
        tn_obs = integer(0),
        remaining_max = integer(0),
        chunk_size = 1L
    )
    out_big_chunk <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 0,
        tn_obs = integer(0),
        remaining_max = integer(0),
        chunk_size = 1000L
    )
    out_small_chunk <- order_tbl(out_small_chunk, "t_j")
    out_big_chunk <- order_tbl(out_big_chunk, "t_j")
    expect_equal(out_small_chunk, out_big_chunk, tolerance = 1e-10)
})

test_that(".convolve_joint_pnuis_vec prunes (Rule A/B) on a nuisance column", {
    d1 <- list(
        t_j = c(0L, 1L),
        cond_1 = c(0L, 1L),
        lw = log(c(1, 1))
    )
    d2 <- list(
        t_j = c(0L, 1L),
        cond_1 = c(0L, 1L),
        lw = log(c(1, 1))
    )
    # tn_obs = 1, remaining_max = 0 (no more strata after this) => only
    # cond_1 == 1 survives (by hand: pairs (0,0)->cond0 dropped by Rule B,
    # (0,1) and (1,0) -> cond 1 survive, (1,1) -> cond 2 dropped by Rule A).
    out <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 1,
        tn_obs = 1L,
        remaining_max = 0L
    )
    expect_identical(.tbl_nrow(out), 1L)
    expect_identical(out$t_j, 1L)
    expect_identical(out$cond_1, 1L)
    expect_equal(exp(out$lw), 2, tolerance = 1e-10)
})

test_that(".convolve_joint_pnuis_vec's pruning is chunk_size-invariant", {
    d1 <- list(
        t_j = c(0L, 1L),
        cond_1 = c(0L, 1L),
        lw = log(c(1, 1))
    )
    d2 <- list(
        t_j = c(0L, 1L),
        cond_1 = c(0L, 1L),
        lw = log(c(1, 1))
    )
    out1 <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 1,
        tn_obs = 1L,
        remaining_max = 0L,
        chunk_size = 1L
    )
    out2 <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 1,
        tn_obs = 1L,
        remaining_max = 0L,
        chunk_size = 100L
    )
    expect_equal(out1, out2, tolerance = 1e-10)
})

test_that(".convolve_joint_pnuis_vec returns empty table if all pairs pruned", {
    d1 <- list(
        t_j = c(0L, 1L),
        cond_1 = c(0L, 1L),
        lw = log(c(1, 1))
    )
    d2 <- list(
        t_j = c(0L, 1L),
        cond_1 = c(0L, 1L),
        lw = log(c(1, 1))
    )
    out <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 1,
        tn_obs = 10L,
        remaining_max = 0L
    )
    expect_identical(.tbl_nrow(out), 0L)
})

test_that(".convolve_joint_pnuis_vec is invariant to d1/d2 argument order", {
    d1 <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    d2 <- list(t_j = 0:1, lw = log(c(1, 3)))
    out_ab <- .convolve_joint_pnuis_vec(
        d1,
        d2,
        p_cond = 0,
        tn_obs = integer(0),
        remaining_max = integer(0)
    )
    out_ba <- .convolve_joint_pnuis_vec(
        d2,
        d1,
        p_cond = 0,
        tn_obs = integer(0),
        remaining_max = integer(0)
    )
    out_ab <- order_tbl(out_ab, "t_j")
    out_ba <- order_tbl(out_ba, "t_j")
    expect_equal(out_ab, out_ba, tolerance = 1e-10)
})

# ------------------------------------------------------------------------
# ep_multi_shift(): end-to-end, anchored on toy2 ----
# ------------------------------------------------------------------------

test_that("ep_multi_shift reproduces every hand-derived toy2 quantity", {
    fit <- ep_multi_shift(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        alpha = 0.10,
        verbose = FALSE
    )
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

test_that("ep_multi_shift is chunk_size-invariant for a multi-stratum case", {
    fit_default <- ep_multi_shift(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        verbose = FALSE
    )
    fit_chunk1 <- ep_multi_shift(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        chunk_size = 1L,
        verbose = FALSE
    )
    expect_equal(
        fit_default$exactparmest,
        fit_chunk1$exactparmest,
        tolerance = 1e-10
    )
})

# ------------------------------------------------------------------------
# ep_multi_shift(): end-to-end with a genuine nuisance covariate ----
# to condition on, anchored on a hand-derived combinatorial argument
# ------------------------------------------------------------------------
#
# multi_df: 4 observations, a 2x2 (group x region) design, unit exposure,
# 2 total cases. Conditional on the total (M = 2) and the observed
# `region1` nuisance total (n_region1 = 1), the exact conditional
# distribution of the group sufficient statistic reduces
# to exactly the same (1, 2, 1) shape used throughout this suite, with
# the observed t_obs = 1 landing on the symmetric mean.

multi_df <- data.frame(
    group = c(1L, 0L, 1L, 0L),
    region = c(1L, 1L, 0L, 0L),
    cases = c(1L, 0L, 0L, 1L)
)

# anchor is combinatorial
test_that("ep_multi_shift matches the hand-derived anchor, nuisance cov", {
    fit <- ep_multi_shift(
        cases ~ group + region,
        data = multi_df,
        exact_vars = "group",
        verbose = FALSE
    )
    row <- fit$exactparmest
    expect_identical(row$Parameter, "group1")
    expect_identical(row$SuffStat, 1L)
    expect_equal(
        as.list(row[c("Estimate", "StdError", "PValue")]),
        list(Estimate = 0, StdError = sqrt(2), PValue = 1),
        tolerance = 1e-6
    )

    od <- fit$OUTDIST[["group1"]]
    expect_identical(od$SuffStat, 0:2)
    expect_equal(od$prob, c(0.25, 0.5, 0.25), tolerance = 1e-10)
})

test_that("ep_multi_shift's condition_on = 'exact' runs for multiple vars", {
    fit <- ep_multi_shift(
        cases ~ group + region,
        data = multi_df,
        exact_vars = c("group", "region"),
        condition_on = "exact",
        verbose = FALSE
    )
    expect_identical(fit$exactparmest$Parameter, c("group1", "region1"))
    expect_true(
        all(
            fit$exactparmest$PValue >= 0 & fit$exactparmest$PValue <= 1
        )
    )
    expect_identical(fit$exactparmest$SuffStat, c(1L, 1L))
})

test_that("ep_multi_shift's condition_on = 'model' returns valid results", {
    fit <- ep_multi_shift(
        cases ~ group + region,
        data = multi_df,
        exact_vars = c("group", "region"),
        condition_on = "model",
        verbose = FALSE
    )
    expect_identical(fit$exactparmest$Parameter, c("group1", "region1"))
    expect_true(
        all(
            fit$exactparmest$PValue >= 0 & fit$exactparmest$PValue <= 1
        )
    )
})

# ------------------------------------------------------------------------
# ep_multi_shift(): an actual interaction term as the exact effect ----
# hand-derived
#
# Testing exact_vars = "group:region" on multi_df (formula
# cases ~ group * region) conditions on BOTH main effects (group1,
# region1) as nuisance columns (p_cond = 2, condition_on = "model"
# default). Re-enumerating the same 10 configurations used for the
# main-effect anchor (see the multi_df comment above), but now tracking
# all three quantities (group total g, region total r, interaction total
# i) per configuration and filtering to the observed (g = 1, r = 1)
# leaves exactly two configurations:
#   - both units split across observations 1 and 4 (the two "diagonal"
#     cells: group=1/region=1, and group=0/region=0): g=1, r=1, i=1, w=1
#   - both units split across observations 2 and 3 (group=0/region=1,
#     and group=1/region=0): g=1, r=1, i=0, w=1
# so the interaction's conditional distribution is exactly (i=0: weight
# 1, i=1: weight 1), i.e. log_ck = log(c(1, 1)), NOT the (1,2,1)
# shape of the main-effect anchors. The observed interaction total is
# t_obs = 1 (only observation 1, group=1 & region=1, has cases > 0
# among the two "both green" rows), which sits at the *upper* boundary
# of this 2-state support, so the CMLE is +Inf and the upper CI limit is
# Inf by definition (not by root-finding).
# ------------------------------------------------------------------------

test_that("ep_multi_shift matches the interaction anchor (group:region)", {
    fit <- ep_multi_shift(
        cases ~ group * region,
        data = multi_df,
        exact_vars = "group:region",
        verbose = FALSE
    )
    row <- fit$exactparmest
    expect_identical(row$Parameter, "group1:region1")
    expect_identical(row$SuffStat, 1L)
    expect_identical(row$Estimate, Inf) # the upper bound of a 2-state support
    expect_true(is.na(row$StdError))
    expect_equal(row$PValue, 1, tolerance = 1e-8)
    expect_identical(row$UpperCL, Inf)
    expect_true(is.finite(row$LowerCL))

    od <- fit$OUTDIST[["group1:region1"]]
    expect_identical(od$SuffStat, 0:1)
    expect_equal(od$prob, c(0.5, 0.5), tolerance = 1e-10)
})

test_that("ep_multi_shift resolves 'group*region' and 'group:region' equally", {
    fit_colon <- ep_multi_shift(
        cases ~ group * region,
        data = multi_df,
        exact_vars = "group:region",
        verbose = FALSE
    )
    fit_star <- ep_multi_shift(
        cases ~ group * region,
        data = multi_df,
        exact_vars = "region*group",
        verbose = FALSE
    )
    expect_exactparmest_equal(
        fit_colon$exactparmest,
        fit_star$exactparmest,
        tolerance = 1e-10
    )
})

# only when main effect alone is requested
test_that("ep_multi_shift never implicitly tests interaction terms", {
    fit <- ep_multi_shift(
        cases ~ group * region,
        data = multi_df,
        exact_vars = "group",
        verbose = FALSE
    )
    expect_identical(fit$exactparmest$Parameter, "group1")
})

# ------------------------------------------------------------------------
# ep_multi_shift(): EXACTMAX handling, verbose, and error propagation ----
# ------------------------------------------------------------------------

test_that("ep_multi_shift honours EXACTMAX=var (rule 1)", {
    fit <- ep_multi_shift(
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

test_that("ep_multi_shift's verbose = TRUE prints per-param/stratum progress", {
    expect_snapshot(
        ep_multi_shift(
            cases ~ group,
            data = toy2,
            strata = "strat",
            exact_vars = "group",
            verbose = TRUE
        )
    )
})

test_that("ep_multi_shift propagates argument-validation errors", {
    expect_error(
        ep_multi_shift(
            cases ~ group,
            data = toy2,
            strata = "strat",
            exact_vars = character(0),
            verbose = FALSE
        ),
        "non-empty character vector"
    )
})

test_that("ep_multi_shift errors when EXACTMAX is below an observed count", {
    # Caught by .stratum_slice()'s own feasibility check now (shared with
    # ep_network()), earlier and with a more specific message
    # than the internal .joint_dist_stratum_dp_vec() check this used to
    # reach; see the dedicated test below and the matching test in
    # test-ep-network.R.
    bad <- data.frame(group = c(1, 0), cases = c(5, 1), n = c(2, 2))
    expect_error(
        ep_multi_shift(
            cases ~ group,
            data = bad,
            exact_vars = "group",
            exactmax = "n",
            verbose = FALSE
        ),
        "exceed the per-observation cap"
    )
})

# not just ep_network()
test_that("ep_multi_shift errors when a count exceeds its SAS-rule-3 cap", {
    # Same reproduction as the ep_network() test in
    # test-ep-network.R: stratum 2's control-arm row has
    # futime = 1, giving a SAS-rule-3 cap of S_i = 1, but an observed count
    # of 3. Both algorithms now raise the same .stratum_slice()-level error.
    bad <- data.frame(
        group = c(1, 0, 1, 0),
        strat = c(1, 1, 2, 2),
        cases = c(2, 1, 0, 3),
        fu = c(5, 5, 5, 1)
    )
    expect_error(
        ep_multi_shift(
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

# unexpected/defensive case
test_that("ep_multi_shift errors if first stratum table is pruned to zero", {
    # Defensive-only: given valid input, the first stratum's own table
    # always contains the state matching the true observed data, and
    # tn_obs is itself derived from that same data, so this Rule B filter
    # can never empty it out in practice. Verifying the guard fires
    # correctly requires forcing .joint_dist_stratum_dp_vec() to return a
    # state that could never actually satisfy the nuisance target, via a mock.
    testthat::local_mocked_bindings(
        .joint_dist_stratum_dp_vec = function(...) {
            list(t_j = 0L, cond_1 = -1000000L, lw = 0.0)
        }
    )
    d <- data.frame(
        group = c(1L, 0L, 1L, 0L),
        strat = c(1L, 1L, 2L, 2L),
        z = c(1L, 0L, 1L, 0L),
        cases = c(2L, 1L, 0L, 3L)
    )
    expect_error(
        ep_multi_shift(
            cases ~ group + z,
            data = d,
            strata = "strat",
            exact_vars = "group",
            verbose = FALSE
        ),
        "All states pruned after first stratum"
    )
})

# unexpected/defensive case
test_that("ep_multi_shift errors if a convolution empties the table", {
    # Same reasoning as the test above, for the *other* defensive
    # "all states pruned" guard: this one fires after combining strata 2+,
    # rather than after the first stratum's own Rule B filter. Unreachable
    # with valid input for the same reason (the true observed allocation
    # always survives a correct bound); forced here by mocking
    # .convolve_joint_pnuis_vec() itself to return an empty table, leaving
    # the (real, correct) per-stratum tables it's given untouched.
    testthat::local_mocked_bindings(
        .convolve_joint_pnuis_vec = function(...) {
            list(t_j = integer(0), lw = numeric(0))
        }
    )
    expect_error(
        ep_multi_shift(
            cases ~ group,
            data = test1,
            strata = "age",
            exact_vars = "group",
            verbose = FALSE
        ),
        "All states pruned after combining \\d+ strata"
    )
})

test_that("ep_multi_shift's print method identifies the algorithm", {
    fit <- ep_multi_shift(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        verbose = FALSE
    )
    expect_snapshot(print(fit))
})

# ------------------------------------------------------------------------
# .multishift_one_parameter() directly; ----
# the per-parameter worker that ep_multi_shift() now delegates its loop
# body to. Its logic (building each stratum's table, pruning, convolving,
# and assembling the final row) is already exercised indirectly through every
# ep_multi_shift() test above, including the verbose-output and
# mock-based pruning-error tests (testthat::local_mocked_bindings patches
# .joint_dist_stratum_dp_vec/.convolve_joint_pnuis_vec at the package
# namespace level, so those mocks still take effect however deep the call
# stack is). These tests instead call the helper directly, against a model
# built by hand with .setup_exact_model(), to pin down its own contract,
# the exact shape of what it returns, and that its result agrees with the
# full ep_multi_shift() fit on the same data, independently of
# how the wrapper happens to invoke it.
# ------------------------------------------------------------------------

test_that(".multishift_one_parameter agrees with ep_multi_shift (p_cond = 0)", {
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
    jcol <- model$exact_cols[1]

    one <- .multishift_one_parameter(
        model,
        jcol,
        condition_on = "model",
        chunk_size = 4000L,
        alpha = 0.05,
        ci_type = "two-sided",
        verbose = FALSE
    )

    expect_type(one, "list")
    expect_named(one, c("parm", "outdist_entry", "result_row"))
    expect_identical(one$parm, "group1")

    fit <- ep_multi_shift(
        cases ~ group,
        data = toy2,
        strata = "strat",
        exact_vars = "group",
        verbose = FALSE
    )
    expect_identical(one$outdist_entry, fit$OUTDIST[["group1"]])
    expect_exactparmest_equal(
        as.data.frame(one$result_row),
        fit$exactparmest[fit$exactparmest$Parameter == "group1", ]
    )

    # Hand-derived anchor (see helper-data.R for the full derivation): at
    # beta = 0 the two-sided p-value is exactly 1 and the MLE is exactly 0.
    expect_equal(
        as.list(one$result_row[c("SuffStat", "Estimate", "PValue")]),
        list(SuffStat = 1, Estimate = 0, PValue = 1),
        tolerance = 1e-8
    )
})

# conditions on a nuisance covariate
test_that(".multishift_one_parameter agrees with ep_multi_shift, p_cond > 0", {
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
        jcol <- model$exact_cols[which(model$cn[model$exact_cols] == "group1")]

        one <- .multishift_one_parameter(
            model,
            jcol,
            condition_on = cond,
            chunk_size = 4000L,
            alpha = 0.05,
            ci_type = "two-sided",
            verbose = FALSE
        )
        expect_identical(one$parm, "group1")

        fit <- ep_multi_shift(
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
        expect_equal(one$result_row$SuffStat, 1)
    }
})

test_that(".multishift_one_parameter prints parameter header when verbose", {
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
    jcol <- model$exact_cols[1]

    expect_snapshot(
        .multishift_one_parameter(
            model,
            jcol,
            condition_on = "model",
            chunk_size = 4000L,
            alpha = 0.05,
            ci_type = "two-sided",
            verbose = TRUE
        )
    )
})

# after-first-stratum case
test_that(".multishift_one_parameter raises pruned-states error when mocked", {
    # Same setup as ep_multi_shift's own "first stratum pruned"
    # test: a nuisance covariate ("z") is required so p_cond > 0, since the
    # Rule B filter this guards is only applied in that branch. The true
    # observed allocation can never actually fail this filter, so it is
    # forced here by mocking .joint_dist_stratum_dp_vec() to return a
    # nuisance total that can never satisfy the observed target.
    d <- data.frame(
        group = c(1L, 0L, 1L, 0L),
        strat = c(1L, 1L, 2L, 2L),
        z = c(1L, 0L, 1L, 0L),
        cases = c(2L, 1L, 0L, 3L)
    )
    model <- .setup_exact_model(
        cases ~ group + z,
        data = d,
        offset = NULL,
        offset_log = TRUE,
        exactmax = NULL,
        strata = "strat",
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    jcol <- model$exact_cols[which(model$cn[model$exact_cols] == "group1")]

    testthat::local_mocked_bindings(
        .joint_dist_stratum_dp_vec = function(...) {
            list(t_j = 0L, cond_1 = -1000000L, lw = 0.0)
        }
    )
    expect_error(
        .multishift_one_parameter(
            model,
            jcol,
            condition_on = "model",
            chunk_size = 4000L,
            alpha = 0.05,
            ci_type = "two-sided",
            verbose = FALSE
        ),
        "All states pruned after first stratum"
    )
})

test_that(".multishift_one_parameter raises the combining-error when mocked", {
    model <- .setup_exact_model(
        cases ~ group,
        data = test1,
        offset = NULL,
        offset_log = TRUE,
        exactmax = NULL,
        strata = "age",
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    jcol <- model$exact_cols[1]

    testthat::local_mocked_bindings(
        .convolve_joint_pnuis_vec = function(...) {
            list(SuffStat = integer(0), LogWeight = numeric(0))
        }
    )
    expect_error(
        .multishift_one_parameter(
            model,
            jcol,
            condition_on = "model",
            chunk_size = 4000L,
            alpha = 0.05,
            ci_type = "two-sided",
            verbose = FALSE
        ),
        "All states pruned after combining \\d+ strata"
    )
})
