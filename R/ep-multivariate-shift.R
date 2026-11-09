# ==============================================================================
# The per-stratum dynamic program and explicit cross-stratum convolution
# at the heart of the multivariate-shift algorithm (Section 5 of the
# theory vignette): a direct, vectorized reference implementation,
# generalizing the multivariate shift algorithm of Hirji, Mehta, and
# Patel (1987) for exact logistic regression to bounded Poisson counts,
# following the EXACTMAX/OFFSET rules documented for SAS PROC GENMOD's
# EXACT statement.
# ==============================================================================

# ------------------------------------------------------------------------
# Vectorized per-stratum dynamic program (log-scale) ----
#
# Per-stratum dynamic program (log-scale) - the multivariate shift
# algorithm, adapted from binary responses to bounded Poisson counts
# ------------------------------------------------------------------------

#' Validate a stratum's raw DP inputs
#'
#' The NA/finiteness/non-negativity checks `.joint_dist_stratum_dp_vec()`
#' runs before doing anything else with `y`/`log_offset`/`si`, split out
#' on its own so the main function's body reads as a sequence of named
#' steps instead of one long block (this is a pure extraction: same
#' checks, same error messages, same order).
#'
#' @param y Integer vector of observed counts (already coerced).
#' @param log_offset Numeric vector of log-offsets.
#' @param si Numeric or integer vector of per-observation caps, still in
#'   whatever type the caller passed (checked for `NA` here, coerced to
#'   integer later by the caller).
#' @return `invisible(NULL)`; called for its `stop()`s only.
#' @noRd
.validate_stratum_dp_inputs <- function(y, log_offset, si) {
    if (anyNA(y) || anyNA(log_offset) || anyNA(si)) {
        stop("NA in y/log_offset/si.", call. = FALSE) # nolint: nonportable_path_linter
    }
    if (!all(is.finite(log_offset))) {
        stop("Non-finite log_offset in stratum.", call. = FALSE)
    }
    if (any(y < 0L)) {
        stop("Negative counts in y.", call. = FALSE)
    }
    invisible(NULL)
}

#' Build the degenerate all-zero stratum table (`sum(y) == 0` shortcut)
#'
#' When a stratum's observed total is zero, the only feasible state is
#' "every observation contributed 0", so the dynamic program's usual
#' per-observation loop is skipped entirely and this single-row table is
#' returned directly. Split out of `.joint_dist_stratum_dp_vec()` as its
#' own function purely to keep that function's main body to the
#' non-degenerate case.
#'
#' @param xc The stratum's raw nuisance-column argument, exactly as
#'   passed to `.joint_dist_stratum_dp_vec()` (i.e. *not* yet coerced to
#'   a matrix by `.prepare_stratum_dp_inputs()`); only its column count
#'   is needed here, and using the raw argument keeps this shortcut
#'   independent of the later coercion step.
#' @return A named list with one row: `t_j = 0L`, one `cond_1`, ...,
#'   `cond_p` element (if any) each `0L`, and `lw = 0.0`.
#' @noRd
.zero_total_stratum_table <- function(xc) {
    p_cond <- if (is.null(xc)) 0L else (ncol(as.matrix(xc)) %||% 0L)
    cond_cols <- make_cond_cols(p_cond)
    out <- list(t_j = 0L)
    for (cc in cond_cols) {
        out[[cc]] <- 0L
    }
    out[["lw"]] <- 0.0
    out
}

#' Coerce and feasibility-check a stratum's DP inputs (non-degenerate case)
#'
#' Once `sum(y) > 0`, `.joint_dist_stratum_dp_vec()` needs `si` capped at
#' the stratum total, `xj`/`xc` rounded to integers, and `xc` normalized
#' to a plain integer matrix with a known column count, all mechanical
#' coercion plus the two feasibility checks (`y_i` never exceeds its own
#' cap; the stratum total is attainable given every cap). Split out on
#' its own so `.joint_dist_stratum_dp_vec()`'s body doesn't interleave
#' this setup with the dynamic program itself.
#'
#' @param y Integer vector of observed counts (already coerced, `sum(y) >
#'   0`).
#' @param xj Numeric vector of the interest column's raw values.
#' @param xc Numeric matrix (or `NULL`) of nuisance-column raw values.
#' @param si Numeric or integer vector of per-observation caps, not yet
#'   coerced.
#' @param ms `sum(y)`, this stratum's observed total.
#' @return A list with `si` (integer, capped at `ms`), `xj` (integer),
#'   `xc` (integer matrix, `nrow(xc) == length(y)`, possibly zero
#'   columns), `p_cond` (`ncol(xc)`), and `cond_cols`
#'   (`make_cond_cols(p_cond)`).
#' @noRd
.prepare_stratum_dp_inputs <- function(y, xj, xc, si, ms) {
    si <- as.integer(si)
    if (any(y > si)) {
        stop(
            "Observed y_i exceeds S_i for at least one row.",
            call. = FALSE
        )
    }
    if (ms > sum(si)) {
        stop(
            "Stratum total M_s exceeds sum(S_i): infeasible under bounds.",
            call. = FALSE
        )
    }
    si <- pmin(si, ms)

    n <- length(y)
    xj <- as.integer(round(xj))
    if (is.null(xc)) {
        xc <- matrix(integer(0L), nrow = n, ncol = 0L)
    } else {
        xc <- as.matrix(xc)
        if (ncol(xc) > 0L) {
            xc <- apply(xc, 2L, function(z) as.integer(round(z)))
            xc <- as.matrix(xc)
        } else {
            xc <- matrix(integer(0L), nrow = n, ncol = 0L)
        }
    }

    p_cond <- as.integer(ncol(xc) %||% 0L)
    cond_cols <- make_cond_cols(p_cond)

    list(
        si = si,
        xj = xj,
        xc = xc,
        p_cond = p_cond,
        cond_cols = cond_cols
    )
}

#' Expand the DP frontier by one observation, prune, and re-merge
#'
#' The body of `.joint_dist_stratum_dp_vec()`'s per-observation loop: cross
#' every currently live state with every feasible value the new
#' observation's own count could take (one vectorized row-replication,
#' not a nested loop), apply Rule A pruning (drop states that have
#' already overshot the stratum total or, if `tn_obs` is supplied, any
#' nuisance coordinate's observed total, safe because every increment
#' is non-negative, so an overshoot can never be corrected by a later
#' observation), and merge states landing on the same key by log-sum-exp
#' (one grouping call for this observation, not one per candidate value).
#' Split out of `.joint_dist_stratum_dp_vec()` as its own function so that
#' function's `for` loop reads as "expand, then move on" instead of
#' inlining the whole cross-product-and-prune step; this is also the one
#' piece of the dynamic program most worth unit-testing in isolation
#' (Rule A's boundary behaviour in particular).
#'
#' @param state The current frontier: a named list of equal-length
#'   vectors with elements `m_j`, `t_j`, one `cond_1`, ..., `cond_p`
#'   element per nuisance coordinate, and `lw`.
#' @param yy_vec Integer vector `0:si[i]`: every feasible value this
#'   observation's own count could take.
#' @param xj_i `xj[i]`, the interest column's (rounded, integer) value
#'   for this observation.
#' @param log_offset_i `log_offset[i]` for this observation.
#' @param xc_i Numeric vector of length `p_cond`: this observation's own
#'   nuisance-column values (`xc[i, ]`), or `numeric(0)` if `p_cond == 0`.
#' @param lgam Precomputed `lgamma(0:ms + 1)`, indexed as
#'   `lgam[yy_rep + 1L]` to get `lgamma(yy + 1)` for every candidate `yy`
#'   at once.
#' @param cond_cols Character vector of nuisance-column names in `state`
#'   (`make_cond_cols(p_cond)`).
#' @inheritParams .cumulative_from_end
#' @param ms This stratum's observed total (`sum(y)`); states with
#'   `m_j > ms` are pruned.
#' @param tn_obs Integer vector of observed nuisance totals (length
#'   `p_cond`), or `NULL` to skip nuisance-coordinate pruning.
#' @return The new frontier after this observation: a named list in the
#'   same shape as `state`, one row per surviving, merged state.
#' @noRd
.dp_expand_observation <- function(
    state,
    yy_vec,
    xj_i,
    log_offset_i,
    xc_i,
    lgam,
    cond_cols,
    p_cond,
    ms,
    tn_obs
) {
    n_yy <- length(yy_vec)
    n_state <- length(state$m_j)

    # Every existing state crossed with every candidate value of y_i, in
    # one shot: `idx_rep` repeats each existing state n_yy times
    # consecutively, `yy_rep` cycles the candidate values once per
    # original state, so position-wise they enumerate the full (state,
    # candidate) cross product, exactly what a nested (m0, yy) loop
    # would do, just built as one vector replication instead of many.
    idx_rep <- rep(seq_len(n_state), each = n_yy)
    yy_rep <- rep.int(yy_vec, n_state)

    new_state <- list(
        m_j = state$m_j[idx_rep] + yy_rep,
        t_j = state$t_j[idx_rep] + yy_rep * xj_i,
        lw = state$lw[idx_rep] + yy_rep * log_offset_i - lgam[yy_rep + 1L]
    )
    if (p_cond > 0L) {
        for (k in seq_len(p_cond)) {
            xik <- xc_i[k]
            new_state[[cond_cols[k]]] <- if (xik != 0L) {
                state[[cond_cols[k]]][idx_rep] + yy_rep * xik
            } else {
                state[[cond_cols[k]]][idx_rep]
            }
        }
    }

    # Rule A (drop states that have already overshot either the stratum
    # total or the observed nuisance total. It is safe because every
    # increment is non-negative, so overshooting can never be corrected
    # by a later observation).
    keep <- new_state$m_j <= ms
    if (p_cond > 0L && !is.null(tn_obs)) {
        for (k in seq_len(p_cond)) {
            keep <- keep & (new_state[[cond_cols[k]]] <= tn_obs[k])
        }
    }
    if (!any(keep)) {
        stop(
            "All states pruned mid-stratum (check tn_obs feasibility).",
            call. = FALSE
        )
    }
    if (!all(keep)) {
        new_state <- .tbl_subset(new_state, keep)
    }

    # One grouping call per observation (not one per (m0, yy) pair).
    log_sum_exp_by(
        new_state,
        by_cols = c("m_j", "t_j", cond_cols),
        lw_col = "lw"
    )
}

#' Extract the final joint table from a completed DP frontier
#'
#' The last step of `.joint_dist_stratum_dp_vec()`: after every
#' observation has been folded in, only the states whose running total
#' lands exactly on the stratum's observed total are feasible.
#'
#' @param state The frontier after the last observation.
#' @param ms This stratum's observed total (`sum(y)`).
#' @param cond_cols Character vector of nuisance-column names in `state`.
#' @return A named list of equal-length vectors with elements `t_j`, one
#'   `cond_1`, ..., `cond_p` element per nuisance coordinate, and `lw`.
#' @noRd
.finalize_stratum_dp <- function(state, ms, cond_cols) {
    keep_final <- state$m_j == ms
    if (!any(keep_final)) {
        stop("No feasible states at m == ms.", call. = FALSE)
    }
    .tbl_subset(state[c("t_j", cond_cols, "lw")], keep_final)
}

#' Per-stratum dynamic program for the multivariate-shift algorithm
#'
#' Builds one stratum's own joint coefficient table `C_s(t_s, c_s |
#' M_s^obs)` (Section 5.3 of the theory vignette) by processing that
#' stratum's observations one at a time: at each observation, every
#' currently live state is crossed with every feasible value the new
#' observation's own count could take (vectorized as a single
#' row-replication instead of a nested loop), the interest statistic
#' `t_j`, each nuisance statistic (`cond_1`, ..., `cond_p`), and the
#' running within-stratum total `m_j` are updated accordingly, and the
#' base-measure log-weight is incremented by `yy * log_offset[i] -
#' lgamma(yy + 1)`. States sharing the same key are merged by
#' log-sum-exp after each observation (one grouping call per
#' observation, not one per candidate value). Rule A pruning (dropping
#' states that have already overshot the stratum total or, if `tn_obs`
#' is supplied, any nuisance coordinate's observed total) is applied
#' after every observation, since non-negative increments mean an
#' overshoot can never be corrected by a later observation. Only the
#' states that land exactly on the stratum's own observed total
#' (`m_j == sum(y)`) are kept in the returned table.
#'
#' This function is a thin orchestrator over five helpers
#' (`.validate_stratum_dp_inputs()`, `.zero_total_stratum_table()`,
#' `.prepare_stratum_dp_inputs()`, `.dp_expand_observation()`,
#' `.finalize_stratum_dp()`).
#'
#' @param y Integer vector of observed counts for this stratum's
#'   observations.
#' @param xj Numeric vector of the interest column's values (the model
#'   effect under test) for this stratum's observations; rounded to the
#'   nearest integer internally.
#' @param xc Numeric matrix of nuisance-column values for this stratum's
#'   observations (`nrow(xc) == length(y)`), or `NULL` if there are no
#'   nuisance columns to condition on (`p_cond == 0`).
#' @param log_offset Numeric vector of log-offsets (log exposure) for
#'   this stratum's observations.
#' @param si Integer vector of per-observation caps `S_i` for this
#'   stratum's observations (already resolved to whichever SAS cap rule
#'   applies; see `.stratum_slice()`).
#' @param tn_obs Integer vector of the observed nuisance totals (length
#'   `p_cond`), used for Rule A pruning against the global nuisance
#'   target as states are generated; or `NULL` to skip that pruning
#'   (still pruning against the stratum's own total).
#' @return A named list of equal-length vectors with elements `t_j`, one
#'   `cond_1`, ..., `cond_p` element per nuisance coordinate (if
#'   `p_cond > 0`), and `lw` (log-weight): this stratum's own joint
#'   coefficient table, restricted
#'   to the states whose running total equals `sum(y)`.
#' @noRd
.joint_dist_stratum_dp_vec <- function(
    y,
    xj,
    xc,
    log_offset,
    si,
    tn_obs = NULL
) {
    y <- as.integer(y)
    ms <- as.integer(sum(y))
    n <- length(y)

    .validate_stratum_dp_inputs(y, log_offset, si)

    if (ms == 0L) {
        return(.zero_total_stratum_table(xc))
    }

    prepped <- .prepare_stratum_dp_inputs(y, xj, xc, si, ms)
    si <- prepped$si
    xj <- prepped$xj
    xc <- prepped$xc
    p_cond <- prepped$p_cond
    cond_cols <- prepped$cond_cols

    if (!is.null(tn_obs)) {
        stopifnot(length(tn_obs) == p_cond)
        tn_obs <- as.integer(tn_obs)
    }

    lgam <- lgamma(0L:ms + 1L)

    # One combined frontier for the whole stratum, instead of a list of
    # per-bucket tables; `m_j` (the running within-stratum total) is just
    # another element, filtered to `m_j == ms` at the very end. Represented
    # as a plain named list of atomic vectors (not a data.frame): every
    # operation below (row replication via `[idx_rep]`, in-place-style
    # updates, logical filtering) is then a base-vector operation.
    state <- list(m_j = 0L, t_j = 0L)
    for (cc in cond_cols) {
        state[[cc]] <- 0L
    }
    state[["lw"]] <- 0.0

    for (i in seq_len(n)) {
        yy_vec <- 0L:si[i]
        xc_i <- if (p_cond > 0L) xc[i, ] else numeric(0L)
        state <- .dp_expand_observation(
            state = state,
            yy_vec = yy_vec,
            xj_i = xj[i],
            log_offset_i = log_offset[i],
            xc_i = xc_i,
            lgam = lgam,
            cond_cols = cond_cols,
            p_cond = p_cond,
            ms = ms,
            tn_obs = tn_obs
        )
    }

    .finalize_stratum_dp(state, ms, cond_cols)
}

#' Order already-built per-stratum tables by ascending row count
#'
#' Used by `ep_multi_shift()` which build every
#' stratum's own coefficient table before combining any of them,
#' regardless of processing order. So, unlike
#' `.order_strata_by_size()` above, an exact measurement of each
#' stratum's contribution to combination cost (its table's own row
#' count) is available at no extra cost, and is a strictly better
#' ordering signal than the raw stratum total: two strata with the same
#' observed total can still produce very differently sized tables,
#' depending on the number of nuisance columns, the number of
#' observations in the stratum, and the pattern of per-observation caps.
#' Combining the smallest table into the running accumulator first keeps
#' the accumulator itself as small as possible for as long as possible,
#' and makes the Rule B bound tightest early.
#'
#' @param joint_list A list of already-built per-stratum tables.
#' @param max_add_list The parallel list of each stratum's Rule B bound
#'   vector (Section 5.5 of the theory document), reordered alongside
#'   `joint_list`.
#' @return A list with `joint_list` and `max_add_list`, both reordered by
#'   ascending `nrow()` of `joint_list`'s original entries.
#' @noRd
.order_joint_tables_by_size <- function(joint_list, max_add_list) {
    if (length(joint_list) <= 1L) {
        return(
            list(
                joint_list = joint_list,
                max_add_list = max_add_list
            )
        )
    }
    ord_s <- order(vapply(joint_list, .tbl_nrow, integer(1L)))
    list(
        joint_list = joint_list[ord_s],
        max_add_list = max_add_list[ord_s]
    )
}

#' Extract one stratum's data for one requested parameter
#'
#' Slices out the observations belonging to one stratum, resolves their
#' per-observation caps under whichever SAS rule applies, and computes
#' the stratum's own forced-maximum contribution to each nuisance
#' coordinate (used by every algorithm's Rule B bound). This covers
#' everything ep_multi_shift() need per stratum;
#' ep_network() extends this with its own additional
#' backward-induction bounds (see `.network_stratum_slice()` in
#' ep-network.R).
#'
#' @param model The list returned by `.setup_exact_model()`.
#' @param level One value from `model$slevels0`/`model$sid`.
#' @param jcol Model-matrix column index of the parameter under test.
#' @param cond_cols_idx Model-matrix column indices being conditioned on.
#' @param p_cond `length(cond_cols_idx)`.
#' @param verbose If `TRUE`, prints this stratum's size and observed total.
#' @return A list with `idx`, `y_s`, `xj_s`, `off_s`, `si_s`, `xc_s`,
#'   `m_s_obs`, `n_s`, and `max_add_s`.
#' @noRd
.stratum_slice <- function(model, level, jcol, cond_cols_idx, p_cond, verbose) {
    idx <- which(model$sid == level)
    y_s <- model$y[idx]
    xj_s <- as.integer(round(model$x[idx, jcol]))
    off_s <- model$log_off[idx]

    si_s <- if (model$use_rule1) {
        model$exactmax_vec[idx]
    } else if (model$use_rule2) {
        si_from_stratum_total(y_s)
    } else {
        si_from_offset(off_s)
    }
    si_s <- as.integer(si_s)

    # Feasibility check, shared by every algorithm that calls this function
    # (ep_multi_shift() directly, ep_network() via
    # .network_stratum_slice()): every observation must satisfy 0 <= y_i <=
    # S_i under whichever SAS cap rule applies (Section 5.2 of the theory
    # vignette). This used to be checked only inside
    # .joint_dist_stratum_dp_vec(), which ep_multi_shift() calls
    # but ep_network() does not, so a stratum whose observed
    # count already exceeds its own cap (e.g. a near-zero offset under SAS
    # rule 3) silently produced a bogus fit instead of an error under
    # ep_network(). Checking here, in the one function both
    # algorithms call before doing anything else with y_s/si_s, catches it
    # for both, with an error that names the offending row(s) instead of the
    # low-level "Observed y_i exceeds S_i" check inside the dynamic program.
    if (any(y_s > si_s)) {
        bad <- which(y_s > si_s)
        stop(
            sprintf(
                "Observed count(s) exceed the per-observation cap S_i in stratum
                '%s': data row(s) %s have y_i > S_i (observed %s vs. cap %s).
                Every observation must satisfy 0 <= y_i <= S_i under whichever SAS
                EXACTMAX/OFFSET rule applies (Section 5.2 of the theory vignette);
                check `offset`/`exactmax` against the observed counts for this
                stratum.",
                level,
                toString(idx[bad]),
                toString(y_s[bad]),
                toString(si_s[bad])
            ),
            call. = FALSE
        )
    }

    xc_s <- if (p_cond > 0L) {
        as.matrix(model$x[idx, cond_cols_idx, drop = FALSE])
    } else {
        NULL
    }
    m_s_obs <- as.integer(sum(y_s))

    if (verbose) {
        message(
            sprintf(
                "  stratum %s: n=%d, M=%d",
                level,
                length(idx),
                m_s_obs
            )
        )
    }

    max_add_s <- if (p_cond > 0L) {
        compute_stratum_max_add(y_s, xc_s, si_s)
    } else {
        integer(0L)
    }

    list(
        idx = idx,
        y_s = y_s,
        xj_s = xj_s,
        off_s = off_s,
        si_s = si_s,
        xc_s = xc_s,
        m_s_obs = m_s_obs,
        n_s = length(idx),
        max_add_s = max_add_s
    )
}

# ------------------------------------------------------------------------
# Vectorized chunked combination across strata (log-scale) ----
#
# Same recursion as convolve_joint_pnuis_chunked() in
# ep-multivariate-shift.R, including the same chunking
# strategy (chunk the larger table by `chunk_size` to bound peak memory)
# and the same Rule A/B pruning. The difference: within each chunk, the
# cross join against the smaller table (`d1`) is built as one vectorized
# row-replication instead of a loop over `seq_len(nrow(d1))`.
# ------------------------------------------------------------------------

#' Combine two strata's joint tables by chunked, vectorized convolution
#'
#' Computes the pairwise convolution of two joint coefficient tables
#' (Section 5.4 of the theory vignette: statistics add and log-weights
#' add across independent strata), processing the larger of the two
#' tables in row-chunks of at most `chunk_size` to bound peak memory.
#' Within each chunk, the cross join against the smaller table is built
#' as a single vectorized row-replication (rather than a loop over the
#' smaller table's rows), duplicate resulting keys are merged by
#' log-sum-exp, and Rule A (cannot exceed the observed nuisance total)
#' and Rule B (cannot possibly still reach it, given `remaining_max`)
#' pruning are applied to the whole chunk's cross join at once. A chunk
#' that happens to produce no surviving rows after pruning is skipped
#' without discarding results already accumulated from other chunks.
#'
#' @param d1,d2 The two joint tables to combine (named lists of
#'   equal-length vectors with elements `t_j`, `cond_1`, ..., `cond_p`,
#'   and `lw`); the smaller of the two (`d1`, `d2`, whichever has fewer
#'   rows) is swapped in first for the outer loop.
#' @param p_cond Number of nuisance coordinates (`cond_*` elements).
#' @param tn_obs Integer vector of the observed nuisance totals, length
#'   `p_cond`, used by Rule A pruning.
#' @param remaining_max Integer vector, length `p_cond`: the Rule B bound,
#'   the most every stratum still to be combined after this step
#'   could still contribute to each nuisance coordinate (see
#'   `.cumulative_from_end()`).
#' @param chunk_size Maximum number of rows of the larger table processed
#'   in one chunk, bounding how many `(d1-row, chunk-row)` pairs are
#'   materialized at once.
#' @return A named list of equal-length vectors with elements `t_j`,
#'   `cond_1`, ..., `cond_p`, and `lw`: the combined, pruned, and merged
#'   joint table, ordered by `t_j` (and the `cond_*` elements); a table
#'   with zero-length vectors if every
#'   chunk was pruned to zero rows.
#' @noRd
.convolve_joint_pnuis_vec <- function(
    d1,
    d2,
    p_cond,
    tn_obs,
    remaining_max,
    chunk_size = 4000L
) {
    cond_cols <- make_cond_cols(p_cond)

    if (.tbl_nrow(d2) < .tbl_nrow(d1)) {
        tmp <- d1
        d1 <- d2
        d2 <- tmp
    }
    n1 <- .tbl_nrow(d1)
    n2 <- .tbl_nrow(d2)

    out_all <- NULL
    for (j0 in seq(1L, n2, by = chunk_size)) {
        j1 <- min(n2, j0 + chunk_size - 1L)
        b <- .tbl_subset(d2, j0:j1)
        nb <- .tbl_nrow(b)

        # Full cross join of this chunk against every row of d1, in one shot:
        # idx1 repeats each d1 row nb times consecutively, idxb cycles the
        # chunk's own row order once per d1 row, together they enumerate
        # every (d1-row, chunk-row) pair exactly once, matching what a
        # `for (i in seq_len(nrow(d1)))` loop would do one row at a time.
        idx1 <- rep(seq_len(n1), each = nb)
        idxb <- rep(seq_len(nb), times = n1)

        tmp <- list(
            t_j = d1$t_j[idx1] + b$t_j[idxb],
            lw = d1$lw[idx1] + b$lw[idxb]
        )

        if (p_cond > 0L) {
            for (k in seq_len(p_cond)) {
                cc <- cond_cols[k]
                tmp[[cc]] <- d1[[cc]][idx1] + b[[cc]][idxb]
            }
            keep <- rep(TRUE, .tbl_nrow(tmp))
            for (k in seq_len(p_cond)) {
                cc <- cond_cols[k]
                # Rule A (cannot exceed the observed nuisance total) and Rule B
                # (cannot possibly still reach it, given what remains), applied
                # as one vectorized filter across the whole chunk's cross join.
                keep <- keep &
                    (tmp[[cc]] <= tn_obs[k]) &
                    ((tmp[[cc]] + remaining_max[k]) >= tn_obs[k])
            }
            if (!all(keep)) {
                tmp <- .tbl_subset(tmp, keep)
            }
        }

        if (.tbl_nrow(tmp) == 0L) {
            # A chunk of d2 can legitimately produce no surviving pairs after
            # Rule A/B pruning without that saying anything about the other
            # chunks; skip it and keep accumulating rather than discarding
            # out_all or stopping early.
            next
        }

        tmp <- log_sum_exp_by(tmp, by_cols = c("t_j", cond_cols), lw_col = "lw")
        out_all <- if (is.null(out_all)) {
            tmp
        } else {
            log_sum_exp_by(
                .tbl_rbind(out_all, tmp),
                by_cols = c("t_j", cond_cols),
                lw_col = "lw"
            )
        }
    }

    if (is.null(out_all) || .tbl_nrow(out_all) == 0L) {
        empty <- list(t_j = integer(0L))
        for (cc in cond_cols) {
            empty[[cc]] <- integer(0L)
        }
        empty[["lw"]] <- numeric(0.0)
        return(empty)
    }
    ord <- do.call(order, out_all[c("t_j", cond_cols)])
    .tbl_subset(out_all, ord)
}

# ==============================================================================
# Per-parameter worker ----
# ==============================================================================

#' One exact-effect column's full multivariate-shift computation
#'
#' Everything `ep_multi_shift()`'s per-parameter loop does for
#' one requested column `jcol`: resolve which columns are conditioned on,
#' build every stratum's own table, order and combine them into the
#' joint one-dimensional distribution, and assemble the final output row.
#' Split out of the main function so its own body reads as "for each
#' requested column, compute its result," with the actual computation
#' isolated and independently readable (and testable) here.
#'
#' @inheritParams ep_multi_shift
#' @inheritParams .stratum_slice
#' @return A list with `parm` (the parameter name), `outdist_entry`, and
#'   `result_row`, ready to be stored under `parm` in the main function's
#'   `outdist`/`results` lists.
#' @noRd
.multishift_one_parameter <- function(
    model,
    jcol,
    condition_on,
    chunk_size,
    alpha,
    ci_type,
    verbose
) {
    parm <- model$cn[jcol]
    if (verbose) {
        message(
            sprintf(
                "\n=== EXACT for parameter column: %s ===",
                parm
            )
        )
    }

    cond_cols_idx <- if (condition_on == "model") {
        setdiff(seq_len(ncol(model$x)), jcol)
    } else {
        setdiff(model$exact_cols, jcol)
    }
    p_cond <- length(cond_cols_idx)
    tn_obs <- if (p_cond > 0L) {
        as.integer(crossprod(model$x[, cond_cols_idx, drop = FALSE], model$y))
    } else {
        integer(0L)
    }

    # Build every stratum's own table first, in whatever order the strata
    # naturally come in, this cost is paid regardless of processing
    # order, so it is not worth reordering before doing it.
    joint_list <- vector("list", length(model$slevels0))
    max_add_list <- vector("list", length(model$slevels0))

    for (s in seq_along(model$slevels0)) {
        info <- .stratum_slice(
            model,
            model$slevels0[s],
            jcol,
            cond_cols_idx,
            p_cond,
            verbose
        )
        max_add_list[[s]] <- info$max_add_s

        joint_list[[s]] <- .joint_dist_stratum_dp_vec(
            y = info$y_s,
            xj = info$xj_s,
            xc = info$xc_s,
            log_offset = info$off_s,
            si = info$si_s,
            tn_obs = if (p_cond > 0L) tn_obs else NULL
        )
    }

    # Now that every table is built, its exact row count is known for
    # free; combine the smallest tables first (see
    # .order_joint_tables_by_size()'s own documentation for why this is a
    # strictly better signal than any proxy computed before building).
    ord <- .order_joint_tables_by_size(joint_list, max_add_list)
    joint_list <- ord$joint_list
    max_add_list <- ord$max_add_list

    rem_after <- .cumulative_from_end(max_add_list, p_cond)

    joint <- joint_list[[1L]]
    if (p_cond > 0L) {
        cond_names <- make_cond_cols(p_cond)
        keep <- rep(TRUE, .tbl_nrow(joint))
        for (k in seq_len(p_cond)) {
            keep <- keep &
                ((joint[[cond_names[k]]] + rem_after[[1L]][k]) >= tn_obs[k])
        }
        if (!all(keep)) {
            joint <- .tbl_subset(joint, keep)
        }
        if (.tbl_nrow(joint) == 0L) {
            stop(
                "All states pruned after first stratum (unexpected).",
                call. = FALSE
            )
        }
    }

    if (length(joint_list) > 1L) {
        for (s in 2L:length(joint_list)) {
            joint <- .convolve_joint_pnuis_vec(
                d1 = joint,
                d2 = joint_list[[s]],
                p_cond = p_cond,
                tn_obs = tn_obs,
                remaining_max = rem_after[[s]],
                chunk_size = chunk_size
            )
            if (verbose) {
                message(
                    sprintf(
                        "  after combining %d/%d strata: joint rows = %d",
                        s,
                        length(joint_list),
                        .tbl_nrow(joint)
                    )
                )
            }
            if (.tbl_nrow(joint) == 0L) {
                stop(
                    sprintf(
                        "All states pruned after combining %d strata for %s.",
                        s,
                        parm
                    ),
                    call. = FALSE
                )
            }
        }
    }

    dist1d <- filter_and_marginalize_pn_lw(joint, tn_obs)
    out <- .assemble_parameter_output(
        dist1d,
        parm,
        model$y,
        model$x,
        jcol,
        alpha,
        ci_type
    )
    list(
        parm = parm,
        outdist_entry = out$outdist_entry,
        result_row = out$result_row
    )
}

# ==============================================================================
# Main function ----
# ==============================================================================

#' Exact conditional Poisson regression (multivariate-shift algorithm)
#'
#' Exact conditional inference for one or more categorical effects in a
#' Poisson regression with a log link, an optional offset, and optional
#' stratification. Nuisance parameters (stratum intercepts and every model
#' effect not listed in `exact_vars`) are eliminated by conditioning on
#' their sufficient statistics, and the resulting exact conditional
#' distribution of each requested effect's sufficient statistic is
#' obtained by the multivariate-shift algorithm: a per-stratum dynamic
#' program combined with explicit cross-stratum convolution, generalizing
#' the multivariate shift algorithm of Hirji, Mehta, and Patel (1987) for
#' exact logistic regression to bounded Poisson counts, following the
#' `EXACTMAX`/`OFFSET` rules documented for SAS `PROC GENMOD`'s `EXACT`
#' statement.
#'
#' `formula` may contain interaction terms (`a*b`, `a:b`, or higher-order
#' interactions) freely, in two distinct roles. An interaction can be a
#' plain adjustment covariate, left out of `exact_vars` entirely, in
#' which case it is conditioned on like any other nuisance term (under
#' the default `condition_on = "model"`) but never itself tested. Or an
#' interaction can be the effect of interest, including an interaction
#' of the effect of interest with another covariate, by listing it as
#' its own entry in `exact_vars`; see that argument's documentation for
#' the exact naming rules. The two roles are independent: listing a main
#' effect in `exact_vars` never implicitly tests its interactions, and
#' listing an interaction never implicitly tests the main effects it is
#' built from, each term to be tested must be named explicitly.
#'
#' @param formula A model formula.
#' @param data A data frame.
#' @param offset Name of an offset variable in `data`, or `NULL` (default).
#'   If `NULL`, every observation is treated as contributing one unit of
#'   exposure (offset = 0 on the log scale).
#' @param offset_log If `TRUE` (default), `offset` already holds
#'   log(exposure); if `FALSE`, it holds raw exposure.
#' @param exactmax `NULL` (not specified), `TRUE` (specified without a
#'   variable), or a character variable name (SAS `EXACTMAX=`).
#' @param strata Optional character vector of stratification variable
#'   names in `data`.
#' @param exact_vars Character vector of model *terms* to compute exact
#'   results for. Each entry must name a term that appears in
#'   `formula`, exactly as SAS's `EXACT` statement requires. An entry can
#'   be a main effect (`"region"`) or an interaction of any order, named
#'   by joining every variable in it with `*` or `:` (`"region:gender"`
#'   and `"region*gender"` both work, and the variables may be listed in
#'   either order). Interaction columns are matched via the model's own
#'   `terms()`/`assign` structure, not by name-guessing, so a term is
#'   never confused with a different term that happens to share a text
#'   prefix (e.g. requesting `"region"` in a model that also contains
#'   `region:gender` selects only region's own main-effect columns).
#'   Requesting a term does *not* automatically request any other term
#'   built from the same variable(s): a main effect and its interactions
#'   are entirely independent entries, so testing `region` while leaving
#'   `region:gender` untested keeps the interaction as an ordinary
#'   nuisance covariate (conditioned on under `condition_on = "model"`,
#'   never tested); to test the interaction itself, add it as its own
#'   entry, e.g. `exact_vars = c("region", "region:gender")`.
#' @param ref_levels Optional named list of reference levels for
#'   categorical variables.
#' @param alpha Significance level for confidence intervals.
#' @param ci_type `"two-sided"` (default) or `"one-sided"`, matching
#'   [ep_closed_form()]'s argument of the same name. For `"one-sided"`,
#'   only the confidence limit on the informative side is computed for
#'   each parameter (the other is left at its boundary, `-Inf` or `Inf`);
#'   which side is informative is determined per parameter from the two
#'   one-sided tail probabilities at the null (see
#'   `two_sided_from_tails()`'s `direction`), exactly as in
#'   [ep_closed_form()].
#' @param chunk_size Maximum number of rows processed at once while
#'   combining two strata's coefficient tables. Combining strata requires
#'   pairing every row of one stratum's table with every row of the
#'   other's; `chunk_size` bounds how many such pairs are materialized
#'   simultaneously, trading memory for speed. The default is generous
#'   for typical problems; lower it if combining strata runs out of
#'   memory, or raise it if memory is plentiful and speed matters more.
#' @param condition_on `"model"` (default; condition on every other model
#'   effect, matching SAS) or `"exact"` (condition only on the other
#'   requested `exact_vars`).
#' @param verbose If `TRUE` (default), prints progress and which SAS rule
#'   determined the per-observation caps.
#'
#' @return An object of class `c("ep_multi_shift", "ep_fit")`: a list with
#'   elements `exactparmest` (a one-row data frame with `Parameter`,
#'   `Estimate` (the coefficient beta = log rate ratio, not the rate
#'   ratio), `StdError`, `Alpha`, `LowerCL`, `UpperCL` (beta scale),
#'   `PValue`, and `SuffStat`), `OUTDIST` (the exact distribution of the
#'   sufficient statistic), and `call`.
#'
#' @examples
#' df <- data.frame(
#'   group = c(1L, 1L, 1L, 0L, 0L, 0L),
#'   age = c(50L, 60L, 70L, 50L, 60L, 70L),
#'   subjects = c(3492, 2141, 1711, 3525, 2166, 1724),
#'   cases = c(3, 2, 1, 10, 7, 5),
#'   futime = c(50, 60, 40, 51, 62, 42)
#' )
#'
#' df$lnfutime <- log(df$futime)
#'
#' ep_multi_shift(
#'   formula = cases ~ group,
#'   data = df,
#'   offset = "lnfutime",
#'   offset_log = TRUE,
#'   strata = "age",
#'   exact_vars = "group",
#'   ref_levels = list(group="0", age="50"),
#'   alpha = 0.05,
#'   verbose = TRUE,
#'   ci_type = "two-sided"
#' )
#'
#' @references Hirji KF, Mehta CR, Patel NR. Computing Distributions for
#'   Exact Logistic Regression. JASA, 82(400):1110-1117, 1987.
#' @references Mehta CR, Patel NR. Exact logistic regression: theory and
#'   examples. Statistics in Medicine, 14(19):2143-2160, 1995.
#' @export
ep_multi_shift <- function(
    formula,
    data,
    exact_vars,
    offset = NULL,
    offset_log = TRUE,
    exactmax = NULL,
    strata = NULL,
    ref_levels = list(),
    alpha = 0.05,
    ci_type = c("two-sided", "one-sided"),
    chunk_size = 4000L,
    condition_on = c("model", "exact"),
    verbose = TRUE
) {
    condition_on <- rlang::arg_match(condition_on)
    ci_type <- rlang::arg_match(ci_type)
    .validate_exact_poisson_args(
        formula = formula,
        data = data,
        offset = offset,
        offset_log = offset_log,
        exactmax = exactmax,
        strata = strata,
        exact_vars = exact_vars,
        ref_levels = ref_levels,
        alpha = alpha,
        verbose = verbose,
        chunk_size = chunk_size
    )

    model <- .setup_exact_model(
        formula = formula,
        data = data,
        offset = offset,
        offset_log = offset_log,
        exactmax = exactmax,
        strata = strata,
        exact_vars = exact_vars,
        ref_levels = ref_levels,
        verbose = verbose
    )

    results <- list()
    outdist <- list()

    for (jcol in model$exact_cols) {
        one <- .multishift_one_parameter(
            model,
            jcol,
            condition_on,
            chunk_size,
            alpha,
            ci_type,
            verbose
        )
        outdist[[one$parm]] <- one$outdist_entry
        results[[one$parm]] <- one$result_row
    }

    exactparmest <- do.call(rbind, results)
    rownames(exactparmest) <- NULL
    ep_fit(
        exactparmest = exactparmest,
        outdist = outdist,
        call = match.call(),
        subclass = "ep_multi_shift",
        algorithm_label = paste(
            "Exact conditional Poisson (categorical + STRATA;",
            "multivariate-shift algorithm)"
        )
    )
}
