# ==============================================================================
# ep_network(): a generalization, to Poisson counts, of the
# network algorithm of Mehta, Patel, and Senchaudhuri (1992), "Exact
# Stratified Linear Rank Tests for Ordered Categorical and Binary Data,"
# Journal of Computational and Graphical Statistics, 1(1):21-40.
# This method builds a network for each parameter
# that you are conditioning out, combines the networks, then uses the
# multivariate shift algorithm to create the exact distribution.").
# Hirji, Mehta, and Patel (1987) is a different paper, for the
# multivariate shift algorithm itself (SAS's own EXACTOPTIONS
# METHOD=DIRECT), see ep-multivariate-shift.R, which generalize
# that algorithm to Poisson counts. The network method reuses the
# multivariate shift algorithm as an internal subroutine but organizes
# the computation completely differently, as a network (graph) of
# nodes and arcs per stratum, with a backward induction pass over that
# network computed *before* the forward pass, so that infeasible partial
# states can be discarded the instant they are generated, mid-stratum,
# rather than only at a stratum boundary or after an explicit table
# cross-join. This file implements that; see below for exactly how the
# paper's formalism (nodes, arcs, LP/SP, backward induction, forward
# record elimination) maps onto Poisson regression with categorical
# covariates.
#
# THE PAPER'S NETWORK, TRANSLATED TO THIS PROBLEM
# --------------------------------------------------
# The paper builds one network per stratum k, with stages 0,1,...,c
# (c = number of response categories) and nodes (j, m_jk) = "the first j
# columns of the table have used up m_jk of the stratum's row total."
# Arcs from (j, m_jk) to (j+1, m_{j+1,k}) correspond to every feasible
# value the (j+1)-th column's count could take. A path from the start
# node to the terminal node is exactly one contingency table, and the
# node/arc structure lets every feasible table be enumerated *implicitly*
# via a forward recursion, instead of one at a time.
#
# For Poisson regression with per-observation caps S_i (Section 5.2 of
# the theory vignette), the direct analog is: stage i = "after processing
# observation i of the stratum" (i = 0,...,n_s), node (i, m) = "the first
# i observations have contributed a total of m toward the stratum's
# required total M_s^obs," and the arc out of (i, m) for observation
# i+1's candidate value yy in {0,...,S_{i+1}} leads to node (i+1, m+yy).
# This is exactly the per-observation forward recursion already used by
# ep-multivariate-shift.R, the
# network's nodes are precisely that recursion's DP buckets. The piece
# this file adds is the paper's OTHER half: a **backward induction** pass
# over the same network (Section 3.1.1, "Backward Induction on Gamma_k"),
# computed once per stratum before the forward pass, giving, at every
# node, the exact longest-path (LP) and shortest-path (SP) values, the
# maximum and minimum a given quantity could still become, using only the
# observations not yet processed. The paper's forward algorithm (its
# 5-step procedure) then discards a partial record the instant its key,
# plus the node's own LP, cannot reach the target ("Step 4.a": "If
# u* + LP(j+1,m_{j+1,k}) < a_k, drop this record"). That is a continuous,
# per-observation, per-node prune.
#
# THE CLOSED FORM FOR LP/SP, AND WHY IT IS EXACT (NOT AN APPROXIMATION)
# --------------------------------------------------------------------
# A full backward induction, as in the paper, would tabulate LP/SP by
# dynamic programming, one entry per node. That is unnecessary extra
# work here because every nuisance (and interest) column in this package
# is a 0/1 dummy (Section 2.3 of the theory vignette: "categorical
# effects only"). For a 0/1 column k, splitting the not-yet-processed
# observations of a stratum into an "active" group (x_{i,k} = 1) with
# total remaining capacity suffix1_k, and an "inactive" group
# (x_{i,k} = 0) with total remaining capacity suffix0_k, allocating an
# exact remaining budget of r more Poisson counts among them (r = how
# much of the stratum's M_s^obs is still unaccounted for) can push
# coordinate k's own total anywhere in the closed interval
#     [ max(0, r - suffix0_k),  min(r, suffix1_k) ],
# and every value in that interval is achievable (load the "active"
# group first to hit the upper end, the "inactive" group first to hit
# the lower end; anything in between follows from swapping one unit of
# budget at a time between an active and an inactive observation with
# unused capacity on each side, which always exists strictly inside the
# interval). This is the exact solution to the paper's backward
# recursion (3.7) for this restricted (0/1-coefficient) case, not a
# looser stand-in for it.
#
# WHAT THIS FILE DOES, END TO END
# -----------------------------------
# For each stratum, in increasing order of its own observed total (the
# ordering that keeps the number of live states, and hence the value of
# the bound below, as small as possible for as long as possible):
# precompute, per nuisance dimension k, the suffix capacities
# suffix1_k(i)/suffix0_k(i) for every position i in the stratum (a plain
# backward cumulative sum, the "backward induction" pass, done once,
# before that stratum's own forward pass), and the whole-stratum forced
# maximum/minimum contribution to k (`compute_stratum_max_add()`, shared
# via ep-common.R, and this file's own
# `.compute_stratum_min_add()`). Then run one continuous forward pass
# over every observation, in every stratum, on a single running frontier
# of live states (no separate per-stratum table, no explicit cross-join
# against another stratum's table): at each observation, expand the live
# states by that observation's own candidate values (vectorized), and
# immediately drop any resulting state whose nuisance keys either (a)
# already exceed the observed total (Rule A), (b) cannot possibly reach
# it even in the best remaining case, current value plus this
# stratum's own remaining LP plus every future stratum's maximum
# (Rule B, computed from the exact closed form), or (c) cannot possibly
# avoid overshooting it, current value plus this stratum's own
# remaining SP plus every future stratum's forced minimum. At the end of
# each stratum's own observations, states are filtered to the ones that
# hit that stratum's required total exactly (eliminating its intercept,
# matching the paper's own per-stratum network), and the within-stratum
# counter resets before moving into the next stratum's observations on
# the same frontier.
#
# This file does not materialize an explicit graph object (nodes and
# arcs as R objects); the "network" here is realized implicitly, the
# same way joint_dist_stratum_dp_vec() already realizes it, via a plain
# named list of equal-length vectors holding live (key, log-weight)
# records (see ep-utils.R's .tbl_nrow()/.tbl_subset()).
# Building literal node/arc objects would cost more in R than it could
# ever save; what matters for both correctness and speed is the
# *backward-induction-derived bound* and *continuous forward pruning*,
# both of which are fully present below.
#
# This file adds its own stratum-level and per-observation helpers, and
# reuses the shared setup/ordering/output-assembly helpers defined in
# ep-common.R (see that file's "Shared orchestration helpers"
# section) so that its main function expresses only what is actually
# specific to the network algorithm: the LP/SP backward induction and the
# continuous forward pass.
# ==============================================================================

# ------------------------------------------------------------------------
# Per-stratum backward induction (closed form) ----
# for each position r in a
# stratum (1-indexed, with an extra row for "after the last observation"),
# the exact longest-path (LP-supporting) and shortest-path (SP-supporting)
# remaining capacities for each nuisance coordinate, i.e. splitting
# observations r, r+1, ..., n_s by whether their own coefficient for that
# coordinate is 0 or 1, and summing each group's own S_i capacity
# (`suffix1[r, k]`, `suffix0[r, k]`). `suffix1[n_s + 1, ]` and
# `suffix0[n_s + 1, ]` are defined as all zero (nothing left once every
# observation has been processed). See the file header for why these two
# suffix sums are the *exact* solution to the paper's backward-induction
# recursion for 0/1-coefficient columns, not merely a bound on it.
# ------------------------------------------------------------------------

#' Per-position backward-induction suffix capacities (LP/SP bounds)
#'
#' For each position `r` in a stratum (1-indexed, with an extra row for
#' "after the last observation"), computes the exact longest-path- and
#' shortest-path-supporting remaining capacities for each nuisance
#' coordinate: splitting the not-yet-processed observations
#' `r, r+1, ..., n_s` by whether their own coefficient for that
#' coordinate is 0 or 1, and summing each group's own `S_i` capacity.
#' As the file header explains in detail, this backward cumulative sum
#' is the *exact* solution to the paper's backward-induction recursion
#' for 0/1-coefficient columns, computed once per stratum, not merely a
#' bound on it.
#'
#' @param xc_s Numeric matrix of nuisance-column values for this
#'   stratum's observations (`nrow(xc_s) == length(si_s)`).
#' @param si_s Integer vector of per-observation caps `S_i` for this
#'   stratum's observations.
#' @param p_cond Number of nuisance coordinates (`ncol(xc_s)`).
#' @return A list with `suffix1` and `suffix0`, each an
#'   `(n_s + 1) x p_cond` integer matrix: `suffix1[r, k]`/`suffix0[r, k]`
#'   give the total remaining capacity, from position `r` onward, of this
#'   stratum's coordinate-`k`-active/-inactive observations respectively;
#'   row `n_s + 1` of each matrix is all zero, since nothing remains once
#'   every observation has been processed.
#' @noRd
.stratum_suffix_bounds <- function(xc_s, si_s, p_cond) {
    n_s <- length(si_s)
    suffix1 <- matrix(0L, nrow = n_s + 1L, ncol = p_cond)
    suffix0 <- matrix(0L, nrow = n_s + 1L, ncol = p_cond)
    if (p_cond > 0L && n_s > 0L) {
        xc_s <- as.matrix(xc_s)
        si_int <- as.integer(si_s)
        for (k in seq_len(p_cond)) {
            is_active <- xc_s[, k] != 0L
            contrib1 <- ifelse(is_active, si_int, 0L)
            contrib0 <- ifelse(is_active, 0L, si_int)
            suffix1[seq_len(n_s), k] <- rev(cumsum(rev(contrib1)))
            suffix0[seq_len(n_s), k] <- rev(cumsum(rev(contrib0)))
        }
    }
    list(suffix1 = suffix1, suffix0 = suffix0)
}

#' Order strata by ascending observed total (used only where no built
#' table is available to measure directly)
#'
#' `ep_network()` never materializes a separate per-stratum
#' table, avoiding exactly that is the point of its continuous-frontier
#' design (Section 7 of the theory document), so it must decide its
#' stratum processing order *before* doing any per-stratum work, from
#' whatever is cheaply available in the raw data. The observed stratum
#' total `M_s^obs` is that cheap proxy: smaller-total strata tend to
#' produce smaller live-state frontiers, so processing them first keeps
#' the frontier small for as long as possible and makes the Rule B bound
#' tightest early. This is only ever a proxy for the quantity that
#' actually matters (how many live states a stratum contributes); where
#' an exact measurement of that quantity is available for free instead
#' (`ep_multi_shift()`, which build every
#' per-stratum table regardless of order, see
#' `.order_joint_tables_by_size()` below), that exact measurement is used
#' instead of this proxy.
#'
#' @param slevels0 Distinct stratum IDs.
#' @param sid Stratum ID per observation.
#' @param y Response per observation.
#' @return `slevels0`, reordered by ascending stratum total.
#' @noRd
.order_strata_by_size <- function(slevels0, sid, y) {
    if (length(slevels0) <= 1L) {
        return(slevels0)
    }
    totals <- vapply(slevels0, function(lv) sum(y[sid == lv]), numeric(1L))
    slevels0[order(totals)]
}

# ------------------------------------------------------------------------
# Whole-stratum forced *minimum* contribution to each nuisance coordinate,
# the SP-side counterpart to `compute_stratum_max_add()` (defined in
# ep-multivariate-shift.R and reused here unchanged for the LP/
# forced-maximum side). If a stratum's own required total M_s^obs exceeds
# the total capacity of its coordinate-k-inactive observations, the
# residual is unavoidably forced onto coordinate-k-active observations,
# no matter how the stratum's total is ultimately allocated.
# ------------------------------------------------------------------------

#' One stratum's own forced-minimum contribution to each nuisance coordinate
#'
#' The SP-side counterpart to `compute_stratum_max_add()` (the LP-side,
#' forced-*maximum* bound defined in `ep-common.R` and reused
#' unchanged by this algorithm). Computes, for each nuisance column, the
#' least amount this stratum could possibly contribute: if the stratum's
#' own required total `M_s^obs` exceeds the total capacity of its
#' coordinate-`k`-inactive observations, the residual is unavoidably
#' forced onto its coordinate-`k`-active observations, no matter how the
#' stratum's total is ultimately allocated among its observations.
#'
#' @inheritParams compute_stratum_max_add
#' @return An integer vector of length `ncol(xc_s)`: this stratum's own
#'   forced-minimum contribution to each nuisance coordinate.
#' @noRd
.compute_stratum_min_add <- function(y_s, xc_s, si_s) {
    ms <- as.integer(sum(y_s))
    p_cond <- ncol(xc_s)
    out <- integer(p_cond)
    for (k in seq_len(p_cond)) {
        idx0 <- which(xc_s[, k] == 0.0)
        cap0 <- if (length(idx0) == 0L) 0L else sum(si_s[idx0])
        out[k] <- as.integer(max(0L, ms - cap0))
    }
    out
}

# ------------------------------------------------------------------------
# Network-specific per-stratum precomputation ----
#
# Extends the shared `.stratum_slice()` (ep-multivariate-
# shift.R) with the extra backward-induction bounds this algorithm needs
# on top of what ep_multi_shift() require:
# the forced-minimum contribution (SP side) and the per-position suffix
# capacities used to tighten the bound continuously within a stratum.
# ------------------------------------------------------------------------

#' Per-stratum precomputation for the network algorithm
#'
#' Extends the shared `.stratum_slice()` (`ep-multivariate-
#' shift.R`) with the extra backward-induction bounds the network
#' algorithm needs on top of what `ep_multi_shift()` requires:
#' the forced-minimum contribution (SP side, via
#' `.compute_stratum_min_add()`) and the per-position suffix capacities
#' (via `.stratum_suffix_bounds()`) used to tighten the pruning bound
#' continuously within a stratum's own forward pass.
#'
#' @inheritParams .stratum_slice
#' @return A list with every element of `.stratum_slice()`'s return value
#'   (`idx`, `y_s`, `xj_s`, `off_s`, `si_s`, `xc_s`, `m_s_obs`, `n_s`,
#'   `max_add_s`), plus `min_add_s` (this stratum's forced-minimum
#'   contribution to each nuisance coordinate) and `suffix1`/`suffix0`
#'   (the per-position backward-induction suffix capacities from
#'   `.stratum_suffix_bounds()`).
#' @noRd
.network_stratum_slice <- function(
    model,
    level,
    jcol,
    cond_cols_idx,
    p_cond,
    verbose
) {
    info <- .stratum_slice(model, level, jcol, cond_cols_idx, p_cond, verbose)
    info$min_add_s <- if (p_cond > 0L) {
        .compute_stratum_min_add(info$y_s, info$xc_s, info$si_s)
    } else {
        integer(0L)
    }
    suf <- .stratum_suffix_bounds(info$xc_s, info$si_s, p_cond)
    info$suffix1 <- suf$suffix1
    info$suffix0 <- suf$suffix0
    info
}

#' One observation's forward-pass update (network algorithm)
#'
#' Expands the current frontier `dt` by every candidate value of
#' observation `r` in stratum `info` (vectorized cross join, as in
#' `.joint_dist_stratum_dp_vec()`), then immediately applies the
#' continuous LP/SP-derived pruning described in the file header: a
#' resulting state survives only if it has not yet overshot the
#' stratum's own running total, cannot be ruled out from reaching
#' `tn_obs` even in the best remaining case (LP side), and cannot be
#' ruled out from avoiding an overshoot of `tn_obs` even in the best
#' remaining case (SP side).
#'
#' @param state The current frontier (a named list of equal-length
#'   vectors, see `.tbl_nrow()`/`.tbl_subset()` in
#'   `ep-utils.R`, of live states) before processing
#'   observation `r`.
#' @param r Index, within this stratum, of the observation to process
#'   next.
#' @param info This stratum's precomputed slice, as returned by
#'   `.network_stratum_slice()`.
#' @inheritParams .cumulative_from_end
#' @param cond_cols Character vector of nuisance-column names (from
#'   `make_cond_cols()`).
#' @param tn_obs Integer vector of the observed nuisance totals, length
#'   `p_cond`.
#' @param ms_here This stratum's own required total (`M_s^obs`).
#' @param rem_future Integer vector, length `p_cond`: the sum of
#'   `max_add_s` (LP side) over every stratum *after* this one in the
#'   processing order.
#' @param min_future Integer vector, length `p_cond`: the sum of
#'   `min_add_s` (SP side) over every stratum *after* this one in the
#'   processing order.
#' @return The filtered, expanded frontier (a named list of equal-length
#'   vectors), possibly with zero-length vectors, the caller checks for
#'   and reports that.
#' @noRd
.network_expand_observation <- function(
    state,
    r,
    info,
    p_cond,
    cond_cols,
    tn_obs,
    ms_here,
    rem_future,
    min_future
) {
    yy_vec <- 0L:min(info$si_s[r], ms_here)
    n_yy <- length(yy_vec)
    n_state <- length(state$m_j)

    idx_rep <- rep(seq_len(n_state), each = n_yy)
    yy_rep <- rep.int(yy_vec, n_state)

    new_state <- list(
        m_j = state$m_j[idx_rep] + yy_rep,
        t_j = state$t_j[idx_rep] + yy_rep * info$xj_s[r],
        lw = state$lw[idx_rep] + yy_rep * info$off_s[r] - lgamma(yy_rep + 1L)
    )

    if (p_cond > 0L) {
        for (k in seq_len(p_cond)) {
            xik <- info$xc_s[r, k]
            new_state[[cond_cols[k]]] <- if (xik != 0L) {
                state[[cond_cols[k]]][idx_rep] + yy_rep * xik
            } else {
                state[[cond_cols[k]]][idx_rep]
            }
        }
    }

    # Rule "count budget": cannot overshoot this stratum's own total.
    keep <- new_state$m_j <= ms_here

    if (p_cond > 0L) {
        suf1_next <- info$suffix1[r + 1L, ]
        suf0_next <- info$suffix0[r + 1L, ]
        budget <- ms_here - new_state$m_j

        for (k in seq_len(p_cond)) {
            cc <- cond_cols[k]
            # LP side (reachability from below): the most coordinate k can
            # still receive is bounded by (i) this state's own remaining count
            # budget and the raw capacity of this stratum's remaining
            # coordinate-k-active observations, whichever is smaller, plus
            # (ii) every future stratum's own forced maximum. If even this
            # best case falls short of tn_obs[k], the state is dead.
            bound_max_k <- pmin(budget, suf1_next[k]) + rem_future[k]
            # SP side (cannot avoid overshooting): the least coordinate k can
            # still receive is bounded below by (i) whatever of this stratum's
            # remaining budget cannot be absorbed by its remaining
            # coordinate-k-inactive observations, plus (ii) every future
            # stratum's own forced minimum. If even this smallest case already
            # exceeds tn_obs[k], the state is dead.
            bound_min_k <- pmax(0L, budget - suf0_next[k]) + min_future[k]

            keep <- keep &
                (new_state[[cc]] <= tn_obs[k]) &
                ((new_state[[cc]] + bound_max_k) >= tn_obs[k]) &
                ((new_state[[cc]] + bound_min_k) <= tn_obs[k])
        }
    }

    if (!all(keep)) {
        new_state <- .tbl_subset(new_state, keep)
    }
    new_state
}

#' End-of-stratum finalization (network algorithm)
#'
#' Filters the frontier to the states that hit this stratum's required
#' total exactly (eliminating its own stratum intercept, matching the
#' paper's per-stratum network), then resets the within-stratum counter
#' so the same frontier can continue directly into the next stratum's
#' observations.
#'
#' @param state The frontier (a named list of equal-length vectors) after
#'   this stratum's last observation has been processed.
#' @inheritParams .network_expand_observation
#' @param cond_cols Character vector of nuisance-column names.
#' @return The filtered, reset frontier (a named list of equal-length
#'   vectors), possibly with zero-length vectors, the caller checks for
#'   and reports that.
#' @noRd
.network_finish_stratum <- function(state, ms_here, cond_cols) {
    keep <- state$m_j == ms_here
    if (!any(keep)) {
        return(.tbl_subset(state, keep))
    }
    state <- .tbl_subset(state, keep)
    state$m_j <- rep(0L, length(state$m_j))
    log_sum_exp_by(state, by_cols = c("m_j", "t_j", cond_cols), lw_col = "lw")
}

# ==============================================================================
# Per-parameter worker
# ==============================================================================

#' One exact-effect column's full network-algorithm computation
#'
#' Everything `ep_network()`'s per-parameter loop does for one
#' requested column `jcol`: per-stratum precomputation, the single
#' continuous forward pass over every observation in every stratum
#' (Rule A/B pruning via `.network_expand_observation()`, per-stratum
#' finalization via `.network_finish_stratum()`), and assembling the
#' final output row. Split out for the same reason as
#' `.multishift_one_parameter()` in `ep-multivariate-shift.R`:
#' so the main function's own body reads as "for each requested column,
#' compute its result."
#'
#' @inheritParams ep_multi_shift
#' @param slevels_ord Stratum levels in processing order, as returned by
#'   `.order_strata_by_size()`.
#' @inheritParams .stratum_slice
#' @return A list with `parm` (the parameter name), `outdist_entry`, and
#'   `result_row`, ready to be stored under `parm` in the main function's
#'   `outdist`/`results` lists.
#' @noRd
.network_one_parameter <- function(
    model,
    jcol,
    slevels_ord,
    condition_on,
    alpha,
    ci_type,
    verbose
) {
    parm <- model$cn[jcol]
    if (verbose) {
        message(sprintf(
            "\n=== EXACT for parameter column: %s (network algorithm) ===",
            parm
        ))
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
    cond_cols <- make_cond_cols(p_cond)

    # Per-stratum precomputation (cheap; no per-stratum table built)
    strata_info <- vector("list", length(slevels_ord))
    for (s in seq_along(slevels_ord)) {
        strata_info[[s]] <- .network_stratum_slice(
            model,
            slevels_ord[s],
            jcol,
            cond_cols_idx,
            p_cond,
            verbose
        )
    }

    # rem_after_strata[[s]] = sum of max_add_s (LP side) over every stratum
    # AFTER s in the chosen processing order; min_after_strata[[s]] = the
    # same, for min_add_s (SP side).
    max_add_list <- lapply(strata_info, `[[`, "max_add_s")
    min_add_list <- lapply(strata_info, `[[`, "min_add_s")
    rem_after_strata <- .cumulative_from_end(max_add_list, p_cond)
    min_after_strata <- .cumulative_from_end(min_add_list, p_cond)

    state <- list(m_j = 0L, t_j = 0L)
    for (cc in cond_cols) {
        state[[cc]] <- 0L
    }
    state[["lw"]] <- 0.0

    for (s in seq_along(strata_info)) {
        info <- strata_info[[s]]
        ms_here <- info$m_s_obs
        rem_future <- rem_after_strata[[s]]
        min_future <- min_after_strata[[s]]

        for (r in seq_len(info$n_s)) {
            new_state <- .network_expand_observation(
                state,
                r,
                info,
                p_cond,
                cond_cols,
                tn_obs,
                ms_here,
                rem_future,
                min_future
            )
            if (.tbl_nrow(new_state) == 0L) {
                stop(
                    sprintf(
                        "All states pruned mid-stratum %s at observation %d of %d
            for parameter %s.",
                        slevels_ord[s],
                        r,
                        info$n_s,
                        parm
                    ),
                    call. = FALSE
                )
            }
            state <- log_sum_exp_by(
                new_state,
                by_cols = c("m_j", "t_j", cond_cols),
                lw_col = "lw"
            )
        }

        state <- .network_finish_stratum(state, ms_here, cond_cols)
        if (.tbl_nrow(state) == 0L) {
            stop(
                sprintf(
                    "All states pruned after finishing stratum %s for parameter %s.",
                    slevels_ord[s],
                    parm
                ),
                call. = FALSE
            )
        }

        if (verbose) {
            message(sprintf(
                "  after finishing stratum %d/%d (%s): live states = %d",
                s,
                length(strata_info),
                slevels_ord[s],
                .tbl_nrow(state)
            ))
        }
    }

    state$m_j <- NULL
    dist1d <- filter_and_marginalize_pn_lw(state, tn_obs)
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

#' Exact conditional Poisson regression (network algorithm)
#'
#' A structurally different, more aggressively pruned alternative to
#' [ep_multi_shift()], generalizing the network algorithm of Mehta, Patel,
#' and Senchaudhuri (1992), SAS's own `EXACTOPTIONS METHOD=NETWORK`,
#' to Poisson counts. This algorithm avoids constructing separate per-stratum
#' tables and cross-joining them. Instead, it processes all observations across
#' every stratum in a single, continuous forward pass. To optimize performance,
#' it instantly discards any partial state the moment it becomes mathematically
#' impossible to hit the observed nuisance total (either because it cannot
#' reach the target or will inevitably overshoot it). These mid-stratum pruning
#' boundaries are calculated just once per stratum using a closed-form backward
#' induction, bypassing the need to wait for stratum boundaries or explicit
#' cross-joins.
#' Every argument, cap rule, interaction-naming convention, and output
#' format is identical to [ep_multi_shift()], this function reuses the
#' same data-preparation, cap-rule, stratum-ordering,
#' `compute_stratum_max_add()`, CI, p-value, and CMLE helpers; only the
#' table-construction strategy differs. See [ep_multi_shift()]'s
#' documentation for the full description of every argument. Unlike
#' [ep_multi_shift()], this function has no `chunk_size` argument: it
#' never materializes a cross-join between two strata's tables for
#' `chunk_size` to bound.
#'
#' @inheritParams ep_multi_shift
#'
#' @return An object of class `c("ep_network", "ep_multi_shift", "ep_fit")`,
#'   in the same format as [ep_multi_shift()].
#'
#' @examples
#' df <- data.frame(
#'   group = c( 1L,  1L,  1L,  0L,  0L,  0L),
#'   age      = c(50L, 60L, 70L, 50L, 60L, 70L),
#'   subjects = c(3492, 2141, 1711, 3525, 2166, 1724),
#'   cases    = c(   3,    2,    1,   10,   7,   5),
#'   futime   = c(50, 60, 40, 51, 62, 42)
#' )
#'
#' df$lnfutime <- log(df$futime)
#'
#' ep_network(
#'   formula    = cases ~ group,
#'   data       = df,
#'   offset     = "lnfutime",
#'   offset_log = TRUE,
#'   strata     = "age",
#'   exact_vars = "group",
#'   ref_levels = list(group="0", age="50"),
#'   alpha      = 0.05,
#'   verbose    = TRUE,
#'   ci_type = "two-sided"
#' )
#'
#' @references Mehta CR, Patel NR, Senchaudhuri P. Exact Stratified Linear
#'   Rank Tests for Ordered Categorical and Binary Data. Journal of
#'   Computational and Graphical Statistics, 1(1):21-40, 1992.
#' @export
ep_network <- function(
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
        chunk_size = NULL
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
    slevels_ord <- .order_strata_by_size(model$slevels0, model$sid, model$y)

    results <- list()
    outdist <- list()

    for (jcol in model$exact_cols) {
        one <- .network_one_parameter(
            model,
            jcol,
            slevels_ord,
            condition_on,
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
        subclass = c("ep_network", "ep_multi_shift"),
        algorithm_label = paste(
            "Exact conditional Poisson",
            "(categorical + STRATA; network algorithm)"
        )
    )
}
