################################################################################
# Shared internal helpers used by every algorithm in this package,
# ep_closed_form(), ep_multi_shift() and
# ep_network(). Nothing in this file is exported; it exists
# so that no algorithm re-implements the same numerics, data preparation,
# or reporting code.
#
# What is unified here, and why:
#
#   * Low-level numerics (log_sum_exp) and data preparation (stratum-ID
#     construction, offset resolution, reference-level handling) are
#     identical work across algorithms and are shared outright.
#
#   * The exact-inference primitives (tilt_probs, solve_cmle,
#     solve_tail_root, two_sided_from_tails) share the SAME formulas
#     (every algorithm produces a one-dimensional exponential-tilt family
#     Pr_beta(T=t) proportional to C(t) exp(beta t)), and the conditional
#     MLE genuinely uses one shared function everywhere, because every
#     algorithm already needs the *complete* coefficient table for that
#     calculation.
#
#     The exact confidence interval and p-value, by contrast, share only
#     the root-finding scaffolding here (solve_tail_root,
#     two_sided_from_tails), not a single table-based evaluator. That is
#     deliberate: ep_closed_form()'s speed advantage for the
#     test and the confidence interval comes specifically from evaluating
#     tail probabilities from a *truncated* table via an analytic
#     normalizing factor, instead of summing over an explicit, complete
#     table the way the multivariate-shift and network algorithms do.
#     Forcing both through one function that consumes "a table" would
#     silently remove that optimization from the closed-form algorithm.
#     See the accompanying theory document for the full argument. Each
#     algorithm therefore supplies its own small callback to
#     solve_tail_root() that knows how to evaluate a tail probability at
#     a given beta in the way that suits its own representation of the
#     distribution.
################################################################################

# ------------------------------------------------------------------------
# Low-level numerics ----
# ------------------------------------------------------------------------

#' Numerically stable log-sum-exp
#'
#' Computes `log(sum(exp(v)))` by subtracting off the running maximum
#' before exponentiating and adding it back afterward, avoiding the
#' overflow/underflow that evaluating `exp()` on raw log-weights would
#' otherwise risk. Used throughout the package (`tilt_probs()`, every
#' `log_sum_exp_by()` grouping step in the dynamic program and
#' convolution) any time several log-scale contributions to the same
#' state need to be combined into one.
#'
#' @param v Numeric vector of log-values to combine.
#' @return The single value `log(sum(exp(v)))`, computed in a
#'   numerically stable way; `-Inf` if `v` has length zero or every
#'   element of `v` is `-Inf`.
#' @noRd
log_sum_exp <- function(v) {
    if (!length(v)) {
        return(-Inf)
    }
    m <- max(v)
    if (!is.finite(m)) {
        return(m)
    } # all -Inf (or +Inf, in which case the sum is +Inf too)
    m + log(sum(exp(v - m)))
}

#' Number of rows in a "table" represented as a plain named list of
#' equal-length vectors, used throughout the general algorithms
#'
#' The general algorithms' hot loops, the per-stratum dynamic program in
#' `.joint_dist_stratum_dp_vec()` and the cross-stratum convolution in
#' `.convolve_joint_pnuis_vec()`, build and discard many small
#' intermediate "joint tables" (one set of parallel vectors `m_j`, `t_j`,
#' `cond_1`, ..., `lw` per live state) far more often than they build one
#' large table. A plain named list of atomic vectors, subset directly
#' with base `[`, avoids the per-call data.frame-method dispatch overhead
#' that dominates when a grouping/filtering step like this runs once per
#' observation on a small working set. Converting to an actual `data.frame`
#' happens only once, at
#' the very end of the whole pipeline, inside `.assemble_parameter_output()`.
#'
#' @param tbl A named list of equal-length vectors, or `NULL`.
#' @return The common length of every element of `tbl` (`0L` if `tbl` is
#'   `NULL` or empty).
#' @noRd
.tbl_nrow <- function(tbl) {
    if (is.null(tbl) || length(tbl) == 0L) {
        return(0L)
    }
    length(tbl[[1L]])
}

#' Subset every column of a list-of-vectors "table" by the same index
#'
#' @param tbl A named list of equal-length vectors.
#' @param idx A logical or integer index, applied identically to every
#'   element of `tbl` (integer indices may repeat, e.g. to replicate
#'   rows for a cross-join-style expansion).
#' @return A named list with the same names as `tbl`, each element
#'   indexed by `idx`.
#' @noRd
.tbl_subset <- function(tbl, idx) {
    lapply(tbl, function(col) col[idx])
}

#' Row-bind two list-of-vectors "tables" with identical column sets
#'
#' @param tbl1,tbl2 Named lists of equal-length vectors, with identical
#'   `names()` (not necessarily in the same order).
#' @return A named list, in `tbl1`'s column order, with each column the
#'   concatenation of the corresponding columns of `tbl1` and `tbl2`.
#' @noRd
.tbl_rbind <- function(tbl1, tbl2) {
    tbl2 <- tbl2[names(tbl1)]
    Map(c, tbl1, tbl2)
}

#' Group-wise log_sum_exp() for a list-of-vectors "table", used by
#' ep_multi_shift() and ep_network()
#'
#' Merges every row sharing the same combination of `by_cols` into a
#' single row by combining their `lw_col` values via `log_sum_exp()`,
#' the step that collapses distinct (state, candidate-value) expansions
#' of the per-observation dynamic program, or distinct pairs from the
#' cross-stratum convolution, that happen to land on the same key back
#' down to one log-weight per key.
#'
#' Groups are identified by encoding `by_cols` into a single numeric key,
#' one column at a time, using each column's own observed range
#' (`max(v) + 1`) as its radix. Using each column's own range, rather
#' than a fixed-width multiplier, guarantees distinct combinations can
#' never collide onto the same key regardless of how large any one
#' coordinate happens to be. The key is sorted once (`order()`, radix
#' sort for numeric input), and each group's combined log-weight is
#' obtained from the numerically stable maximum (`tapply()`) and a
#' single vectorized, C-level grouped sum (`rowsum()`), rather than one
#' `log_sum_exp()` call per group.
#'
#' @param tbl A named list containing at least the columns named in
#'   `by_cols` and `lw_col`, each an equal-length vector.
#' @param by_cols Character vector of column names identifying a state
#'   (e.g. `c("m_j", "t_j", cond_cols)`).
#' @param lw_col Name of the log-weight column to combine within each
#'   group (default `"lw"`).
#' @return A named list with one row per distinct combination of
#'   `by_cols` and a single `lw` element holding the log-sum-exp of
#'   `lw_col` within that group.
#' @noRd
log_sum_exp_by <- function(tbl, by_cols, lw_col = "lw") {
    n <- .tbl_nrow(tbl)
    if (n == 0L) {
        out <- tbl[by_cols]
        out[[lw_col]] <- numeric(0.0)
        return(out)
    }

    key <- rep(0L, n)
    mult <- 1L
    for (col_name in by_cols) {
        v <- as.double(tbl[[col_name]])
        key <- key + mult * v
        mult <- mult * (max(v) + 1L)
    }

    ord <- order(key)
    key_ord <- key[ord]
    new_grp <- c(TRUE, key_ord[-1L] != key_ord[-n])
    grp_id <- cumsum(new_grp)

    lw_ord <- tbl[[lw_col]][ord]
    grp_max <- as.numeric(tapply(lw_ord, grp_id, max))
    shifted <- exp(lw_ord - grp_max[grp_id])
    grp_sum <- as.numeric(rowsum(shifted, grp_id, reorder = FALSE))
    combined_lw <- grp_max + log(grp_sum)

    first_idx <- ord[!duplicated(grp_id)]
    out <- .tbl_subset(tbl[by_cols], first_idx)
    out[[lw_col]] <- combined_lw
    out
}

#' Column names for `p_cond` nuisance-statistic axes, used by
#' ep_multi_shift()'s dynamic program and convolution.
#'
#' Generates the joint-table element names (`"cond_1"`, `"cond_2"`, ...)
#' used throughout the general algorithms to carry each nuisance
#' coordinate's running value, so that every function that builds or
#' consumes a joint table agrees on the same naming convention.
#'
#' @param p_cond Number of nuisance coordinates; coerced to a
#'   non-negative integer, treated as `0L` if `NULL`, non-scalar, `NA`,
#'   or negative.
#' @return A character vector `c("cond_1", ..., "cond_<p_cond>")`, or
#'   `character(0)` if `p_cond == 0`.
#' @noRd
make_cond_cols <- function(p_cond) {
    p_cond <- as.integer(p_cond %||% 0L)
    if (length(p_cond) != 1L || is.na(p_cond) || p_cond < 0L) {
        p_cond <- 0L
    }
    if (p_cond == 0L) character(0L) else paste0("cond_", seq_len(p_cond))
}

# ------------------------------------------------------------------------
# Shared data preparation ----
# ------------------------------------------------------------------------

#' Build a single stratum-ID vector from zero or more stratification columns
#'
#' Both algorithms treat "no strata supplied" as a single pseudo-stratum
#' called `"ALL"`, so that every statement about "strata" in either
#' algorithm's documentation also covers the unstratified case.
#'
#' @inheritParams ep_multi_shift
#' @param strata Character vector of column names in `data`, or `NULL`.
#' @param row_idx Integer row indices (post `model.frame()`/`na.omit()`
#'   alignment) to subset `data` by.
#' @return A character vector, one entry per element of `row_idx`: the
#'   values of every `strata` column pasted together with `"__"` as a
#'   single stratum ID per observation, or `"ALL"` for every observation
#'   if `strata` is `NULL`/empty.
#' @noRd
build_stratum_id <- function(data, strata, row_idx) {
    if (is.null(strata) || !length(strata)) {
        return(rep("ALL", length(row_idx)))
    }
    strat_cols <- lapply(strata, function(s) as.character(data[[s]][row_idx]))
    do.call(paste, c(strat_cols, sep = "__"))
}

#' Apply a user-supplied reference level to one factor column
#'
#' Used both by ep_multi_shift()'s multi-variable
#' factor preparation and directly by ep_closed_form()'s
#' single-variable group factor.
#'
#' @param f A vector or factor.
#' @param var The name of `f`, used only for error messages.
#' @param ref_levels A named list; if it has an entry named `var`, that
#'   level is moved to the front via `stats::relevel()`.
#' @return `f`, coerced to a factor if it was not already one, with its
#'   reference level moved to the front if `ref_levels` has an entry
#'   named `var`; unchanged otherwise.
#' @noRd
apply_ref_level <- function(f, var, ref_levels) {
    if (!is.factor(f)) {
        f <- factor(f)
    }
    if (var %in% names(ref_levels)) {
        ref <- as.character(ref_levels[[var]])
        if (!ref %in% levels(f)) {
            stop(
                sprintf("ref_levels '%s' is not a level of '%s'.", ref, var),
                call. = FALSE
            )
        }
        f <- stats::relevel(f, ref = ref)
    }
    f
}

#' Is an `offset` argument actually specified?
#'
#' Centralizes the package's one notion of "no offset supplied", so that
#' every function deciding whether to treat an observation as having an
#' offset (as opposed to zero log-offset / unit exposure) agrees with
#' every other.
#'
#' @param offset The `offset` argument as supplied to a main function:
#'   `NULL`, `""`, or the name of a column in `data`.
#' @return `TRUE` if `offset` names an actual column (i.e. is not `NULL`
#'   and not a zero-length or empty-string character value), `FALSE`
#'   otherwise.
#' @noRd
has_offset_spec <- function(offset) {
    !(is.null(offset) ||
        (is.character(offset) && length(offset) == 1L && !nzchar(offset)))
}

#' Resolve an offset variable to both the log-offset and raw-exposure scales
#'
#' Both exported functions accept the same `offset`/`offset_log` pair of
#' arguments. `ep_multi_shift()` needs the *log*-offset
#' (it enters the dynamic program's log-weight additively, as `y * o_i`).
#' `ep_closed_form()` needs the raw exposure (it enters the
#' binomial success probability as a ratio, `n_s / m_s`). Resolving both
#' here in one place keeps the meaning of `offset_log` identical between
#' the two functions.
#'
#' @inheritParams ep_multi_shift
#' @param offset Name of the offset column in `data`, or `NULL`/`""`.
#' @param offset_log If `TRUE`, `offset` already holds log(exposure); if
#'   `FALSE`, `offset` holds raw exposure.
#' @param row_idx Integer row indices to subset `data` by.
#' @param required If `TRUE`, an absent offset is an error rather than
#'   being treated as "no offset" (o_i = 0, exposure = 1).
#' @return A list with `log_offset` and `exposure` (`exposure =
#'   exp(log_offset)`), each aligned to `row_idx`.*
#' @noRd
resolve_offset <- function(
    data,
    offset,
    offset_log,
    row_idx,
    required = FALSE
) {
    if (!has_offset_spec(offset)) {
        if (required) {
            stop(
                "`offset` (name of the person-time / exposure variable) is required.", # nolint: nonportable_path_linter
                call. = FALSE
            )
        }
        return(
            list(
                log_offset = rep(0.0, length(row_idx)),
                exposure = rep(1.0, length(row_idx))
            )
        )
    }
    if (!is.character(offset) || length(offset) != 1L) {
        stop("offset must be a single variable name or NULL.", call. = FALSE)
    }
    if (!offset %in% names(data)) {
        stop(sprintf("offset '%s' not found in data.", offset), call. = FALSE)
    }
    off_raw <- as.numeric(data[[offset]][row_idx])
    if (anyNA(off_raw)) {
        stop("offset contains NA after alignment.", call. = FALSE)
    }

    # Validate off_raw *before* calling log()/exp() below, rather than
    # after, so a bad value produces this function's own clean error
    # instead of a base-R "NaNs produced" warning (log() of a non-positive
    # value) on the way to it. What "valid" means depends on offset_log:
    # when offset_log = FALSE, off_raw is the raw exposure itself, so it
    # must be strictly positive (and finite) for log() to be well-defined;
    # when offset_log = TRUE, off_raw is already a log-exposure, which is
    # legitimately allowed to be negative or zero (log(exposure) < 0
    # whenever exposure < 1), only finiteness is required in that case,
    # since exp() of any finite number is already strictly positive.
    off_raw_valid <- if (offset_log) {
        all(is.finite(off_raw))
    } else {
        all(is.finite(off_raw)) && all(off_raw > 0.0)
    }
    if (!off_raw_valid) {
        stop(
            "offset/exposure must be strictly positive and finite after alignment.", # nolint: nonportable_path_linter
            call. = FALSE
        )
    }

    log_offset <- if (offset_log) off_raw else log(off_raw)
    exposure <- if (offset_log) exp(off_raw) else off_raw

    if (!all(is.finite(log_offset)) || !all(is.finite(exposure))) {
        stop(
            "offset/exposure must be strictly positive and finite after alignment.", # nolint: nonportable_path_linter
            call. = FALSE
        )
    }
    list(
        log_offset = log_offset,
        exposure = exposure
    )
}

# ------------------------------------------------------------------------
# Shared exact-inference primitives (exponential-tilt family) ----
# ------------------------------------------------------------------------
# Both algorithms ultimately produce a distribution of a scalar sufficient
# statistic T of the form  Pr_beta(T=t) proportional to C(t) * exp(beta*t).

#' Tilted probabilities from a complete (t, log-weight) table
#'
#' @param beta The natural parameter (log rate ratio).
#' @param t_vals Attainable values of the sufficient statistic.
#' @param log_weight The corresponding natural-log coefficients, log C(t).
#' @return A numeric vector of probabilities, same length/order as `t_vals`.
#' @noRd
tilt_probs <- function(beta, t_vals, log_weight) {
    lw <- log_weight + beta * t_vals
    lw <- lw - max(lw)
    p <- exp(lw)
    p / sum(p)
}

#' Conditional maximum likelihood estimate from a *complete* coefficient table
#'
#' Shared by all algorithms: `ep_multi_shift()` and
#' `ep_network()` passes its `dist1d` table directly
#' (`t_vals = dist1d$t_j`, `log_weight = dist1d$lw`);
#' `ep_closed_form()` passes its complete `c_k` table.
#' All algorithms already need the complete table for this step -
#' the truncation trick used elsewhere in `ep_closed_form()`
#' only benefits the confidence interval and the p-value (see
#' `solve_tail_root()` below and the theory vignette).
#'
#' @param t_vals Attainable values of the sufficient statistic (the full
#'   support, not truncated).
#' @param log_weight The corresponding natural-log coefficients.
#' @param t_obs The observed sufficient statistic.
#' @return A list with `estimate` (possibly `+-Inf` at the boundary of
#'   `t_vals`) and `se` (`NA` at the boundary).
#' @noRd
solve_cmle <- function(t_vals, log_weight, t_obs) {
    tmin <- min(t_vals)
    tmax <- max(t_vals)
    if (t_obs <= tmin) {
        return(list(estimate = -Inf, se = NA_real_))
    }
    if (t_obs >= tmax) {
        return(list(estimate = Inf, se = NA_real_))
    }

    mean_t <- function(b) sum(t_vals * tilt_probs(b, t_vals, log_weight))
    root <- solve_tail_root(mean_t, t_obs)

    p <- tilt_probs(root, t_vals, log_weight)
    mu <- sum(t_vals * p)
    v <- sum(t_vals^2L * p) - mu^2L
    se <- if (v > 0.0) 1.0 / sqrt(v) else NA_real_
    list(estimate = root, se = se)
}

#' Generic bracket-and-root-find for one exact confidence limit (or the CMLE)
#'
#' Finds `beta` such that `stat_fun(beta) == target`, where `stat_fun` is a
#' caller-supplied function of `beta` alone. This is the shared control
#' flow (progressively widening bracket, then `uniroot()`) behind the
#' confidence interval and CMLE calculations of both algorithms; the two
#' algorithms differ only in how `stat_fun` is computed (see the file
#' header), which is why it is passed in rather than hard-coded here.
#'
#' The caller is expected to have already handled the boundary case where
#' no finite root exists (i.e. the observed statistic sits at the edge of
#' its attainable range) *before* calling this function, since detecting
#' that requires knowing the support of the statistic, which this function
#' does not see.
#'
#' @param stat_fun A function of `beta` returning the quantity to match to
#'   `target` (a tail probability, or the conditional mean of the
#'   sufficient statistic).
#' @param target The target value `stat_fun(beta)` should equal at the
#'   root.
#' @param lo,hi Starting bracket, in `beta` (log rate ratio) units.
#' @param expand_factor Multiplicative factor used to widen the bracket.
#' @param max_expand Maximum number of widening attempts.
#' @return The root, `beta_hat`.
#' @noRd
solve_tail_root <- function(
    stat_fun,
    target,
    lo = -50.0,
    hi = 50.0,
    expand_factor = 1.4,
    max_expand = 1000L
) {
    f <- function(b) stat_fun(b) - target
    for (i in seq_len(max_expand)) {
        s_lo <- tryCatch(f(lo), error = function(e) NA_real_)
        s_hi <- tryCatch(f(hi), error = function(e) NA_real_)
        if (is.finite(s_lo) && is.finite(s_hi) && sign(s_lo) != sign(s_hi)) {
            break
        }
        lo <- lo * expand_factor
        hi <- hi * expand_factor
    }
    stats::uniroot(f, c(lo, hi), tol = 1e-12, maxiter = 5000L)$root
}

#' Two-sided exact p-value from two one-sided tail probabilities
#'
#' Both algorithms define the two-sided exact p-value the same way (and so
#' does SAS's `PROC GENMOD` `EXACT` statement): twice the smaller one-sided
#' tail probability at the null, capped at 1.
#'
#' @param p_lower `Pr_(beta=0)(T <= t_obs)`.
#' @param p_upper `Pr_(beta=0)(T >= t_obs)`.
#' @return A list with `pvalue` (two-sided), `pvalue_onesided`, and
#'   `direction`: `"upper"` when the lower tail is the smaller (or equal)
#'   one, i.e. the observed statistic sits low enough that only an upper
#'   confidence limit is informative (the lower limit is the boundary),
#'   and `"lower"` when the upper tail is smaller, i.e. only a lower
#'   confidence limit is informative. This matches the convention of the
#'   original `exactpoisson_test_log()`/`direction` logic: `direction =
#'   "upper"` iff `p_lower <= p_upper`.
#' @noRd
two_sided_from_tails <- function(p_lower, p_upper) {
    onesided <- min(p_lower, p_upper)
    list(
        pvalue = min(1.0, 2L * onesided),
        pvalue_onesided = onesided,
        direction = if (p_lower <= p_upper) "upper" else "lower"
    )
}

#' Shared exact-confidence-interval branching logic
#'
#' All algorithms' confidence intervals reduce to the same four cases
#' (two-sided; one-sided reporting only the upper limit; one-sided
#' reporting only the lower limit; and, within each, the boundary case
#' where `t_obs` sits at the edge of its support and the corresponding
#' limit is infinite by definition rather than by root-finding). This
#' function is that shared branching, expressed in terms of two
#' caller-supplied tail-probability callbacks rather than a table, so it
#' can sit in front of either algorithm's own way of evaluating a tail
#' probability:
#' `ep_multi_shift()` and `ep_network()`
#' evaluates `tail_le`/`tail_ge` from its *complete* one-dimensional
#' table via `tilt_probs()`. `ep_closed_form()` evaluates
#' them from a *truncated* table plus an analytic normalizing prefactor (see
#' `cumulative_dist_natlog()` in `ep-closed-form.R` and the
#' theory vignette).
#'
#' @param tail_le A function of `beta` returning `Pr_beta(T <= t_obs)`.
#' @param tail_ge A function of `beta` returning `Pr_beta(T >= t_obs)`.
#' @param alpha Significance level (halved internally for `"two-sided"`).
#' @param ci_type `"two-sided"` or `"one-sided"`.
#' @param direction For `"one-sided"`: `"upper"` reports only the upper
#'   limit (lower limit fixed at `-Inf`); `"lower"` reports only the lower
#'   limit (upper limit fixed at `Inf`). Ignored for `"two-sided"`.
#' @param at_lower_boundary `TRUE` if `t_obs` sits at the minimum of its
#'   support (the lower limit is `-Inf` by definition, not by root-finding).
#' @param at_upper_boundary `TRUE` if `t_obs` sits at the maximum of its
#'   support (the upper limit is `Inf` by definition).
#' @return `c(lo, hi)`, on the beta scale.
#' @noRd
exact_ci_from_tails <- function(
    tail_le,
    tail_ge,
    alpha,
    ci_type = c("two-sided", "one-sided"),
    direction = c("upper", "lower"),
    at_lower_boundary = FALSE,
    at_upper_boundary = FALSE
) {
    ci_type <- rlang::arg_match(ci_type)
    direction <- rlang::arg_match(direction)

    if (ci_type == "two-sided") {
        lo <- if (at_lower_boundary) {
            -Inf
        } else {
            solve_tail_root(tail_ge, alpha / 2L)
        }
        hi <- if (at_upper_boundary) {
            Inf
        } else {
            solve_tail_root(tail_le, alpha / 2L)
        }
    } else if (direction == "upper") {
        lo <- -Inf
        hi <- if (at_upper_boundary) Inf else solve_tail_root(tail_le, alpha)
    } else {
        lo <- if (at_lower_boundary) -Inf else solve_tail_root(tail_ge, alpha)
        hi <- Inf
    }
    c(lo, hi)
}
