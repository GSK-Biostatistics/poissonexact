################################################################################
# Shared infrastructure for the general (arbitrary categorical covariates
# + optional strata) family of exact-inference algorithms: the model/data
# preparation, SAS EXACTMAX/OFFSET cap-rule resolution, per-stratum
# slicing, stratum-ordering heuristics, and output assembly used by all
# of:
#
#   * ep_network()             (ep-network.R)
#   * ep_multi_shift()       (ep-multivariate-shift.R).
#
# These two functions take an identical argument list, prepare the
# model in an identical way, loop over the same requested parameters, and
# turn each parameter's final one-dimensional table into the same
# `exactparmest`/`OUTDIST` output. They differ only in how each stratum's
# own table is built and how strata are combined; that algorithm-
# specific logic lives in each function's own file. Everything here is
# what all three share.
#
# This is deliberately a different file from ep-utils.R.
# ep-utils.R holds low-level numerics and exact-inference
# primitives (log-sum-exp, the exponential-tilt family, root-finding)
# shared by *every* algorithm in the package, including
# ep_closed_form(), which does not use anything in this file,
# since it has its own, much simpler two-group/strata-only data
# preparation (Section 6 of the theory vignette). This file, by contrast,
# is specific to the general categorical-covariate model shared only by
# the multivariate-shift and network algorithms.
#
# Nothing in this file is exported.
################################################################################

# ------------------------------------------------------------------------
# SAS EXACTMAX / OFFSET option parsing ----
# (specific to this family of algorithms: the closed-form algorithm
# always requires no per-observation caps, since it has no individual
# counts to bound, only stratum arm totals)
# ------------------------------------------------------------------------

#' Parse the user's `exactmax` argument into a resolved spec
#'
#' Determines which of the three SAS `EXACTMAX` forms was supplied,
#' not specified, specified without a variable, or `EXACTMAX=varname`,
#' and returns a small structure the rest of `.setup_exact_model()` uses
#' to pick between the three per-observation cap rules of Section 5.2 of
#' the theory vignette. This is parsing only; it does not itself resolve
#' any per-observation cap values (see `si_from_exactmax_var()`,
#' `si_from_stratum_total()`, and `si_from_offset()` for that).
#'
#' @param exactmax `NULL` (not specified), `TRUE` (specified without a
#'   variable), or a single non-empty character variable name (SAS
#'   `EXACTMAX=`), exactly as passed to `ep_multi_shift()`/
#'   `ep_network()`.
#' @return A list with `mode` (one of `"none"`, `"novar"`, `"var"`) and
#'   `var` (the variable name if `mode == "var"`, otherwise
#'   `character(0)`).
#' @noRd
parse_exactmax_spec <- function(exactmax) {
    # Supported forms:
    # - exactmax = NULL         => EXACTMAX not specified
    # - exactmax = TRUE         => EXACTMAX specified without variable
    # - exactmax = "varname"    => EXACTMAX=varname
    if (is.null(exactmax)) {
        return(list(mode = "none", var = character(0L)))
    }
    if (isTRUE(exactmax)) {
        return(list(mode = "novar", var = character(0L)))
    }
    if (is.character(exactmax) && length(exactmax) == 1L && nzchar(exactmax)) {
        return(list(mode = "var", var = exactmax))
    }
    stop(
        "exactmax must be NULL, TRUE, or a single non-empty character
    variable name.",
        call. = FALSE
    )
}

# SAS rule (2): every S_i = the intercept's sufficient statistic (the
# stratum total), computed within each stratum.

#' Per-observation caps under SAS rule 2 (stratum total)
#'
#' Implements SAS cap rule 2 (Section 5.2 of the theory vignette): when
#' `EXACTMAX` is specified without a variable, or neither `OFFSET` nor
#' `EXACTMAX` is specified, every observation's cap is set to its own
#' stratum's observed total. This is not an extra modelling assumption,
#' it simply restates the bound that conditioning on the stratum total
#' already implies for every non-negative count in that stratum.
#'
#' @inheritParams compute_stratum_max_add
#' @return An integer vector, the same length as `y_s`, with every entry
#'   equal to `sum(y_s)`.
#' @noRd
si_from_stratum_total <- function(y_s) {
    ms <- as.integer(sum(as.integer(y_s)))
    rep.int(ms, length(y_s))
}

# SAS rule (1): EXACTMAX=var, so S_i := e_i (per row).

#' Per-observation caps under SAS rule 1 (`EXACTMAX=var`)
#'
#' Implements SAS cap rule 1 (Section 5.2 of the theory vignette): when
#' `EXACTMAX=var` is specified, each observation's cap `S_i` is simply
#' the observed value of that variable for that observation, validated
#' to be a non-negative, finite number.
#'
#' @param exactmax_values Numeric vector of the supplied `EXACTMAX`
#'   variable's values, one per observation.
#' @return An integer vector of per-observation caps `S_i`
#'   (`floor()`-truncated).
#' @noRd
si_from_exactmax_var <- function(exactmax_values) {
    ex <- as.numeric(exactmax_values)
    if (anyNA(ex) || !all(is.finite(ex))) {
        stop("EXACTMAX variable contains NA/NaN/Inf.", call. = FALSE) # nolint: nonportable_path_linter
    }
    exi <- as.integer(floor(ex))
    if (any(exi < 0L)) {
        stop(
            "EXACTMAX variable has negative values (not allowed).",
            call. = FALSE
        )
    }
    exi
}

# SAS rule (3): offset specified, EXACTMAX not specified =>
#  S_i = max(ceil(exp(o_i)), 1).

#' Per-observation caps under SAS rule 3 (`OFFSET` only)
#'
#' Implements SAS cap rule 3 (Section 5.2 of the theory vignette): when
#' `OFFSET` is specified and `EXACTMAX` is not, each observation's cap is
#' derived from its own exposure, `S_i = max(ceil(exp(o_i)), 1)`.
#'
#' @param log_offset Numeric vector of log-offsets (log exposure), one
#'   per observation.
#' @return An integer vector of per-observation caps `S_i`.
#' @noRd
si_from_offset <- function(log_offset) {
    pmax(as.integer(ceiling(exp(log_offset))), 1L)
}

# ------------------------------------------------------------------------
# Design / factor helpers ----
# (specific to this family of algorithms: each
# needs a full model matrix with contrasts for potentially many
# categorical effects, not just one two-level factor)
# ------------------------------------------------------------------------

#' Apply SAS-style `PARAM=REF` treatment coding to factor variables
#'
#' For every variable named in `class_vars`, applies the user-requested
#' reference level (via `apply_ref_level()`) and then sets treatment
#' contrasts with that (first) level as the reference, matching SAS's
#' `PARAM=REF` convention for the `CLASS` statement. Variables not
#' present in `data`, or not coercible to a factor, are left untouched.
#'
#' @param class_vars Character vector of variable names to treat as
#'   SAS-style reference-coded factors (every categorical effect named
#'   in `strata`, `exact_vars`, or `ref_levels`).
#' @inheritParams ep_multi_shift
#' @return `data`, with every variable in `class_vars` converted (or
#'   confirmed) to be a factor with reference-level treatment contrasts.
#' @noRd
prepare_sas_factors <- function(data, class_vars, ref_levels = list()) {
    d <- data
    for (v in class_vars) {
        if (v %in% names(d)) d[[v]] <- apply_ref_level(d[[v]], v, ref_levels)
    }
    # SAS-like PARAM=REF: treatment coding; first level is reference
    for (v in class_vars) {
        if (!v %in% names(d) || !is.factor(d[[v]])) {
            next
        }
        lev <- levels(d[[v]])
        base <- 1L
        contr_mat <- stats::contr.treatment(n = length(lev), base = base)
        colnames(contr_mat) <- lev[-base]
        stats::contrasts(d[[v]]) <- contr_mat
    }
    d
}

#' Build a no-intercept model matrix, preserving the `assign` attribute
#'
#' Builds the model frame and a model matrix with the intercept column
#' removed, by hand carrying the `assign` attribute (mapping each
#' surviving column back to its model term) across the intercept-removal
#' subsetting step, since ordinary `[` subsetting on a matrix does not
#' preserve it. This `assign` attribute is exactly what
#' `exact_term_columns()`/`match_term_label()` rely on to resolve
#' `exact_vars` entries to columns, including interaction terms.
#'
#' @inheritParams ep_multi_shift
#' @param data A data frame, already SAS-style factor-coded (e.g. by
#'   `prepare_sas_factors()`).
#' @return A list with `mf` (the model frame), `x` (the no-intercept
#'   model matrix, with its `assign` attribute preserved), and `terms`
#'   (`stats::terms(mf)`).
#' @noRd
build_mm_no_intercept <- function(formula, data) {
    mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
    x <- stats::model.matrix(
        stats::delete.response(stats::terms(mf)),
        data = mf
    )
    # `model.matrix()` attaches an `assign` attribute (one entry per column,
    # mapping it back to its term, 0 for the intercept, an index into
    # `term.labels` for everything else); this is exactly what
    # `exact_term_columns()`/`match_term_label()` rely on to find the columns
    # for a requested term, including interaction terms. Ordinary `[`
    # subsetting on a matrix does not preserve this (or any other
    # non-standard) attribute.
    asg <- attr(x, "assign")
    int <- which(colnames(x) == "(Intercept)")
    if (length(int)) {
        x <- x[, -int, drop = FALSE]
        asg <- asg[-int]
    }
    attr(x, "assign") <- asg
    list(mf = mf, x = x, terms = stats::terms(mf))
}

#' Generate every permutation of a short vector
#'
#' Used only by `match_term_label()` to make interaction-term matching
#' order-insensitive (R's own `term.labels` order follows the order
#' variables first appear in the formula, e.g. `gender + region +
#' region:gender` labels the interaction `"region:gender"`, never
#' `"gender:region"`, so a user-supplied `exact_vars` entry has to be
#' checked against every ordering to be matched reliably). Guarded to
#' vectors of length <= 6 (720 permutations) since this is only ever
#' applied to the handful of variable names in one interaction term.
#'
#' @param x A short character vector (a term's constituent variable names).
#' @return A list of every permutation of `x`, each itself a character
#'   vector of the same length; `list(x)` unchanged if `length(x) <= 1`
#'   or `length(x) > 6` (the latter to avoid combinatorial blow-up, since
#'   the caller falls back to its own handling in that case).
#' @noRd
.permutations <- function(x) {
    if (length(x) <= 1L) {
        return(list(x))
    }
    if (length(x) > 6L) {
        return(list(x))
    } # refuse to blow up; caller falls back
    out <- list()
    for (i in seq_along(x)) {
        for (p in .permutations(x[-i])) {
            out[[length(out) + 1L]] <- c(x[i], p)
        }
    }
    out
}

#' Match one `exact_vars` entry to its term label
#'
#' The single, authoritative mechanism for turning a user-supplied
#' `exact_vars` entry into the model-matrix columns for that model
#' *term*, a main effect (`"region"`), an interaction of any order
#' (`"region:gender"` or `"region*gender"`), or a `factor()`-wrapped
#' variable (`"factor(region)"`), using R's own `terms()`/`assign`
#' machinery, which already parses interactions correctly. There is no
#' separate, looser name-matching fallback based on column-name prefixes:
#' matching only through `terms()`/`assign` ensures a term is never
#' confused with a different term that happens to share a text prefix,
#' e.g. requesting `exact_vars = "region"` in a model that also contains
#' a `region:gender` interaction resolves only to region's own
#' main-effect columns.
#'
#' @param v One `exact_vars` entry, as supplied by the user.
#' @param term_labels `attr(terms(formula), "term.labels")`.
#' @return The matching index into `term_labels`, or `NA_integer_` if
#'   no candidate spelling of `v` matches any term.
#' @noRd
match_term_label <- function(v, term_labels) {
    candidates <- v
    if (grepl("[*:]", v)) {
        parts <- trimws(strsplit(v, "[*:]")[[1L]])
        parts <- parts[nzchar(parts)]
        for (p in .permutations(parts)) {
            candidates <- c(candidates, paste(p, collapse = ":"))
        }
    }
    candidates <- c(candidates, paste0("factor(", v, ")"))
    for (cand in candidates) {
        tidx <- match(cand, term_labels)
        if (!is.na(tidx)) return(tidx)
    }
    NA_integer_
}

#' Resolve every `exact_vars` entry to its model-matrix columns
#'
#' Applies `match_term_label()` to each entry of `exact_vars` in turn,
#' translating a requested term (main effect or interaction) into the
#' set of model-matrix columns that belong to it, via the model's own
#' `term.labels`/`assign` structure. Errors, with a message listing every
#' term actually available in the model, if any entry does not match a
#' term in `formula`.
#'
#' @param mm The model matrix (as returned by `build_mm_no_intercept()`),
#'   with its `assign` attribute intact.
#' @param terms_obj The model's `terms()` object.
#' @param exact_vars Character vector of requested term names (main
#'   effects or interactions), exactly as the main functions' own
#'   `exact_vars` argument.
#' @return A sorted, deduplicated integer vector of model-matrix column
#'   indices covering every requested term.
#' @noRd
exact_term_columns <- function(mm, terms_obj, exact_vars) {
    cn <- colnames(mm)
    if (is.null(cn)) {
        stop("Model matrix has no column names.", call. = FALSE)
    }
    term_labels <- attr(terms_obj, "term.labels")
    a <- attr(mm, "assign")
    if (is.null(term_labels) || is.null(a)) {
        stop(
            "Model matrix is missing its `term.labels`/`assign` attributes; ", # nolint: nonportable_path_linter
            "cannot map exact_vars to columns.",
            call. = FALSE
        )
    }
    cols <- integer(0L)
    for (v in exact_vars) {
        tidx <- match_term_label(v, term_labels)
        if (is.na(tidx)) {
            stop(
                sprintf(
                    paste0(
                        "exact_vars entry '%s' does not match any term in ",
                        "`formula`.\n",
                        "Terms available in this model are: %s.\n",
                        "To request an interaction, name every variable in it, ",
                        "joined by '*' or ':' (e.g. \"region:gender\" or ",
                        "\"region*gender\"); it must be a term that actually ",
                        "appears in `formula` (requesting a main effect does ",
                        "not automatically include its interactions, or vice ",
                        "versa, each term to be tested must be listed ",
                        "explicitly)."
                    ),
                    v,
                    toString(term_labels)
                ),
                call. = FALSE
            )
        }
        cols <- c(cols, which(a == tidx))
    }
    sort(unique(cols))
}

# ------------------------------------------------------------------------
# Pruning helper ----
# ------------------------------------------------------------------------

#' One stratum's own forced-maximum contribution to each nuisance coordinate
#'
#' Computes, for one stratum and every nuisance column, the largest total
#' that column could possibly reach from this stratum alone: the sum of
#' the per-observation caps `si_s` over the rows where that column is
#' non-zero, capped by the stratum's own observed total. This is the
#' per-stratum bound `U_{s,k}` of Section 5.5 of the theory vignette,
#' used by Rule B pruning during cross-stratum combination.
#'
#' @param y_s Integer vector of observed counts for one stratum.
#' @param xc_s Numeric matrix of nuisance-column values for that
#'   stratum's observations (`nrow(xc_s) == length(y_s)`).
#' @param si_s Integer vector of per-observation caps `S_i` for that
#'   stratum.
#' @return An integer vector of length `ncol(xc_s)`, one bound per
#'   nuisance column.
#' @noRd
compute_stratum_max_add <- function(y_s, xc_s, si_s) {
    ms <- as.integer(sum(y_s))
    p_cond <- ncol(xc_s)
    out <- integer(p_cond)
    for (k in seq_len(p_cond)) {
        idx1 <- which(xc_s[, k] != 0.0)
        max_k <- if (length(idx1) == 0L) 0L else sum(si_s[idx1])
        out[k] <- as.integer(min(ms, max_k))
    }
    out
}

# ------------------------------------------------------------------------
# Global filtering and marginalization ----
#
# Shared as-is by every algorithm that builds a joint (t_j, cond_*, lw) table,
# ep_multi_shift() and ep_network() alike,
# since this step only depends on the *shape* of that table,
# never on how it was built.
# ------------------------------------------------------------------------

#' Filter a joint table to the observed nuisance value and marginalize
#'
#' Performs the final two steps common to every general algorithm once
#' its complete joint table `C(t, c)` has been built (Section 5.6 of the
#' theory vignette): filter to the rows whose nuisance columns equal the
#' observed nuisance totals `tn_obs`, then marginalize by dropping the
#' now-constant nuisance columns and merging any resulting duplicate `t_j`
#' values by log-sum-exp. When there are no nuisance columns at all
#' (`length(tn_obs) == 0`), the filtering step is skipped and only the
#' marginalization is performed.
#'
#' @param joint A named list (see `.tbl_nrow()`/`.tbl_subset()` in
#'   `ep-utils.R`) with element `t_j`, one `cond_1`, ...,
#'   `cond_p` element per nuisance coordinate, and `lw` (log-weight), each
#'   an equal-length vector.
#' @param tn_obs Numeric vector of the observed nuisance totals, of
#'   length `p_cond` (the number of `cond_*` elements in `joint`).
#' @return A named list with elements `t_j` and `lw`, ordered by `t_j`:
#'   the one-dimensional table `{(t, log C(t))}` used by every downstream
#'   exact-inference computation.
#' @noRd
filter_and_marginalize_pn_lw <- function(joint, tn_obs) {
    p_cond <- length(tn_obs)

    if (p_cond == 0L) {
        out <- log_sum_exp_by(joint, by_cols = "t_j", lw_col = "lw")
        return(.tbl_subset(out, order(out$t_j)))
    }

    cond_cols <- paste0("cond_", seq_len(p_cond))
    keep <- rep(TRUE, .tbl_nrow(joint))
    for (k in seq_len(p_cond)) {
        keep <- keep & (joint[[cond_cols[k]]] == tn_obs[k])
    }

    jf <- .tbl_subset(joint, keep)
    if (.tbl_nrow(jf) == 0L) {
        stop("No states remain after global nuisance filtering.", call. = FALSE)
    }

    out <- log_sum_exp_by(jf, by_cols = "t_j", lw_col = "lw")
    ord <- order(out$t_j)
    list(t_j = out$t_j[ord], lw = out$lw[ord])
}

#' Exact confidence interval for the general algorithms' 1-D table
#'
#' Builds the two tail-probability callbacks from the complete `dist1d`
#' table via `tilt_probs()`, then defers all `ci_type`/`direction`/boundary
#' branching to the shared `exact_ci_from_tails()` (see that function's
#' documentation in `ep-utils.R` for why the branching is
#' shared but this table-based evaluation is not: `ep_closed_form()`
#' evaluates the same two callbacks from a truncated table instead, which
#' is precisely its speed advantage for the CI).
#'
#' @param dist1d The 1-D coefficient table (`t_j`, `lw`) for this parameter.
#' @param t_obs Observed value of the sufficient statistic.
#' @inheritParams exact_ci_from_tails
#' @return `c(lo, hi)`.
#' @noRd
exact_ci_from_dist1d <- function(
    dist1d,
    t_obs,
    alpha,
    ci_type = c("two-sided", "one-sided"),
    direction = c("upper", "lower")
) {
    ci_type <- rlang::arg_match(ci_type)
    direction <- rlang::arg_match(direction)

    tail_ge <- function(b) {
        sum(tilt_probs(b, dist1d$t_j, dist1d$lw)[dist1d$t_j >= t_obs])
    }
    tail_le <- function(b) {
        sum(tilt_probs(b, dist1d$t_j, dist1d$lw)[dist1d$t_j <= t_obs])
    }

    exact_ci_from_tails(
        tail_le,
        tail_ge,
        alpha,
        ci_type,
        direction,
        at_lower_boundary = (t_obs <= min(dist1d$t_j)),
        at_upper_boundary = (t_obs >= max(dist1d$t_j))
    )
}

# ==============================================================================
# Shared orchestration helpers ----
#
# ep_multi_shift() and ep_network()
# take an identical argument list, prepare
# the model in an identical way, loop over the same requested parameters,
# and turn each parameter's final one-dimensional table into the same
# exactparmest/OUTDIST output. They differ only in how each stratum's
# table is built and how strata are combined. The helpers below factor
# out everything the three functions share, so that each main function's
# own body only ever expresses what is actually specific to it.
# ==============================================================================

#' Prepare a model for exact conditional inference
#'
#' Performs every step common to all three main functions before their
#' per-parameter loop begins: resolving `strata`/`exactmax`, preparing
#' SAS-style reference-coded factors, building the no-intercept model
#' matrix, resolving the offset, choosing which of the three SAS cap
#' rules applies (Section 5.2 of the theory vignette), building the
#' stratum-ID vector, and resolving `exact_vars` to model-matrix columns.
#'
#' @inheritParams ep_multi_shift
#' @return A list with `x` (model matrix), `y` (response), `log_off`
#'   (log-offset per observation), `exactmax_vec` (per-observation
#'   EXACTMAX values, or `NULL` if rule 1 does not apply), `use_rule1`/
#'   `use_rule2`/`use_rule3` (which SAS cap rule applies), `sid` (stratum
#'   ID per observation), `slevels0` (distinct stratum IDs), `exact_cols`
#'   (model-matrix columns to compute exact results for), and `cn`
#'   (`colnames(x)`).
#' @noRd
.setup_exact_model <- function(
    formula,
    data,
    offset,
    offset_log,
    exactmax,
    strata,
    exact_vars,
    ref_levels,
    verbose
) {
    if (is.null(strata)) {
        strata <- character(0L)
    }
    exm <- parse_exactmax_spec(exactmax)

    # exact_vars entries may name an interaction (e.g. "region:gender");
    # every individual variable inside such a term.
    atomic_exact_vars <- unique(trimws(unlist(strsplit(exact_vars, "[*:]"))))
    class_vars <- unique(c(
        strata,
        exact_vars,
        atomic_exact_vars,
        names(ref_levels)
    ))
    d <- prepare_sas_factors(data, class_vars, ref_levels)

    built <- build_mm_no_intercept(formula, d)
    x <- built$x
    y <- as.integer(stats::model.response(built$mf))
    row_idx <- as.integer(rownames(built$mf))

    use_offset <- has_offset_spec(offset)
    off <- resolve_offset(d, offset, offset_log, row_idx, required = FALSE)
    log_off <- off$log_offset

    exactmax_vec <- NULL
    if (exm$mode == "var") {
        v <- exm$var
        if (!v %in% names(d)) {
            stop(
                sprintf("EXACTMAX variable '%s' not found in data.", v),
                call. = FALSE
            )
        }
        exactmax_vec <- si_from_exactmax_var(d[[v]][row_idx])
    }

    use_rule1 <- (exm$mode == "var")
    use_rule2 <- (!use_rule1) && (exm$mode == "novar" || !use_offset)
    use_rule3 <- (!use_rule1) && (!use_rule2)

    if (verbose) {
        if (use_rule1) {
            message(
                "\nSAS rule (1): EXACTMAX=var specified => S_i = exactmax_var."
            )
        }
        if (use_rule2) {
            message(
                "\nSAS rule (2): EXACTMAX(no var) or no OFFSET/EXACTMAX =>
                S_i = intercept sufficient statistic (within stratum)."
            )
        }
        if (use_rule3) {
            message(
                "\nSAS rule (3): OFFSET specified and EXACTMAX not specified =>
                S_i = max(ceil(exp(o_i)),1)."
            )
        }
    }

    if (length(strata)) {
        sid <- build_stratum_id(d, strata, row_idx)
        slevels0 <- sort(unique(sid))
    } else {
        sid <- rep("ALL", length(y))
        slevels0 <- "ALL"
    }

    exact_cols <- exact_term_columns(x, built$terms, exact_vars)
    if (length(exact_cols) == 0L) {
        stop("No exact_vars columns found in model matrix.", call. = FALSE)
    }

    list(
        x = x,
        y = y,
        log_off = log_off,
        exactmax_vec = exactmax_vec,
        use_rule1 = use_rule1,
        use_rule2 = use_rule2,
        use_rule3 = use_rule3,
        sid = sid,
        slevels0 = slevels0,
        exact_cols = exact_cols,
        cn = colnames(x)
    )
}


#' Running "sum over strata not yet reached" for a Rule B-style bound
#'
#' Given one bound vector per stratum (e.g. each stratum's own forced
#' maximum contribution to every nuisance coordinate), returns, for each
#' stratum `s`, the sum of that bound over every stratum *after* `s` in
#' the processing order; the quantity needed to bound what the strata
#' still to come could contribute, on top of whatever the current
#' stratum contributes.
#'
#' @param vec_list A list of integer vectors, one per stratum, each of
#'   length `p_cond`.
#' @param p_cond Number of nuisance coordinates.
#' @return A list of the same length as `vec_list`.
#' @noRd
.cumulative_from_end <- function(vec_list, p_cond) {
    n <- length(vec_list)
    out <- vector("list", n)
    if (p_cond == 0L) {
        for (s in seq_len(n)) {
            out[[s]] <- integer(0L)
        }
        return(out)
    }
    running <- rep(0L, p_cond)
    if (n >= 1L) {
        out[[n]] <- running
    }
    if (n >= 2L) {
        for (s in (n - 1L):1L) {
            running <- running + vec_list[[s + 1L]]
            out[[s]] <- running
        }
    }
    out
}

#' Turn one parameter's final 1-D table into its output rows
#'
#' Computes the normalized probability table (`OUTDIST` entry), the exact
#' p-value, the conditional MLE, and the exact confidence interval from a
#' parameter's complete `(t_j, lw)` table, identical across all three
#' main functions once that table has been obtained.
#'
#' @param dist1d The complete 1-D coefficient table for this parameter.
#' @param parm The parameter's (model-matrix column) name.
#' @param y,x,jcol Used to recover the observed sufficient statistic.
#' @param alpha,ci_type As in the calling main function.
#' @return A list with `outdist_entry` (data frame for `OUTDIST[[parm]]`)
#'   and `result_row` (one row of `exactparmest`).
#' @noRd
.assemble_parameter_output <- function(
    dist1d,
    parm,
    y,
    x,
    jcol,
    alpha,
    ci_type
) {
    prob <- exp(dist1d$lw - log_sum_exp(dist1d$lw))
    outdist_entry <- data.frame(
        SuffStat = dist1d$t_j,
        prob = prob,
        LogWeight = dist1d$lw
    )

    t_obs <- as.integer(sum(y * as.integer(round(x[, jcol]))))

    tails <- two_sided_from_tails(
        p_lower = sum(prob[dist1d$t_j <= t_obs]),
        p_upper = sum(prob[dist1d$t_j >= t_obs])
    )
    pvalue <- if (ci_type == "two-sided") {
        tails$pvalue
    } else {
        tails$pvalue_onesided
    }

    cmle <- solve_cmle(dist1d$t_j, dist1d$lw, t_obs)
    ci <- exact_ci_from_dist1d(
        dist1d,
        t_obs,
        alpha = alpha,
        ci_type = ci_type,
        direction = tails$direction
    )

    result_row <- data.frame(
        Parameter = parm,
        Estimate = cmle$estimate,
        StdError = cmle$se,
        Alpha = alpha,
        LowerCL = ci[1L],
        UpperCL = ci[2L],
        PValue = pvalue,
        SuffStat = t_obs
    )

    list(
        outdist_entry = outdist_entry,
        result_row = result_row
    )
}
