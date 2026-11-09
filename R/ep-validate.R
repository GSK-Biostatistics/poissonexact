# ==============================================================================
# Shared argument validation, used by every exported main function
# (ep_multi_shift(), ep_closed_form(), ep_network()).
#
# None of the functions in this file are exported. They exist so that
# each main function's own body can start with one short call to a
# single per-function validator (e.g. `.validate_exact_poisson_args(...)`)
# instead of repeating the same handful of type/shape checks inline;
# each validator is built out of these small, single-purpose checkers,
# one per kind of argument, so that a check only needs to be written
# (and reviewed) once even though several main functions share the same
# argument. Every check here is written to match what that argument's
# `@param` documentation promises in the function that uses it; e.g.
# "Optional character vector of stratification variable names in `data`"
# is checked as exactly that: NULL, or a character vector whose entries
# are all non-missing, non-empty, and present as columns of `data`.
# ==============================================================================

# ------------------------------------------------------------------------
# Single-purpose checkers, one per kind of argument ----
# ------------------------------------------------------------------------

#' Validate a scalar `TRUE`/`FALSE` argument
#'
#' Shared by every main function for its own logical flag arguments
#' (e.g. `offset_log`, `verbose`), so that all of them report the same
#' error message for the same kind of mistake.
#'
#' @param x The value to check.
#' @param name The argument's name, used in the error message.
#' @return `invisible(TRUE)` if `x` is a single, non-`NA` logical value;
#'   otherwise raises an error.
#' @noRd
.check_scalar_flag <- function(x, name) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
        stop(
            sprintf(
                "`%s` must be a single `TRUE` or `FALSE` (not `NA`).",
                name
            ),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Validate the `formula` argument
#'
#' No `missing(formula)` special-case here: `formula` has no default in
#' every main function that calls this, so simply touching it below
#' (via `inherits()`) already forces R's own clear "argument ... is
#' missing, with no default" error if it was never supplied and that is
#' enough, without relying on `missing()` propagating correctly through
#' this validator's own extra layer of indirection.
#'
#' @param formula The `formula` argument as supplied by the caller.
#' @return `invisible(TRUE)` if `formula` is a model formula; otherwise
#'   raises an error.
#' @noRd
.check_formula_arg <- function(formula) {
    # No `missing(formula)` special-case here: `formula` has no default in
    # every main function that calls this, so simply touching it below
    # (via `inherits()`) already forces R's own clear "argument ... is
    # missing, with no default" error if it was never supplied and that is
    # enough, without relying on `missing()` propagating correctly through
    # this validator's own extra layer of indirection.
    if (!inherits(formula, "formula")) {
        stop(
            "`formula` must be a model formula, e.g. `cases ~ group`.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Validate the `data` argument
#'
#' @param data The `data` argument as supplied by the caller.
#' @return `invisible(TRUE)` if `data` is a data frame with at least one
#'   row; otherwise raises an error.
#' @noRd
.check_data_arg <- function(data) {
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.", call. = FALSE)
    }
    if (nrow(data) == 0L) {
        stop("`data` has zero rows.", call. = FALSE)
    }
    invisible(TRUE)
}

#' Validate an `offset` argument
#'
#' Matches `has_offset_spec()`'s own notion of "no offset supplied"
#' (`NULL`, or a zero-length/empty-string character) so this check never
#' rejects an offset value that the rest of the package treats as valid.
#'
#' @param offset The `offset` argument as supplied by the caller: `NULL`,
#'   `""`, or the name of a column in `data`.
#' @param data A data frame, checked for a column named `offset` when one
#'   is supplied.
#' @param required If `TRUE`, an unspecified offset is an error rather
#'   than being treated as "no offset".
#' @return `invisible(TRUE)` if valid; otherwise raises an error.
#' @noRd
.check_offset_name_arg <- function(offset, data, required) {
    no_spec <- is.null(offset) ||
        (is.character(offset) && length(offset) == 1L && !nzchar(offset))
    if (no_spec) {
        if (required) {
            stop(
                "`offset` (name of the person-time / exposure variable) is required.", # nolint: nonportable_path_linter
                call. = FALSE
            )
        }
        return(invisible(TRUE))
    }
    if (
        !is.character(offset) ||
            length(offset) != 1L ||
            is.na(offset) ||
            !nzchar(offset)
    ) {
        stop(
            "`offset` must be a single non-empty character string naming
            a column in `data`, or `NULL`.",
            call. = FALSE
        )
    }
    if (!offset %in% names(data)) {
        stop(
            sprintf("`offset` = '%s' not found in `data`.", offset),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Validate the `strata` argument
#'
#' @param strata The `strata` argument as supplied by the caller: `NULL`
#'   or a character vector of column names.
#' @param data A data frame, checked for every column named in `strata`.
#' @return `invisible(TRUE)` if valid; otherwise raises an error.
#' @noRd
.check_strata_arg <- function(strata, data) {
    if (is.null(strata)) {
        return(invisible(TRUE))
    }
    if (!is.character(strata) || anyNA(strata) || !all(nzchar(strata))) {
        stop(
            "`strata` must be `NULL` or a character vector of non-empty
            column names.",
            call. = FALSE
        )
    }
    missing_cols <- setdiff(strata, names(data))
    if (length(missing_cols)) {
        stop(
            sprintf(
                "`strata` names column(s) not found in `data`: %s.",
                toString(missing_cols)
            ),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Validate the `exact_vars` argument
#'
#' @param exact_vars The `exact_vars` argument as supplied by the caller.
#' @return `invisible(TRUE)` if `exact_vars` is a non-empty character
#'   vector of non-missing, non-empty term names; otherwise raises an
#'   error.
#' @noRd
.check_exact_vars_arg <- function(exact_vars) {
    if (
        !length(exact_vars) ||
            !is.character(exact_vars) ||
            anyNA(exact_vars) ||
            !all(nzchar(exact_vars))
    ) {
        stop(
            "`exact_vars` must be a non-empty character vector of model term names.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Validate the `ref_levels` argument
#'
#' @param ref_levels The `ref_levels` argument as supplied by the caller:
#'   a possibly-empty named list.
#' @param data A data frame, checked for every variable named in
#'   `ref_levels`.
#' @return `invisible(TRUE)` if valid; otherwise raises an error.
#' @noRd
.check_ref_levels_arg <- function(ref_levels, data) {
    if (!is.list(ref_levels)) {
        stop(
            "`ref_levels` must be a (possibly empty) named list.",
            call. = FALSE
        )
    }
    if (!length(ref_levels)) {
        return(invisible(TRUE))
    }
    nm <- names(ref_levels)
    if (is.null(nm) || anyNA(nm) || !all(nzchar(nm))) {
        stop(
            "Every element of `ref_levels` must be named for the variable it
            sets a reference level for.",
            call. = FALSE
        )
    }
    bad_len <- lengths(ref_levels) != 1L
    if (any(bad_len)) {
        stop(
            sprintf(
                "Each `ref_levels` entry must be a single value; check: %s.",
                toString(nm[bad_len])
            ),
            call. = FALSE
        )
    }
    missing_vars <- setdiff(nm, names(data))
    if (length(missing_vars)) {
        stop(
            sprintf(
                "`ref_levels` names variable(s) not found in `data`: %s.",
                toString(missing_vars)
            ),
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Validate the `alpha` argument
#'
#' @param alpha The `alpha` argument as supplied by the caller.
#' @return `invisible(TRUE)` if `alpha` is a single numeric value
#'   strictly between 0 and 1; otherwise raises an error.
#' @noRd
.check_alpha_arg <- function(alpha) {
    if (
        !is.numeric(alpha) ||
            length(alpha) != 1L ||
            is.na(alpha) ||
            alpha <= 0.0 ||
            alpha >= 1.0
    ) {
        stop(
            "`alpha` must be a single numeric value strictly between 0 and 1.",
            call. = FALSE
        )
    }
    invisible(TRUE)
}

#' Validate the `chunk_size` argument
#'
#' @param chunk_size The `chunk_size` argument as supplied by the caller.
#' @return `invisible(TRUE)` if `chunk_size` is a single positive
#'   integer; otherwise raises an error.
#' @noRd
.check_chunk_size_arg <- function(chunk_size) {
    if (
        !is.numeric(chunk_size) ||
            length(chunk_size) != 1L ||
            is.na(chunk_size) ||
            chunk_size != as.integer(chunk_size) ||
            chunk_size < 1L
    ) {
        stop("`chunk_size` must be a single positive integer.", call. = FALSE)
    }
    invisible(TRUE)
}

# ------------------------------------------------------------------------
# Per-function validators ----
# one call, at the top of each main function's body, in place of
# repeating the checks above inline.
# ------------------------------------------------------------------------

#' Validate arguments shared by ep_multi_shift() and ep_network()
#'
#' All three functions take the same core argument list (see
#' [ep_multi_shift()]'s documentation for what each
#' argument means), so they share this one validator.
#' `ep_network()` has no `chunk_size` argument (it never
#' materializes a cross-join between two strata's tables for that argument
#' to bound); pass `chunk_size = NULL` to skip that one check.
#'
#' @inheritParams ep_multi_shift
#' @param exactmax As in [ep_multi_shift()]; not
#'   itself re-checked here (validated where it is parsed, by
#'   `parse_exactmax_spec()`).
#' @param chunk_size As in [ep_multi_shift()], or
#'   `NULL` to skip this one check (used by [ep_network()],
#'   which has no `chunk_size` argument of its own).
#' @return `invisible(TRUE)` if every argument is valid; otherwise raises
#'   an error from whichever single-purpose checker fails first.
#' @noRd
.validate_exact_poisson_args <- function(
    formula,
    data,
    offset,
    offset_log,
    exactmax,
    strata,
    exact_vars,
    ref_levels,
    alpha,
    verbose,
    chunk_size = NULL
) {
    .check_formula_arg(formula)
    .check_data_arg(data)
    .check_offset_name_arg(offset, data, required = FALSE)
    .check_scalar_flag(offset_log, "offset_log")
    .check_strata_arg(strata, data)
    .check_exact_vars_arg(exact_vars)
    .check_ref_levels_arg(ref_levels, data)
    .check_alpha_arg(alpha)
    if (!is.null(chunk_size)) {
        .check_chunk_size_arg(chunk_size)
    }
    .check_scalar_flag(verbose, "verbose")
    invisible(TRUE)
}

#' Validate arguments for ep_closed_form()
#'
#' Shares most of its checks with `.validate_exact_poisson_args()`.
#' `ep_closed_form()` has no `exact_vars` argument of its own (the single
#' effect to test is read directly off `formula`'s one right-hand-side
#' term), so unlike `.validate_exact_poisson_args()` this validator does
#' not call `.check_exact_vars_arg()`; the "formula has exactly one term"
#' requirement (Section 6 of the theory vignette) is instead checked by
#' `ep_closed_form()` itself, once `formula`'s terms are available.
#'
#' @inheritParams ep_closed_form
#' @return `invisible(TRUE)` if every argument is valid; otherwise raises
#'   an error.
#' @noRd
.validate_closedform_args <- function(
    formula,
    data,
    offset,
    offset_log,
    strata,
    ref_levels,
    alpha,
    verbose
) {
    .check_formula_arg(formula)
    .check_data_arg(data)
    .check_offset_name_arg(offset, data, required = FALSE)
    .check_scalar_flag(offset_log, "offset_log")
    .check_strata_arg(strata, data)
    .check_ref_levels_arg(ref_levels, data)
    .check_alpha_arg(alpha)
    .check_scalar_flag(verbose, "verbose")
    invisible(TRUE)
}

#' Validate arguments for recompute_exact_inference()
#'
#' Checks the four arguments that must be valid before re-tilting an
#' already-built `OUTDIST` table: that `fit` is a fit this package
#' produced, that `parm` names a parameter actually present in it (with
#' exactly one matching `exactparmest` row), and that `beta0`/`alpha` are
#' well-formed. Shares its `alpha` check with every main function via
#' `.check_alpha_arg()`, rather than repeating that logic a fourth time.
#'
#' @param fit The `fit` argument as supplied by the caller.
#' @param parm The `parm` argument as supplied by the caller.
#' @param beta0 The `beta0` argument as supplied by the caller.
#' @inheritParams recompute_exact_inference
#' @return `invisible(NULL)` if every argument is valid; otherwise raises
#'   an error from whichever check fails first. Does not extract or
#'   otherwise process `fit`; see `.recompute_lookup()` for that.
#' @noRd
.validate_recompute_args <- function(fit, parm, beta0, alpha) {
    if (!inherits(fit, "ep_fit")) {
        stop(
            "`fit` must be an object returned by ep_multi_shift(), ",
            "ep_network(), or ep_closed_form().",
            call. = FALSE
        )
    }
    if (!is.character(parm) || length(parm) != 1L || !nzchar(parm)) {
        stop("`parm` must be a single non-empty parameter name.", call. = FALSE)
    }
    if (is.null(fit$OUTDIST[[parm]])) {
        stop(
            sprintf(
                "'%s' is not a parameter in fit$OUTDIST. Available: %s",
                parm,
                toString(names(fit$OUTDIST))
            ),
            call. = FALSE
        )
    }
    if (!is.numeric(beta0) || length(beta0) != 1L || !is.finite(beta0)) {
        stop("`beta0` must be a single finite numeric value.", call. = FALSE)
    }
    .check_alpha_arg(alpha)

    n_rows <- sum(fit$exactparmest$Parameter == parm)
    if (n_rows != 1L) {
        stop(
            sprintf(
                "Expected exactly one exactparmest row for '%s'; found %d.",
                parm,
                n_rows
            ),
            call. = FALSE
        )
    }
    invisible(NULL)
}

#' Look up a parameter's one-dimensional exact distribution and observed
#' sufficient statistic from a fitted `ep_fit` object
#'
#' Pure extraction, no validation: callers must run
#' `.validate_recompute_args()` first, which is what guarantees `parm`
#' names exactly one row of `fit$exactparmest` and one entry of
#' `fit$OUTDIST` for this function to pull out.
#'
#' @param fit A validated `ep_fit` object.
#' @param parm A validated parameter name.
#' @return A list with `dist1d` (the `fit$OUTDIST[[parm]]` table) and
#'   `t_obs` (the observed sufficient statistic for `parm`).
#' @noRd
.recompute_lookup <- function(fit, parm) {
    list(
        dist1d = fit$OUTDIST[[parm]],
        t_obs = fit$exactparmest$SuffStat[fit$exactparmest$Parameter == parm]
    )
}
