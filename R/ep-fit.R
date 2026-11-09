#' Construct a fitted-object of (base) class `"ep_fit"`, without checking it
#'
#' The bare constructor: builds the list and attaches the class/attribute,
#' nothing else. Does not check that `exactparmest` has the required
#' columns, that `outdist` matches it, or that `call` is really a call;
#' see `ep_fit()` for the checked entry point every algorithm in this
#' package actually uses. `new_ep_fit()` itself exists for callers that
#' already know their input is well-formed and want to build (or
#' deliberately build a malformed) `"ep_fit"` object without paying for,
#' or being blocked by, those checks - in this package, that is the test
#' suite's own fixtures for exercising `print.ep_fit()`.
#'
#' @param exactparmest A data frame, intended to have one row per
#'   requested parameter and, at least, the columns `Parameter`,
#'   `Estimate`, `StdError`, `Alpha`, `LowerCL`, `UpperCL`, `PValue`,
#'   `SuffStat`, but not checked here.
#' @param outdist Intended to be a named list, one entry per parameter
#'   named in `exactparmest$Parameter`, each a data frame with that
#'   parameter's exact distribution, but not checked here.
#' @param call Intended to be the calling function's `match.call()`, but
#'   not checked here.
#' @param subclass Character vector of one or more algorithm-specific
#'   classes, most specific first (e.g. `"ep_closed_form"`, or
#'   `c("ep_network", "ep_multi_shift")`); the shared base class
#'   `"ep_fit"` is appended automatically and must not be included here.
#' @param algorithm_label A one-line label used by `print.ep_fit()`.
#'
#' @return An object of class `c(subclass, "ep_fit")`: a list with
#'   `exactparmest`, `OUTDIST`, and `call`, carrying `algorithm_label` as
#'   an attribute for `print.ep_fit()`.
#'
#' @noRd
new_ep_fit <- function(exactparmest, outdist, call, subclass, algorithm_label) {
    structure(
        list(
            exactparmest = exactparmest,
            OUTDIST = outdist,
            call = call
        ),
        class = c(subclass, "ep_fit"),
        algorithm_label = algorithm_label
    )
}

#' Check that an about-to-be-built `"ep_fit"` object's pieces are well-formed
#'
#' Validation only: does not build or return the fit object itself, does
#' not modify its arguments, and has no side effects beyond raising an
#' error from whichever check fails first. Shared by `ep_fit()` so the
#' checks are written once.
#'
#' A failure here always indicates a bug inside this package (a caller
#' assembling `exactparmest`/`outdist` incorrectly), never a user input
#' mistake; user input is already validated earlier, by
#' `.validate_exact_poisson_args()`/`.validate_closedform_args()` and by
#' `exact_term_columns()`.
#'
#' @inheritParams new_ep_fit
#' @return `invisible(NULL)` if every piece is well-formed; otherwise
#'   raises an error from whichever check fails first.
#' @noRd
validate_ep_fit <- function(exactparmest, outdist, call) {
    if (!is.data.frame(exactparmest)) {
        stop("`exactparmest` must be a data frame.", call. = FALSE)
    }
    required_cols <- c(
        "Parameter",
        "Estimate",
        "StdError",
        "Alpha",
        "LowerCL",
        "UpperCL",
        "PValue",
        "SuffStat"
    )
    missing_cols <- setdiff(required_cols, names(exactparmest))
    if (length(missing_cols)) {
        stop(
            sprintf(
                "`exactparmest` is missing column(s): %s.",
                toString(missing_cols)
            ),
            call. = FALSE
        )
    }
    if (
        !is.list(outdist) ||
            !identical(sort(names(outdist)), sort(exactparmest$Parameter))
    ) {
        stop(
            "`outdist` must be a named list with exactly one ",
            "entry per parameter in `exactparmest$Parameter`.",
            call. = FALSE
        )
    }
    if (!is.call(call)) {
        stop("`call` must be a call object.", call. = FALSE)
    }
    invisible(NULL)
}

#' Construct and check a fitted-object of (base) class `"ep_fit"`
#'
#' The checked entry point behind every algorithm's returned fit object:
#' [ep_closed_form()], [ep_multi_shift()], and [ep_network()] all end by
#' calling this instead of building their result with `structure()` (or
#' the unchecked `new_ep_fit()`) by hand, so that the fitted-object's
#' shape - `exactparmest`'s required columns, and `outdist` having
#' exactly one entry per parameter - is validated in exactly one place,
#' and every algorithm-specific class sits under the same shared base
#' class that `print.ep_fit()` and `recompute_exact_inference()`
#' dispatch/check against.
#'
#' @inheritParams new_ep_fit
#' @return An object of class `c(subclass, "ep_fit")`; see `new_ep_fit()`.
#' @noRd
ep_fit <- function(exactparmest, outdist, call, subclass, algorithm_label) {
    validate_ep_fit(exactparmest, outdist, call)
    new_ep_fit(exactparmest, outdist, call, subclass, algorithm_label)
}

#' @exportS3Method
print.ep_fit <- function(x, ...) {
    label <- attr(x, "algorithm_label") %||%
        "Exact conditional Poisson regression"
    cat("\n", label, "\n", sep = "")
    cat(strrep("-", nchar(label)), "\n", sep = "")
    print(x$exactparmest, row.names = FALSE)
    invisible(x)
}
