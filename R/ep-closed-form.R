################################################################################
# Exact conditional inference for a single two-level ("treatment" vs.
# "control") factor in a stratified Poisson regression, with strata as
# the *only* nuisance structure (no other categorical covariates). This
# is the special case of ep_multi_shift()
# (R/ep-multivariate-shift.R) for which the exact distribution
# of the sufficient statistic can be obtained in closed form, stratum by
# stratum, via the Poisson-thinning / binomial-expansion argument of
# Guess, Lydick, Small and Miller (1987) and Guess and Thomas (1990),
# instead of by the general per-observation dynamic program. See the
# accompanying theory vignette for exactly which properties of this
# special case make it faster (it is not simply "aggregation":
# ep_multi_shift() can also aggregate identical rows
# before running its dynamic program).
#
# All coefficient tables and tilted probabilities in this file are on
# the natural log scale, matching ep-utils.R's shared
# tilt_probs()/solve_cmle()/solve_tail_root() primitives.
################################################################################

# ------------------------------------------------------------------------
# Core recursive coefficient computation ----
# ------------------------------------------------------------------------

#' Build the (natural-log) coefficient table of the sufficient statistic
#'
#' The `log_ck` values are the natural-log coefficients required to build
#' the distribution of the group (treatment vs. control) sufficient
#' statistic, see Guess and Thomas (1990). Each stratum contributes a
#' closed-form binomial polynomial in the rate ratio theta = exp(beta),
#' and the coefficient sequence for the whole sample is obtained by
#' multiplying (convolving) these polynomials one stratum at a time.
#'
#' @param x A vector containing the number of cases in the treatment group
#'   for each stratum.
#' @param y A vector containing the number of cases in the control group
#'   for each stratum.
#' @param n A vector containing the person-time (exposure) in the treatment
#'   group for each stratum.
#' @param m A vector containing the person-time (exposure) in the control
#'   group for each stratum.
#' @param max_k The largest value of `k` (treatment-group case count summed
#'   over strata) to include in the table. Defaults to `sum(x)`, i.e. the
#'   *observed* value: this is enough to obtain every quantity needed for
#'   the exact test and the exact confidence interval, because their tail
#'   probabilities beyond the observed value can be recovered by
#'   complementation (the full distribution sums to one, at every value of
#'   beta, see `cumulative_dist_natlog()` below) without ever building
#'   the rest of the table. Pass `max_k = sum(x) + sum(y)` to obtain the
#'   complete table, which is what the conditional maximum likelihood
#'   estimate requires.
#'
#' @return A data frame with columns `k` and `log_ck` (natural log scale).
#'
#' @references Guess HA, Thomas JE. A rapidly converging algorithm
#' for exact binomial confidence intervals about the relative risk
#' in follow-up studies with stratified incidence-density data. Epidemiology.
#' 1990 Jan;1(1):75-7. doi: 10.1097/00001648-199001000-00016. PMID: 2081244.
#'
#' @noRd
compute_ck_natlog <- function(x, y, n, m, max_k = sum(x)) {
    ss <- length(x) # Number of strata
    tot <- x + y # Total cases per stratum
    # Cumulative total cases per stratum, capped at max_k
    max_cases <- as.integer(max_k)
    vvi <- pmin(max_cases, cumsum(tot))
    # Coefficient matrix, natural-log scale throughout
    coef_matrix <- matrix(NA_real_, nrow = max_cases + 1L, ncol = ss + 2L)
    coef_matrix[1L, 1L] <- 0.0 # Starting point
    for (s in 1L:ss) {
        v_i <- ifelse(s == 1L, 0L, vvi[s - 1L])
        for (k in 0L:vvi[s]) {
            if (s == ss) {
                coef_matrix[k + 1L, s + 2L] <- k
            }
            l_ik <- max(0L, k - v_i)
            h_ik <- min(k, tot[s])
            for (r in l_ik:h_ik) {
                log_multiplier <- coef_matrix[k - r + 1L, s] +
                    lchoose(tot[s], r) +
                    r * log(n[s] / m[s])
                delta <- max(
                    coef_matrix[k + 1L, s + 1L],
                    log_multiplier,
                    na.rm = TRUE
                )
                coef_matrix[k + 1L, s + 1L] <-
                    log(sum(
                        exp(coef_matrix[k + 1L, s + 1L] - delta),
                        exp(log_multiplier - delta),
                        na.rm = TRUE
                    )) +
                    delta
            }
        }
    }
    matrix(coef_matrix[, c(ss + 2L, ss + 1L)], ncol = 2L) |>
        data.frame() |>
        stats::setNames(c("k", "log_ck"))
}

#' Cumulative distribution of the sufficient statistic at a given beta
#'
#' See equation 1 in Guess and Thomas (1990). Reparametrized directly on
#' `beta` (the natural parameter / log rate ratio), matching the shared
#' exponential-tilt convention used by `ep_multi_shift()`
#' and by `ep-utils.R`'s `tilt_probs()`; the rate ratio itself is
#' `theta = exp(beta)`.
#'
#' Because `Pr_beta(T=t)` always sums to one over the *complete* support
#' of `T`, regardless of `beta`, this normalizing total can be obtained
#' analytically (`gg`, below) without ever building the complete table:
#' it is exactly the product, over strata, of `Pr_beta(no case falls in
#' the treatment group | stratum total)`. This is what lets
#' `ck_tab` be the truncated table (`max_k = sum(x)`, the default from
#' `compute_ck_natlog()`) for every quantity this function is used for.
#'
#' @param ck_tab The `k`/`log_ck` table from `compute_ck_natlog()`.
#' @param beta The natural parameter (log rate ratio) at which the
#'   distribution is evaluated.
#' @param equal `TRUE` for `P(K <= k_obs)`, `FALSE` for `P(K < k_obs)`,
#'   where `k_obs` is the largest `k` present in `ck_tab` (i.e. its
#'   truncation point).
#'
#' @return The requested cumulative probability.
#'
#' @noRd
cumulative_dist_natlog <- function(
    x,
    y,
    n,
    m,
    ck_tab,
    beta = 0.0,
    equal = TRUE
) {
    theta <- exp(beta)
    log_1mp <- log(m) - log(theta * n + m) # log(1 - p_theta), per stratum
    gg <- sum((x + y) * log_1mp)
    idx <- if (equal) seq_len(nrow(ck_tab)) else seq_len(nrow(ck_tab) - 1L)
    sum(exp(ck_tab$log_ck[idx] + beta * ck_tab$k[idx] + gg))
}

# ------------------------------------------------------------------------
# User-facing function ----
# ------------------------------------------------------------------------

#' Exact conditional Poisson regression (fast, closed-form algorithm)
#'
#' Computes the exact conditional maximum likelihood estimate, exact
#' confidence interval, and exact p-value for the coefficient of a single
#' two-level (treatment vs. control) factor in a Poisson regression with
#' log offset and, optionally, one or more stratification variables acting
#' purely as nuisance intercepts. This is the special case handled by
#' [ep_multi_shift()] when the model contains exactly one exact effect
#' with two levels and no other categorical covariates; in that special
#' case, the exact distribution can be obtained in closed form, stratum by
#' stratum, following Guess, Lydick, Small and Miller (1987) and Guess and
#' Thomas (1990), instead of by the general per-observation network
#' algorithm. See the accompanying theory note for why this is faster in
#' this scenario, and use [ep_multi_shift()] whenever the model contains
#' additional categorical covariates to adjust for, or more than one exact
#' effect.
#'
#' Unlike [ep_multi_shift()] and [ep_network()], this function has no
#' `exactmax` argument and never imposes a per-observation cap beyond what
#' conditioning on the stratum total already implies (SAS's own "Rule 2";
#' see the theory vignette, Sections 5.2 and 10). In particular, supplying
#' `offset` here does *not* trigger the offset-implied per-observation cap
#' ("Rule 3": \eqn{S_i = \max(\lceil e^{o_i} \rceil, 1)}) that the general
#' algorithms apply automatically whenever `offset` is supplied without
#' `exactmax`. When that cap would bind, small exposures relative to the
#' observed counts, typically in unstratified or lightly-stratified
#' two-group data, `ep_closed_form()` and the general algorithms are, by
#' design, answering subtly different questions, and can legitimately
#' disagree on the estimate, confidence interval, and p-value even given
#' identical call arguments. This is expected, documented behaviour, not a
#' bug in either algorithm; see the theory vignette for a worked numerical
#' example and `tests/testthat/test-cross-algorithm-agreement.R` for a
#' regression test that locks in the expected divergence.
#'
#' @param formula A model formula of the form `response ~ group`, where
#'   `group` is the two-level factor of interest and is the formula's only
#'   right-hand-side term. Strata and the offset are supplied separately
#'   (as in [ep_multi_shift()]), not as terms in this formula.
#' @param data A data frame. One row per observation; rows belonging to the
#'   same stratum and the same level of `group` are summed internally, so
#'   `data` may be individual-level or already aggregated.
#' @param offset Name of the person-time / exposure variable in `data`, or
#'   `NULL` (default). If `NULL`, every observation is treated as
#'   contributing one unit of exposure (offset = 0 on the log scale),
#'   matching [ep_multi_shift()]. Note that, unlike the general algorithms,
#'   supplying `offset` here never imposes a per-observation cap on any
#'   observation's count (see Details).
#' @param offset_log If `TRUE` (default, matching [ep_multi_shift()]),
#'   `offset` already holds log(exposure). If `FALSE`, `offset` holds raw
#'   exposure and is logged internally.
#' @param strata Optional character vector naming one or more
#'   stratification variables in `data`. If `NULL` (default), all
#'   observations form a single stratum, i.e. an unstratified two-sample
#'   exact Poisson test.
#' @param ref_levels Optional named list giving the reference (control)
#'   level of the formula's group variable, e.g. `list(trt = "0")`; the
#'   other level is treated as the treatment group. Defaults to R's usual
#'   first level.
#' @param alpha Type I error rate for the confidence interval (default
#'   0.05).
#' @param ci_type `"two-sided"` or `"one-sided"`.
#' @param verbose If `TRUE` (default), prints a short note identifying this
#'   as the closed-form two-group/strata-only algorithm.
#'
#' @return An object of class `c("ep_closed_form", "ep_fit")`: a list with
#'   elements `exactparmest` (a one-row data frame with `Parameter`,
#'   `Estimate` (the coefficient beta = log rate ratio, not the rate
#'   ratio), `StdError`, `Alpha`, `LowerCL`, `UpperCL` (beta scale),
#'   `PValue`, and `SuffStat`), `OUTDIST` (the exact distribution of the
#'   sufficient statistic), and `call`.
#'
#' @references Guess HA, Lydick EG, Small RD, Miller LP,
#'  Exact binomial confidence intervals for the relative risk in follow-up
#'  studies with sparsely stratified incidence density data, American Journal
#'  of Epidemiology, Volume 125, Issue 2, February 1987, Pages 340-347,
#'  https://doi.org/10.1093/oxfordjournals.aje.a114535
#' @references Guess HA, Thomas JE. A rapidly converging algorithm
#'  for exact binomial confidence intervals about the relative risk
#'  in follow-up studies with stratified incidence-density data. Epidemiology.
#'  1990 Jan;1(1):75-7. doi: 10.1097/00001648-199001000-00016. PMID: 2081244.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   cases = c(3, 32, 1, 19, 2, 24, 4, 40),
#'   fu = c(100, 100, 120, 120, 150, 150, 170, 170),
#'   trt = c(1, 0, 1, 0, 1, 0, 1, 0),
#'   strata = c(1, 1, 2, 2, 3, 3, 4, 4)
#' )
#' ep_closed_form(
#'   formula = cases ~ trt,
#'   data = df,
#'   offset = "fu",
#'   offset_log = FALSE,
#'   strata = "strata"
#' )
ep_closed_form <- function(
    formula,
    data,
    offset = NULL,
    offset_log = TRUE,
    strata = NULL,
    ref_levels = list(),
    alpha = 0.05,
    ci_type = c("two-sided", "one-sided"),
    verbose = TRUE
) {
    ci_type <- rlang::arg_match(ci_type)
    .validate_closedform_args(
        formula,
        data,
        offset,
        offset_log,
        strata,
        ref_levels,
        alpha,
        verbose
    )

    mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
    row_idx <- as.integer(rownames(mf))
    resp <- as.numeric(stats::model.response(mf))

    term_labels <- attr(stats::terms(mf), "term.labels")
    if (length(term_labels) != 1L) {
        stop(
            sprintf(
                "`formula` must contain exactly one right-hand-side term (the
                effect to test), with strata and the offset supplied separately
                via the `strata`/`offset` arguments, not as additional terms in
                `formula`; found term(s): %s. Use `ep_multi_shift()` or
                `ep_network()` for models with additional categorical
                covariates.",
                if (length(term_labels)) {
                    toString(term_labels)
                } else {
                    "(none)"
                }
            ),
            call. = FALSE
        )
    }
    grp_name <- term_labels[1L]
    grp <- apply_ref_level(data[[grp_name]][row_idx], grp_name, ref_levels)
    lev <- levels(droplevels(grp))
    if (length(lev) != 2L) {
        stop(
            sprintf(
                "`%s` must have exactly two levels (found %d).",
                grp_name,
                length(lev)
            ),
            call. = FALSE
        )
    }
    z <- as.integer(grp == lev[2L]) # 1 = treatment (non-ref), 0 = control (ref)

    off <- resolve_offset(data, offset, offset_log, row_idx, required = FALSE)
    expo <- off$exposure

    sid <- build_stratum_id(data, strata, row_idx)
    slevels <- sort(unique(sid))
    ss <- length(slevels)

    x <- y <- n <- m <- numeric(ss)
    for (i in seq_len(ss)) {
        idx <- which(sid == slevels[i])
        x[i] <- sum(resp[idx][z[idx] == 1L])
        y[i] <- sum(resp[idx][z[idx] == 0L])
        n[i] <- sum(expo[idx][z[idx] == 1L])
        m[i] <- sum(expo[idx][z[idx] == 0L])
    }
    if (any(n <= 0.0) || any(m <= 0.0)) {
        stop(
            "Every stratum must contribute positive exposure in both levels ",
            "of the formula's group variable; a stratum with no exposure in ",
            "one arm carries no information about the rate ratio and cannot ",
            "be used by this method.",
            call. = FALSE
        )
    }

    if (verbose) {
        message(
            sprintf(
                "\nClosed-form two-group/strata-only exact Poisson algorithm
                (Guess & Thomas, 1990): %d stratum(strata), %d treatment cases,
                %d control cases.",
                ss,
                sum(x),
                sum(y)
            )
        )
    }

    k_obs <- sum(x)
    k_max <- sum(x) + sum(y)

    ck_trunc <- compute_ck_natlog(x, y, n, m) # fast, truncated at k_obs

    # Tail probabilities at an arbitrary beta, using only the truncated
    # table plus the analytic normalizing total (see cumulative_dist_natlog()).
    tail_le <- function(beta) {
        cumulative_dist_natlog(x, y, n, m, ck_trunc, beta, equal = TRUE)
    }
    tail_ge <- function(beta) {
        1.0 - cumulative_dist_natlog(x, y, n, m, ck_trunc, beta, equal = FALSE)
    }

    p_lower <- tail_le(0.0)
    p_upper <- tail_ge(0.0)
    tails <- two_sided_from_tails(p_lower, p_upper)
    pvalue <- if (ci_type == "two-sided") {
        tails$pvalue
    } else {
        tails$pvalue_onesided
    }

    ci <- exact_ci_from_tails(
        tail_le,
        tail_ge,
        alpha,
        ci_type,
        tails$direction,
        at_lower_boundary = (k_obs <= 0.0),
        at_upper_boundary = (k_obs >= k_max)
    )
    beta_ll <- ci[1L]
    beta_ul <- ci[2L]

    ck_full <- compute_ck_natlog(x, y, n, m, max_k = k_max)
    cmle <- solve_cmle(ck_full$k, ck_full$log_ck, k_obs)

    parm_label <- paste0(grp_name, lev[2L])
    exactparmest <- data.frame(
        Parameter = parm_label,
        Estimate = cmle$estimate,
        StdError = cmle$se,
        Alpha = alpha,
        LowerCL = beta_ll,
        UpperCL = beta_ul,
        PValue = pvalue,
        SuffStat = k_obs
    )

    prob0 <- tilt_probs(0.0, ck_full$k, ck_full$log_ck)
    outdist <- list()
    outdist[[parm_label]] <- data.frame(
        SuffStat = ck_full$k,
        prob = prob0,
        LogWeight = ck_full$log_ck
    )

    ep_fit(
        exactparmest = exactparmest,
        outdist = outdist,
        call = match.call(),
        subclass = "ep_closed_form",
        algorithm_label = paste(
            "Exact conditional Poisson, two groups + strata",
            "(closed-form Guess & Thomas algorithm)"
        )
    )
}
