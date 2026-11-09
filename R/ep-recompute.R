################################################################################
# recompute_exact_inference(): a small, exported convenience wrapper around
# the exponential-tilt machinery already shared by every algorithm in this
# package (tilt_probs(), two_sided_from_tails(), exact_ci_from_dist1d(), all
# defined in ep-utils.R/ep-common.R). It exists so
# that a *different* null hypothesis or a *different* confidence level can
# be answered directly from an already-fitted object's `OUTDIST` table,
# without repeating the combinatorial work (the per-stratum dynamic
# program and cross-stratum convolution, or the closed-form enumeration)
# that built it in the first place.
#
# Why this is correct without re-enumerating anything: every algorithm's
# reference distribution is a one-dimensional exponential-tilt family,
# Pr_beta(T = t) proportional to C(t) * exp(beta * t), and that family
# sums to one at *every* value of beta, not just at the null the fitting
# function happened to test. Once C(t) (the `LogWeight` column of an
# `OUTDIST` entry) has been built once, re-tilting it at a different beta0
# (tilt_probs()) or re-deriving a confidence limit at a different alpha
# (exact_ci_from_dist1d(), via exact_ci_from_tails()'s root-finding) reuses
# that same table; nothing about either calculation depends on which beta0
# or alpha the original fit used. See the theory vignette's discussion of
# reusing an already-built distribution for a different null/confidence
# level for the full derivation.
################################################################################

#' Recompute a p-value and/or confidence interval from an already-fitted
#' exact distribution
#'
#' Answers a *different* question about a parameter already fit by
#' [ep_multi_shift()], [ep_network()], or [ep_closed_form()]; a different
#' null value `beta0`, a different significance level `alpha`, or a
#' different `ci_type`/`direction`, directly from that fit's `OUTDIST`
#' table, without rebuilding the joint distribution. This is possible, and
#' exact (not an approximation), because every algorithm's reference
#' distribution is the same one-dimensional exponential-tilt family,
#' \eqn{\Pr_\beta(T = t) \propto C(t)\, e^{\beta t}}, which sums to one at
#' every value of \eqn{\beta}: re-tilting an already-built `C(t)` table at
#' a new `beta0`, or re-deriving a confidence limit at a new `alpha`, only
#' ever needs that one table, not a new enumeration over the model's
#' stratification and nuisance structure. Building the joint table in the
#' first place is the expensive part of every algorithm in this package;
#' this function is the cheap part, made reusable.
#'
#' @param fit A fitted object returned by [ep_multi_shift()],
#'   [ep_network()], or [ep_closed_form()] (any object of class
#'   `"ep_fit"`), i.e. anything with an `OUTDIST` list and an
#'   `exactparmest` data frame in the format those three functions return.
#' @param parm The name of one parameter already present in
#'   `fit$OUTDIST`/`fit$exactparmest$Parameter` (e.g. `"armvaccinated"`);
#'   matches the column name `exactparmest` reports for that term.
#' @param beta0 The null value of `beta` (the log rate ratio) to test
#'   against; `0` (the default) reproduces the same null the original fit
#'   itself tested. Supply a non-zero value to test a different margin,
#'   e.g. `log(0.5)` for a 50% rate-ratio margin.
#' @param alpha Significance level for the confidence interval; need not
#'   match the `alpha` the original fit used.
#' @param ci_type `"two-sided"` (default) or `"one-sided"`; need not match
#'   the original fit's own `ci_type`.
#' @param direction For `ci_type = "one-sided"`: `"upper"` reports only the
#'   upper confidence limit (lower limit fixed at `-Inf`); `"lower"`
#'   reports only the lower limit (upper limit fixed at `Inf`). Ignored
#'   for `ci_type = "two-sided"`.
#'
#' @return A list with:
#' * `beta0`: the null value tested, echoed back.
#' * `p_lower`: the lower tail probability at `beta0`,
#'   \eqn{\Pr_{\beta_0}(T \le t_{obs})}.
#' * `p_upper`: the upper tail probability at `beta0`,
#'   \eqn{\Pr_{\beta_0}(T \ge t_{obs})}.
#' * `pvalue_two_sided`: twice the smaller one-sided tail probability,
#'   capped at 1 (the package's two-sided exact p-value convention,
#'   matching SAS's `PROC GENMOD` `EXACT` statement).
#' * `pvalue_one_sided`: `min(p_lower, p_upper)`.
#' * `alpha`: the significance level used, echoed back.
#' * `LowerCL`, `UpperCL`: the recomputed confidence limits, on the same
#'   scale (`beta`, the log rate ratio) as `fit$exactparmest$Estimate`.
#'
#' @seealso [ep_multi_shift()], [ep_network()], [ep_closed_form()], and
#'   `vignette("Theory-Implementation")`'s discussion of reusing an already-built
#'   exact distribution.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   group = c(1L, 1L, 1L, 0L, 0L, 0L),
#'   age = c(50L, 60L, 70L, 50L, 60L, 70L),
#'   subjects = c(3492, 2141, 1711, 3525, 2166, 1724),
#'   cases = c(3, 2, 1, 10, 7, 5),
#'   futime = c(50, 60, 40, 51, 62, 42)
#' )
#' df$lnfutime <- log(df$futime)
#'
#' fit <- ep_multi_shift(
#'   formula = cases ~ group,
#'   data = df,
#'   offset = "lnfutime",
#'   offset_log = TRUE,
#'   strata = "age",
#'   exact_vars = "group",
#'   ref_levels = list(group = "0", age = "50"),
#'   verbose = FALSE
#' )
#'
#' # Same question the fit itself already answered (beta0 = 0, alpha = 0.05):
#' # pvalue_two_sided/LowerCL/UpperCL should match fit$exactparmest exactly.
#' recompute_exact_inference(fit, "group1")
#'
#' # A different confidence level, reusing the same OUTDIST table:
#' recompute_exact_inference(fit, "group1", alpha = 0.10)
recompute_exact_inference <- function(
    fit,
    parm,
    beta0 = 0.0,
    alpha = 0.05,
    ci_type = c("two-sided", "one-sided"),
    direction = c("upper", "lower")
) {
    ci_type <- rlang::arg_match(ci_type)
    direction <- rlang::arg_match(direction)

    .validate_recompute_args(fit, parm, beta0, alpha)
    looked_up <- .recompute_lookup(fit, parm)
    dist1d <- looked_up$dist1d
    t_vals <- dist1d$SuffStat
    lw <- dist1d$LogWeight
    t_obs <- looked_up$t_obs

    # Re-tilt the *same* table at beta0; nothing here re-enumerates the
    # reference set (see the file header, and vignette("Theory-Implementation")'s
    # discussion of Sections 6.5/8.4).
    p0 <- tilt_probs(beta0, t_vals, lw)
    p_lower <- sum(p0[t_vals <= t_obs])
    p_upper <- sum(p0[t_vals >= t_obs])
    tails <- two_sided_from_tails(p_lower, p_upper)

    ci <- exact_ci_from_dist1d(
        list(t_j = t_vals, lw = lw),
        t_obs,
        alpha = alpha,
        ci_type = ci_type,
        direction = direction
    )

    list(
        beta0 = beta0,
        p_lower = p_lower,
        p_upper = p_upper,
        pvalue_two_sided = tails$pvalue,
        pvalue_one_sided = tails$pvalue_onesided,
        alpha = alpha,
        LowerCL = ci[1L],
        UpperCL = ci[2L]
    )
}
