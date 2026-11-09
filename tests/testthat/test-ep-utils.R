# ==============================================================================
# Tests for exact-poisson-utils.R; shared internal helpers used by all
# three algorithms (log_sum_exp, log_sum_exp_by, make_cond_cols,
# build_stratum_id, apply_ref_level, has_offset_spec, resolve_offset,
# tilt_probs, solve_cmle, solve_tail_root, two_sided_from_tails,
# exact_ci_from_tails). print.ep_fit()'s own tests live in
# tests/testthat/test-ep-fit.R, next to R/ep-fit.R.
#
# These are internal (non-exported) functions, reachable directly by name
# inside the package's test environment (see tests/testthat.R).
# ==============================================================================

# ------------------------------------------------------------------------
# log_sum_exp ----
# ------------------------------------------------------------------------

test_that("log_sum_exp matches direct computation for ordinary values", {
    v <- log(c(2, 3))
    expect_equal(log_sum_exp(v), log(5))
})

test_that("log_sum_exp handles a single value", {
    expect_equal(log_sum_exp(3.7), 3.7)
})

test_that("log_sum_exp returns -Inf for a zero-length vector", {
    expect_equal(log_sum_exp(numeric(0)), -Inf)
})

test_that("log_sum_exp returns -Inf when every element is -Inf", {
    expect_equal(log_sum_exp(c(-Inf, -Inf, -Inf)), -Inf)
})

test_that("log_sum_exp is numerically stable for very large log-values", {
    # A naive log(sum(exp(v))) would overflow to Inf for v this large.
    v <- c(1000, 1000.5, 999)
    expect_true(is.finite(log_sum_exp(v)))
    # Shift-invariance check: log_sum_exp(v) - k should equal
    # log_sum_exp(v - k) for any constant k.
    expect_equal(log_sum_exp(v) - 1000, log_sum_exp(v - 1000))
})

test_that("log_sum_exp propagates +Inf", {
    expect_equal(log_sum_exp(c(1, Inf, 2)), Inf)
})

# ------------------------------------------------------------------------
# log_sum_exp_by ----
# ------------------------------------------------------------------------

test_that("log_sum_exp_by combines rows sharing a single grouping column", {
    dt <- list(
        t_j = c(1L, 1L, 2L),
        lw = log(c(2, 3, 5))
    )
    out <- log_sum_exp_by(dt, "t_j")
    out <- order_tbl(out, "t_j")
    expect_equal(out$t_j, c(1L, 2L))
    expect_equal(out$lw, log(c(5, 5)))
})

test_that("log_sum_exp_by combines rows sharing a multi-column key", {
    dt <- list(
        m_j = c(1L, 1L, 1L, 2L),
        t_j = c(0L, 0L, 1L, 0L),
        lw = log(c(1, 1, 4, 7))
    )
    out <- log_sum_exp_by(dt, c("m_j", "t_j"))
    out <- order_tbl(out, c("m_j", "t_j"))
    expect_identical(.tbl_nrow(out), 3L)
    expect_equal(
        list(
            lw_1_0 = out$lw[out$m_j == 1L & out$t_j == 0L],
            lw_1_1 = out$lw[out$m_j == 1L & out$t_j == 1L],
            lw_2_0 = out$lw[out$m_j == 2L & out$t_j == 0L]
        ),
        list(lw_1_0 = log(2), lw_1_1 = log(4), lw_2_0 = log(7))
    )
})

test_that("log_sum_exp_by respects a custom lw_col name", {
    dt <- list(t_j = c(1L, 1L), weight = log(c(2, 2)))
    out <- log_sum_exp_by(dt, "t_j", lw_col = "weight")
    expect_equal(out$weight, log(4))
})

test_that("log_sum_exp_by leaves singleton groups unchanged", {
    dt <- list(t_j = c(1L, 2L, 3L), lw = log(c(1, 2, 3)))
    out <- log_sum_exp_by(dt, "t_j")
    out <- order_tbl(out, "t_j")
    expect_equal(out$lw, log(c(1, 2, 3)))
})

test_that("log_sum_exp_by returns an empty table unchanged for zero-row input", {
    # covr::zero_coverage() flagged this early-return branch (n == 0L) as
    # never hit: every other test above feeds it a non-empty table, but
    # both .joint_dist_stratum_dp_vec() (after Rule A/tn_obs pruning empties
    # a mid-stratum frontier) and .convolve_joint_pnuis_vec() (a chunk with
    # no surviving pairs) can and do call this with a zero-row table in
    # normal operation, so it needs its own direct test.
    dt <- list(m_j = integer(0), t_j = integer(0), lw = numeric(0))
    out <- log_sum_exp_by(dt, c("m_j", "t_j"))
    expect_identical(out$m_j, integer(0))
    expect_identical(out$t_j, integer(0))
    expect_identical(out$lw, numeric(0))
})

test_that("log_sum_exp_by returns an empty table unchanged with no nuisance cols", {
    dt <- list(t_j = integer(0), lw = numeric(0))
    out <- log_sum_exp_by(dt, "t_j")
    expect_identical(out$t_j, integer(0))
    expect_identical(out$lw, numeric(0))
})

# ------------------------------------------------------------------------
# .tbl_nrow / .tbl_subset / .tbl_rbind ----
# ------------------------------------------------------------------------

test_that(".tbl_nrow returns 0L for NULL", {
    # Flagged as never hit: every other caller in the package always passes
    # an actual list (possibly with zero-length elements, which is a
    # different case, see below), never a literal NULL. .tbl_nrow()'s own
    # documented contract explicitly covers NULL, so it gets its own test.
    expect_identical(.tbl_nrow(NULL), 0L)
})

test_that(".tbl_nrow returns 0L for an empty list", {
    expect_identical(.tbl_nrow(list()), 0L)
})

test_that(".tbl_nrow returns the common column length for a non-empty table", {
    expect_identical(.tbl_nrow(list(t_j = c(1L, 2L, 3L), lw = c(0, 0, 0))), 3L)
})

test_that(".tbl_nrow returns 0L for a table whose columns are zero-length vectors", {
    # Distinct from the NULL/empty-list cases above: here `tbl` itself is a
    # non-empty *list* (length(tbl) > 0), but every element it contains has
    # length zero, the `length(tbl[[1L]])` branch, not the
    # `is.null(tbl) || length(tbl) == 0L` shortcut.
    expect_identical(.tbl_nrow(list(t_j = integer(0), lw = numeric(0))), 0L)
})

test_that(".tbl_subset indexes every column identically, including repeats", {
    tbl <- list(t_j = c(10L, 20L, 30L), lw = c(0.1, 0.2, 0.3))
    out <- .tbl_subset(tbl, c(2L, 2L, 1L))
    expect_identical(out$t_j, c(20L, 20L, 10L))
    expect_identical(out$lw, c(0.2, 0.2, 0.1))
})

test_that(".tbl_subset works with a logical index", {
    tbl <- list(t_j = c(10L, 20L, 30L), lw = c(0.1, 0.2, 0.3))
    out <- .tbl_subset(tbl, c(TRUE, FALSE, TRUE))
    expect_identical(out$t_j, c(10L, 30L))
})

test_that(".tbl_rbind concatenates matching columns in tbl1's column order", {
    tbl1 <- list(t_j = 1L, lw = 0.1)
    tbl2 <- list(lw = 0.2, t_j = 2L) # deliberately different column order
    out <- .tbl_rbind(tbl1, tbl2)
    expect_named(out, c("t_j", "lw"))
    expect_identical(out$t_j, c(1L, 2L))
    expect_identical(out$lw, c(0.1, 0.2))
})

# ------------------------------------------------------------------------
# make_cond_cols ----
# ------------------------------------------------------------------------

test_that("make_cond_cols builds the expected names for a positive count", {
    expect_identical(make_cond_cols(3), c("cond_1", "cond_2", "cond_3"))
    expect_identical(make_cond_cols(1), "cond_1")
})

test_that("make_cond_cols returns character(0) for zero nuisance coordinates", {
    expect_identical(make_cond_cols(0), character(0))
})

test_that("make_cond_cols treats NULL as zero", {
    expect_identical(make_cond_cols(NULL), character(0))
})

test_that("make_cond_cols falls back to zero for invalid input", {
    expect_identical(make_cond_cols(NA), character(0))
    expect_identical(make_cond_cols(-1), character(0))
    expect_identical(make_cond_cols(c(1, 2)), character(0))
})

# ------------------------------------------------------------------------
# build_stratum_id ----
# ------------------------------------------------------------------------

test_that("build_stratum_id returns 'ALL' for every row when strata is NULL", {
    d <- data.frame(x = 1:4)
    expect_identical(build_stratum_id(d, NULL, 1:4), rep("ALL", 4))
})

test_that("build_stratum_id returns 'ALL' for every row when strata is empty", {
    d <- data.frame(x = 1:3)
    expect_identical(build_stratum_id(d, character(0), 1:3), rep("ALL", 3))
})

test_that("build_stratum_id pastes a single stratification column", {
    d <- data.frame(age = c(50L, 60L, 50L))
    expect_identical(build_stratum_id(d, "age", 1:3), c("50", "60", "50"))
})

test_that("build_stratum_id pastes multiple stratification columns with '__'", {
    d <- data.frame(group = c(1L, 1L, 0L), age = c(50L, 60L, 50L))
    out <- build_stratum_id(d, c("group", "age"), 1:3)
    expect_identical(out, c("1__50", "1__60", "0__50"))
})

test_that("build_stratum_id respects row_idx subsetting/reordering", {
    d <- data.frame(age = c(50L, 60L, 70L))
    expect_identical(build_stratum_id(d, "age", c(3L, 1L)), c("70", "50"))
})

# ------------------------------------------------------------------------
# apply_ref_level ----
# ------------------------------------------------------------------------

test_that("apply_ref_level coerces a non-factor to a factor when no ref given", {
    out <- apply_ref_level(c("a", "b", "a"), "x", list())
    expect_s3_class(out, "factor")
    expect_identical(levels(out), c("a", "b"))
})

test_that("apply_ref_level moves the requested level to the front", {
    f <- factor(c("a", "b", "c"))
    out <- apply_ref_level(f, "grp", list(grp = "c"))
    expect_identical(levels(out)[1], "c")
    expect_identical(as.character(out), c("a", "b", "c"))
})

test_that("apply_ref_level is a no-op when var has no entry in ref_levels", {
    f <- factor(c("a", "b"))
    out <- apply_ref_level(f, "grp", list(other = "b"))
    expect_identical(levels(out), levels(f))
})

test_that("apply_ref_level coerces the ref level to character before matching", {
    f <- factor(c("1", "0", "1"))
    out <- apply_ref_level(f, "grp", list(grp = 0))
    expect_identical(levels(out)[1], "0")
})

test_that("apply_ref_level errors when the requested level does not exist", {
    f <- factor(c("a", "b"))
    expect_error(
        apply_ref_level(f, "grp", list(grp = "z")),
        "not a level"
    )
})

# ------------------------------------------------------------------------
# has_offset_spec ----
# ------------------------------------------------------------------------

test_that("has_offset_spec is FALSE for NULL", {
    expect_false(has_offset_spec(NULL))
})

test_that("has_offset_spec is FALSE for an empty string", {
    expect_false(has_offset_spec(""))
})

test_that("has_offset_spec is TRUE for a column name", {
    expect_true(has_offset_spec("futime"))
})

# ------------------------------------------------------------------------
# resolve_offset ----
# ------------------------------------------------------------------------

test_that("resolve_offset returns 0 log-offset/unit exposure when absent", {
    d <- data.frame(x = 1:3)
    out <- resolve_offset(d, NULL, offset_log = FALSE, row_idx = 1:3)
    expect_equal(out$log_offset, rep(0, 3))
    expect_equal(out$exposure, rep(1, 3))
})

test_that("resolve_offset errors when required but absent", {
    d <- data.frame(x = 1:3)
    expect_error(
        resolve_offset(
            d,
            NULL,
            offset_log = FALSE,
            row_idx = 1:3,
            required = TRUE
        ),
        "required"
    )
    expect_error(
        resolve_offset(
            d,
            "",
            offset_log = FALSE,
            row_idx = 1:3,
            required = TRUE
        ),
        "required"
    )
})

test_that("resolve_offset converts raw exposure to log-offset correctly", {
    d <- data.frame(futime = c(1, exp(1), exp(2)))
    out <- resolve_offset(d, "futime", offset_log = FALSE, row_idx = 1:3)
    expect_equal(out$log_offset, c(0, 1, 2))
    expect_equal(out$exposure, c(1, exp(1), exp(2)))
})

test_that("resolve_offset passes through an already-logged offset unchanged", {
    d <- data.frame(lnfutime = c(0, 1, 2))
    out <- resolve_offset(d, "lnfutime", offset_log = TRUE, row_idx = 1:3)
    expect_equal(out$log_offset, c(0, 1, 2))
    expect_equal(out$exposure, exp(c(0, 1, 2)))
})

test_that("resolve_offset respects row_idx subsetting", {
    d <- data.frame(futime = c(1, 2, 4))
    out <- resolve_offset(d, "futime", offset_log = FALSE, row_idx = c(3L, 1L))
    expect_equal(out$exposure, c(4, 1))
})

test_that("resolve_offset errors for a non-scalar-character offset argument", {
    d <- data.frame(futime = 1:3)
    expect_error(
        resolve_offset(
            d,
            c("futime", "futime"),
            offset_log = FALSE,
            row_idx = 1:3
        ),
        "single variable name"
    )
})

test_that("resolve_offset errors when the named column is absent from data", {
    d <- data.frame(x = 1:3)
    expect_error(
        resolve_offset(d, "futime", offset_log = FALSE, row_idx = 1:3),
        "not found in data"
    )
})

test_that("resolve_offset errors on NA in the offset column", {
    d <- data.frame(futime = c(1, NA, 3))
    expect_error(
        resolve_offset(d, "futime", offset_log = FALSE, row_idx = 1:3),
        "NA"
    )
})

test_that("resolve_offset errors on non-positive raw exposure", {
    d <- data.frame(futime = c(1, 0, -3))
    expect_error(
        resolve_offset(d, "futime", offset_log = FALSE, row_idx = 1:3),
        "strictly positive"
    )
})

test_that("resolve_offset errors on non-finite log-offset/exposure", {
    d <- data.frame(lnfutime = c(0, Inf, 1))
    expect_error(
        resolve_offset(d, "lnfutime", offset_log = TRUE, row_idx = 1:3),
        "strictly positive and finite"
    )
})

test_that("resolve_offset errors when a finite log-offset overflows exp() to Inf", {
    # A distinct code path from the two tests above: off_raw itself is
    # finite (so the offset_log = TRUE finiteness check earlier in the
    # function passes), but exp(off_raw) overflows to +Inf for a large
    # enough value.
    d <- data.frame(lnfutime = c(0, log(.Machine$double.xmax) + 10, 1))
    expect_error(
        resolve_offset(d, "lnfutime", offset_log = TRUE, row_idx = 1:3),
        "strictly positive and finite"
    )
})

# ------------------------------------------------------------------------
# tilt_probs ----
# ------------------------------------------------------------------------
# Anchor: t_vals = 0:2, log_weight = log(c(1, 2, 1)) (i.e. C(t) = 1, 2, 1,
# the Binomial(2, 0.5)-shaped coefficient sequence used throughout this
# suite as `toy2`'s exact distribution). At beta = 0, probabilities are
# exactly (1, 2, 1) / 4. At beta = log(2), C(t) * exp(beta * t) becomes
# (1*1, 2*2, 1*4) = (1, 4, 4), so probabilities are (1, 4, 4) / 9.

test_that("tilt_probs at beta = 0 reproduces the un-tilted coefficient shape", {
    p <- tilt_probs(0, t_vals = 0:2, log_weight = log(c(1, 2, 1)))
    expect_equal(p, c(0.25, 0.5, 0.25))
    expect_equal(sum(p), 1)
})

test_that("tilt_probs applies the exponential tilt correctly at beta != 0", {
    p <- tilt_probs(log(2), t_vals = 0:2, log_weight = log(c(1, 2, 1)))
    expect_equal(p, c(1, 4, 4) / 9)
})

test_that("tilt_probs always sums to 1 regardless of beta magnitude", {
    for (b in c(-30, -5, 0, 5, 30)) {
        p <- tilt_probs(b, t_vals = 0:2, log_weight = log(c(1, 2, 1)))
        expect_equal(sum(p), 1, tolerance = 1e-10)
    }
})

test_that("tilt_probs is stable for extreme beta (no NaN from overflow)", {
    p <- tilt_probs(500, t_vals = 0:2, log_weight = log(c(1, 2, 1)))
    expect_false(anyNA(p))
    # All mass should concentrate on the largest attainable t.
    expect_equal(p[3], 1, tolerance = 1e-10)
})

# ------------------------------------------------------------------------
# solve_cmle ----
# ------------------------------------------------------------------------
# Same anchor table as above: t_vals = 0:2, log_weight = log(c(1, 2, 1)).
# At t_obs = 1 (the symmetric mean), the CMLE is exactly beta = 0, and the
# variance of T at beta = 0 is E[T^2] - E[T]^2 = 1.5 - 1 = 0.5, so
# se = 1 / sqrt(0.5) = sqrt(2).

test_that("solve_cmle recovers beta = 0 at the symmetric mean of the anchor table", {
    out <- solve_cmle(t_vals = 0:2, log_weight = log(c(1, 2, 1)), t_obs = 1)
    expect_equal(out$estimate, 0, tolerance = 1e-8)
    expect_equal(out$se, sqrt(2), tolerance = 1e-6)
})

test_that("solve_cmle returns -Inf/NA at the lower boundary of the support", {
    out <- solve_cmle(t_vals = 0:2, log_weight = log(c(1, 2, 1)), t_obs = 0)
    expect_identical(out$estimate, -Inf)
    expect_true(is.na(out$se))
})

test_that("solve_cmle returns +Inf/NA at the upper boundary of the support", {
    out <- solve_cmle(t_vals = 0:2, log_weight = log(c(1, 2, 1)), t_obs = 2)
    expect_identical(out$estimate, Inf)
    expect_true(is.na(out$se))
})

test_that("solve_cmle finds a nonzero root when t_obs is off-center", {
    # Skewed table: C(t) = (1, 1, 4) for t = 0, 1, 2; mean_t(0) = (0+1+8)/6 = 1.5.
    # Requesting t_obs = 1 (below the un-tilted mean) should require beta < 0.
    out <- solve_cmle(t_vals = 0:2, log_weight = log(c(1, 1, 4)), t_obs = 1)
    expect_lt(out$estimate, 0)
    expect_true(is.finite(out$se))
    # Sanity check via tilt_probs: the returned root should actually reproduce
    # a conditional mean equal to t_obs.
    p <- tilt_probs(out$estimate, 0:2, log(c(1, 1, 4)))
    expect_equal(sum((0:2) * p), 1, tolerance = 1e-6)
})

# ------------------------------------------------------------------------
# solve_tail_root ----
# ------------------------------------------------------------------------

test_that("solve_tail_root finds a simple linear root within the default bracket", {
    root <- solve_tail_root(function(b) 2 * b + 3, target = 7)
    expect_equal(root, 2, tolerance = 1e-8)
})

test_that("solve_tail_root finds the root of a logistic-shaped function", {
    stat_fun <- function(b) 1 / (1 + exp(-b))
    root <- solve_tail_root(stat_fun, target = 0.05)
    expect_equal(root, -log(19), tolerance = 1e-8)
})

test_that("solve_tail_root expands bracket when root lies outside default range", {
    # Root at b = 80, well outside the default (-50, 50) starting bracket;
    # this only succeeds if the widening loop actually runs.
    root <- solve_tail_root(function(b) b - 80, target = 0)
    expect_equal(root, 80, tolerance = 1e-6)
})

# ------------------------------------------------------------------------
# two_sided_from_tails ----
# ------------------------------------------------------------------------

test_that("two_sided_from_tails doubles the smaller tail and caps at 1", {
    out <- two_sided_from_tails(0.3, 0.9)
    expect_equal(out$pvalue, 0.6)
    expect_equal(out$pvalue_onesided, 0.3)
    expect_identical(out$direction, "upper")
})

test_that("two_sided_from_tails picks 'lower' when the upper tail is smaller", {
    out <- two_sided_from_tails(0.9, 0.3)
    expect_equal(out$pvalue_onesided, 0.3)
    expect_identical(out$direction, "lower")
})

test_that("two_sided_from_tails caps the two-sided p-value at 1", {
    out <- two_sided_from_tails(0.6, 0.7)
    expect_equal(out$pvalue, 1)
})

test_that("two_sided_from_tails resolves ties toward 'upper'", {
    out <- two_sided_from_tails(0.5, 0.5)
    expect_identical(out$direction, "upper")
    expect_equal(out$pvalue, 1)
})

# ------------------------------------------------------------------------
# exact_ci_from_tails ----
# ------------------------------------------------------------------------
# Anchor: tail_ge(b) = 1 / (1 + exp(-b)) (increasing logistic, 0 to 1) and
# tail_le(b) = 1 / (1 + exp(b)) (decreasing logistic, 1 to 0). At
# alpha = 0.10 (so alpha/2 = 0.05 for the two-sided case):
#   solve_tail_root(tail_ge, 0.05): 1/(1+exp(-b)) = 0.05 => b = -log(19)
#   solve_tail_root(tail_le, 0.05): 1/(1+exp(b))  = 0.05 => b =  log(19)

test_that("exact_ci_from_tails reproduces hand-derived 2sided logistic anchor", {
    tail_ge <- function(b) 1 / (1 + exp(-b))
    tail_le <- function(b) 1 / (1 + exp(b))
    ci <- exact_ci_from_tails(
        tail_le,
        tail_ge,
        alpha = 0.10,
        ci_type = "two-sided"
    )
    expect_equal(ci, c(-log(19), log(19)), tolerance = 1e-8)
})

test_that("exact_ci_from_tails one-sided 'upper' fixes the lower limit at -Inf", {
    tail_ge <- function(b) 1 / (1 + exp(-b))
    tail_le <- function(b) 1 / (1 + exp(b))
    ci <- exact_ci_from_tails(
        tail_le,
        tail_ge,
        alpha = 0.10,
        ci_type = "one-sided",
        direction = "upper"
    )
    expect_identical(ci[1], -Inf)
    # solve_tail_root(tail_le, alpha = 0.10): 1/(1+exp(b)) = 0.10 => b = log(9)
    expect_equal(ci[2], log(9), tolerance = 1e-8)
})

test_that("exact_ci_from_tails one-sided 'lower' fixes the upper limit at Inf", {
    tail_ge <- function(b) 1 / (1 + exp(-b))
    tail_le <- function(b) 1 / (1 + exp(b))
    ci <- exact_ci_from_tails(
        tail_le,
        tail_ge,
        alpha = 0.10,
        ci_type = "one-sided",
        direction = "lower"
    )
    expect_identical(ci[2], Inf)
    # solve_tail_root(tail_ge, alpha = 0.10):
    # 1/(1+exp(-b)) = 0.10 => b = -log(9)
    expect_equal(ci[1], -log(9), tolerance = 1e-8)
})

test_that("exact_ci_from_tails honours at_lower_boundary/at_upper_boundary flags", {
    tail_ge <- function(b) 1 / (1 + exp(-b))
    tail_le <- function(b) 1 / (1 + exp(b))
    ci <- exact_ci_from_tails(
        tail_le,
        tail_ge,
        alpha = 0.10,
        ci_type = "two-sided",
        at_lower_boundary = TRUE,
        at_upper_boundary = TRUE
    )
    expect_identical(ci, c(-Inf, Inf))
})

test_that("exact_ci_from_tails defaults to two-sided/upper when args are omitted", {
    tail_ge <- function(b) 1 / (1 + exp(-b))
    tail_le <- function(b) 1 / (1 + exp(b))
    ci <- exact_ci_from_tails(tail_le, tail_ge, alpha = 0.10)
    expect_equal(ci, c(-log(19), log(19)), tolerance = 1e-8)
})

test_that("exact_ci_from_tails rejects an invalid ci_type/direction", {
    tail_ge <- function(b) 1 / (1 + exp(-b))
    tail_le <- function(b) 1 / (1 + exp(b))
    expect_error(
        exact_ci_from_tails(tail_le, tail_ge, alpha = 0.10, ci_type = "bogus")
    )
    expect_error(
        exact_ci_from_tails(
            tail_le,
            tail_ge,
            alpha = 0.10,
            ci_type = "one-sided",
            direction = "bogus"
        )
    )
})

# print.ep_fit()'s tests live in tests/testthat/test-ep-fit.R, alongside
# new_ep_fit()/validate_ep_fit()/ep_fit() -- one test file per R file
# (see R/ep-fit.R), matching the reviewer's naming-convention comment.
