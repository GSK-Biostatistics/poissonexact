# ==============================================================================
# Tests for ep-common.R, shared infrastructure specific to the
# general (arbitrary categorical covariates + optional strata) family of
# algorithms: exactmax/offset cap-rule resolution, SAS-style factor/model-
# matrix preparation, exact_vars-to-column resolution, pruning bounds,
# global filtering/marginalization, the shared CI wrapper, and the
# .setup_exact_model()/.assemble_parameter_output() orchestration helpers
# used by both ep_multi_shift() and ep_network(). new_ep_fit()'s own
# tests live in tests/testthat/test-ep-fit.R, next to R/ep-fit.R.
# ==============================================================================

# ------------------------------------------------------------------------
# parse_exactmax_spec ----
# ------------------------------------------------------------------------

test_that("parse_exactmax_spec resolves NULL to mode 'none'", {
    out <- parse_exactmax_spec(NULL)
    expect_identical(out$mode, "none")
    expect_identical(out$var, character(0))
})

test_that("parse_exactmax_spec resolves TRUE to mode 'novar'", {
    out <- parse_exactmax_spec(TRUE)
    expect_identical(out$mode, "novar")
    expect_identical(out$var, character(0))
})

test_that("parse_exactmax_spec resolves a variable name to mode 'var'", {
    out <- parse_exactmax_spec("subjects")
    expect_identical(out$mode, "var")
    expect_identical(out$var, "subjects")
})

test_that("parse_exactmax_spec rejects invalid specifications", {
    expect_error(parse_exactmax_spec(5), "exactmax must be")
    expect_error(parse_exactmax_spec(c("a", "b")), "exactmax must be")
    expect_error(parse_exactmax_spec(""), "exactmax must be")
    expect_error(parse_exactmax_spec(FALSE), "exactmax must be")
})

# ------------------------------------------------------------------------
# si_from_stratum_total / si_from_exactmax_var / si_from_offset ----
# ------------------------------------------------------------------------

test_that("si_from_stratum_total caps every observation at the stratum total", {
    expect_identical(si_from_stratum_total(c(3L, 4L, 5L)), rep(12L, 3))
})

test_that("si_from_stratum_total handles a single-row stratum", {
    expect_identical(si_from_stratum_total(7L), 7L)
})

test_that("si_from_exactmax_var floors and returns the supplied values", {
    expect_identical(si_from_exactmax_var(c(2.7, 3.2, 0)), c(2L, 3L, 0L))
})

test_that("si_from_exactmax_var errors on NA/NaN/Inf", {
    expect_error(si_from_exactmax_var(c(1, NA, 3)), "NA/NaN/Inf")
    expect_error(si_from_exactmax_var(c(1, NaN, 3)), "NA/NaN/Inf")
    expect_error(si_from_exactmax_var(c(1, Inf, 3)), "NA/NaN/Inf")
})

test_that("si_from_exactmax_var errors on negative values", {
    expect_error(si_from_exactmax_var(c(1, -2, 3)), "negative values")
})

test_that("si_from_offset applies ceil-and-floor-at-1 to exposures", {
    log_offset <- c(0, log(2), log(0.5))
    expect_identical(si_from_offset(log_offset), c(1L, 2L, 1L))
})

test_that("si_from_offset floors nothing below 1 even for tiny exposures", {
    expect_identical(si_from_offset(log(0.001)), 1L)
})

# ------------------------------------------------------------------------
# prepare_sas_factors ----
# ------------------------------------------------------------------------

test_that("prepare_sas_factors relevels and sets treatment contrasts", {
    d <- data.frame(grp = c("b", "a", "b"), x = 1:3)
    out <- prepare_sas_factors(
        d,
        class_vars = "grp",
        ref_levels = list(grp = "b")
    )
    expect_s3_class(out$grp, "factor")
    expect_identical(levels(out$grp)[1], "b")
    contr_mat <- attr(out$grp, "contrasts")
    expect_identical(colnames(contr_mat), setdiff(levels(out$grp), "b"))
})

test_that("prepare_sas_factors ignores class_vars not present in data", {
    d <- data.frame(x = 1:3)
    out <- prepare_sas_factors(
        d,
        class_vars = "nonexistent",
        ref_levels = list()
    )
    expect_named(out, names(d))
})

test_that("prepare_sas_factors coerces a character/numeric column to factor", {
    d <- data.frame(grp = c(1, 0, 1))
    out <- prepare_sas_factors(d, class_vars = "grp")
    expect_s3_class(out$grp, "factor")
})

# ------------------------------------------------------------------------
# build_mm_no_intercept ----
# ------------------------------------------------------------------------

test_that("build_mm_no_intercept removes the intercept and preserves assign", {
    d <- data.frame(y = 1:3, g = factor(c("a", "b", "a")))
    built <- build_mm_no_intercept(y ~ g, d)
    expect_identical(colnames(built$x), "gb")
    expect_equal(as.numeric(built$x[, "gb"]), c(0, 1, 0))
    expect_identical(attr(built$x, "assign"), 1L)
    expect_equal(nrow(built$mf), 3L)
})

# uses a model with multiple terms
test_that("build_mm_no_intercept assigns each column to its own term", {
    d <- data.frame(
        y = 1:6,
        b = factor(c(0, 1, 0, 1, 0, 1)),
        a = factor(c("x", "x", "y", "y", "z", "z"))
    )
    built <- build_mm_no_intercept(y ~ b + a, d)
    expect_identical(colnames(built$x), c("b1", "ay", "az"))
    expect_identical(attr(built$x, "assign"), c(1L, 2L, 2L))
    expect_identical(attr(stats::terms(built$mf), "term.labels"), c("b", "a"))
})

# ------------------------------------------------------------------------
# .permutations ----
# ------------------------------------------------------------------------

test_that(".permutations enumerates all orderings of a short vector", {
    out <- .permutations(c("a", "b"))
    expect_length(out, 2)
    expect_true(any(vapply(out, identical, logical(1), c("a", "b"))))
    expect_true(any(vapply(out, identical, logical(1), c("b", "a"))))
})

test_that(".permutations returns the vector unchanged for length <= 1", {
    expect_identical(.permutations(character(0)), list(character(0)))
    expect_identical(.permutations("a"), list("a"))
})

test_that(".permutations enumerates 3! = 6 orderings for a 3-element vector", {
    out <- .permutations(c("a", "b", "c"))
    expect_length(out, 6)
    expect_length(unique(lapply(out, sort)), 1) # All has the same set
})

test_that(".permutations refuses to expand vectors longer than 6", {
    x <- letters[1:7]
    out <- .permutations(x)
    expect_identical(out, list(x))
})

# ------------------------------------------------------------------------
# match_term_label ----
# ------------------------------------------------------------------------

test_that("match_term_label matches a plain main-effect term", {
    expect_identical(match_term_label("a", c("a", "b", "a:b")), 1L)
})

test_that("match_term_label matches an interaction regardless of order", {
    expect_identical(match_term_label("b:a", c("a", "b", "a:b")), 3L)
    expect_identical(match_term_label("a*b", c("a", "b", "a:b")), 3L)
})

test_that("match_term_label matches a factor()-wrapped variable", {
    expect_identical(match_term_label("g", c("factor(g)", "x")), 1L)
})

test_that("match_term_label returns NA for an unmatched term", {
    expect_true(is.na(match_term_label("zzz", c("a", "b", "a:b"))))
})

# ------------------------------------------------------------------------
# exact_term_columns ----
# ------------------------------------------------------------------------

test_that("exact_term_columns resolves a main effect and an interaction", {
    d <- data.frame(
        y = 1:6,
        b = factor(c(0, 1, 0, 1, 0, 1)),
        a = factor(c("x", "x", "y", "y", "z", "z"))
    )
    built <- build_mm_no_intercept(y ~ b + a, d)

    expect_identical(
        exact_term_columns(built$x, built$terms, "a"),
        c(2L, 3L)
    )
    expect_identical(
        exact_term_columns(built$x, built$terms, "b"),
        1L
    )
    expect_identical(
        exact_term_columns(built$x, built$terms, c("a", "b")),
        c(1L, 2L, 3L)
    )
})

test_that("exact_term_columns errors for a term absent from the model", {
    d <- data.frame(y = 1:6, b = factor(c(0, 1, 0, 1, 0, 1)))
    built <- build_mm_no_intercept(y ~ b, d)
    expect_error(
        exact_term_columns(built$x, built$terms, "nonexistent"),
        "does not match any term"
    )
})

test_that("exact_term_columns errors when model matrix has no column names", {
    mm <- matrix(1, nrow = 2, ncol = 1)
    attr(mm, "assign") <- 1L
    terms_obj <- structure(list(), term.labels = "x")
    expect_error(
        exact_term_columns(mm, terms_obj, "x"),
        "no column names"
    )
})

test_that("exact_term_columns errors when assign/term.labels attrs missing", {
    mm <- matrix(1, nrow = 2, ncol = 1, dimnames = list(NULL, "x"))
    terms_obj <- structure(list(), term.labels = "x")
    expect_error(
        exact_term_columns(mm, terms_obj, "x"),
        "missing its"
    )
})

# ------------------------------------------------------------------------
# Interaction terms, resolved from a genuine `a * b` formula ----
# (not fabricated term_labels): d has a (2 levels: x, y; ref x) and
# b (2 levels: 0, 1; ref 0), so model.matrix(~a*b) yields columns
# "ay", "b1", "ay:b1" with assign = c(1, 2, 3) and
# term.labels = c("a", "b", "a:b").
# ------------------------------------------------------------------------

test_that("exact_term_columns resolves interaction to its own column", {
    d <- data.frame(
        y = 1:4,
        a = factor(c("x", "y", "x", "y")),
        b = factor(c(0, 0, 1, 1))
    )
    built <- build_mm_no_intercept(y ~ a * b, d)
    expect_identical(colnames(built$x), c("ay", "b1", "ay:b1"))
    expect_identical(attr(built$x, "assign"), c(1L, 2L, 3L))
    expect_identical(
        attr(stats::terms(built$mf), "term.labels"),
        c("a", "b", "a:b")
    )

    expect_identical(exact_term_columns(built$x, built$terms, "a:b"), 3L)
})

test_that("exact_term_columns excludes interaction when main effect only", {
    d <- data.frame(
        y = 1:4,
        a = factor(c("x", "y", "x", "y")),
        b = factor(c(0, 0, 1, 1))
    )
    built <- build_mm_no_intercept(y ~ a * b, d)

    # Requesting "a" alone must resolve only to column 1 ("ay"), never to
    # column 3 (the interaction), matching ep_multi_shift()'s
    # documented contract that a main effect never implicitly pulls in its
    # interactions.
    expect_identical(exact_term_columns(built$x, built$terms, "a"), 1L)
    expect_identical(exact_term_columns(built$x, built$terms, "b"), 2L)

    # Requesting all three explicitly covers every column.
    expect_identical(
        exact_term_columns(built$x, built$terms, c("a", "b", "a:b")),
        c(1L, 2L, 3L)
    )
})

test_that("exact_term_columns resolves interaction via '*' syntax, any order", {
    d <- data.frame(
        y = 1:4,
        a = factor(c("x", "y", "x", "y")),
        b = factor(c(0, 0, 1, 1))
    )
    built <- build_mm_no_intercept(y ~ a * b, d)

    # "b*a" (reversed order, '*' instead of ':') must still resolve to the
    # same interaction column as "a:b", since R's own term.labels always
    # records it as "a:b" (first-appearance order) regardless of how the
    # user spells the request.
    expect_identical(exact_term_columns(built$x, built$terms, "b*a"), 3L)
    expect_identical(exact_term_columns(built$x, built$terms, "b:a"), 3L)
})

test_that("exact_term_columns resolves a three-way interaction in any order", {
    d <- data.frame(
        y = 1:8,
        a = factor(rep(c("x", "y"), 4)),
        b = factor(rep(c(0, 1), each = 2, times = 2)),
        cc = factor(rep(c(0, 1), each = 4))
    )
    built <- build_mm_no_intercept(y ~ a * b * cc, d)
    term_labels <- attr(stats::terms(built$mf), "term.labels")
    expect_true("a:b:cc" %in% term_labels)

    three_way_col <- which(
        attr(built$x, "assign") == match("a:b:cc", term_labels)
    )
    expect_identical(
        exact_term_columns(built$x, built$terms, "cc:a:b"),
        three_way_col
    )
    expect_identical(
        exact_term_columns(built$x, built$terms, "b*cc*a"),
        three_way_col
    )
})

# ------------------------------------------------------------------------
# compute_stratum_max_add ----
# ------------------------------------------------------------------------

test_that("compute_stratum_max_add matches the hand-derived example", {
    y_s <- c(3, 4)
    xc_s <- matrix(c(1, 0), ncol = 1)
    si_s <- c(2, 5)
    expect_identical(compute_stratum_max_add(y_s, xc_s, si_s), 2L)
})

test_that("compute_stratum_max_add returns 0 for a column that is all-zero", {
    y_s <- c(3, 4)
    xc_s <- matrix(c(0, 0), ncol = 1)
    si_s <- c(2, 5)
    expect_identical(compute_stratum_max_add(y_s, xc_s, si_s), 0L)
})

test_that("compute_stratum_max_add caps the bound at the stratum total", {
    # Both rows contribute to the column, so the naive sum(si_s) = 20 would
    # exceed the stratum's own observed total (ms = 7); the bound must be
    # capped at 7.
    y_s <- c(3, 4)
    xc_s <- matrix(c(1, 1), ncol = 1)
    si_s <- c(10, 10)
    expect_identical(compute_stratum_max_add(y_s, xc_s, si_s), 7L)
})

# nuisance columns are handled independently
test_that("compute_stratum_max_add handles multiple nuisance columns", {
    y_s <- c(3, 4)
    xc_s <- matrix(c(1, 0, 0, 1), ncol = 2)
    si_s <- c(2, 5)
    expect_identical(compute_stratum_max_add(y_s, xc_s, si_s), c(2L, 5L))
})

# ------------------------------------------------------------------------
# filter_and_marginalize_pn_lw ----
# ------------------------------------------------------------------------

test_that("filter_and_marginalize_pn_lw marginalizes with no nuisance cols", {
    joint <- list(t_j = c(1L, 1L, 2L), lw = log(c(2, 3, 5)))
    out <- filter_and_marginalize_pn_lw(joint, tn_obs = numeric(0))
    expect_identical(out$t_j, c(1L, 2L))
    expect_equal(out$lw, log(c(5, 5)))
})

test_that("filter_and_marginalize_pn_lw filters nuisance total, marginalizes", {
    joint <- list(
        t_j = c(0L, 0L, 1L, 1L),
        cond_1 = c(5, 5, 5, 7),
        lw = log(c(1, 2, 3, 4))
    )
    out <- filter_and_marginalize_pn_lw(joint, tn_obs = 5)
    expect_identical(out$t_j, c(0L, 1L))
    expect_equal(out$lw, log(c(3, 3)))
})

test_that("filter_and_marginalize_pn_lw errors when no rows match total", {
    joint <- list(
        t_j = c(0L, 1L),
        cond_1 = c(5, 5),
        lw = log(c(1, 2))
    )
    expect_error(
        filter_and_marginalize_pn_lw(joint, tn_obs = 999),
        "No states remain"
    )
})

# columns are filtered simultaneously
test_that("filter_and_marginalize_pn_lw filters on multiple nuisance columns", {
    joint <- list(
        t_j = c(0L, 0L, 1L),
        cond_1 = c(5, 5, 5),
        cond_2 = c(1, 2, 1),
        lw = log(c(1, 100, 3))
    )
    out <- filter_and_marginalize_pn_lw(joint, tn_obs = c(5, 1))
    expect_identical(out$t_j, c(0L, 1L))
    expect_equal(out$lw, log(c(1, 3)))
})

# ------------------------------------------------------------------------
# exact_ci_from_dist1d ----
# ------------------------------------------------------------------------
# Anchor table: t_j = 0:2, lw = log(c(1, 2, 1)), the same Binomial(2,0.5)-
# shaped coefficients used throughout this suite. Since the coefficient
# sequence is palindromic (c_0 == c_2), the two-sided CI at any t_obs = 1
# must be symmetric about 0 by the theta -> 1/theta relabeling argument
# documented in helper-data.R.

# anchor is the palindromic point
test_that("exact_ci_from_dist1d produces symmetric interval at anchor", {
    dist1d <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    ci <- exact_ci_from_dist1d(
        dist1d,
        t_obs = 1,
        alpha = 0.10,
        ci_type = "two-sided"
    )
    expect_equal(ci[1], -ci[2], tolerance = 1e-8)

    # Cross-check the lower limit actually solves tail_ge(beta) = alpha/2 by
    # recomputing tail_ge at the returned root via tilt_probs() directly.
    p_at_lo <- tilt_probs(ci[1], dist1d$t_j, dist1d$lw)
    expect_equal(sum(p_at_lo[dist1d$t_j >= 1]), 0.05, tolerance = 1e-8)

    p_at_hi <- tilt_probs(ci[2], dist1d$t_j, dist1d$lw)
    expect_equal(sum(p_at_hi[dist1d$t_j <= 1]), 0.05, tolerance = 1e-8)
})

test_that("exact_ci_from_dist1d fixes lower limit at -Inf at minimum t_obs", {
    dist1d <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    ci <- exact_ci_from_dist1d(
        dist1d,
        t_obs = 0,
        alpha = 0.10,
        ci_type = "two-sided"
    )
    expect_identical(ci[1], -Inf)
    expect_true(is.finite(ci[2]))
})

test_that("exact_ci_from_dist1d fixes upper limit at Inf at maximum t_obs", {
    dist1d <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    ci <- exact_ci_from_dist1d(
        dist1d,
        t_obs = 2,
        alpha = 0.10,
        ci_type = "two-sided"
    )
    expect_identical(ci[2], Inf)
    expect_true(is.finite(ci[1]))
})

test_that("exact_ci_from_dist1d one-sided 'upper' reports only upper limit", {
    dist1d <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    ci <- exact_ci_from_dist1d(
        dist1d,
        t_obs = 1,
        alpha = 0.10,
        ci_type = "one-sided",
        direction = "upper"
    )
    expect_identical(ci[1], -Inf)
    expect_true(is.finite(ci[2]))
})

# ------------------------------------------------------------------------
# .cumulative_from_end ----
# ------------------------------------------------------------------------

test_that(".cumulative_from_end matches the hand-derived 3-stratum example", {
    vec_list <- list(c(1, 2), c(3, 4), c(5, 6))
    out <- .cumulative_from_end(vec_list, p_cond = 2)
    expect_equal(out, list(c(8, 10), c(5, 6), c(0, 0)))
})

test_that(".cumulative_from_end handles a single stratum (nothing after it)", {
    out <- .cumulative_from_end(list(c(1, 2)), p_cond = 2)
    expect_equal(out[[1]], c(0, 0))
})

# vectors are integer-typed
test_that(".cumulative_from_end returns empty vectors when p_cond == 0", {
    out <- .cumulative_from_end(list(1, 2, 3), p_cond = 0)
    expect_length(out, 3)
    expect_identical(out[[1]], integer(0))
    expect_identical(out[[3]], integer(0))
})

# ------------------------------------------------------------------------
# .setup_exact_model ----
# ------------------------------------------------------------------------

# verified against a hand-derived fixture
test_that(".setup_exact_model applies SAS rule 3 (offset, no exactmax)", {
    out <- .setup_exact_model(
        formula = cases ~ group,
        data = test1,
        offset = "lnfutime",
        offset_log = TRUE,
        exactmax = NULL,
        strata = "age",
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    expect_true(out$use_rule3)
    expect_false(out$use_rule1)
    expect_false(out$use_rule2)
    expect_identical(out$y, test1$cases)
    expect_equal(out$log_off, test1$lnfutime)
    expect_identical(out$sid, as.character(test1$age))
    expect_identical(out$slevels0, c("50", "60", "70"))
    expect_identical(out$cn, "group1")
    expect_identical(out$exact_cols, 1L)
})

test_that(".setup_exact_model applies SAS rule 2 with no offset or exactmax", {
    out <- .setup_exact_model(
        formula = cases ~ group,
        data = test1,
        offset = NULL,
        offset_log = FALSE,
        exactmax = NULL,
        strata = character(0),
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    expect_true(out$use_rule2)
    expect_false(out$use_rule1)
    expect_false(out$use_rule3)
    expect_identical(out$sid, rep("ALL", nrow(test1)))
    expect_identical(out$slevels0, "ALL")
    expect_equal(out$log_off, rep(0, nrow(test1)))
})

test_that(".setup_exact_model applies SAS rule 1 when EXACTMAX=var is given", {
    out <- .setup_exact_model(
        formula = cases ~ group,
        data = test1,
        offset = NULL,
        offset_log = FALSE,
        exactmax = "subjects",
        strata = character(0),
        exact_vars = "group",
        ref_levels = list(),
        verbose = FALSE
    )
    expect_true(out$use_rule1)
    expect_false(out$use_rule2)
    expect_false(out$use_rule3)
    expect_identical(out$exactmax_vec, as.integer(test1$subjects))
})

test_that(".setup_exact_model errors when EXACTMAX variable is absent", {
    expect_error(
        .setup_exact_model(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = "nonexistent_var",
            strata = character(0),
            exact_vars = "group",
            ref_levels = list(),
            verbose = FALSE
        ),
        "not found in data"
    )
})

test_that(".setup_exact_model verbose = TRUE prints the applicable SAS rule", {
    expect_message(
        .setup_exact_model(
            formula = cases ~ group,
            data = test1,
            offset = "lnfutime",
            offset_log = TRUE,
            exactmax = NULL,
            strata = character(0),
            exact_vars = "group",
            ref_levels = list(),
            verbose = TRUE
        ),
        "SAS rule \\(3\\)"
    )
    expect_message(
        .setup_exact_model(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = NULL,
            strata = character(0),
            exact_vars = "group",
            ref_levels = list(),
            verbose = TRUE
        ),
        "SAS rule \\(2\\)"
    )
    expect_message(
        .setup_exact_model(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = "subjects",
            strata = character(0),
            exact_vars = "group",
            ref_levels = list(),
            verbose = TRUE
        ),
        "SAS rule \\(1\\)"
    )
})

test_that(".setup_exact_model errors when exact_term_columns() finds none", {
    # This branch (`if (length(exact_cols) == 0L) stop("No exact_vars
    # columns found in model matrix.")`) is defensive, not reachable through
    # the public API with well-formed input.
    testthat::local_mocked_bindings(
        exact_term_columns = function(...) integer(0)
    )
    expect_error(
        .setup_exact_model(
            formula = cases ~ group,
            data = test1,
            offset = "lnfutime",
            offset_log = TRUE,
            exactmax = NULL,
            strata = character(0),
            exact_vars = "group",
            ref_levels = list(),
            verbose = FALSE
        ),
        "No exact_vars columns found"
    )
})

# ------------------------------------------------------------------------
# .assemble_parameter_output ----
# ------------------------------------------------------------------------
# Anchor: dist1d = (t_j = 0:2, lw = log(c(1, 2, 1))); y/x chosen so that
# t_obs = sum(y * x[, jcol]) == 1, landing exactly on the symmetric mean
# (same toy2-style derivation as elsewhere in this suite): Estimate ~ 0,
# PValue == 1 (two-sided), SuffStat == 1.

# hand-derived, two-sided case
test_that(".assemble_parameter_output reproduces symmetric anchor", {
    dist1d <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    y <- c(1, 0, 0, 1)
    x <- matrix(c(1, 0, 1, 0), ncol = 1)

    out <- .assemble_parameter_output(
        dist1d,
        parm = "group1",
        y = y,
        x = x,
        jcol = 1,
        alpha = 0.10,
        ci_type = "two-sided"
    )

    expect_equal(
        as.list(out$result_row[c("Estimate", "PValue")]),
        list(Estimate = 0, PValue = 1),
        tolerance = 1e-8
    )
    expect_identical(out$result_row$SuffStat, 1L)
    expect_equal(out$result_row$Parameter, "group1")
    expect_equal(
        as.list(out$outdist_entry[c("prob", "SuffStat")]),
        list(prob = c(0.25, 0.5, 0.25), SuffStat = 0:2),
        tolerance = 1e-10
    )
    # Symmetric interval, same argument as exact_ci_from_dist1d's own test.
    expect_equal(
        out$result_row$LowerCL,
        -out$result_row$UpperCL,
        tolerance = 1e-8
    )
})

test_that(".assemble_parameter_output reports one-sided p-value un-doubled", {
    dist1d <- list(t_j = 0:2, lw = log(c(1, 2, 1)))
    y <- c(1, 0, 0, 1)
    x <- matrix(c(1, 0, 1, 0), ncol = 1)

    out <- .assemble_parameter_output(
        dist1d,
        parm = "group1",
        y = y,
        x = x,
        jcol = 1,
        alpha = 0.10,
        ci_type = "one-sided"
    )
    expect_equal(out$result_row$PValue, 0.75, tolerance = 1e-8)
    expect_identical(out$result_row$LowerCL, -Inf)
})

test_that(".assemble_parameter_output computes SuffStat from asymmetric obs", {
    # Skewed table C(t) = (1, 1, 4); an observation landing at t_obs = 2
    # (the boundary) should report Estimate = Inf.
    dist1d <- list(t_j = 0:2, lw = log(c(1, 1, 4)))
    y <- c(1, 1)
    x <- matrix(c(1, 1), ncol = 1)
    out <- .assemble_parameter_output(
        dist1d,
        parm = "x1",
        y = y,
        x = x,
        jcol = 1,
        alpha = 0.10,
        ci_type = "two-sided"
    )
    expect_identical(out$result_row$SuffStat, 2L)
    expect_identical(out$result_row$Estimate, Inf)
})

# Tests for new_ep_fit()/validate_ep_fit()/ep_fit() (R/ep-fit.R) live in
# tests/testthat/test-ep-fit.R, alongside print.ep_fit()'s tests -- one
# test file per R file, matching the reviewer's naming-convention comment
# on ep-fit.R.
