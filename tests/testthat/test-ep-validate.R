# ==============================================================================
# Tests for ep-validate.R; shared argument validation used by
# every exported main function. Covers each single-purpose `.check_*_arg()`
# checker individually, then the two per-function validators that compose
# them (`.validate_exact_poisson_args()`, `.validate_closedform_args()`).
# ==============================================================================

# ------------------------------------------------------------------------
# .check_scalar_flag ----
# ------------------------------------------------------------------------

test_that(".check_scalar_flag accepts TRUE and FALSE", {
    expect_true(.check_scalar_flag(TRUE, "verbose"))
    expect_true(.check_scalar_flag(FALSE, "verbose"))
})

test_that(".check_scalar_flag rejects NA, non-logical, and non-scalar values", {
    expect_error(
        .check_scalar_flag(NA, "verbose"),
        "verbose.*must be a single"
    )
    expect_error(
        .check_scalar_flag(1, "verbose"),
        "must be a single"
    )
    expect_error(
        .check_scalar_flag("TRUE", "verbose"),
        "must be a single"
    )
    expect_error(
        .check_scalar_flag(c(TRUE, FALSE), "verbose"),
        "must be a single"
    )
})

test_that(".check_scalar_flag includes the arg name in the error message", {
    expect_error(
        .check_scalar_flag(NA, "offset_log"),
        "`offset_log`"
    )
})

# ------------------------------------------------------------------------
# .check_formula_arg ----
# ------------------------------------------------------------------------

test_that(".check_formula_arg accepts a genuine formula", {
    expect_true(.check_formula_arg(cases ~ group))
})

test_that(".check_formula_arg rejects a non-formula object", {
    expect_error(
        .check_formula_arg("cases ~ group"),
        "must be a model formula"
    )
    expect_error(
        .check_formula_arg(list(1, 2)),
        "must be a model formula"
    )
})

# ------------------------------------------------------------------------
# .check_data_arg ----
# ------------------------------------------------------------------------

test_that(".check_data_arg accepts a non-empty data frame", {
    expect_true(.check_data_arg(data.frame(x = 1)))
})

test_that(".check_data_arg rejects a non-data-frame", {
    expect_error(
        .check_data_arg(matrix(1)),
        "must be a data frame"
    )
    expect_error(
        .check_data_arg(list(x = 1)),
        "must be a data frame"
    )
})

test_that(".check_data_arg rejects a zero-row data frame", {
    expect_error(
        .check_data_arg(data.frame(x = numeric(0))),
        "zero rows"
    )
})

# ------------------------------------------------------------------------
# .check_offset_name_arg ----
# ------------------------------------------------------------------------

test_that(".check_offset_name_arg accepts NULL/'' when not required", {
    d <- data.frame(x = 1)
    expect_true(.check_offset_name_arg(NULL, d, required = FALSE))
    expect_true(.check_offset_name_arg("", d, required = FALSE))
})

test_that(".check_offset_name_arg errors when required but absent", {
    d <- data.frame(x = 1)
    expect_error(
        .check_offset_name_arg(NULL, d, required = TRUE),
        "is required"
    )
    expect_error(
        .check_offset_name_arg("", d, required = TRUE),
        "is required"
    )
})

test_that(".check_offset_name_arg accepts a column that exists in data", {
    d <- data.frame(futime = 1)
    expect_true(.check_offset_name_arg("futime", d, required = FALSE))
})

test_that(".check_offset_name_arg errors when the named column is absent", {
    d <- data.frame(x = 1)
    expect_error(
        .check_offset_name_arg("futime", d, required = FALSE),
        "not found in `data`"
    )
})

test_that(".check_offset_name_arg errors for a malformed offset value", {
    d <- data.frame(futime = 1)
    expect_error(
        .check_offset_name_arg(NA_character_, d, FALSE),
        "must be a single non-empty"
    )
    expect_error(
        .check_offset_name_arg(c("a", "b"), d, FALSE),
        "must be a single non-empty"
    )
    expect_error(
        .check_offset_name_arg(123, d, FALSE),
        "must be a single non-empty"
    )
})

# ------------------------------------------------------------------------
# .check_strata_arg ----
# ------------------------------------------------------------------------

test_that(".check_strata_arg accepts NULL", {
    expect_true(.check_strata_arg(NULL, data.frame(x = 1)))
})

test_that(".check_strata_arg accepts existing column names", {
    d <- data.frame(age = 1, group = 1)
    expect_true(.check_strata_arg(c("age", "group"), d))
})

test_that(".check_strata_arg rejects non-character/NA/empty-string entries", {
    d <- data.frame(age = 1)
    expect_error(
        .check_strata_arg(1, d),
        "must be `NULL` or a character"
    )
    expect_error(
        .check_strata_arg(NA_character_, d),
        "must be `NULL` or a character"
    )
    expect_error(
        .check_strata_arg(c("age", ""), d),
        "must be `NULL` or a character"
    )
})

test_that(".check_strata_arg errors listing every missing column", {
    d <- data.frame(age = 1)
    expect_error(
        .check_strata_arg(c("age", "region", "site"), d),
        "region, site"
    )
})

# ------------------------------------------------------------------------
# .check_exact_vars_arg ----
# ------------------------------------------------------------------------

test_that(".check_exact_vars_arg accepts a non-empty character vector", {
    expect_true(.check_exact_vars_arg("group"))
    expect_true(.check_exact_vars_arg(c("group", "region")))
})

test_that(".check_exact_vars_arg rejects empty/NA/non-character/blank", {
    expect_error(
        .check_exact_vars_arg(character(0)),
        "non-empty character vector"
    )
    expect_error(
        .check_exact_vars_arg(NULL),
        "non-empty character vector"
    )
    expect_error(
        .check_exact_vars_arg(NA_character_),
        "non-empty character vector"
    )
    expect_error(
        .check_exact_vars_arg(1),
        "non-empty character vector"
    )
    expect_error(
        .check_exact_vars_arg(c("group", "")),
        "non-empty character vector"
    )
})

# ------------------------------------------------------------------------
# .check_ref_levels_arg ----
# ------------------------------------------------------------------------

test_that(".check_ref_levels_arg accepts an empty list", {
    expect_true(.check_ref_levels_arg(list(), data.frame(x = 1)))
})

test_that(".check_ref_levels_arg accepts valid named list matching data cols", {
    d <- data.frame(group = 1, region = 1)
    expect_true(
        .check_ref_levels_arg(
            list(group = "0", region = "1"),
            d
        )
    )
})

test_that(".check_ref_levels_arg rejects a non-list", {
    expect_error(
        .check_ref_levels_arg("not a list", data.frame(x = 1)),
        "possibly empty. named list"
    )
})

test_that(".check_ref_levels_arg rejects an unnamed or partially-named list", {
    d <- data.frame(group = 1)
    expect_error(.check_ref_levels_arg(list("0"), d), "must be named")
    expect_error(
        .check_ref_levels_arg(list(group = "0", "1"), d),
        "must be named"
    )
})

test_that(".check_ref_levels_arg rejects a multi-value entry", {
    d <- data.frame(group = 1)
    expect_error(
        .check_ref_levels_arg(list(group = c("0", "1")), d),
        "must be a single value"
    )
})

test_that(".check_ref_levels_arg errors listing variables missing from data", {
    d <- data.frame(group = 1)
    expect_error(
        .check_ref_levels_arg(list(group = "0", region = "1"), d),
        "region"
    )
})

# ------------------------------------------------------------------------
# .check_alpha_arg ----
# ------------------------------------------------------------------------

test_that(".check_alpha_arg accepts values strictly between 0 and 1", {
    expect_true(.check_alpha_arg(0.05))
    expect_true(.check_alpha_arg(0.999))
    expect_true(.check_alpha_arg(0.001))
})

test_that(".check_alpha_arg rejects boundary and out-of-range values", {
    expect_error(
        .check_alpha_arg(0),
        "strictly between 0 and 1"
    )
    expect_error(
        .check_alpha_arg(1),
        "strictly between 0 and 1"
    )
    expect_error(
        .check_alpha_arg(-0.1),
        "strictly between 0 and 1"
    )
    expect_error(
        .check_alpha_arg(1.5),
        "strictly between 0 and 1"
    )
})

test_that(".check_alpha_arg rejects NA, non-numeric, and non-scalar values", {
    expect_error(
        .check_alpha_arg(NA_real_),
        "strictly between 0 and 1"
    )
    expect_error(
        .check_alpha_arg("0.05"),
        "strictly between 0 and 1"
    )
    expect_error(
        .check_alpha_arg(c(0.05, 0.1)),
        "strictly between 0 and 1"
    )
})

# ------------------------------------------------------------------------
# .check_chunk_size_arg ----
# ------------------------------------------------------------------------

test_that(".check_chunk_size_arg accepts a positive integer (as a double)", {
    expect_true(.check_chunk_size_arg(1L))
    expect_true(.check_chunk_size_arg(5000))
})

test_that(".check_chunk_size_arg rejects non-int/non-pos/NA/non-numeric", {
    expect_error(
        .check_chunk_size_arg(0),
        "single positive integer"
    )
    expect_error(
        .check_chunk_size_arg(-3),
        "single positive integer"
    )
    expect_error(
        .check_chunk_size_arg(2.5),
        "single positive integer"
    )
    expect_error(
        .check_chunk_size_arg(NA_real_),
        "single positive integer"
    )
    expect_error(
        .check_chunk_size_arg("100"),
        "single positive integer"
    )
    expect_error(
        .check_chunk_size_arg(c(1, 2)),
        "single positive integer"
    )
})

# ------------------------------------------------------------------------
# .validate_exact_poisson_args ----
# ------------------------------------------------------------------------

test_that(".validate_exact_poisson_args accepts a fully valid argument set", {
    expect_true(
        .validate_exact_poisson_args(
            formula = cases ~ group,
            data = test1,
            offset = "lnfutime",
            offset_log = TRUE,
            exactmax = NULL,
            strata = "age",
            exact_vars = "group",
            ref_levels = list(),
            alpha = 0.05,
            chunk_size = 10000,
            verbose = FALSE
        )
    )
})

test_that(".validate_exact_poisson_args accepts chunk_size=NULL", {
    expect_true(
        .validate_exact_poisson_args(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = NULL,
            strata = NULL,
            exact_vars = "group",
            ref_levels = list(),
            alpha = 0.05,
            chunk_size = NULL,
            verbose = FALSE
        )
    )
    # An invalid chunk_size must be silently ignored when NULL is passed.
    expect_true(
        .validate_exact_poisson_args(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = NULL,
            strata = NULL,
            exact_vars = "group",
            ref_levels = list(),
            alpha = 0.05,
            chunk_size = NULL,
            verbose = FALSE
        )
    )
})

test_that(".validate_exact_poisson_args propagates the first failing check", {
    expect_error(
        .validate_exact_poisson_args(
            formula = "not a formula",
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = NULL,
            strata = NULL,
            exact_vars = "group",
            ref_levels = list(),
            alpha = 0.05,
            verbose = FALSE
        ),
        "must be a model formula"
    )
    expect_error(
        .validate_exact_poisson_args(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = NULL,
            strata = NULL,
            exact_vars = character(0),
            ref_levels = list(),
            alpha = 0.05,
            verbose = FALSE
        ),
        "non-empty character vector"
    )
    expect_error(
        .validate_exact_poisson_args(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = NULL,
            strata = NULL,
            exact_vars = "group",
            ref_levels = list(),
            alpha = 1.5,
            verbose = FALSE
        ),
        "strictly between 0 and 1"
    )
    expect_error(
        .validate_exact_poisson_args(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = NULL,
            strata = NULL,
            exact_vars = "group",
            ref_levels = list(),
            alpha = 0.05,
            chunk_size = -1,
            verbose = FALSE
        ),
        "single positive integer"
    )
    expect_error(
        .validate_exact_poisson_args(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            exactmax = NULL,
            strata = NULL,
            exact_vars = "group",
            ref_levels = list(),
            alpha = 0.05,
            verbose = NA
        ),
        "`verbose`"
    )
})

# ------------------------------------------------------------------------
# .validate_closedform_args ----
# ------------------------------------------------------------------------

test_that(".validate_closedform_args accepts a valid single-effect arg set", {
    expect_true(
        .validate_closedform_args(
            formula = cases ~ group,
            data = test1,
            offset = "lnfutime",
            offset_log = TRUE,
            strata = "age",
            ref_levels = list(),
            alpha = 0.05,
            verbose = FALSE
        )
    )
})

test_that(".validate_closedform_args runs shared checks before specific ones", {
    # An invalid `data` should be caught immediately, since
    # .check_data_arg() runs first in the function body.
    expect_error(
        .validate_closedform_args(
            formula = cases ~ group,
            data = list(x = 1),
            offset = NULL,
            offset_log = FALSE,
            strata = NULL,
            ref_levels = list(),
            alpha = 0.05,
            verbose = FALSE
        ),
        "must be a data frame"
    )
})

test_that(".validate_closedform_args propagate errors", {
    expect_error(
        .validate_closedform_args(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            strata = NULL,
            ref_levels = list(nonexistent = "x"),
            alpha = 0.05,
            verbose = FALSE
        ),
        "not found in `data`"
    )
    expect_error(
        .validate_closedform_args(
            formula = cases ~ group,
            data = test1,
            offset = NULL,
            offset_log = FALSE,
            strata = NULL,
            ref_levels = list(),
            alpha = 0,
            verbose = FALSE
        ),
        "strictly between 0 and 1"
    )
})
