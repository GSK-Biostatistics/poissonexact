# ==============================================================================
# Tests for R/ep-fit.R: the "ep_fit" object's bare constructor
# (new_ep_fit()), the validation it deliberately skips (validate_ep_fit()),
# the checked entry point every algorithm actually calls (ep_fit()), and
# the shared print.ep_fit() method. See R/ep-fit.R's own file header for
# why these four live together in one file; this test file mirrors that
# by covering all four in one place too, rather than splitting them
# across the files that happened to hold them before ep-fit.R existed.
# ==============================================================================

# ------------------------------------------------------------------------
# new_ep_fit ----
# ------------------------------------------------------------------------
#
# The bare, unchecked constructor. Every failure branch that used to live
# here now lives in validate_ep_fit()/ep_fit() below; new_ep_fit() itself
# is only ever expected to build whatever it's given, checked or not.

valid_exactparmest <- data.frame(
    Parameter = c("a", "b"),
    Estimate = c(0.1, -0.2),
    StdError = c(0.2, 0.3),
    Alpha = c(0.05, 0.05),
    LowerCL = c(-0.1, -0.6),
    UpperCL = c(0.3, 0.1),
    PValue = c(0.5, 0.2),
    SuffStat = c(3, 5)
)
valid_outdist <- list(
    a = data.frame(SuffStat = 0:1, prob = c(0.5, 0.5)),
    b = data.frame(SuffStat = 0:1, prob = c(0.3, 0.7))
)
valid_call <- quote(ep_multi_shift(cases ~ group, data = test1))

# also preserves exactparmest/outdist/call
test_that("new_ep_fit builds object, subclass plus 'ep_fit' base class", {
    out <- new_ep_fit(
        exactparmest = valid_exactparmest,
        outdist = valid_outdist,
        call = valid_call,
        subclass = "ep_multi_shift",
        algorithm_label = paste(
            "Exact conditional Poisson regression",
            "(multivariate shift)"
        )
    )

    expect_s3_class(out, c("ep_multi_shift", "ep_fit"))
    expect_s3_class(out, "ep_multi_shift")
    expect_s3_class(out, "ep_fit")
    expect_identical(out$exactparmest, valid_exactparmest)
    expect_identical(out$OUTDIST, valid_outdist)
    expect_identical(out$call, valid_call)
    expect_identical(
        attr(out, "algorithm_label"),
        "Exact conditional Poisson regression (multivariate shift)"
    )
})

# most-specific first; mirrors how ep_network() uses it
test_that("new_ep_fit accepts a multi-element subclass vector, 'ep_fit' last", {
    out <- new_ep_fit(
        exactparmest = valid_exactparmest,
        outdist = valid_outdist,
        call = valid_call,
        subclass = c("ep_network", "ep_multi_shift"),
        algorithm_label = "Exact conditional Poisson (network algorithm)"
    )

    expect_s3_class(out, c("ep_network", "ep_multi_shift", "ep_fit"))
    expect_s3_class(out, "ep_network")
    expect_s3_class(out, "ep_multi_shift")
    expect_s3_class(out, "ep_fit")
})

# unlike ep_fit(), which does validate
test_that("new_ep_fit does not validate: builds object from malformed pieces", {
    # No required columns, outdist doesn't match, call isn't a call: none
    # of this is checked by new_ep_fit() itself (see validate_ep_fit()
    # and ep_fit() below for the checked entry point every algorithm
    # actually uses). This is deliberate: it is what lets the test suite
    # build deliberately malformed fixtures for print.ep_fit() without
    # validation getting in the way.
    out <- new_ep_fit(
        exactparmest = data.frame(Parameter = "a"),
        outdist = NULL,
        call = "not a call",
        subclass = "ep_closed_form",
        algorithm_label = "label"
    )
    expect_s3_class(out, "ep_fit")
    expect_null(out$OUTDIST)
})

# ------------------------------------------------------------------------
# validate_ep_fit / ep_fit ----
# ------------------------------------------------------------------------

test_that("validate_ep_fit passes silently on well-formed pieces", {
    expect_silent(
        validate_ep_fit(valid_exactparmest, valid_outdist, valid_call)
    )
})

test_that("ep_fit builds the same object new_ep_fit would, once validated", {
    out <- ep_fit(
        exactparmest = valid_exactparmest,
        outdist = valid_outdist,
        call = valid_call,
        subclass = "ep_multi_shift",
        algorithm_label = "label"
    )
    expect_identical(
        out,
        new_ep_fit(
            exactparmest = valid_exactparmest,
            outdist = valid_outdist,
            call = valid_call,
            subclass = "ep_multi_shift",
            algorithm_label = "label"
        )
    )
})

test_that("ep_fit errors when exactparmest is not a data frame", {
    expect_error(
        ep_fit(
            exactparmest = list(Parameter = "a"),
            outdist = valid_outdist,
            call = valid_call,
            subclass = "ep_closed_form",
            algorithm_label = "label"
        ),
        "`exactparmest` must be a data frame"
    )
})

test_that("ep_fit errors when exactparmest is missing a required column", {
    keep_cols <- setdiff(names(valid_exactparmest), "PValue")
    incomplete <- valid_exactparmest[, keep_cols]
    expect_error(
        ep_fit(
            exactparmest = incomplete,
            outdist = valid_outdist,
            call = valid_call,
            subclass = "ep_closed_form",
            algorithm_label = "label"
        ),
        "missing column\\(s\\): PValue"
    )
})

test_that("ep_fit reports every missing column at once, not just the first", {
    incomplete <- valid_exactparmest[, c("Parameter", "Estimate")]
    expect_error(
        ep_fit(
            exactparmest = incomplete,
            outdist = valid_outdist,
            call = valid_call,
            subclass = "ep_closed_form",
            algorithm_label = "label"
        ),
        "StdError, Alpha, LowerCL, UpperCL, PValue, SuffStat"
    )
})

test_that("ep_fit errors when outdist is not a list", {
    expect_error(
        ep_fit(
            exactparmest = valid_exactparmest,
            outdist = c(a = 1, b = 2),
            call = valid_call,
            subclass = "ep_closed_form",
            algorithm_label = "label"
        ),
        "`outdist` must be a named list"
    )
})

test_that("ep_fit errors when outdist names don't match exactparmest names", {
    bad_outdist <- valid_outdist
    names(bad_outdist) <- c("a", "wrong_name")
    expect_error(
        ep_fit(
            exactparmest = valid_exactparmest,
            outdist = bad_outdist,
            call = valid_call,
            subclass = "ep_closed_form",
            algorithm_label = "label"
        ),
        "`outdist` must be a named list"
    )
})

test_that("ep_fit errors when outdist and exactparmest lengths differ", {
    expect_error(
        ep_fit(
            exactparmest = valid_exactparmest,
            outdist = valid_outdist["a"],
            call = valid_call,
            subclass = "ep_closed_form",
            algorithm_label = "label"
        ),
        "`outdist` must be a named list"
    )
})

test_that("ep_fit errors when call is not a call object", {
    expect_error(
        ep_fit(
            exactparmest = valid_exactparmest,
            outdist = valid_outdist,
            call = "not a call",
            subclass = "ep_closed_form",
            algorithm_label = "label"
        ),
        "`call` must be a call object"
    )
})

# ------------------------------------------------------------------------
# print.ep_fit ----
# ------------------------------------------------------------------------

test_that("print.ep_fit prints the stored algorithm label and a matching rule", {
    fake <- new_ep_fit(
        exactparmest = data.frame(
            Parameter = "group",
            Estimate = 0.5,
            StdError = 0.1,
            Alpha = 0.05,
            LowerCL = 0.3,
            UpperCL = 0.7,
            PValue = 0.01,
            SuffStat = 12
        ),
        outdist = list(
            group = data.frame(
                SuffStat = 12,
                LogWeight = 0
            )
        ),
        call = quote(
            ep_multi_shift(
                cases ~ group,
                data = test1
            )
        ),
        subclass = "ep_multi_shift",
        algorithm_label = "Exact conditional Poisson regression (multivariate shift)"
    )

    expect_snapshot(print(fake))
})

test_that("print.ep_fit falls back to a default label when none is set", {
    # new_ep_fit() always requires an algorithm_label, so to exercise the
    # %||% fallback for a *missing* label we drop the attribute afterward
    # rather than passing one in; this is the one case new_ep_fit() itself
    # cannot construct directly.
    fake <- new_ep_fit(
        exactparmest = data.frame(
            Parameter = "x",
            Estimate = 0,
            StdError = 0,
            Alpha = 0.05,
            LowerCL = 0,
            UpperCL = 0,
            PValue = 1,
            SuffStat = 0
        ),
        outdist = list(x = data.frame(SuffStat = 0, LogWeight = 0)),
        call = quote(f()),
        subclass = character(0),
        algorithm_label = NULL
    )

    expect_snapshot(print(fake))
})

test_that("print.ep_fit returns its argument invisibly", {
    fake <- new_ep_fit(
        exactparmest = data.frame(
            Parameter = "x",
            Estimate = 0,
            StdError = 0,
            Alpha = 0.05,
            LowerCL = 0,
            UpperCL = 0,
            PValue = 1,
            SuffStat = 0
        ),
        outdist = list(
            x = data.frame(
                SuffStat = 0,
                LogWeight = 0
            )
        ),
        call = quote(f()),
        subclass = character(0),
        algorithm_label = "label"
    )
    expect_snapshot({
        ret <- withVisible(print(fake))
    })
    expect_false(ret$visible)
    expect_identical(ret$value, fake)
})
