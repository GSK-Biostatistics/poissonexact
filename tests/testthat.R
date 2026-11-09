# This file is part of the standard testthat setup for an R package. It is
# picked up automatically by `R CMD check` and by `devtools::test()`.
#
# Internal (non-exported, including dot-prefixed) functions such as
# parse_exactmax_spec(), log_sum_exp(), .setup_exact_model(), etc. are
# reachable by name (no `:::` needed) from the test files below, because
# test_check()/devtools::test() run the tests inside the package's own
# namespace environment.

library(testthat)
library(poissonexact)

test_check("poissonexact")
