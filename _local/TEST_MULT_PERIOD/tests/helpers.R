# _local/TEST_MULT_PERIOD/tests/helpers.R
# Common setup for all test files. Source this at top of each test file.

# Ensure we're run from package root
if (!file.exists("data/sru.rda")) {
  stop("Run tests from the package root directory (lpt/).")
}

source("_local/TEST_MULT_PERIOD/source_all.R")
library(testthat)
