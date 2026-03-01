source("_local/TEST_MULT_PERIOD/tests/helpers.R")
library(testthat)
library(ggplot2)

test_files <- list.files("_local/TEST_MULT_PERIOD/tests",
                          pattern = "^test-.*[.]R$", full.names = TRUE)
cat(sprintf("Running %d test files...\n", length(test_files)))

for (f in test_files) {
  cat(sprintf("\n=== %s ===\n", basename(f)))
  source(f)
}
cat("\n=== ALL TESTS COMPLETE ===\n")
