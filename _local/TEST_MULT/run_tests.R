# run_tests.R — Run from package root: source("_local/TEST_MULT/run_tests.R")
source("_local/TEST_MULT/source_all.R")

test_files <- list.files("_local/TEST_MULT/tests", pattern = "^test-.*\\.R$",
                          full.names = TRUE)
for (tf in test_files) {
  cat(sprintf("\n=== Running %s ===\n", basename(tf)))
  tryCatch(source(tf, local = TRUE), error = function(e) {
    cat(sprintf("FAIL: %s\n", conditionMessage(e)))
  })
}
cat("\n=== All test files executed ===\n")
