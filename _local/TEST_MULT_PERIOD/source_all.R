# _local/TEST_MULT_PERIOD/source_all.R
# Run this from the package root to load all modified functions.
# Usage: source("_local/TEST_MULT_PERIOD/source_all.R")

if (!requireNamespace("mgcv", quietly = TRUE)) stop("mgcv required")
if (!requireNamespace("ggplot2", quietly = TRUE)) message("ggplot2 not found - plot tests will be skipped")

# Load unmodified package functions from package R/ directory
source("R/estimate_dose_slope.R")

# Load modified functions from test directory (order matters: dependencies first)
source("_local/TEST_MULT_PERIOD/R/calibrate_B.R")
source("_local/TEST_MULT_PERIOD/R/att_bounds.R")
source("_local/TEST_MULT_PERIOD/R/simulate.R")
source("_local/TEST_MULT_PERIOD/R/lpt.R")
source("_local/TEST_MULT_PERIOD/R/summary.R")
source("_local/TEST_MULT_PERIOD/R/plot.R")

# Load SRU data
load("data/sru.rda")

message("TEST_MULT_PERIOD environment loaded.")
