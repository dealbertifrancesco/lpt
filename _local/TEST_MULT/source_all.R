# source_all.R — Run from package root: source("_local/TEST_MULT/source_all.R")
if (!requireNamespace("mgcv", quietly = TRUE)) stop("mgcv required")

# Load unmodified package functions
source("R/estimate_dose_slope.R")

# Load modified functions (order matters: dependencies first)
source("_local/TEST_MULT/R/calibrate_B.R")
source("_local/TEST_MULT/R/att_bounds.R")
source("_local/TEST_MULT/R/simulate.R")
source("_local/TEST_MULT/R/lpt.R")
source("_local/TEST_MULT/R/summary.R")
source("_local/TEST_MULT/R/plot.R")

# Load data
load("data/sru.rda")

message("TEST_MULT environment loaded.")
