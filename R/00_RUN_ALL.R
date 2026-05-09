###############################################################################
# 00_RUN_ALL.R
# Run all analysis scripts in order (01 -> 10)
# Open this file in RStudio and click Source — all scripts will run in sequence.
#
# BEFORE RUNNING:
#   1. Place nuclear_cpi_panel.xlsx inside the data/ sub-folder
#   2. Run from the R/ folder (or set working directory to R/)
#   3. Install required packages if needed:
#      install.packages(c("readxl","dplyr","tidyr","plm","Matrix",
#                         "MASS","sandwich","car","urca","ggplot2"))
###############################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

run_script <- function(fname) {
  cat(rep("=", 55), "\n", sep = "")
  cat("Running:", fname, "\n")
  cat(rep("=", 55), "\n", sep = "")
  tryCatch(
    source(fname, local = FALSE),
    error = function(e) {
      cat("ERROR in", fname, ":\n")
      cat(conditionMessage(e), "\n")
      cat("Fix the error above then re-run.\n")
      stop(e)
    }
  )
  cat("\nDone:", fname, "\n\n")
}

run_script("01_data_loading.R")
run_script("02_descriptive_statistics.R")
run_script("03_unit_root_tests.R")
run_script("04_crosssection_dependence.R")
run_script("05_cointegration_tests.R")
run_script("06_fmsur_main.R")
run_script("07_groupmean_fmols.R")
run_script("08_pooled_cpr.R")
run_script("09_fmsur_ardl.R")
run_script("10_causality.R")

cat(rep("=", 55), "\n", sep = "")
cat("ALL SCRIPTS COMPLETE.\n")
cat("Outputs saved to the output/ folder:\n")
cat("  table_2_descriptive.csv\n")
cat("  table_3_cd_heterogeneity.csv\n")
cat("  table_4_unit_root.csv\n")
cat("  table_5_cointegration.csv\n")
cat("  table_6_fmsur.csv\n")
cat("  table_7_robustness.csv\n")
cat("  table_8_causality.csv\n")
cat(rep("=", 55), "\n", sep = "")
