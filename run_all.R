###############################################################################
#  run_all.R — MASTER SCRIPT
#
#  Reproduces the entire study from the reference series
#  data/derived/reference_adamaoua_weekly_2018_2025.csv (417 ISO weeks).
#  That file is assembled beforehand, from the extractions in data/raw/, by
#      python data/extend_reference_2018_2025.py
#  which only needs to be rerun if a raw extraction changes.
#  Open malaria-imputation-benchmark.Rproj in RStudio, then:
#      source("run_all.R")
#
#  TOTAL RUNTIME at CFG$n_rep = 20: roughly 7 hours.
#  For a dry run, set CFG$n_rep = 2 in R/functions.R first: roughly 2.5 to
#  3 hours, because the mechanism arms (R/05, R/07) keep CFG_MECH$n_rep = 10;
#  roughly 1.5 hours if CFG_MECH$n_rep is also set to 2 in
#  R/functions_mechanism.R.
#
#  Step 4 (R/10) re-imputes the single example dataset drawn in the
#  illustrative figure. Steps 5 to 8 are the reduced mechanism design
#  (Additional files 2 and 6). Step 9 (R/04) comes after them because its
#  decision guide quotes the reduced mechanism design.
###############################################################################

STEPS <- list(
  list(file = "R/01_build_reference.R",
       what = "Certify the reference dataset (2018-2025)",      minutes = 1),
  list(file = "R/02_run_benchmark.R",
       what = "Run the benchmark (masks, imputations, metrics)", minutes = 300),
  list(file = "R/03_benchmark_analysis.R",
       what = "Analysis, statistical tests, downstream forecasting", minutes = 60),
  list(file = "R/10_example_series.R",
       what = "Example imputed series for the illustrative figure", minutes = 2),
  list(file = "R/05_mechanism_sensitivity.R",
       what = "Mechanism sensitivity: MAR and MNAR arms",      minutes = 75),
  list(file = "R/06_mechanism_diagnostics.R",
       what = "Corrected mechanism diagnostics",              minutes = 1),
  list(file = "R/07_mar_strong_arm.R",
       what = "Strong MAR arm, calibrated to MNAR intensity",  minutes = 35),
  list(file = "R/06_mechanism_diagnostics.R",
       what = "Diagnostics again, now including MAR_strong",   minutes = 1),
  list(file = "R/04_final_analysis.R",
       what = "District analysis, decision guide, figures",   minutes = 3),
  list(file = "R/11_supplementary_statistics.R",
       what = "Supplementary statistics (Additional files 2 and 6)", minutes = 20),
  list(file = "R/99_session_info.R",
       what = "Record the software environment",              minutes = 1)
)

cat("=====================================================================\n")
cat(" MALARIA IMPUTATION BENCHMARK — full reproduction\n")
cat(" Estimated total:", sum(sapply(STEPS, `[[`, "minutes")), "minutes\n")
cat("=====================================================================\n\n")

t0 <- Sys.time()
for (i in seq_along(STEPS)) {
  s <- STEPS[[i]]
  cat(sprintf("\n>>> STEP %d/%d — %s\n    (%s, ~%d min)\n",
              i, length(STEPS), s$what, s$file, s$minutes))
  ti <- Sys.time()
  # Each script runs in the GLOBAL environment: they are written to be
  # source()d, and several of them depend on objects created by R/functions.R.
  source(s$file, echo = FALSE)
  cat(sprintf("    done in %.1f min\n",
              as.numeric(difftime(Sys.time(), ti, units = "mins"))))
}

cat(sprintf("\n=== COMPLETE in %.1f h ===\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
cat("Outputs: output/ (main benchmark), output_mechanism/ (Additional files 2 and 6),\n")
cat("         docs/ENVIRONMENT.md (software record).\n")
