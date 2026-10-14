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
#  For a dry run, set CFG$n_rep = 2 in R/functions.R first (~30 min).
#
#  Steps 6 to 8 are the mechanism sensitivity analysis (supplementary
#  material S4). Step 9 is the optional choropleth map, which needs `sf`
#  and an internet connection; it is skipped if either is unavailable.
###############################################################################

STEPS <- list(
  list(file = "R/01_build_reference.R",
       what = "Certify the reference dataset (2018-2025)",      minutes = 1),
  list(file = "R/02_run_benchmark.R",
       what = "Run the benchmark (masks, imputations, metrics)", minutes = 300),
  list(file = "R/03_benchmark_analysis.R",
       what = "Analysis, statistical tests, downstream forecasting", minutes = 60),
  list(file = "R/04_final_analysis.R",
       what = "District analysis, decision guide, figures",   minutes = 3),
  list(file = "R/05_mechanism_sensitivity.R",
       what = "Mechanism sensitivity: MAR and MNAR arms",      minutes = 75),
  list(file = "R/06_mechanism_diagnostics.R",
       what = "Corrected mechanism diagnostics",              minutes = 1),
  list(file = "R/07_mar_strong_arm.R",
       what = "Strong MAR arm, calibrated to MNAR intensity",  minutes = 35),
  list(file = "R/06_mechanism_diagnostics.R",
       what = "Diagnostics again, now including MAR_strong",   minutes = 1),
  list(file = "R/11_supplementary_statistics.R",
       what = "Calibration curve, paired test, reduced-design CD", minutes = 2),
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
cat("Outputs: output/ (main benchmark), output_mechanism/ (supplement S4),\n")
cat("         docs/ENVIRONMENT.md (software record).\n")
cat("\nOptional, not run automatically:\n")
cat("  R/08_district_map.R   choropleth maps (needs sf + internet)\n")
cat("  R/09_explore.R        interactive helpers for looking at the data\n")
