###############################################################################
#  10 — EXAMPLE IMPUTED SERIES FOR THE ILLUSTRATIVE FIGURE (Additional file 5)
#
#  Script 02 writes the imputed series to disk only for gap_length = 1 (they
#  feed the downstream forecasting validation). The illustrative figure also
#  needs ONE imputed dataset at gap_length = 8, to show what each method does
#  inside a long gap. This script re-imputes that single dataset.
#
#  Nothing new is simulated: the mask is the one saved by script 02
#  (repetition 1, rate 10%, gap 8), and every method is deterministic given
#  the mask (missRanger runs with seed = 1, Kalman and interpolation have no
#  random component). The script checks that the MAE it obtains equals the
#  MAE stored in output/metrics/metrics_long.csv, so the series drawn in the
#  figure are exactly those scored in the benchmark.
#
#  Prerequisite : 02_run_benchmark.R has been run (masks and metrics exist)
#  Usage        : source("R/10_example_series.R")   (project root as wd)
#  Runtime      : about one minute
#  Then         : python3 redaction/figures/make_figures.py
###############################################################################

source("R/functions.R")

EX <- list(rep = 1, rate = 10, gap = 8)
COMBOS <- list(c("Kalman", "NoTime"), c("LinearInterp", "NoTime"),
               c("MissForest_ST", "SeasonalTrend"), c("Mean", "NoTime"))

out_dir <- file.path(CFG$output_dir, "imputed_examples")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mask_file <- file.path(CFG$output_dir, "masks",
                       sprintf("mask_rep%02d_rate%02d_gap%d.rds",
                               EX$rep, EX$rate, EX$gap))
stopifnot(file.exists(mask_file))
M <- readRDS(mask_file)
X <- REF; X[M] <- NA

stored <- read_csv(file.path(CFG$output_dir, "metrics", "metrics_long.csv"),
                   show_col_types = FALSE) |>
  filter(rep == EX$rep, missing_rate == EX$rate, gap_length == EX$gap)

for (cb in COMBOS) {
  nm <- cb[1]; v <- cb[2]
  Z <- METHODS[[nm]](X, period_block(v))
  stopifnot(isTRUE(all.equal(Z[!M], REF[!M], tolerance = 1e-6)), !anyNA(Z))
  Z[Z < 0] <- 0                                  # same step as in script 02

  mae_now    <- compute_metrics(Z, M)$MAE
  mae_stored <- stored |> filter(method == nm, variant == v) |> pull(MAE)
  cat(sprintf("%-13s %-13s MAE = %.4f   stored = %.4f\n",
              nm, v, mae_now, mae_stored))
  if (!isTRUE(all.equal(mae_now, mae_stored, tolerance = 1e-6)))
    stop("The re-imputed series does not reproduce the benchmark for ", nm,
         ". Do not use it for the figure.")

  write_csv(as_tibble(Z) |> mutate(periodid = dates, .before = 1),
            file.path(out_dir, sprintf("imp_rep%02d_rate%02d_gap%d_%s_%s.csv",
                                       EX$rep, EX$rate, EX$gap, nm, v)))
}
cat("Example series written to", out_dir, "\n")
