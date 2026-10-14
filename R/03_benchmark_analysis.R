###############################################################################
#  03 — BENCHMARK ANALYSIS AND VALIDATION
#
#  Run AFTER 02_run_benchmark.R.
#  This script RE-READS the saved result files and recomputes NO imputation.
#  It therefore picks up exactly where 02 stopped, without risking the hours
#  of computation already invested.
#
#  Usage : source("R/03_benchmark_analysis.R")
#
#  SEPARATION OF CONCERNS
#  02 GENERATES (masks, imputations, metrics). 03 ANALYSES.
#  In the French version both scripts performed the analysis, which meant the
#  expensive downstream forecasting ran twice and two scripts wrote the same
#  CSV files. Everything analytical now lives here, and here only.
###############################################################################

# functions.R gives us CFG, the reference data, districts, UNIVARIATE and the
# packages, with exactly the same configuration the benchmark ran under.
source("R/functions.R")

MET <- file.path(CFG$output_dir, "metrics")

## ------------------------------------------------- 0. RE-READ THE OUTPUTS
metrics       <- read_csv(file.path(MET, "metrics_long.csv"),       show_col_types = FALSE)
metrics_distr <- read_csv(file.path(MET, "metrics_by_district.csv"), show_col_types = FALSE)
validation    <- read_csv(file.path(MET, "mcar_validation.csv"),    show_col_types = FALSE)
hyperparams   <- read_csv(file.path(MET, "hyperparameters.csv"),    show_col_types = FALSE)

cat(sprintf("Re-read: %d imputations, %d repetitions, %d method x variant combinations\n",
            nrow(metrics), n_distinct(metrics$rep),
            n_distinct(paste(metrics$method, metrics$variant))))
cat(sprintf("Computation time already invested: %.1f h\n\n", sum(metrics$sec) / 3600))


## ------------------------------ 1. MCAR VALIDATION (CORRECTED CRITERION) ---
#
#  THIS REPLACES AN EARLIER, INCORRECT ASSERTION.
#  The old version required ALL 60 gap_length = 1 masks to pass both of their
#  tests at alpha = 0.05. But under H0 -- that is, when the simulation is
#  PERFECTLY correct -- each test fails 5 times out of 100 by pure chance.
#  Demanding 120 consecutive successes means demanding an event of probability
#  0.95^120 = 0.2%. The assertion was therefore doomed to fail as soon as
#  n_rep grew: it passed at n_rep = 2 (12 tests, 54% chance) and failed at
#  n_rep = 20.
#
#  What must be tested is not "no failures" but "no MORE failures than chance
#  produces". Two checks:
#    (1) binomial test: does the failure rate exceed 5%?
#    (2) uniformity of the p-values, which is the true signature of H0.
mcar_summary <- validation |>
  group_by(gap_length) |>
  summarise(n_masks = n(),
            pct_valid = 100 * mean(valid),
            p_ks_min = min(p_values_ks), .groups = "drop")
print(mcar_summary)

val1    <- validation |> filter(gap_length == 1)
n_tests <- 2 * nrow(val1)
n_fail  <- sum(val1$p_time <= 0.05) + sum(val1$p_values_ks <= 0.05)
bt      <- binom.test(n_fail, n_tests, p = 0.05, alternative = "greater")

cat(sprintf("\ngap 1: %d failures out of %d tests (expected under H0: %.1f)\n",
            n_fail, n_tests, 0.05 * n_tests))
cat(sprintf("Binomial test (are there TOO MANY failures?): p = %.3f\n", bt$p.value))
cat("Uniformity of the p-values (expected if the simulation is correct):\n")
for (cc in c("p_districts", "p_time", "p_values_ks"))
  cat(sprintf("   %-14s KS vs U(0,1): p = %.3f | mean = %.3f\n", cc,
              suppressWarnings(ks.test(val1[[cc]], "punif")$p.value),
              mean(val1[[cc]])))

# We halt only if failures are ABNORMALLY numerous.
stopifnot(bt$p.value > 0.01)
cat("-> MCAR simulation validated: the failure rate is consistent with chance.\n")
write_csv(mcar_summary, file.path(MET, "mcar_validation_summary.csv"))

cat("\n--- k selected by internal re-masking ---\n")
print(hyperparams |> count(missing_rate, value) |>
        pivot_wider(names_from = missing_rate, values_from = n, values_fill = 0))


## ------------------------------------ 2. ACHIEVED PRECISION (n_rep) --------
# Justifies n_rep by a precision calculation rather than by habit.
achieved_precision <- metrics |>
  filter(gap_length == 1, variant == "NoTime",
         method %in% c("LinearInterp", "MissForest")) |>
  select(rep, missing_rate, method, MAE) |>
  pivot_wider(names_from = method, values_from = MAE) |>
  mutate(diff = MissForest - LinearInterp) |>
  summarise(mean_difference = mean(diff), sd_diff = sd(diff),
            half_ci95 = 1.96 * sd(diff) / sqrt(n()), n = n())
cat("\n--- Achieved precision on the paired difference ---\n")
print(achieved_precision)


## --------------------------------------------------- 3. MAIN TABLE --------
#
#  BEST VARIANT OF EACH METHOD
#  Table 1 of the paper reports every method AT ITS BEST. That same selection
#  is then reused for the Friedman test (section 6), so that the statistical
#  ranking and the table describe the same object. Reporting the best variant
#  in the table while testing the baseline variant is a classic inconsistency
#  -- and a reviewer will catch it.
#  with_ties = FALSE is essential: Mean ignores the period covariates, so its
#  three variants are strictly identical and slice_min would return all three,
#  which produces list-columns downstream.
best_variant <- metrics |>
  filter(gap_length == 1) |>
  group_by(method, variant) |>
  summarise(m = mean(MAE), .groups = "drop") |>
  group_by(method) |> slice_min(m, n = 1, with_ties = FALSE) |>
  ungroup() |> select(method, variant)
cat("\n--- Variant retained for each method ---\n"); print(best_variant, n = 20)
write_csv(best_variant, file.path(MET, "selected_variant.csv"))

# Full detail: every variant, every metric (-> supplementary material)
table_full <- metrics |>
  filter(gap_length == 1) |>
  group_by(method, variant, missing_rate) |>
  summarise(across(c(MAE, RMSE, MAPE, bias, var_ratio, youden),
                   list(m = mean, s = sd)),
            sec = mean(sec), .groups = "drop") |>
  mutate(MAE_txt  = sprintf("%.1f (%.1f)", MAE_m, MAE_s),
         RMSE_txt = sprintf("%.1f (%.1f)", RMSE_m, RMSE_s)) |>
  arrange(missing_rate, MAE_m)
write_csv(table_full, file.path(MET, "table_full_supplement.csv"))

#  TABLE 1 OF THE PAPER: one row per method, the three error metrics side by
#  side. The ranking stays based on MAE, but RMSE and MAPE are shown because
#  readers look for them.
table1_paper <- metrics |>
  filter(gap_length == 1) |>
  inner_join(best_variant, by = c("method", "variant")) |>
  group_by(method, variant, missing_rate) |>
  summarise(across(c(MAE, RMSE, MAPE), list(m = mean, s = sd)), .groups = "drop") |>
  mutate(MAE_fmt  = sprintf("%.1f (%.1f)", MAE_m,  MAE_s),
         RMSE_fmt = sprintf("%.1f (%.1f)", RMSE_m, RMSE_s),
         MAPE_fmt = sprintf("%.1f (%.1f)", MAPE_m, MAPE_s)) |>
  select(method, variant, missing_rate, MAE_fmt, RMSE_fmt, MAPE_fmt) |>
  pivot_wider(names_from = missing_rate,
              values_from = c(MAE_fmt, RMSE_fmt, MAPE_fmt), names_sep = "_") |>
  left_join(
    metrics |> filter(gap_length == 1) |>
      inner_join(best_variant, by = c("method", "variant")) |>
      group_by(method) |>
      summarise(MAE_overall = mean(MAE),
                sensitivity = mean(sensitivity, na.rm = TRUE),
                youden      = mean(youden, na.rm = TRUE),
                sec         = mean(sec), .groups = "drop"),
    by = "method") |>
  arrange(MAE_overall)
cat("\n=== TABLE 1 OF THE PAPER (best variant, isolated gaps) ===\n")
print(table1_paper, n = 20, width = Inf)
write_csv(table1_paper, file.path(MET, "table1_paper.csv"))


## --------------------------------- 4. ACCURACY VS PEAK DETECTION ----------
cat("\n--- Average accuracy (MAE) VS peak detection (Youden) ---\n")
peak_detection <- metrics |>
  filter(gap_length == 1) |>
  group_by(method, variant) |>
  summarise(MAE = mean(MAE), sensitivity = mean(sensitivity, na.rm = TRUE),
            specificity = mean(specificity, na.rm = TRUE),
            youden = mean(youden, na.rm = TRUE), .groups = "drop") |>
  mutate(rank_MAE = rank(MAE), rank_youden = rank(-youden),
         discordance = rank_youden - rank_MAE) |>
  arrange(rank_MAE)
print(peak_detection, n = 100)
write_csv(peak_detection, file.path(MET, "peak_detection.csv"))


## ----------------------------------------- 5. EFFECT OF THE PERIOD BLOCK ---
# gap_length IS ESSENTIAL in the pivot: without it every (rep, rate, method,
# variant) combination exists in triplicate (one per missingness structure),
# pivot_wider creates list-columns and mean() returns NA everywhere.
period_effect <- metrics |>
  filter(!method %in% UNIVARIATE) |>
  select(rep, missing_rate, gap_length, method, variant, MAE) |>
  pivot_wider(names_from = variant, values_from = MAE) |>
  group_by(method, gap_length, missing_rate) |>
  summarise(
    MAE_NoTime        = mean(NoTime),
    MAE_Seasonal      = mean(Seasonal),
    MAE_SeasonalTrend = mean(SeasonalTrend),
    # gain from SEASONALITY alone: the epidemiologically interesting result
    gain_seasonal_pct = 100 * (mean(NoTime) - mean(Seasonal)) / mean(NoTime),
    # extra gain contributed by the bare time index
    gain_trend_pct    = 100 * (mean(Seasonal) - mean(SeasonalTrend)) / mean(NoTime),
    gain_total_pct    = 100 * (mean(NoTime) - mean(SeasonalTrend)) / mean(NoTime),
    # paired tests: same mask on both sides, therefore legitimate
    p_seasonal = tryCatch(t.test(NoTime, Seasonal, paired = TRUE)$p.value,
                          error = \(e) NA),
    p_total    = tryCatch(t.test(NoTime, SeasonalTrend, paired = TRUE)$p.value,
                          error = \(e) NA),
    .groups = "drop")
write_csv(period_effect, file.path(MET, "period_effect.csv"))
cat("\n--- Period effect, isolated gaps ---\n")
print(period_effect |> filter(gap_length == 1), n = 100)
cat("\n--- Period effect, 8-week gaps ---\n")
print(period_effect |> filter(gap_length == 8), n = 100)


## ----------------------------------------------- 6. STATISTICAL TESTS -----
#
#  Three analyses:
#   (A) MAIN         — best variant of each method, MAE criterion.
#                      This is the one that matches Table 1.
#   (B) SUPPLEMENT   — NoTime variant only: "methods without any temporal
#                      information". The contrast with (A) is informative.
#   (C) ROBUSTNESS   — best variant, RMSE criterion. Lets us state that the
#                      conclusion does not depend on the chosen metric.

block_matrix <- function(df, column) {
  df |> unite("blk", rep, missing_rate) |>
    select(blk, method, all_of(column)) |>
    pivot_wider(names_from = method, values_from = all_of(column)) |>
    column_to_rownames("blk") |> as.matrix()
}

friedman_analysis <- function(mat, title, nemenyi = TRUE) {
  k <- ncol(mat); N <- nrow(mat)
  CD <- qtukey(0.95, k, Inf) / sqrt(2) * sqrt(k * (k + 1) / (6 * N))
  cat(sprintf("\n=========== %s ===========\n", title))
  cat(sprintf("Friedman on %d blocks x %d methods | critical difference CD = %.2f\n",
              N, k, CD))
  print(friedmanTest(mat))
  if (nemenyi) print(frdAllPairsNemenyiTest(mat))
  r <- tibble(method = colnames(mat),
              mean_rank = colMeans(t(apply(mat, 1, rank)))) |> arrange(mean_rank)
  cat("\nMean ranks (1 = best):\n"); print(r, n = 100)
  attr(r, "CD") <- CD
  r
}

base_best   <- metrics |> filter(gap_length == 1) |>
  inner_join(best_variant, by = c("method", "variant"))
base_notime <- metrics |> filter(gap_length == 1, variant == "NoTime")

# (A) MAIN
ranks_main <- friedman_analysis(block_matrix(base_best, "MAE"),
                                "A. MAIN — best variant, MAE criterion")
write_csv(ranks_main, file.path(MET, "friedman_ranks_main.csv"))

# (B) SUPPLEMENT
ranks_notime <- friedman_analysis(block_matrix(base_notime, "MAE"),
                                  "B. SUPPLEMENT — NoTime variant (no temporal information)")
write_csv(ranks_notime, file.path(MET, "friedman_ranks_notime.csv"))

# (C) ROBUSTNESS
ranks_rmse <- friedman_analysis(block_matrix(base_best, "RMSE"),
                                "C. ROBUSTNESS — best variant, RMSE criterion",
                                nemenyi = FALSE)
write_csv(ranks_rmse, file.path(MET, "friedman_ranks_rmse.csv"))

cat("\n--- Do the three rankings agree? ---\n")
rank_comparison <- ranks_main |> rename(rank_MAE_best = mean_rank) |>
  left_join(ranks_notime |> rename(rank_MAE_notime = mean_rank), by = "method") |>
  left_join(ranks_rmse   |> rename(rank_RMSE_best  = mean_rank), by = "method")
print(rank_comparison, n = 20)
cat(sprintf("\n  Spearman(MAE best variant, MAE NoTime variant) = %+.3f\n",
            cor(rank_comparison$rank_MAE_best, rank_comparison$rank_MAE_notime,
                method = "spearman")))
cat(sprintf("  Spearman(MAE best variant, RMSE)               = %+.3f\n",
            cor(rank_comparison$rank_MAE_best, rank_comparison$rank_RMSE_best,
                method = "spearman")))
write_csv(rank_comparison, file.path(MET, "ranking_comparison.csv"))


## ---------------------------- 6bis. CONCORDANCE BETWEEN METRICS -----------
#
#  Justifies not drowning the reader in six metrics: MAE, RMSE, MAPE and R2
#  measure the same thing. What adds NEW information is the variance ratio
#  and Youden's index.
#  (MSE is not reported: the square root is a monotone transformation, so its
#   ranking is strictly identical to that of RMSE.)
metric_ranks <- base_best |>
  group_by(method) |>
  summarise(MAE = mean(MAE), RMSE = mean(RMSE), MAPE = mean(MAPE),
            R2 = mean(R2), var_ratio = mean(var_ratio),
            youden = mean(youden, na.rm = TRUE), .groups = "drop") |>
  mutate(r_MAE = rank(MAE), r_RMSE = rank(RMSE), r_MAPE = rank(MAPE),
         r_R2 = rank(-R2), r_var_ratio = rank(abs(var_ratio - 1)),
         r_youden = rank(-youden))

rank_cols <- c("r_MAE", "r_RMSE", "r_MAPE", "r_R2", "r_var_ratio", "r_youden")
metric_concordance <- expand_grid(a = rank_cols, b = rank_cols) |>
  filter(a < b) |>
  mutate(rho = map2_dbl(a, b, ~ cor(metric_ranks[[.x]], metric_ranks[[.y]],
                                    method = "spearman")),
         p   = map2_dbl(a, b, ~ suppressWarnings(
           cor.test(metric_ranks[[.x]], metric_ranks[[.y]],
                    method = "spearman")$p.value)),
         a = str_remove(a, "^r_"), b = str_remove(b, "^r_"))
cat("\n=== Concordance between the per-metric rankings ===\n")
print(metric_concordance |> arrange(desc(rho)) |>
        mutate(rho = round(rho, 3), p = signif(p, 2)), n = 30)
write_csv(metric_ranks,       file.path(MET, "ranks_by_metric.csv"))
write_csv(metric_concordance, file.path(MET, "metric_concordance.csv"))


## ------------------------- 7. SENSITIVITY TO THE MISSINGNESS STRUCTURE ----
gap_sensitivity <- metrics |>
  group_by(method, variant, gap_length) |>
  summarise(MAE = mean(MAE), .groups = "drop") |>
  group_by(gap_length) |> mutate(rank = rank(MAE)) |> ungroup() |>
  select(-MAE) |>
  pivot_wider(names_from = gap_length, values_from = rank, names_prefix = "rank_gap")
print(gap_sensitivity, n = 100)
write_csv(gap_sensitivity, file.path(MET, "gap_structure_sensitivity.csv"))


## ---------------------------- 8. DOWNSTREAM VALIDATION: FORECASTING -------
# WEEKLY data -> s = 52. A SARIMA[52] on a few hundred points is ill-posed, so we use
# Fourier terms + ARIMA errors.
#
# THREE CORRECTIONS ARE BAKED INTO THIS SECTION, each found during dry runs:
#   1. Train on the IMPUTED series but evaluate against the REFERENCE.
#      Evaluating against the imputed test window made mean imputation the
#      best forecaster, purely because it flattens the series
#      (rho = 0.86 between smoothing and the fake "improvement").
#   2. Fix the ARIMA order ONCE on the reference and impose it on every
#      imputed series (coefficients are re-estimated). Re-running auto.arima
#      per series measured model selection, not data quality.
#   3. Add a SHORT horizon (4 weeks, rolling origin). At 52 weeks the
#      RMSE-optimal forecast of a noisy weekly series is close to its mean,
#      so any smoothing method wins mechanically: that horizon simply cannot
#      rank imputation methods. The contrast between the two horizons is
#      itself a publishable result.
H_LONG  <- 52
H_SHORT <- 4

REF_ARIMA_ORDERS <- map(seq_len(n_dis), function(j) {
  y  <- ts(REF[, j], frequency = 52)
  tr <- subset(y, end = n_obs - H_LONG)
  arimaorder(auto.arima(tr, xreg = fourier(tr, K = 2), seasonal = FALSE))[1:3]
})
cat("\nARIMA orders fixed on the reference:\n")
walk2(districts, REF_ARIMA_ORDERS,
      ~ cat(sprintf("  %-18s ARIMA(%d,%d,%d)\n", .x, .y[1], .y[2], .y[3])))

forecast_rmse <- function(imputed_series, true_series, j) {
  y   <- ts(imputed_series, frequency = 52)
  n   <- length(y); n_tr <- n - H_LONG
  tr  <- subset(y, end = n_tr)
  xr  <- fourier(tr, K = 2)
  fit <- tryCatch(Arima(tr, order = REF_ARIMA_ORDERS[[j]], xreg = xr),
                  error = function(e) auto.arima(tr, xreg = xr, seasonal = FALSE))

  # long horizon: 52 weeks forecast in one shot
  fc  <- forecast(fit, xreg = fourier(tr, K = 2, h = H_LONG), h = H_LONG)
  h52 <- sqrt(mean((as.numeric(fc$mean) - tail(true_series, H_LONG))^2))

  # short horizon: rolling origin, 4 weeks at a time, coefficients frozen
  err <- numeric(0)
  for (t in seq(n_tr, n - H_SHORT, by = H_SHORT)) {
    hist <- ts(imputed_series[1:t], frequency = 52)
    f <- tryCatch(Arima(hist, model = fit, xreg = fourier(hist, K = 2)),
                  error = function(e) NULL)
    if (is.null(f)) next
    fcc <- forecast(f, xreg = fourier(hist, K = 2, h = H_SHORT), h = H_SHORT)
    err <- c(err, as.numeric(fcc$mean) - true_series[(t + 1):(t + H_SHORT)])
  }
  c(h52 = h52, h4 = if (length(err)) sqrt(mean(err^2)) else NA_real_)
}

# Reference: trained AND evaluated on the complete, never-punched series.
# This is the attainable error floor; any degradation must be >= 0.
ref_forecast <- map(seq_len(n_dis), ~ forecast_rmse(REF[, .x], REF[, .x], .x))
ref_h52 <- map_dbl(ref_forecast, "h52")
ref_h4  <- map_dbl(ref_forecast, "h4")
cat(sprintf("Reference forecast RMSE  h=52: %.1f   h=4: %.1f\n",
            mean(ref_h52), mean(ref_h4)))

# MIND THE COST: each file requires 11 ARIMA fits. With 20 repetitions that
# would be over 10,000 fits (several days). We restrict to the first
# CFG$forecast_reps repetitions: forecast error is far more stable across
# repetitions than imputation error, so 5 are plenty. TO BE DECLARED in the
# Methods section.
files <- list.files(file.path(CFG$output_dir, "imputed"), full.names = TRUE)
files <- files[as.integer(str_match(basename(files), "imp_rep(\\d+)")[, 2])
               <= CFG$forecast_reps]
cat(sprintf("Downstream validation: %d datasets x %d districts = %d ARIMA fits\n",
            length(files), n_dis, length(files) * n_dis))

# NOTE ON THE REGEX: SeasonalTrend must come BEFORE Seasonal in the
# alternation, otherwise "Seasonal" matches first and the trailing "Trend"
# breaks the parse, silently producing NA in the variant column.
FILE_PATTERN <- "imp_rep(\\d+)_rate(\\d+)_(.+)_(SeasonalTrend|Seasonal|NoTime)\\.csv$"
downstream <- map_dfr(files, function(f) {
  meta <- str_match(basename(f), FILE_PATTERN)
  Z  <- as.matrix(read_csv(f, show_col_types = FALSE)[, districts])
  pr <- map(seq_len(n_dis), ~ forecast_rmse(Z[, .x], REF[, .x], .x))
  tibble(rep = as.integer(meta[2]), missing_rate = as.integer(meta[3]),
         method = meta[4], variant = meta[5], district = districts,
         rmse_h52 = map_dbl(pr, "h52"), ref_h52 = ref_h52,
         rmse_h4  = map_dbl(pr, "h4"),  ref_h4  = ref_h4) |>
    mutate(degradation_h52 = rmse_h52 - ref_h52,
           degradation_h4  = rmse_h4  - ref_h4)
})
stopifnot(!any(is.na(downstream$variant)))   # guards the regex above
write_csv(downstream, file.path(MET, "downstream_forecast.csv"))

criteria_comparison <- downstream |>
  group_by(method, variant) |>
  summarise(degradation_h4  = mean(degradation_h4, na.rm = TRUE),
            degradation_h52 = mean(degradation_h52), .groups = "drop") |>
  left_join(metrics |> filter(gap_length == 1) |> group_by(method, variant) |>
              summarise(MAE = mean(MAE), youden = mean(youden, na.rm = TRUE),
                        .groups = "drop"),
            by = c("method", "variant")) |>
  mutate(rank_MAE = rank(MAE), rank_h4 = rank(degradation_h4),
         rank_h52 = rank(degradation_h52)) |>
  arrange(rank_MAE)
print(criteria_comparison, n = 100)
write_csv(criteria_comparison, file.path(MET, "criteria_comparison.csv"))

# READ THIS FIRST: if the h = 4 correlation is POSITIVE, the downstream
# validation is valid (better imputation gives better forecasts). If it stays
# negative, do NOT publish this section as a ranking -- report it as an
# explained null result.
cat("\n--- Downstream validation diagnostics ---\n")
cat(sprintf("  Spearman(MAE, degradation h=4 )  = %+.2f   <- must be POSITIVE\n",
            cor(criteria_comparison$MAE, criteria_comparison$degradation_h4,
                method = "spearman")))
cat(sprintf("  Spearman(MAE, degradation h=52)  = %+.2f   <- expected negative\n",
            cor(criteria_comparison$MAE, criteria_comparison$degradation_h52,
                method = "spearman")))
cat("  At 52 weeks every forecast converges to the mean, so smoothing methods\n",
    "  win mechanically. At 4 weeks the recent dynamics dominate and imputation\n",
    "  quality becomes measurable again.\n")


## --------------------------------- 9. REGIONAL TOTAL FIGURE ---------------
# This is the use case the Palu_Adamaoua column was originally created for.
# That column no longer exists in the files: the total is RECOMPUTED on both
# sides. It is the only way to be sure we compare the same quantity.
plot_regional_total <- function(rep = 1, rate = 20,
                                meths = c("LinearInterp", "MissForest", "Mean")) {
  ref_series <- tibble(date = dates, total = rowSums(REF), source = "Reference")
  imp <- map_dfr(meths, function(m) {
    f <- file.path(CFG$output_dir, "imputed",
                   sprintf("imp_rep%02d_rate%02d_%s_%s.csv", rep, rate, m,
                           if (m %in% UNIVARIATE) "NoTime" else "SeasonalTrend"))
    if (!file.exists(f)) return(NULL)
    Z <- as.matrix(read_csv(f, show_col_types = FALSE)[, districts])
    tibble(date = dates, total = rowSums(Z), source = m)   # <- RECOMPUTED
  })
  ggplot(bind_rows(ref_series, imp), aes(date, total, colour = source)) +
    geom_line(linewidth = .4) +
    scale_colour_manual(values = c(Reference = "black", Mean = "#d95f02",
                                   MissForest = "#1b9e77",
                                   LinearInterp = "#7570b3")) +
    labs(title = sprintf("Regional weekly malaria incidence — %d%% missing data", rate),
         subtitle = "Total recomputed from the 11 districts on both sides",
         x = NULL, y = "Malaria cases", colour = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
}
ggsave(file.path(CFG$output_dir, "figures", "regional_total_20pct.png"),
       plot_regional_total(1, 20), width = 10, height = 4.5, dpi = 300)

cat("\n=== Analysis complete. No imputation was recomputed. ===\n")
