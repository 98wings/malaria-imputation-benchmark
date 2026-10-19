###############################################################################
#  02 — RUN THE IMPUTATION BENCHMARK
#
#  This script GENERATES. It does not analyse.
#  It regenerates everything from the certified reference: masks, imputations
#  and per-imputation metrics. Every table, test and figure lives in
#  03_benchmark_analysis.R, which re-reads the files written here.
#
#  Prerequisite : 01_build_reference.R
#  Usage        : source("R/02_run_benchmark.R")   (project root as wd)
#  Runtime      : ~5 h at n_rep = 20. Set CFG$n_rep = 2 in R/functions.R for
#                 a dry run.
#
#  WHY GENERATION AND ANALYSIS ARE SPLIT
#  The French version performed the analysis inside this script AND again in
#  the next one. The expensive downstream forecasting therefore ran twice, and
#  two scripts wrote the same CSV files. Worse, a failure in a summary table
#  aborted the script after five hours of imputation. Here, once the CSVs are
#  written the expensive work is safe.
#
#  KEY DESIGN POINT: the mask is generated ONLY on the district columns.
#  The period covariates are added AFTERWARDS and never enter the computation
#  of the missingness rate -- while remaining available as predictors.
###############################################################################

source("R/functions.R")

res_global <- list(); res_district <- list(); res_valid <- list()
res_hp <- list(); k <- 0

# gap_length = 1, 4, 8 : design factor (isolated gaps, short and prolonged
#                        interruptions of reporting); 1 week is the reference condition
for (gap in CFG$gap_lengths) {
for (rep in seq_len(CFG$n_rep)) {
  for (rate in CFG$missing_rates) {

    M <- simulate_mask(rate, seed = 1000 * rep + rate + 7 * gap, gap_length = gap)
    saveRDS(M, file.path(CFG$output_dir, "masks",
                         sprintf("mask_rep%02d_rate%02d_gap%d.rds", rep, rate, gap)))
    X <- REF; X[M] <- NA

    # (a) Is the mask really MCAR? The trace is kept and tested in script 03.
    res_valid[[length(res_valid) + 1]] <- validate_mcar(M) |>
      mutate(rep = rep, missing_rate = rate, gap_length = gap, .before = 1)

    # (b) Tuning k WITHOUT ever looking at the truth.
    #     COST: tuning runs B x 6 KNN imputations per configuration.
    #     At 20 reps x 3 rates x 3 structures that would be 32,400 calls.
    #     We therefore tune on the first CFG$tuning_reps repetitions and reuse
    #     the modal k afterwards. This is what a practitioner would do, and
    #     the variability of k across repetitions is low.
    #     TIE RULE: the modal k is names(sort(table(.), decreasing = TRUE))[1].
    #     sort() keeps the original (increasing) order of tied values, so a
    #     tie is broken towards the SMALLEST k. This happened once in the
    #     2018-2025 run (4-week gaps, 5%: tuned values 3, 3, 5, 7, 7 -> k = 3
    #     reused for repetitions 6-20). The mechanism arms (05, 07) apply the
    #     same rule with CFG_MECH$tuning_reps = 3.
    if (rep <= CFG$tuning_reps) {
      k_opt <- tune_k(X, rate, B = CFG$tuning_B)
    } else {
      seen  <- bind_rows(res_hp) |>
        filter(missing_rate == !!rate, gap_length == !!gap)
      k_opt <- as.integer(names(sort(table(seen$value), decreasing = TRUE))[1])
    }
    res_hp[[length(res_hp) + 1]] <- tibble(
      rep = rep, missing_rate = rate, gap_length = gap,
      parameter = "k_KNN", value = k_opt, tuned = rep <= CFG$tuning_reps)
    METHODS$KNN <- function(X, P) imp_knn(X, P, k = k_opt)

    for (nm in names(METHODS)) {
      variants <- if (nm %in% UNIVARIATE) "NoTime" else CFG$variants

      for (v in variants) {
        P <- period_block(v)

        tic(); Z <- METHODS[[nm]](X, P); tps <- toc(quiet = TRUE)
        seconds <- as.numeric(tps$toc - tps$tic)

        # Guard rails: observed values untouched, no NA, no negative counts.
        # isTRUE() is essential: all.equal() returns a CHARACTER STRING when
        # it fails, and stopifnot() would then print something unreadable.
        stopifnot(isTRUE(all.equal(Z[!M], REF[!M], tolerance = 1e-6)), !anyNA(Z))
        n_neg <- sum(Z < 0)          # negative counts are a defect to report
        Z[Z < 0] <- 0

        k <- k + 1
        res_global[[k]] <- compute_metrics(Z, M) |>
          mutate(rep = rep, missing_rate = rate, gap_length = gap,
                 method = nm, variant = v,
                 sec = seconds, n_negative = n_neg, .before = 1)
        res_district[[k]] <- compute_metrics_by_district(Z, M) |>
          mutate(rep = rep, missing_rate = rate, gap_length = gap,
                 method = nm, variant = v, .before = 1)

        # Only the primary scenario is written out as CSV: these files feed
        # the downstream forecasting validation in script 03.
        if (gap == 1)
          write_csv(as_tibble(Z) |> mutate(periodid = dates, .before = 1),
                    file.path(CFG$output_dir, "imputed",
                              sprintf("imp_rep%02d_rate%02d_%s_%s.csv",
                                      rep, rate, nm, v)))
      }
    }
    cat(sprintf("  gap %d / rep %02d / rate %02d%% done\n", gap, rep, rate))
  }
}
}

metrics       <- bind_rows(res_global)
metrics_distr <- bind_rows(res_district)
validation    <- bind_rows(res_valid)
hyperparams   <- bind_rows(res_hp)

write_csv(metrics,       file.path(CFG$output_dir, "metrics", "metrics_long.csv"))
write_csv(metrics_distr, file.path(CFG$output_dir, "metrics", "metrics_by_district.csv"))
write_csv(validation,    file.path(CFG$output_dir, "metrics", "mcar_validation.csv"))
write_csv(hyperparams,   file.path(CFG$output_dir, "metrics", "hyperparameters.csv"))

# Completeness check: n_rep x rates x gaps x (7 multivariate methods x 3
# variants + 3 univariate methods x 1 variant) = n_rep x 3 x 3 x 24.
expected <- CFG$n_rep * length(CFG$missing_rates) * length(CFG$gap_lengths) *
  (sum(!names(METHODS) %in% UNIVARIATE) * length(CFG$variants) +
     sum(names(METHODS) %in% UNIVARIATE))
stopifnot(nrow(metrics) == expected)

if (N_KALMAN_FALLBACKS > 0)
  cat(sprintf("\nNOTE: %d fallbacks on na_kalman(StructTS). If this number is\n",
              N_KALMAN_FALLBACKS),
      "large, report it in the Methods section.\n")

cat(sprintf("\n=== Generation complete: %d imputations, %.1f h of compute ===\n",
            nrow(metrics), sum(metrics$sec) / 3600))
cat("Files written to '", CFG$output_dir, "'.\n", sep = "")
cat("Next step: source(\"R/03_benchmark_analysis.R\")\n")
