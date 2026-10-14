###############################################################################
#  05 — SENSITIVITY TO THE MISSINGNESS MECHANISM (MAR and MNAR arms)
#
#  SCIENTIFIC QUESTION
#  The whole benchmark rests on MCAR (data-entry error). That is a justified
#  assumption for these data, but a reviewer will object:
#    "Do your conclusions hold if missingness is not completely at random?"
#  This script answers by re-running a subset of the benchmark under two
#  other mechanisms, EVERYTHING ELSE BEING IDENTICAL.
#
#  What is tested: the ROBUSTNESS OF THE RANKING, not the level of error.
#  (Under MAR/MNAR the masked cells are systematically higher, so MAE rises
#   mechanically for every method. Comparing LEVELS across mechanisms is
#   meaningless; comparing RANKS is not.)
#
#  Prerequisite : 02_run_benchmark.R (output/ must exist)
#  Runtime      : 60-80 min (against 5.1 h for the main run)
#  Next         : 06_mechanism_diagnostics.R, then 07_mar_strong_arm.R
###############################################################################

source("R/functions_mechanism.R")

MECHANISMS <- c("MAR", "MNAR")
#  MNAR is included by default: it costs only ~30 min more and it is the
#  SECOND question a reviewer systematically asks after MAR. Answering both
#  at once avoids a round of revision. Set MECHANISMS <- "MAR" to skip it.

## -------------------------------------------------------------------------
## 1. LOOP
## -------------------------------------------------------------------------
res_g <- list(); res_d <- list(); res_v <- list(); res_k <- list(); kk <- 0

for (mech in MECHANISMS) {
for (gap in CFG_MECH$gap_lengths) {
for (rep in seq_len(CFG_MECH$n_rep)) {
  rate <- CFG_MECH$missing_rate
  # Seeds disjoint from the MCAR run (1000*rep + rate + 7*gap)
  seed <- switch(mech, MAR = 500000L, MNAR = 900000L) + 1000L*rep + rate + 7L*gap
  M <- simulate_mask_mechanism(rate, seed, gap, mech)
  saveRDS(M, file.path(DOS_MECH, "masks",
                       sprintf("mask_%s_rep%02d_gap%d.rds", mech, rep, gap)))
  X <- REF; X[M] <- NA

  res_v[[length(res_v) + 1]] <- diagnose_mechanism(M, mech) |>
    mutate(rep = rep, gap_length = gap, .before = 1)

  # Tuning k: we keep the MCAR internal re-masking of tune_k(). This is
  # DELIBERATE. A practitioner does not know the true mechanism; they
  # re-mask naively at random. Using the true mechanism here would give KNN
  # information nobody has in practice.
  if (rep <= CFG_MECH$tuning_reps) {
    k_opt <- tune_k(X, rate, B = CFG_MECH$tuning_B)
  } else {
    seen  <- bind_rows(res_k) |> filter(mechanism == mech, gap_length == !!gap)
    k_opt <- as.integer(names(sort(table(seen$value), decreasing = TRUE))[1])
  }
  res_k[[length(res_k) + 1]] <- tibble(
    mechanism = mech, rep = rep, gap_length = gap,
    parameter = "k_KNN", value = k_opt, tuned = rep <= CFG_MECH$tuning_reps)
  METHODS$KNN <- function(X, P) imp_knn(X, P, k = k_opt)

  for (nm in names(METHODS)) {
    variants <- if (nm %in% UNIVARIATE) "NoTime" else CFG$variants
    for (v in variants) {
      P <- period_block(v)
      tic(); Z <- METHODS[[nm]](X, P); tps <- toc(quiet = TRUE)
      stopifnot(isTRUE(all.equal(Z[!M], REF[!M], tolerance = 1e-6)), !anyNA(Z))
      Z[Z < 0] <- 0
      kk <- kk + 1
      res_g[[kk]] <- compute_metrics(Z, M) |>
        mutate(mechanism = mech, rep = rep, missing_rate = rate,
               gap_length = gap, method = nm, variant = v,
               sec = as.numeric(tps$toc - tps$tic), .before = 1)
      res_d[[kk]] <- compute_metrics_by_district(Z, M) |>
        mutate(mechanism = mech, rep = rep, gap_length = gap,
               method = nm, variant = v, .before = 1)
    }
  }
  cat(sprintf("  %s / gap %d / rep %02d done\n", mech, gap, rep))
}}}

metrics_mech <- bind_rows(res_g); distr_mech <- bind_rows(res_d)
valid_mech   <- bind_rows(res_v); hp_mech    <- bind_rows(res_k)
write_csv(metrics_mech, file.path(DOS_MECH, "metrics", "metrics_mechanism.csv"))
write_csv(distr_mech,   file.path(DOS_MECH, "metrics", "district_mechanism.csv"))
write_csv(valid_mech,   file.path(DOS_MECH, "metrics", "validation_mechanism.csv"))
write_csv(hp_mech,      file.path(DOS_MECH, "metrics", "hyperparameters_mechanism.csv"))

# The row cap must never have bound for these two arms; if it did, the masks
# are not identical to an unconstrained draw and it must be reported.
cat(sprintf("\nRow constraint hits during this run: %d (expected 0)\n",
            N_CONSTRAINT_HITS))


## -------------------------------------------------------------------------
## 2. WAS THE MECHANISM ACTUALLY INSTALLED? (positive control)
## -------------------------------------------------------------------------
cat("\n=== POSITIVE CONTROL OF THE MECHANISM ===\n")
ctrl <- valid_mech |> group_by(mechanism, gap_length) |>
  summarise(pct_tests_detect  = 100 * mean(detected),
            pct_reject_season = 100 * mean(p_season < 0.05),
            or_season         = mean(or_season),
            median_p_ks       = median(p_values_ks),
            mean_ratio        = mean(mean_ratio),
            cohen_d           = mean(cohen_d),
            pct_in_rainy      = mean(pct_missing_in_rainy), .groups = "drop")
print(as.data.frame(ctrl), row.names = FALSE, digits = 3)
cat("\nReading: mean_ratio > 1 means the masked cells are higher than the\n",
    "observed ones. pct_in_rainy is to be compared with the share of rainy\n",
    "season weeks in the calendar (",
    sprintf("%.1f %%", 100 * mean(PERIOD_SEASONAL$rainy_season)), ").\n", sep = "")
if (any(ctrl$pct_reject_season < 50))
  cat("\nWARNING: mechanism barely detectable -> increase FORCE_SEASON.\n")


## -------------------------------------------------------------------------
## 3. DOES THE RANKING SURVIVE?  (the real question)
## -------------------------------------------------------------------------
#  We reload the MCAR run and restrict it to EXACTLY the same perimeter:
#  10% rate, gaps 1 and 8, repetitions 1 to 10. A paired comparison.
mcar <- read_csv("output/metrics/metrics_long.csv", show_col_types = FALSE) |>
  filter(missing_rate == CFG_MECH$missing_rate,
         gap_length %in% CFG_MECH$gap_lengths,
         rep <= CFG_MECH$n_rep) |>
  mutate(mechanism = "MCAR")
stopifnot(nrow(mcar) > 0)

all_arms <- bind_rows(mcar |> select(any_of(names(metrics_mech))), metrics_mech)

#  with_ties = FALSE is essential: Mean ignores the period covariates, so its
#  three variants are strictly identical (see 04_final_analysis.R).
best_variant <- all_arms |>
  group_by(mechanism, gap_length, method, variant) |>
  summarise(m = mean(MAE), .groups = "drop") |>
  group_by(mechanism, gap_length, method) |>
  slice_min(m, n = 1, with_ties = FALSE) |> ungroup() |>
  select(mechanism, gap_length, method, variant)

perf <- all_arms |> inner_join(best_variant,
                               by = c("mechanism", "gap_length", "method", "variant"))

ranks <- perf |>
  group_by(mechanism, gap_length, rep) |>
  mutate(rank = rank(MAE)) |> ungroup() |>
  group_by(mechanism, gap_length, method) |>
  summarise(mean_rank = mean(rank), MAE = mean(MAE),
            youden = mean(youden), sensitivity = mean(sensitivity),
            bias = mean(bias), .groups = "drop")

cat("\n=== MEAN RANK BY MECHANISM (best variant of each method) ===\n")
rank_table <- ranks |> select(mechanism, gap_length, method, mean_rank) |>
  pivot_wider(names_from = c(mechanism, gap_length), values_from = mean_rank)
ref_col <- if ("MCAR_1" %in% names(rank_table)) "MCAR_1" else names(rank_table)[2]
rank_table <- rank_table |> arrange(.data[[ref_col]])
print(as.data.frame(rank_table), row.names = FALSE, digits = 3)
write_csv(rank_table, file.path(DOS_MECH, "metrics", "ranks_by_mechanism.csv"))

cat("\n=== AGREEMENT OF THE RANKINGS (Spearman on mean ranks) ===\n")
#  THIS is the number to put in the paper. rho close to 1 means the ranking
#  does not depend on the mechanism -> the conclusions are robust.
pairs <- expand.grid(a = unique(ranks$mechanism), b = unique(ranks$mechanism),
                     gap = CFG_MECH$gap_lengths, stringsAsFactors = FALSE) |>
  filter(a < b)
for (i in seq_len(nrow(pairs))) {
  ra <- ranks |> filter(mechanism == pairs$a[i], gap_length == pairs$gap[i]) |>
    arrange(method) |> pull(mean_rank)
  rb <- ranks |> filter(mechanism == pairs$b[i], gap_length == pairs$gap[i]) |>
    arrange(method) |> pull(mean_rank)
  ct <- suppressWarnings(cor.test(ra, rb, method = "spearman"))
  cat(sprintf("  gap %d: %-5s vs %-5s   rho = %+.3f   p = %.4f\n",
              pairs$gap[i], pairs$a[i], pairs$b[i], ct$estimate, ct$p.value))
}


## -------------------------------------------------------------------------
## 4. DOES THE GAP-LENGTH REVERSAL SURVIVE?
## -------------------------------------------------------------------------
#  The central contribution. If Kalman degrades and MissForest improves as
#  gaps grow from 1 to 8 weeks UNDER ALL MECHANISMS, then the result concerns
#  the STRUCTURE of missingness and not its mechanism: precisely the axis
#  orthogonal to Rubin's taxonomy that this paper argues for.
cat("\n=== REVERSAL WITH GAP LENGTH (rank at gap 1 -> rank at gap 8) ===\n")
inv <- ranks |> select(mechanism, gap_length, method, mean_rank) |>
  pivot_wider(names_from = gap_length, values_from = mean_rank,
              names_prefix = "gap") |>
  mutate(delta = gap8 - gap1) |> arrange(mechanism, gap1)
print(as.data.frame(inv), row.names = FALSE, digits = 3)
write_csv(inv, file.path(DOS_MECH, "metrics", "reversal_by_mechanism.csv"))

cat("\nConsistency of the reversal direction across mechanisms:\n")
d_wide <- inv |> select(mechanism, method, delta) |>
  pivot_wider(names_from = mechanism, values_from = delta)
for (m2 in setdiff(names(d_wide), c("method", "MCAR"))) {
  ct <- suppressWarnings(cor.test(d_wide$MCAR, d_wide[[m2]], method = "spearman"))
  cat(sprintf("  MCAR vs %-5s: rho(delta) = %+.3f   p = %.4f\n",
              m2, ct$estimate, ct$p.value))
}


## -------------------------------------------------------------------------
## 5. FRIEDMAN UNDER EACH MECHANISM
## -------------------------------------------------------------------------
#  CAVEAT TO CARRY INTO THE PAPER: with N = 10 repetitions the Nemenyi
#  critical difference is about 4.28 ranks, against 1.75 in the main analysis
#  (N = 180). This section measures RANK AGREEMENT (Spearman above), not the
#  significance of pairwise differences. Do not present these Friedman
#  outputs as established rankings.
friedman_analysis <- function(d, title) {
  Mx <- d |> select(rep, method, MAE) |>
    pivot_wider(names_from = method, values_from = MAE) |>
    select(-rep) |> as.matrix()
  stopifnot(!anyNA(Mx))
  ft <- friedman.test(Mx)
  k <- ncol(Mx); N <- nrow(Mx)
  CD <- qtukey(0.95, k, Inf) / sqrt(2) * sqrt(k * (k + 1) / (6 * N))
  cat(sprintf("\n--- %s ---\n  chi2 = %.1f, df = %d, p = %.3g | k = %d, N = %d, CD = %.2f\n",
              title, ft$statistic, ft$parameter, ft$p.value, k, N, CD))
  rg <- sort(colMeans(t(apply(Mx, 1, rank))))
  for (i in seq_along(rg)) cat(sprintf("    %-20s %.2f\n", names(rg)[i], rg[i]))
  invisible(list(test = ft, CD = CD, ranks = rg))
}
for (mech in unique(perf$mechanism)) for (g in CFG_MECH$gap_lengths)
  friedman_analysis(perf |> filter(mechanism == mech, gap_length == g),
                    sprintf("Friedman — %s, %d-week gaps", mech, g))


## -------------------------------------------------------------------------
## 6. WHAT THE MECHANISM REALLY CHANGES: bias and peak detection
## -------------------------------------------------------------------------
#  The ranking can be robust while the BIAS is not at all. Under MAR/MNAR the
#  masked cells are the highest ones, so any method that pulls towards the
#  mean systematically underestimates the peaks. This is the strongest public
#  health argument in the paper.
cat("\n=== BIAS AND PEAK DETECTION BY MECHANISM (isolated gaps) ===\n")
bias_table <- ranks |> filter(gap_length == 1) |>
  select(mechanism, method, bias, sensitivity, youden) |>
  pivot_wider(names_from = mechanism,
              values_from = c(bias, sensitivity, youden))
print(as.data.frame(bias_table), row.names = FALSE, digits = 3)
write_csv(bias_table, file.path(DOS_MECH, "metrics", "bias_by_mechanism.csv"))


## -------------------------------------------------------------------------
## 7. FIGURE
## -------------------------------------------------------------------------
p <- ranks |>
  mutate(gap_lab = factor(paste0(gap_length, "-week gaps"),
                          levels = paste0(CFG_MECH$gap_lengths, "-week gaps")),
         method = reorder(method, mean_rank)) |>
  ggplot(aes(mean_rank, method, colour = mechanism, group = mechanism)) +
  geom_line(aes(group = method), colour = "grey80", linewidth = 3, alpha = .5) +
  geom_point(size = 3) +
  facet_wrap(~ gap_lab) +
  scale_colour_manual(values = c(MCAR = "#1b9e77", MAR = "#d95f02",
                                 MNAR = "#7570b3"), name = "Mechanism") +
  labs(x = "Mean rank (1 = best)", y = NULL,
       title = "Does the ranking of methods withstand the missingness mechanism?",
       subtitle = sprintf("%d%% missing, %d repetitions, best variant of each method",
                          CFG_MECH$missing_rate, CFG_MECH$n_rep)) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(DOS_MECH, "figures", "fig9_mechanism_sensitivity.png"), p,
       width = 10, height = 6, dpi = 300)

cat("\n=== Done ===\n")
cat("Outputs in ", DOS_MECH, "/ : ranks_by_mechanism.csv, reversal_by_mechanism.csv,\n",
    "bias_by_mechanism.csv, validation_mechanism.csv, fig9_mechanism_sensitivity.png\n",
    sep = "")
cat("Next: source(\"R/06_mechanism_diagnostics.R\")\n")
