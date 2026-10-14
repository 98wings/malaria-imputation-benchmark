###############################################################################
#  11 — SUPPLEMENTARY STATISTICS ADDED DURING THE 2018-2025 REVISION
#
#  Run AFTER 07 and the second run of 06. Recomputes NO imputation: it
#  re-reads result files, re-simulates masks and (section 10) refits the
#  forecasting models of 03 on archived imputed datasets.
#  Runtime: sections 1-9 a few minutes (bootstrap of section 8); section 10
#  about 10-20 minutes (set RUN_A3 <- FALSE in section 10 to skip it).
#
#  WHY THIS SCRIPT EXISTS
#  Three numbers cited in the paper were either printed to the console only
#  (07) or never computed by the pipeline. They are now written to files so
#  that every number in the text can be traced to a CSV.
#
#  Produces
#    output_mechanism/metrics/mar_strong_calibration.csv
#        Calibration curve of the strong-MAR arm (beta -> Cohen's d),
#        recomputed EXACTLY as in section 3 of 07: same seeds, same grid,
#        same row constraint, same retention rule. 07 printed it only.
#    output_mechanism/metrics/mar_strong_contamination.csv
#        Share of neighbouring values unavailable at conditioning time and
#        correlation between true and observable neighbour means (07,
#        section 4, same seed). Also console-only until now.
#    output_mechanism/metrics/nemenyi_cd_reduced.csv
#        Nemenyi critical difference of the reduced mechanism design
#        (k methods, N repetitions per mechanism and gap length).
#    output/metrics/paired_gap8.csv
#        Paired comparison with 8-week gaps, over the 60 blocks
#        (repetition x missingness rate): MissForest-ST (SeasonalTrend)
#        against Kalman smoothing (NoTime).
#    output/metrics/achieved_precision.csv
#        Achieved precision of the design (Methods): half-width of the 95%
#        CI of the paired MAE difference MissForest - linear interpolation,
#        isolated gaps, NoTime, 60 blocks. Same computation as 03, section 2,
#        which printed it only.
#    output/metrics/mcar_uniformity_tests.csv
#        Binomial test on the number of failures and KS uniformity tests of
#        the p-values, 1-week gaps (03, section 1, console only until now).
#    output/metrics/friedman_tests.csv
#        Friedman chi-squared, df, p and Nemenyi CD of the three analyses
#        (03, section 6: main, NoTime, RMSE).
#    output/metrics/district_spearman_tests.csv
#        The eight exploratory Spearman tests at district level (04).
#    output_mechanism/metrics/mechanism_agreement_tests.csv
#        Spearman agreement of the mean rankings between mechanisms, per gap
#        length, and consistency of the reversal direction (07).
#
#  Added after the devil's-advocate review (step 8, October 2026):
#    output/metrics/friedman_by_gap.csv, friedman_ranks_by_gap.csv,
#    nemenyi_pairs_by_gap.csv                                   (7, A2)
#        Friedman and Nemenyi tests within each gap length (1, 4, 8 weeks).
#    output/metrics/rate_vs_gap_agreement.csv, rate_vs_gap_top3.csv,
#    rate_vs_gap_interaction.csv                                (8, A1)
#        Does gap length change the recommended method more than the rate?
#        Three criteria, all reported (see the section header).
#    output/metrics/mcar_pvalue_null.csv                       (9, A8)
#        Null distribution of the MCAR-validation p-values under the exact
#        mask generator (fixed-size sampling without replacement).
#    output/metrics/forecast_origin_errors.csv, forecast_leakage_diagnostic.csv,
#    forecast_leakage_by_rate.csv, forecast_leakage_summary.csv (10, A3 b)
#        Leakage diagnostic of the 4-week forecast of 03, section 8.
#    output/metrics/gap_length_distribution.csv, gap_length_distribution_long.csv
#    output_mechanism/metrics/mechanism_ranks_long.csv           (11, A6)
#        Effective gap lengths, block overlaps, realised vs nominal rate and
#        last masked week of the 180 masks; mean rank, MAE and bias by
#        mechanism and gap length (07, section 6, recomputed and checked).
#    output/metrics/kalman_fallback_check.csv                    (12)
#        How often the first level of the Kalman cascade (StructTS) fails on
#        the archived masks of the main benchmark and the mechanism arms.
###############################################################################

source("R/functions_mechanism.R")
MET  <- file.path("output", "metrics")
MMET <- file.path(DOS_MECH, "metrics")

## -------------------------------------------------------------------------
## 1. STRONG-MAR CALIBRATION CURVE (reproduces 07, section 3)
## -------------------------------------------------------------------------
#  Settings copied from 07: mechanism, row constraint (at least 3 districts
#  observed per week), target, grid, number of masks per beta, seeds, and
#  tolerance. If any of them changes in 07, change it here too: the check at
#  the end of this section compares the retained beta with the one stored by
#  07 in validation_mar_strong.csv.
MECH <- "MAR_strong"
MAX_MISSING_PER_WEEK <- n_dis - 3L
TARGET <- read_csv(file.path(MMET, "mechanism_diagnostics.csv"),
                   show_col_types = FALSE) |>
  filter(mechanism == "MNAR", gap_length == 1) |> pull(cohen_d) |> mean()

observed_d <- function(beta, n_tries = 6) {
  FORCE_NEIGHBOUR <<- beta
  Z <- scale(REF)
  mean(sapply(seq_len(n_tries), function(b) {
    M  <- simulate_mask_mechanism(CFG_MECH$missing_rate, 700000L + b,
                                  gap_length = 1, mechanism = MECH)
    Mv <- M & !EXCLUDE
    (mean(Z[Mv]) - mean(Z[!Mv])) / sd(Z[!Mv])
  }))
}

grid <- seq(0.5, 6.0, by = 0.25)
TOL  <- 0.03
ds   <- sapply(grid, observed_d)
within <- which(abs(ds - TARGET) <= TOL)
retained_beta <- if (length(within)) grid[min(within)] else grid[which.max(ds)]

calib <- tibble(beta = grid, cohen_d = ds, target_mnar = TARGET,
                within_tolerance = abs(ds - TARGET) <= TOL,
                retained = grid == retained_beta,
                target_reached = length(within) > 0)
write_csv(calib, file.path(MMET, "mar_strong_calibration.csv"))
cat(sprintf("Calibration: max d = %.3f at beta = %.2f; MNAR target = %.3f; %s\n",
            max(ds), grid[which.max(ds)], TARGET,
            if (length(within)) "target reached" else "TARGET NOT REACHABLE"))

stored_beta <- read_csv(file.path(MMET, "validation_mar_strong.csv"),
                        show_col_types = FALSE) |> pull(force_neighbour) |> unique()
if (!isTRUE(all.equal(stored_beta, retained_beta)))
  stop(sprintf("Retained beta %.2f differs from the beta used by 07 (%s): ",
               retained_beta, paste(stored_beta, collapse = ", ")),
       "the settings of this section no longer match 07.")
cat(sprintf("Check passed: retained beta %.2f is the one used by 07.\n", retained_beta))

## -------------------------------------------------------------------------
## 2. CONTAMINATION OF THE CONDITIONING VARIABLE (reproduces 07, section 4)
## -------------------------------------------------------------------------
FORCE_NEIGHBOUR <- retained_beta
M_test <- simulate_mask_mechanism(CFG_MECH$missing_rate, 700099L,
                                  gap_length = 1, mechanism = MECH)
X_test <- REF; X_test[M_test] <- NA
NEIGH_obs <- t(sapply(seq_len(n_obs), function(i)
  sapply(seq_len(n_dis), function(j) mean(X_test[i, -j], na.rm = TRUE))))
contam <- tibble(
  beta = retained_beta, seed = 700099L,
  cor_true_vs_observable_neighbours =
    cor(as.numeric(NEIGH), as.numeric(NEIGH_obs), use = "complete.obs"),
  # same quantity as printed by 07 ("share of unavailable neighbours per cell")
  pct_unavailable_neighbours = 100 * mean(is.na(X_test[, -1])))
write_csv(contam, file.path(MMET, "mar_strong_contamination.csv"))
print(contam)

## -------------------------------------------------------------------------
## 3. NEMENYI CRITICAL DIFFERENCE OF THE REDUCED MECHANISM DESIGN
## -------------------------------------------------------------------------
#  Same formula as the main analysis (04): qtukey(0.95, k, Inf) / sqrt(2) *
#  sqrt(k (k + 1) / (6 N)), with N = number of repetitions (blocks) of the
#  reduced design, read from the metric files rather than assumed.
mm <- bind_rows(
  read_csv(file.path(MMET, "metrics_mechanism.csv"),  show_col_types = FALSE),
  read_csv(file.path(MMET, "metrics_mar_strong.csv"), show_col_types = FALSE))
k_meth <- read_csv(file.path(MMET, "ranks_four_mechanisms.csv"),
                   show_col_types = FALSE) |> nrow()
cd_red <- mm |>
  distinct(mechanism, gap_length, rep) |>
  count(mechanism, gap_length, name = "n_blocks") |>
  mutate(k_methods = k_meth,
         cd = qtukey(0.95, k_methods, Inf) / sqrt(2) *
              sqrt(k_methods * (k_methods + 1) / (6 * n_blocks)))
write_csv(cd_red, file.path(MMET, "nemenyi_cd_reduced.csv"))
print(cd_red)

## -------------------------------------------------------------------------
## 4. PAIRED COMPARISON, 8-WEEK GAPS: MissForest-ST vs Kalman
## -------------------------------------------------------------------------
metrics <- read_csv(file.path(MET, "metrics_long.csv"), show_col_types = FALSE)
pair <- metrics |>
  filter(gap_length == 8,
         (method == "MissForest_ST" & variant == "SeasonalTrend") |
         (method == "Kalman"        & variant == "NoTime")) |>
  select(rep, missing_rate, method, MAE) |>
  pivot_wider(names_from = method, values_from = MAE) |>
  mutate(diff = MissForest_ST - Kalman)
stopifnot(nrow(pair) == 60, !anyNA(pair$diff))
wt <- wilcox.test(pair$MissForest_ST, pair$Kalman, paired = TRUE)
bt <- binom.test(sum(pair$diff < 0), nrow(pair), p = 0.5)
paired <- tibble(
  gap_length = 8, comparison = "MissForest_ST SeasonalTrend vs Kalman NoTime",
  n_blocks = nrow(pair),
  mae_missforest_st = mean(pair$MissForest_ST), mae_kalman = mean(pair$Kalman),
  mean_diff = mean(pair$diff),
  pct_diff_vs_kalman = 100 * mean(pair$diff) / mean(pair$Kalman),
  blocks_won_by_missforest_st = sum(pair$diff < 0),
  wilcoxon_V = unname(wt$statistic), wilcoxon_p = wt$p.value,
  sign_test_p = bt$p.value)
write_csv(paired, file.path(MET, "paired_gap8.csv"))
print(paired, width = Inf)

## -------------------------------------------------------------------------
## 5. ACHIEVED PRECISION OF THE DESIGN (reproduces 03, section 2)
## -------------------------------------------------------------------------
achieved_precision <- metrics |>
  filter(gap_length == 1, variant == "NoTime",
         method %in% c("LinearInterp", "MissForest")) |>
  select(rep, missing_rate, method, MAE) |>
  pivot_wider(names_from = method, values_from = MAE) |>
  mutate(diff = MissForest - LinearInterp) |>
  summarise(gap_length = 1, variant = "NoTime",
            comparison = "MissForest - LinearInterp",
            mean_difference = mean(diff), sd_diff = sd(diff),
            half_ci95 = 1.96 * sd(diff) / sqrt(n()), n = n())
stopifnot(achieved_precision$n == 60)
write_csv(achieved_precision, file.path(MET, "achieved_precision.csv"))
print(achieved_precision, width = Inf)

## -------------------------------------------------------------------------
## 6. STATISTICS PRINTED ONLY TO THE CONSOLE BY 03, 04 AND 07
## -------------------------------------------------------------------------
#  Recomputed from the result files with the SAME calls as the original
#  scripts, so that every number cited in the paper can be traced to a CSV.

# 6a. MCAR validation (03, section 1)
val1    <- read_csv(file.path(MET, "mcar_validation.csv"), show_col_types = FALSE) |>
  filter(gap_length == 1)
n_tests <- 2 * nrow(val1)
n_fail  <- sum(val1$p_time <= 0.05) + sum(val1$p_values_ks <= 0.05)
bt1     <- binom.test(n_fail, n_tests, p = 0.05, alternative = "greater")
mcar_tests <- bind_rows(
  tibble(test = "binomial_failures", variable = "p_time + p_values_ks",
         n = n_tests, n_fail = n_fail, expected = 0.05 * n_tests,
         statistic = NA_real_, p_value = bt1$p.value, mean_p = NA_real_),
  lapply(c("p_districts", "p_time", "p_values_ks"), function(cc) {
    kt <- suppressWarnings(ks.test(val1[[cc]], "punif"))
    tibble(test = "ks_uniformity", variable = cc, n = nrow(val1),
           n_fail = NA_integer_, expected = NA_real_,
           statistic = unname(kt$statistic), p_value = kt$p.value,
           mean_p = mean(val1[[cc]]))
  }) |> bind_rows())
write_csv(mcar_tests, file.path(MET, "mcar_uniformity_tests.csv"))
print(mcar_tests, width = Inf)

# 6b. Friedman tests (03, section 6)
sel <- read_csv(file.path(MET, "selected_variant.csv"), show_col_types = FALSE)
g1  <- metrics |> filter(gap_length == 1)
bm  <- function(df, column) {
  df |> unite("blk", rep, missing_rate) |>
    select(blk, method, all_of(column)) |>
    pivot_wider(names_from = method, values_from = all_of(column)) |>
    column_to_rownames("blk") |> as.matrix()
}
fr_one <- function(mat, analysis) {
  ft <- friedman.test(mat); k <- ncol(mat); N <- nrow(mat)
  tibble(analysis = analysis, n_blocks = N, k_methods = k,
         chi_squared = unname(ft$statistic), df = unname(ft$parameter),
         p_value = ft$p.value,
         cd = qtukey(0.95, k, Inf) / sqrt(2) * sqrt(k * (k + 1) / (6 * N)))
}
fried <- bind_rows(
  fr_one(bm(inner_join(g1, sel, by = c("method", "variant")), "MAE"), "main_best_variant_MAE"),
  fr_one(bm(filter(g1, variant == "NoTime"), "MAE"), "NoTime_MAE"),
  fr_one(bm(inner_join(g1, sel, by = c("method", "variant")), "RMSE"), "best_variant_RMSE"))
write_csv(fried, file.path(MET, "friedman_tests.csv"))
print(fried, width = Inf)

# 6c. District-level Spearman tests (04)
gbd   <- read_csv(file.path(MET, "gain_by_district.csv"), show_col_types = FALSE)
CHARS <- c("cv", "mean_corr_other_districts", "seasonal_r2", "autocorr_lag1")
dist_tests <- bind_rows(lapply(c("MissForest", "gain_MissForest_pct"), function(y)
  bind_rows(lapply(CHARS, function(v) {
    ct <- suppressWarnings(cor.test(gbd[[y]], gbd[[v]], method = "spearman"))
    tibble(outcome = y, characteristic = v, n = nrow(gbd),
           rho = unname(ct$estimate), p_value = ct$p.value,
           bonferroni_threshold = 0.05 / 8)
  }))))
stopifnot(nrow(dist_tests) == 8)
write_csv(dist_tests, file.path(MET, "district_spearman_tests.csv"))
print(dist_tests, width = Inf)

# 6d. Agreement between mechanisms (07)
r4   <- read_csv(file.path(MMET, "ranks_four_mechanisms.csv"), show_col_types = FALSE) |>
  arrange(method)
rv4  <- read_csv(file.path(MMET, "reversal_four_mechanisms.csv"), show_col_types = FALSE)
ARMS <- c("MCAR", "MAR", "MAR_strong", "MNAR")
agree <- list()
for (gap in c(1, 8)) for (i in 1:3) for (j in (i + 1):4) {
  a <- paste0(ARMS[i], "_", gap); b <- paste0(ARMS[j], "_", gap)
  ct <- suppressWarnings(cor.test(r4[[a]], r4[[b]], method = "spearman"))
  agree[[length(agree) + 1]] <- tibble(
    analysis = "ranking_agreement", gap_length = gap,
    mechanism_a = ARMS[i], mechanism_b = ARMS[j],
    rho = unname(ct$estimate), p_value = ct$p.value)
}
dw <- rv4 |> select(mechanism, method, delta) |>
  pivot_wider(names_from = mechanism, values_from = delta)
for (m2 in setdiff(ARMS, "MCAR")) {
  ct <- suppressWarnings(cor.test(dw$MCAR, dw[[m2]], method = "spearman"))
  agree[[length(agree) + 1]] <- tibble(
    analysis = "reversal_direction", gap_length = NA_real_,
    mechanism_a = "MCAR", mechanism_b = m2,
    rho = unname(ct$estimate), p_value = ct$p.value)
}
agree <- bind_rows(agree)
write_csv(agree, file.path(MMET, "mechanism_agreement_tests.csv"))
print(agree, n = 20, width = Inf)


## =========================================================================
## SECTIONS 7-10: ADDED AFTER THE DEVIL'S-ADVOCATE REVIEW (step 8, Oct 2026)
## =========================================================================
#  Report: claude/rapport_etape8.md (analyses A1, A2, A3 (b), A8).
#  Decisions: claude/decisions_etape8.md.
#  NO imputation is recomputed. Sections 7-9 re-read result files and
#  re-simulate masks (seconds). Section 10 refits the forecasting models of
#  03 (section 8) on the archived 1-week imputed datasets; it is the only slow
#  part (see RUN_A3 below).

library(PMCMRplus)   # already loaded by functions.R; explicit for clarity

sel <- read_csv(file.path(MET, "selected_variant.csv"), show_col_types = FALSE)

## -------------------------------------------------------------------------
## 7. A2 - FRIEDMAN AND NEMENYI TESTS WITHIN EACH GAP LENGTH
## -------------------------------------------------------------------------
#  Methods (section 2.10) state that the tests are run WITHIN each gap length,
#  but 03 only ran them for 1-week gaps. Same design as 03, section 6 (A):
#  10 methods, variant retained in selected_variant.csv (chosen on 1-week
#  gaps), blocks = repetition x missingness rate (60 blocks), MAE criterion.
#  As a check, the variant that would be best at each gap length is also
#  recorded; if it differs, the column `variant_best_at_gap` shows it.

fr_matrix <- function(df) {
  df |> unite("blk", rep, missing_rate) |>
    select(blk, method, MAE) |>
    pivot_wider(names_from = method, values_from = MAE) |>
    column_to_rownames("blk") |> as.matrix()
}

fried_gap <- list(); ranks_gap <- list(); pairs_gap <- list()
for (gap in CFG$gap_lengths) {
  dg   <- metrics |> filter(gap_length == gap)
  mat  <- fr_matrix(inner_join(dg, sel, by = c("method", "variant")))
  k    <- ncol(mat); N <- nrow(mat)
  cd   <- qtukey(0.95, k, Inf) / sqrt(2) * sqrt(k * (k + 1) / (6 * N))
  ft   <- friedman.test(mat)
  fried_gap[[length(fried_gap) + 1]] <- tibble(
    gap_length = gap, n_blocks = N, k_methods = k,
    chi_squared = unname(ft$statistic), df = unname(ft$parameter),
    p_value = ft$p.value, cd = cd)

  best_here <- dg |> group_by(method, variant) |>
    summarise(m = mean(MAE), .groups = "drop") |>
    group_by(method) |> slice_min(m, n = 1, with_ties = FALSE) |>
    ungroup() |> select(method, variant_best_at_gap = variant)
  mr <- colMeans(t(apply(mat, 1, rank)))
  ranks_gap[[length(ranks_gap) + 1]] <- tibble(
    gap_length = gap, method = names(mr), mean_rank = unname(mr)) |>
    left_join(sel, by = "method") |>
    left_join(best_here, by = "method") |>
    arrange(mean_rank)

  nt <- frdAllPairsNemenyiTest(mat)
  pv <- nt$p.value
  for (a in rownames(pv)) for (b in colnames(pv)) if (!is.na(pv[a, b]))
    pairs_gap[[length(pairs_gap) + 1]] <- tibble(
      gap_length = gap, method_a = a, method_b = b,
      rank_diff = abs(mr[[a]] - mr[[b]]), p_nemenyi = pv[a, b],
      exceeds_cd = abs(mr[[a]] - mr[[b]]) > cd)
}
fried_gap <- bind_rows(fried_gap)
ranks_gap <- bind_rows(ranks_gap)
pairs_gap <- bind_rows(pairs_gap)
# Consistency with 03 for 1-week gaps (same ranks as friedman_ranks_main.csv)
chk <- read_csv(file.path(MET, "friedman_ranks_main.csv"), show_col_types = FALSE) |>
  inner_join(filter(ranks_gap, gap_length == 1), by = "method")
stopifnot(nrow(chk) == 10, max(abs(chk$mean_rank.x - chk$mean_rank.y)) < 1e-9)
write_csv(fried_gap, file.path(MET, "friedman_by_gap.csv"))
write_csv(ranks_gap, file.path(MET, "friedman_ranks_by_gap.csv"))
write_csv(pairs_gap, file.path(MET, "nemenyi_pairs_by_gap.csv"))
print(fried_gap, width = Inf)
print(ranks_gap, n = 30, width = Inf)


## -------------------------------------------------------------------------
## 8. A1 - DOES GAP LENGTH CHANGE THE RECOMMENDATION MORE THAN THE RATE?
## -------------------------------------------------------------------------
#  WHY THE PRIMARY CRITERION IS THE IDENTITY OF THE BEST METHODS
#  The question the paper asks is practical: which method should a programme
#  use? The answer to that question is the identity of the best method(s),
#  so that is the primary criterion. It was chosen because it follows from
#  the question, NOT because it supports the title. The two other criteria
#  are reported as well, including the one that does not favour gap length:
#    C1  global agreement of the ranking (Spearman, 24 combinations);
#    C2  identity of the top places:
#        C2a  overlap of the top 3 combinations (mean MAE over repetitions);
#        C2b  probability that the best method differs between two
#             conditions (repetition level), with a bootstrap 95% CI;
#    C3  size of the method x factor interaction on log(MAE).
#  If C2 and C3 do not show a larger effect of gap length than of rate, the
#  title must be reformulated (decision of the first author, 12/10/2026).
#
#  COMPARISONS: gap lengths 1-4, 1-8, 4-8 at each rate; rates 5-10, 5-20,
#  10-20 at each gap length. Note that the two factors span different ranges
#  (x4 for the rate, x8 for the length): any conclusion holds for the ranges
#  studied only.
#
#  UNIT OF ANALYSIS
#  C1, C2a, C3: the 24 method x predictor-block combinations (as Table 2).
#  C2b: (i) the 10 methods with their retained variant (the recommendation a
#  programme would follow) and (ii) the 24 combinations.
#
#  C2b, DEFINITION. Repetitions are NOT paired across conditions (each
#  condition has its own masks and seeds), so pairing repetition r of one
#  condition with repetition r of another would be arbitrary. We therefore
#  use the probability that the winner differs for one repetition drawn at
#  random in each condition:  p_cross = 1 - sum_i pA_i * pB_i,  where pA_i is
#  the share of repetitions of condition A won by method i. Because winners
#  also vary between repetitions of the SAME condition, p_cross is compared
#  with the within-condition disagreement (two distinct repetitions of A):
#  w_A = 1 - sum_i cA_i (cA_i - 1) / (n (n - 1)). The quantity of interest is
#  excess = p_cross - (w_A + w_B) / 2  (0 = no more change than noise).
#  For transparency the naive paired version (same repetition index) is also
#  given, with an exact (Clopper-Pearson) CI.
#
#  C3, NORMALISATION (fixed before computing; not tuned to match any value):
#  y[c, r, g] = mean over the 20 repetitions of log(MAE) for combination c,
#  rate r and gap g (24 x 3 x 3 table, one value per cell). For a factor F
#  (rate or gap), with O the other factor (3 levels):
#     SS(method x F) = 3 * sum_{c,f} (ybar[c,f] - ybar[c,.] - ybar[.,f] + ybar)^2
#  where ybar[c,f] averages y over O. Both interactions have 23 x 2 = 46 df,
#  so their ratio is directly interpretable. Each SS is also given as a share
#  of the total SS of the table.
#
#  UNCERTAINTY: nonparametric bootstrap of the repetitions, B = 2000. Within
#  each (rate, gap) cell, 20 repetition indices are drawn with replacement,
#  the SAME draw for all methods of that cell (a mask is shared by all
#  methods); cells are resampled independently. Percentile 95% intervals.

B_BOOT <- 2000
RATES  <- CFG$missing_rates; GAPS <- CFG$gap_lengths; NREP <- CFG$n_rep

combo_df <- metrics |> distinct(method, variant) |>
  arrange(method, variant) |> mutate(combo = paste(method, variant, sep = "|"))
COMBOS   <- combo_df$combo
METHODS  <- sort(unique(metrics$method))
sel_key  <- paste(sel$method, sel$variant, sep = "|")

# A[rep, combo, rate, gap] = MAE
A <- array(NA_real_, c(NREP, length(COMBOS), length(RATES), length(GAPS)),
           dimnames = list(NULL, COMBOS, as.character(RATES), as.character(GAPS)))
mm <- metrics |> mutate(combo = paste(method, variant, sep = "|"))
A[cbind(mm$rep, match(mm$combo, COMBOS), match(mm$missing_rate, RATES),
        match(mm$gap_length, GAPS))] <- mm$MAE
stopifnot(!anyNA(A))
meth_idx <- match(sel_key[match(METHODS, sel$method)], COMBOS)  # 10 columns

# --- statistics computed on one (possibly resampled) array ---------------
pairs_spec <- bind_rows(
  expand_grid(factor = "gap",  fixed = RATES,
              tibble(a = c(1, 1, 4), b = c(4, 8, 8))),
  expand_grid(factor = "rate", fixed = GAPS,
              tibble(a = c(5, 5, 10), b = c(10, 20, 20))))
cell <- function(fac, fixed, lvl) {          # returns c(rate, gap) labels
  if (fac == "gap") c(as.character(fixed), as.character(lvl))
  else              c(as.character(lvl),   as.character(fixed))
}

winners <- function(Mx) apply(Mx, 1, which.min)   # rows = repetitions
change_stats <- function(wa, wb, K) {
  pa <- tabulate(wa, K) / length(wa); pb <- tabulate(wb, K) / length(wb)
  n  <- length(wa)
  w_in <- function(w) { cnt <- tabulate(w, K); 1 - sum(cnt * (cnt - 1)) / (n * (n - 1)) }
  p_cross <- 1 - sum(pa * pb)
  c(p_cross = p_cross, within_a = w_in(wa), within_b = w_in(wb),
    excess = p_cross - (w_in(wa) + w_in(wb)) / 2)
}
ss_inter <- function(Y, along) {          # Y: combo x rate x gap
  other <- if (along == 3) 2 else 3
  Yf <- apply(Y, c(1, along), mean)      # combo x level of F
  r  <- sweep(sweep(Yf, 1, rowMeans(Yf)), 2, colMeans(Yf)) + mean(Yf)
  dim(Y)[other] * sum(r^2)
}

compute_all <- function(Ar) {
  mean_mae <- apply(Ar, c(2, 3, 4), mean)                     # combo x rate x gap
  out <- pmap_dfr(pairs_spec, function(factor, fixed, a, b) {
    ca <- cell(factor, fixed, a); cb <- cell(factor, fixed, b)
    xa <- mean_mae[, ca[1], ca[2]]; xb <- mean_mae[, cb[1], cb[2]]
    top_a <- names(sort(xa))[1:3]; top_b <- names(sort(xb))[1:3]
    wm_a <- winners(Ar[, meth_idx, ca[1], ca[2]]); wm_b <- winners(Ar[, meth_idx, cb[1], cb[2]])
    wc_a <- winners(Ar[, , ca[1], ca[2]]);         wc_b <- winners(Ar[, , cb[1], cb[2]])
    sm <- change_stats(wm_a, wm_b, length(meth_idx))
    sc <- change_stats(wc_a, wc_b, length(COMBOS))
    tibble(factor = factor, fixed_level = fixed, level_a = a, level_b = b,
           spearman_24 = cor(xa, xb, method = "spearman"),
           top3_overlap = length(intersect(top_a, top_b)),
           top1_same = top_a[1] == top_b[1],
           method_p_cross = sm[["p_cross"]], method_within_a = sm[["within_a"]],
           method_within_b = sm[["within_b"]], method_excess = sm[["excess"]],
           combo_p_cross = sc[["p_cross"]], combo_excess = sc[["excess"]])
  })
  Y <- apply(log(Ar), c(2, 3, 4), mean)
  ss_gap <- ss_inter(Y, 3); ss_rate <- ss_inter(Y, 2)
  attr(out, "inter") <- c(ss_gap = ss_gap, ss_rate = ss_rate,
                          ss_total = sum((Y - mean(Y))^2),
                          ratio = ss_gap / ss_rate)
  out
}

# --- point estimates -----------------------------------------------------
est <- compute_all(A)
inter_est <- attr(est, "inter")

# top-3 lists, for reading (point estimate only)
mean_mae_all <- apply(A, c(2, 3, 4), mean)
top3_lists <- expand_grid(missing_rate = RATES, gap_length = GAPS) |>
  mutate(top3 = map2_chr(missing_rate, gap_length, ~ paste(
    names(sort(mean_mae_all[, as.character(.x), as.character(.y)]))[1:3],
    collapse = " ; ")))

# --- naive paired version (same repetition index), exact CI ---------------
naive <- pmap_dfr(pairs_spec, function(factor, fixed, a, b) {
  ca <- cell(factor, fixed, a); cb <- cell(factor, fixed, b)
  ch <- winners(A[, meth_idx, ca[1], ca[2]]) != winners(A[, meth_idx, cb[1], cb[2]])
  bt <- binom.test(sum(ch), length(ch))
  tibble(factor = factor, fixed_level = fixed, level_a = a, level_b = b,
         method_changed_same_index = mean(ch),
         same_index_lo = bt$conf.int[1], same_index_hi = bt$conf.int[2])
})

# --- bootstrap -------------------------------------------------------------
set.seed(20261012)
boot <- vector("list", B_BOOT); boot_inter <- matrix(NA_real_, B_BOOT, 4)
for (bb in seq_len(B_BOOT)) {
  Ab <- A
  for (r in seq_along(RATES)) for (g in seq_along(GAPS))
    Ab[, , r, g] <- A[sample.int(NREP, NREP, replace = TRUE), , r, g]
  o <- compute_all(Ab)
  boot[[bb]] <- o; boot_inter[bb, ] <- attr(o, "inter")
}
boot <- bind_rows(boot, .id = "b")
q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
ci <- boot |>
  group_by(factor, fixed_level, level_a, level_b) |>
  summarise(across(c(spearman_24, method_p_cross, method_excess,
                     combo_p_cross, combo_excess),
                   list(lo = ~ q(.x, 0.025), hi = ~ q(.x, 0.975))),
            top1_same_boot_share = mean(top1_same), .groups = "drop")

rank_agreement <- est |>
  left_join(ci, by = c("factor", "fixed_level", "level_a", "level_b")) |>
  left_join(naive, by = c("factor", "fixed_level", "level_a", "level_b")) |>
  arrange(factor, fixed_level, level_a, level_b)

colnames(boot_inter) <- names(inter_est)
interaction_tbl <- tibble(
  term = c("method_x_gap", "method_x_rate", "ratio_gap_over_rate"),
  df = c(46, 46, NA),
  ss = c(inter_est[["ss_gap"]], inter_est[["ss_rate"]], inter_est[["ratio"]]),
  share_of_total_ss = c(inter_est[["ss_gap"]], inter_est[["ss_rate"]], NA) /
    inter_est[["ss_total"]],
  ci_lo = c(q(boot_inter[, "ss_gap"], .025), q(boot_inter[, "ss_rate"], .025),
            q(boot_inter[, "ratio"], .025)),
  ci_hi = c(q(boot_inter[, "ss_gap"], .975), q(boot_inter[, "ss_rate"], .975),
            q(boot_inter[, "ratio"], .975)),
  normalisation = "24x3x3 table of mean log(MAE) over 20 reps; SS = 3*sum(residual^2); bootstrap of reps within cells, B = 2000")

write_csv(rank_agreement,  file.path(MET, "rate_vs_gap_agreement.csv"))
write_csv(top3_lists,      file.path(MET, "rate_vs_gap_top3.csv"))
write_csv(interaction_tbl, file.path(MET, "rate_vs_gap_interaction.csv"))
print(rank_agreement |> select(factor, fixed_level, level_a, level_b, spearman_24,
                               top3_overlap, method_p_cross, method_excess,
                               method_excess_lo, method_excess_hi), n = 30, width = Inf)
print(interaction_tbl, width = Inf)


## -------------------------------------------------------------------------
## 9. A8 - REFERENCE DISTRIBUTION OF THE MCAR-VALIDATION P-VALUES
## -------------------------------------------------------------------------
#  The 1-week masks draw a FIXED number of cells WITHOUT replacement. The
#  chi-squared tests of 02 (uniformity across districts, and over years)
#  compare the counts with a multinomial reference, which is more variable;
#  their p-values are therefore not expected to be uniform. The text used to
#  attribute the departure to the row constraint, but that constraint is only
#  a stopifnot() check in simulate_mask() and never acts.
#  Here we simulate the null distribution of the two p-values under the exact
#  mask generator (simulate_mask, gap 1), with seeds disjoint from those of
#  the benchmark (benchmark: 1000 * rep + rate + 7 * gap < 21 000), and
#  compare the 60 observed p-values with it (two-sample KS test).

N_SIM <- 2000
p_two <- function(M) {                       # same formulas as validate_mcar()
  yr  <- lubridate::year(dates)
  p1  <- chisq.test(colSums(M))$p.value
  obs <- tapply(rowSums(M), yr, sum)
  ex  <- as.numeric(table(yr)) * n_dis; ex <- ex / sum(ex) * sum(obs)
  p2  <- chisq.test(as.numeric(obs), p = ex / sum(ex))$p.value
  c(p_districts = p1, p_time = p2)
}
sim_p <- map_dfr(RATES, function(rate)
  map_dfr(seq_len(N_SIM), function(i) {
    M <- simulate_mask(rate, seed = 5e6 + 1e4 * rate + i, gap_length = 1)
    as_tibble_row(p_two(M)) |> mutate(missing_rate = rate)
  }))
val_obs <- read_csv(file.path(MET, "mcar_validation.csv"), show_col_types = FALSE) |>
  filter(gap_length == 1)
mcar_null <- map_dfr(c("p_districts", "p_time"), function(v)
  map_dfr(c(as.list(RATES), list(RATES)), function(rr) {
    s <- sim_p |> filter(missing_rate %in% rr) |> pull(all_of(v))
    o <- val_obs |> filter(missing_rate %in% rr) |> pull(all_of(v))
    tibble(variable = v,
           missing_rate = if (length(rr) == 1) as.character(rr) else "all",
           n_sim = length(s), sim_mean_p = mean(s),
           sim_mean_p_se_for_n_obs = sd(s) / sqrt(length(o)),
           n_obs = length(o), obs_mean_p = mean(o),
           ks_obs_vs_sim_p = suppressWarnings(ks.test(o, s)$p.value),
           ks_obs_vs_uniform_p = suppressWarnings(ks.test(o, "punif")$p.value))
  }))
write_csv(mcar_null, file.path(MET, "mcar_pvalue_null.csv"))
print(mcar_null, width = Inf)


## -------------------------------------------------------------------------
## 10. A3 (b) - LEAKAGE DIAGNOSTIC OF THE 4-WEEK FORECAST
## -------------------------------------------------------------------------
#  03 (section 8) forecasts 4 weeks ahead from rolling origins t in the last
#  year, using the IMPUTED history 1..t. Masks also cover the last year, and
#  the bilateral imputers (Kalman smoothing, linear interpolation, the lead
#  terms of MissForest-ST) used weeks t+1..t+4, i.e. the forecast targets,
#  to impute masked weeks just before t. Diagnostic: split the origins of
#  each (dataset, district) by whether any of the LOOKBACK weeks ending at t
#  was masked, and recompute the degradation (forecast RMSE on the imputed
#  history minus RMSE on the complete reference, same origins) separately.
#  "clean" origins = none of the last LOOKBACK weeks masked: the targets
#  cannot have entered the imputation of the weeks the forecast starts from.
#  Residual (global) leakage remains and is not removed here: mean
#  imputation, seasonal decomposition and missRanger were fitted on the whole
#  series, last year included.
#  The forecasting code is IDENTICAL to 03 (same orders, same Fourier terms,
#  same fallback); the overall 4-week degradation is recomputed and checked
#  against downstream_forecast.csv before anything else is written.
#
#  COST: ~360 datasets x 11 districts x (1 fit + 13 refilters). The script
#  prints an estimate of the remaining time after the first 10 datasets.
#  Set RUN_A3 <- FALSE to skip it (e.g. when rerunning sections 1-9 only).

RUN_A3   <- TRUE
LOOKBACK <- 4

if (RUN_A3) {
  H_LONG <- 52; H_SHORT <- 4
  n_tr    <- n_obs - H_LONG
  origins <- seq(n_tr, n_obs - H_SHORT, by = H_SHORT)

  REF_ARIMA_ORDERS <- map(seq_len(n_dis), function(j) {      # as in 03
    y  <- ts(REF[, j], frequency = 52)
    tr <- subset(y, end = n_obs - H_LONG)
    arimaorder(auto.arima(tr, xreg = fourier(tr, K = 2), seasonal = FALSE))[1:3]
  })

  origin_sq_errors <- function(series, j) {
    y   <- ts(series, frequency = 52)
    tr  <- subset(y, end = n_tr)
    xr  <- fourier(tr, K = 2)
    fit <- tryCatch(Arima(tr, order = REF_ARIMA_ORDERS[[j]], xreg = xr),
                    error = function(e) auto.arima(tr, xreg = xr, seasonal = FALSE))
    vapply(origins, function(t) {
      hist <- ts(series[1:t], frequency = 52)
      f <- tryCatch(Arima(hist, model = fit, xreg = fourier(hist, K = 2)),
                    error = function(e) NULL)
      if (is.null(f)) return(NA_real_)
      fc <- forecast(f, xreg = fourier(hist, K = 2, h = H_SHORT), h = H_SHORT)
      sum((as.numeric(fc$mean) - REF[(t + 1):(t + H_SHORT), j])^2)
    }, numeric(1))
  }

  ref_sq <- sapply(seq_len(n_dis), function(j) origin_sq_errors(REF[, j], j))
  # ref_sq: origins x districts

  files <- list.files(file.path(CFG$output_dir, "imputed"), full.names = TRUE)
  files <- files[as.integer(str_match(basename(files), "imp_rep(\\d+)")[, 2])
                 <= CFG$forecast_reps]
  FILE_PATTERN <- "imp_rep(\\d+)_rate(\\d+)_(.+)_(SeasonalTrend|Seasonal|NoTime)\\.csv$"
  cat(sprintf("\nA3 (b): %d datasets x %d districts\n", length(files), n_dis))

  t0 <- Sys.time()
  origin_err <- imap_dfr(files, function(f, i) {
    meta <- str_match(basename(f), FILE_PATTERN)
    rep_i <- as.integer(meta[2]); rate_i <- as.integer(meta[3])
    M <- readRDS(file.path(CFG$output_dir, "masks",
                           sprintf("mask_rep%02d_rate%02d_gap1.rds", rep_i, rate_i)))
    Z <- as.matrix(read_csv(f, show_col_types = FALSE)[, districts])
    out <- map_dfr(seq_len(n_dis), function(j) {
      tibble(rep = rep_i, missing_rate = rate_i, method = meta[4], variant = meta[5],
             district = districts[j], origin = origins,
             sq_err = origin_sq_errors(Z[, j], j), sq_err_ref = ref_sq[, j],
             n_masked_lookback = vapply(origins, function(t)
               sum(M[(t - LOOKBACK + 1):t, j]), numeric(1)),
             last_week_masked = M[origins, j])
    })
    if (i == 10) {
      el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("  10 datasets in %.1f min -> about %.0f min remaining\n",
                  el, el / 10 * (length(files) - 10)))
    }
    out
  })
  cat(sprintf("A3 (b) forecasting done in %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))

  # RMSE over a set of origins (4 errors per origin), as in 03
  rmse_o <- function(se) sqrt(sum(se) / (H_SHORT * length(se)))
  deg <- origin_err |>
    group_by(rep, missing_rate, method, variant, district) |>
    summarise(
      deg_all     = rmse_o(sq_err[!is.na(sq_err)]) - rmse_o(sq_err_ref),   # as in 03
      n_clean     = sum(n_masked_lookback == 0 & !is.na(sq_err)),
      deg_clean   = if (n_clean > 0) {
        k <- n_masked_lookback == 0 & !is.na(sq_err)
        rmse_o(sq_err[k]) - rmse_o(sq_err_ref[k]) } else NA_real_,
      n_exposed   = sum(n_masked_lookback > 0 & !is.na(sq_err)),
      deg_exposed = if (n_exposed > 0) {
        k <- n_masked_lookback > 0 & !is.na(sq_err)
        rmse_o(sq_err[k]) - rmse_o(sq_err_ref[k]) } else NA_real_,
      .groups = "drop")

  # Reproducibility check against 03 before writing anything
  ds <- read_csv(file.path(MET, "downstream_forecast.csv"), show_col_types = FALSE) |>
    select(rep, missing_rate, method, variant, district, degradation_h4)
  chk3 <- inner_join(deg, ds, by = c("rep", "missing_rate", "method", "variant", "district"))
  # Arima() maximises the likelihood numerically: on another R / forecast
  # version the coefficients, hence single forecasts, may differ slightly.
  # The check is therefore made at the level used in the paper (the 24
  # combination means), with the per-series discrepancy printed.
  stopifnot(nrow(chk3) == nrow(deg))
  cat(sprintf("A3 (b): max per-series |difference| with 03 = %.2g\n",
              max(abs(chk3$deg_all - chk3$degradation_h4), na.rm = TRUE)))
  chk_combo <- chk3 |> group_by(method, variant) |>
    summarise(d = mean(deg_all) - mean(degradation_h4), .groups = "drop")
  stopifnot(max(abs(chk_combo$d)) < 0.01)
  cat(sprintf("A3 (b): 24 combination means reproduce 03 (max |diff| = %.2g)\n",
              max(abs(chk_combo$d))))

  crit <- read_csv(file.path(MET, "criteria_comparison.csv"), show_col_types = FALSE) |>
    select(method, variant, MAE)
  leak <- deg |>
    group_by(method, variant) |>
    summarise(degradation_h4_all     = mean(deg_all, na.rm = TRUE),
              degradation_h4_clean   = mean(deg_clean, na.rm = TRUE),
              degradation_h4_exposed = mean(deg_exposed, na.rm = TRUE),
              share_origins_exposed  = sum(n_exposed) / (sum(n_exposed) + sum(n_clean)),
              .groups = "drop") |>
    left_join(crit, by = c("method", "variant")) |>
    arrange(MAE)
  leak_by_rate <- deg |>
    group_by(missing_rate, method, variant) |>
    summarise(degradation_h4_all   = mean(deg_all, na.rm = TRUE),
              degradation_h4_clean = mean(deg_clean, na.rm = TRUE),
              degradation_h4_exposed = mean(deg_exposed, na.rm = TRUE),
              .groups = "drop")
  sp <- function(y) {
    ct <- suppressWarnings(cor.test(leak$MAE, y, method = "spearman"))
    c(rho = unname(ct$estimate), p = ct$p.value)
  }
  leak_summary <- tibble(
    criterion = c("degradation_h4_all", "degradation_h4_clean", "degradation_h4_exposed"),
    spearman_with_MAE = c(sp(leak$degradation_h4_all)[1],
                          sp(leak$degradation_h4_clean)[1],
                          sp(leak$degradation_h4_exposed)[1]),
    p_value = c(sp(leak$degradation_h4_all)[2],
                sp(leak$degradation_h4_clean)[2],
                sp(leak$degradation_h4_exposed)[2]),
    n_negative_of_24 = c(sum(leak$degradation_h4_all < 0),
                         sum(leak$degradation_h4_clean < 0),
                         sum(leak$degradation_h4_exposed < 0)),
    lookback_weeks = LOOKBACK)

  write_csv(origin_err,   file.path(MET, "forecast_origin_errors.csv"))
  write_csv(leak,         file.path(MET, "forecast_leakage_diagnostic.csv"))
  write_csv(leak_by_rate, file.path(MET, "forecast_leakage_by_rate.csv"))
  write_csv(leak_summary, file.path(MET, "forecast_leakage_summary.csv"))
  print(leak, n = 30, width = Inf)
  print(leak_summary, width = Inf)
}


## -------------------------------------------------------------------------
## 11. A6 - EFFECTIVE GAP LENGTHS OF THE MASKS, AND 8-WEEK MECHANISM SUMMARY
## -------------------------------------------------------------------------
#  (a) The nominal gap length (1, 4 or 8 weeks) is the length of the blocks
#      that simulate_mask() draws. Blocks are drawn independently, so they
#      can touch or overlap: isolated 1-week cells form runs of 2, 3 or more
#      weeks by chance, and more often when the rate is high; 4- and 8-week
#      blocks can overlap or merge into longer runs, and the realised rate can
#      exceed the target. Blocks start at i <= n_obs - gap, so for 4- and
#      8-week gaps the last week of the series is never masked.
#      This section measures, on every archived mask:
#        - the realised number of masked cells against the nominal target;
#        - the distribution of the lengths of the runs of consecutive missing
#          weeks (per district column);
#        - by REPLAYING simulate_mask() with the seed used by 02
#          (1000 * rep + rate + 7 * gap), the number of blocks drawn and how
#          many overlapped or touched an earlier block; the replay is checked
#          cell by cell against the archived mask;
#        - whether the last week of the series is ever masked.
#      It answers: does a higher missingness rate also lengthen the effective
#      gaps, i.e. are rate and effective length confounded in the design?
#  (b) The mechanism arms (07) wrote mean ranks only (ranks_four_mechanisms
#      .csv) and bias for 1-week gaps only (bias_four_mechanisms.csv). The
#      text needs, for 8-week gaps, which method leads under each mechanism
#      and how large the bias is. The table of 07 (section 6) is recomputed
#      here with the SAME code and checked against ranks_four_mechanisms.csv.

# (a) run lengths, block replay and last week of the masks
mask_files <- list.files(file.path(CFG$output_dir, "masks"),
                         "^mask_rep\\d+_rate\\d+_gap\\d+\\.rds$", full.names = TRUE)
runs_of <- function(M) unlist(lapply(seq_len(ncol(M)), function(j) {
  r <- rle(M[, j]); r$lengths[r$values] }))

# Same random calls, in the same order, as simulate_mask() (functions.R):
# the archived mask must be reproduced exactly (checked below).
replay_blocks <- function(rate, seed, gap) {
  set.seed(seed)
  M <- matrix(FALSE, n_obs, n_dis)
  n_na <- round(n_obs * n_dis * rate / 100)
  if (gap == 1) {
    M[sample.int(n_obs * n_dis, n_na)] <- TRUE
    return(list(M = M, n_na = n_na, n_blocks = n_na, n_overlap = 0L,
                n_touch = NA_integer_, cells_overlapped = 0L))
  }
  n_blocks <- 0L; n_overlap <- 0L; n_touch <- 0L; cells_overlapped <- 0L
  while (sum(M) < n_na) {
    j <- sample.int(n_dis, 1); i <- sample.int(n_obs - gap, 1)
    idx <- i:(i + gap - 1)
    ov  <- sum(M[idx, j])
    touch <- ov == 0 && ((i > 1 && M[i - 1, j]) || (i + gap <= n_obs && M[i + gap, j]))
    n_blocks <- n_blocks + 1L
    if (ov > 0) n_overlap <- n_overlap + 1L
    if (touch)  n_touch   <- n_touch + 1L
    cells_overlapped <- cells_overlapped + ov
    M[idx, j] <- TRUE
  }
  list(M = M, n_na = n_na, n_blocks = n_blocks, n_overlap = n_overlap,
       n_touch = n_touch, cells_overlapped = cells_overlapped)
}

mask_info <- list()
mask_runs <- map_dfr(mask_files, function(f) {
  meta <- as.integer(str_match(basename(f), "mask_rep(\\d+)_rate(\\d+)_gap(\\d+)")[, 2:4])
  M  <- readRDS(f)
  rp <- replay_blocks(meta[2], 1000 * meta[1] + meta[2] + 7 * meta[3], meta[3])
  stopifnot(all(unname(M) == rp$M))          # replay reproduces the archived mask
  mask_info[[length(mask_info) + 1]] <<- tibble(
    rep = meta[1], missing_rate = meta[2], gap_length = meta[3],
    target_cells = rp$n_na, realised_cells = sum(M),
    blocks_drawn = rp$n_blocks, blocks_overlapping = rp$n_overlap,
    blocks_touching = rp$n_touch, cells_overlapped = rp$cells_overlapped,
    last_week_masked = any(M[n_obs, ]), last_masked_row = max(which(rowSums(M) > 0)))
  tibble(rep = meta[1], missing_rate = meta[2], gap_length = meta[3],
         realised_rate = 100 * mean(M), run_length = runs_of(M))
})
mask_info <- bind_rows(mask_info)
gap_dist_summary <- mask_runs |>
  group_by(missing_rate, gap_length) |>
  summarise(
    n_masks                 = n_distinct(rep),
    realised_rate_mean      = mean(realised_rate[!duplicated(rep)]),
    realised_rate_max       = max(realised_rate),
    runs_per_mask           = n() / n_distinct(rep),
    mean_run_length         = mean(run_length),
    max_run_length          = max(run_length),
    # share of masked CELLS lying in runs longer than the nominal length
    share_cells_longer_than_nominal = sum(run_length[run_length > gap_length[1]]) /
                                      sum(run_length),
    share_cells_run_ge_2    = sum(run_length[run_length >= 2]) / sum(run_length),
    share_cells_run_ge_3    = sum(run_length[run_length >= 3]) / sum(run_length),
    share_cells_run_ge_2x_nominal = sum(run_length[run_length >= 2 * gap_length[1]]) /
                                    sum(run_length),
    .groups = "drop") |>
  left_join(mask_info |>
              group_by(missing_rate, gap_length) |>
              summarise(
                target_cells            = first(target_cells),
                realised_cells_mean     = mean(realised_cells),
                overshoot_cells_max     = max(realised_cells - target_cells),
                blocks_drawn_mean       = mean(blocks_drawn),
                share_blocks_overlapping = sum(blocks_overlapping) / sum(blocks_drawn),
                share_blocks_touching   = sum(blocks_touching) / sum(blocks_drawn),
                cells_overlapped_mean   = mean(cells_overlapped),
                masks_with_last_week_masked = sum(last_week_masked),
                last_masked_row_max     = max(last_masked_row),
                n_obs                   = n_obs,
                .groups = "drop"),
            by = c("missing_rate", "gap_length"))
n_masks_tbl <- mask_runs |> distinct(missing_rate, gap_length, rep) |>
  count(missing_rate, gap_length, name = "n_masks")
gap_dist_long <- mask_runs |>
  count(missing_rate, gap_length, run_length, name = "n_runs") |>
  left_join(n_masks_tbl, by = c("missing_rate", "gap_length")) |>
  group_by(missing_rate, gap_length) |>
  mutate(runs_per_mask = n_runs / n_masks,
         share_of_masked_cells = n_runs * run_length / sum(n_runs * run_length)) |>
  ungroup()
write_csv(gap_dist_summary, file.path(MET, "gap_length_distribution.csv"))
write_csv(gap_dist_long,    file.path(MET, "gap_length_distribution_long.csv"))
print(gap_dist_summary, width = Inf)

# (b) mechanism arms: mean rank, MAE, bias by mechanism and gap (07, section 6)
mcar_m <- read_csv(file.path(MET, "metrics_long.csv"), show_col_types = FALSE) |>
  filter(missing_rate == CFG_MECH$missing_rate,
         gap_length %in% CFG_MECH$gap_lengths, rep <= CFG_MECH$n_rep) |>
  mutate(mechanism = "MCAR")
others_m <- read_csv(file.path(MMET, "metrics_mechanism.csv"), show_col_types = FALSE)
strong_m <- read_csv(file.path(MMET, "metrics_mar_strong.csv"), show_col_types = FALSE)
ARM_ORDER <- c("MCAR", "MAR", "MAR_strong", "MNAR")
all_arms <- bind_rows(mcar_m |> select(any_of(names(strong_m))), others_m, strong_m) |>
  mutate(mechanism = factor(mechanism, levels = ARM_ORDER))
bv <- all_arms |>
  group_by(mechanism, gap_length, method, variant) |>
  summarise(m = mean(MAE), .groups = "drop") |>
  group_by(mechanism, gap_length, method) |>
  slice_min(m, n = 1, with_ties = FALSE) |> ungroup() |>
  select(mechanism, gap_length, method, variant)
mech_ranks <- all_arms |>
  inner_join(bv, by = c("mechanism", "gap_length", "method", "variant")) |>
  group_by(mechanism, gap_length, rep) |> mutate(rank = rank(MAE)) |> ungroup() |>
  group_by(mechanism, gap_length, method, variant) |>
  summarise(mean_rank = mean(rank), MAE = mean(MAE), bias = mean(bias),
            sensitivity = mean(sensitivity), youden = mean(youden),
            n_reps = n(), .groups = "drop") |>
  arrange(mechanism, gap_length, mean_rank)
# check against 07
r4 <- read_csv(file.path(MMET, "ranks_four_mechanisms.csv"), show_col_types = FALSE) |>
  pivot_longer(-method, names_to = "arm", values_to = "mean_rank_07") |>
  mutate(gap_length = as.numeric(sub(".*_", "", arm)),
         mechanism  = sub("_[0-9]+$", "", arm))
chk_m <- mech_ranks |> mutate(mechanism = as.character(mechanism)) |>
  inner_join(r4, by = c("mechanism", "gap_length", "method"))
stopifnot(nrow(chk_m) == 80, max(abs(chk_m$mean_rank - chk_m$mean_rank_07)) < 1e-9)
cd_red <- read_csv(file.path(MMET, "nemenyi_cd_reduced.csv"), show_col_types = FALSE)
mech_ranks <- mech_ranks |>
  group_by(mechanism, gap_length) |>
  mutate(rank_diff_to_first = mean_rank - min(mean_rank)) |> ungroup() |>
  left_join(cd_red |> distinct(gap_length, cd_reduced = cd), by = "gap_length")
write_csv(mech_ranks, file.path(MMET, "mechanism_ranks_long.csv"))
print(mech_ranks |> filter(gap_length == 8) |> group_by(mechanism) |>
        slice_min(mean_rank, n = 4) |> ungroup(), n = 20, width = Inf)


## -------------------------------------------------------------------------
## 12. KALMAN FALLBACK CHECK (step 8, October 2026)
## -------------------------------------------------------------------------
#  imp_kalman_vector() (functions.R) tries na_kalman(model = "StructTS")
#  first and falls back to auto.arima, then to linear interpolation. The
#  counter N_KALMAN_FALLBACKS is printed by 02 only when it is > 0 and is
#  never printed by 05 or 07, so no output file records it. This section
#  replays the FIRST level of the cascade, with the same call and the same
#  success test, on every archived mask of the main benchmark (180 masks)
#  and of the mechanism arms (60 masks), 11 districts each, and counts the
#  series on which StructTS would have failed. It imputes nothing that is
#  scored or written; it only records success or failure. Runtime: about
#  one minute.
kalman_try <- function(v) {
  r <- try(as.numeric(imputeTS::na_kalman(v, model = "StructTS")), silent = TRUE)
  !inherits(r, "try-error") && !anyNA(r)          # same test as functions.R
}
kalman_masks <- c(
  list.files(file.path(CFG$output_dir, "masks"), "\\.rds$", full.names = TRUE),
  list.files(file.path(DOS_MECH, "masks"), "\\.rds$", full.names = TRUE))
kalman_check <- map_dfr(kalman_masks, function(f) {
  M <- readRDS(f); X <- REF; X[M] <- NA
  ok <- vapply(seq_len(n_dis), function(j) kalman_try(X[, j]), logical(1))
  tibble(arm = if (grepl(DOS_MECH, f, fixed = TRUE)) "mechanism" else "main",
         mask = basename(f), n_series = n_dis, n_fallback = sum(!ok),
         districts_fallback = paste(districts[!ok], collapse = ";"))
})
kalman_summary <- kalman_check |>
  group_by(arm) |>
  summarise(n_masks = n(), n_series = sum(n_series),
            n_fallback = sum(n_fallback), .groups = "drop")
write_csv(bind_rows(kalman_summary |> mutate(mask = "ALL"),
                    kalman_check |> filter(n_fallback > 0)),
          file.path(MET, "kalman_fallback_check.csv"))
print(kalman_summary)

cat(sprintf("\n11 done: %d files written.\n", 9 + 3 + 3 + 1 + 3 + 1 + if (RUN_A3) 4 else 0))
