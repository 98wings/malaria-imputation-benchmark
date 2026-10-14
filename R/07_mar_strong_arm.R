###############################################################################
#  07 — THE "STRONG MAR" ARM: removing the confound left by 05
#
#  THE PROBLEM THIS SCRIPT SOLVES
#  The corrected diagnostics quantified the selection strength on the value
#  (Cohen's d, isolated gaps):
#        MCAR  d = 0.01      (no selection)
#        MAR   d = 0.14      (NEGLIGIBLE selection)
#        MNAR  d = 0.95      (STRONG selection)
#  Two things therefore vary at once between MAR and MNAR: the NATURE of the
#  mechanism and its INTENSITY. When MNAR behaves differently, we cannot tell
#  which of the two caused it. That is a confound, and it is visible in the
#  table: a reviewer will see it.
#
#  Moreover, d = 0.14 is "negligible" by Cohen's convention: concluding "the
#  ranking withstands MAR" from a mechanism that could barely have broken it
#  would be an empty conclusion.
#
#  THE SOLUTION
#  A fourth arm, MAR_strong, whose propensity depends on the mean of the
#  OTHER districts in the same week. Those values are observed, so the
#  construction remains MAR. Because the districts are strongly correlated,
#  the induced selection on the value becomes strong. The intensity is
#  CALIBRATED AUTOMATICALLY to match that of MNAR.
#
#  This yields the clean contrast:
#        MNAR       : propensity depends on the MISSING value
#        MAR_strong : propensity depends on OBSERVED values
#        ... at EQUAL selection strength.
#  Any remaining difference is attributable to the nature of the mechanism,
#  and to that alone.
#
#  Prerequisites : 02, 05, 06.   Runtime: ~35 min.
###############################################################################

source("R/functions_mechanism.R")

MECH <- "MAR_strong"

## -------------------------------------------------------------------------
## 1. HONESTY ABOUT THE LIMIT — to be written into the paper
## -------------------------------------------------------------------------
#  Exact MAR requires the conditioning variable to be ENTIRELY observed. Here
#  part of the neighbouring districts is itself masked in the same
#  replication. With 11 districts and 10% missingness, roughly 1 neighbour in
#  10 is unavailable: the mean remains computable and close to the true one.
#  We MEASURE that gap below rather than passing over it in silence.
#  (The "exact MAR" alternative would require permanently excluding some
#  districts from masking, which would distort the experimental design.)

## -------------------------------------------------------------------------
## 2. A DIFFICULTY THAT IS ALSO A RESULT
## -------------------------------------------------------------------------
#  First attempt: the draw failed from beta = 1.25 onwards on
#      all(rowSums(!M) > 0) is not TRUE
#  that is, ENTIRELY MISSING WEEKS (11 districts out of 11).
#
#  This is not a programming bug, it is a property of the mechanism.
#  NEIGH[i,j] is the mean of the other ten districts in week i: that quantity
#  is dominated by the WEEK'S LEVEL and barely varies from one district to
#  another. It is therefore very nearly a ROW variable. Inevitable
#  consequence: the more one strengthens a MAR built on spatial covariation,
#  the more missingness concentrates into COMPLETE WEEKS.
#
#  => In a spatial panel, mechanism and missingness STRUCTURE are not
#     independent. One cannot strengthen MAR selection without distorting the
#     structure. This belongs in the Discussion: it is precisely why this
#     paper treats structure as an axis in its own right rather than as a
#     detail of the mechanism.
#
#  Handling: cap missingness at MAX_MISSING_PER_WEEK districts per week (see
#  R/functions_mechanism.R for why this is a scope definition, not a hack).
MAX_MISSING_PER_WEEK <- n_dis - 3L    # at least 3 districts observed per week
N_CONSTRAINT_HITS    <- 0L

## -------------------------------------------------------------------------
## 3. AUTOMATIC CALIBRATION OF THE INTENSITY
## -------------------------------------------------------------------------
#  We do not guess FORCE_NEIGHBOUR: we search for the value that reproduces
#  MNAR's Cohen's d. Only masks are simulated here (no imputation), so the
#  calibration costs seconds.
f_diag <- file.path(DOS_MECH, "metrics", "mechanism_diagnostics.csv")
stopifnot(file.exists(f_diag))
TARGET <- read_csv(f_diag, show_col_types = FALSE) |>
  filter(mechanism == "MNAR", gap_length == 1) |> pull(cohen_d) |> mean()
cat(sprintf("Target selection strength (MNAR Cohen's d, isolated gaps): %.3f\n",
            TARGET))

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
cat("Calibrating (d rises then PLATEAUS -- see the note below):\n")
ds <- sapply(grid, function(b) {
  d <- observed_d(b)
  cat(sprintf("   beta = %.2f  ->  d = %.3f\n", b, d)); d
})

#  PICKING THE BETA. The d curve plateaus well before the end of the grid.
#  Taking the beta closest to the target picks an arbitrarily large value on
#  the plateau, which needlessly pressures the row constraint and distorts
#  the seasonal component at gap 8. We therefore take the SMALLEST beta whose
#  d is within tolerance of the target.
TOL <- 0.03
within <- which(abs(ds - TARGET) <= TOL)
if (length(within)) {
  FORCE_NEIGHBOUR <- grid[min(within)]
  reached <- ds[min(within)]
} else {
  FORCE_NEIGHBOUR <- grid[which.max(ds)]
  reached <- max(ds)
  cat("\n>>> THE TARGET IS NOT REACHABLE. THIS IS A RESULT, NOT A FAILURE.\n")
  cat("A MAR mechanism can select high values only INDIRECTLY, through the\n")
  cat("correlation between districts. That correlation therefore bounds the\n")
  cat("achievable selection. In other words: in these data, no MAR can be as\n")
  cat("selective as an MNAR. Wording for the paper:\n")
  cat(sprintf("  \"The maximum selection achievable under MAR (d = %.2f) remains\n", max(ds)))
  cat(sprintf("   below that of the simulated MNAR (d = %.2f), which structurally\n", TARGET))
  cat("   bounds the severity of MAR missingness in a spatial panel.\"\n")
}
cat(sprintf("\nFORCE_NEIGHBOUR retained: %.2f  (d reached = %.3f, target %.3f)\n",
            FORCE_NEIGHBOUR, reached, TARGET))
cat(sprintf("Maximum d achievable on the grid: %.3f (beta = %.2f)\n",
            max(ds), grid[which.max(ds)]))
cat(sprintf("Row constraint hits during calibration: %d\n", N_CONSTRAINT_HITS))
cat("A high count confirms that strong spatial MAR concentrates missingness\n",
    "into complete weeks (see section 2).\n", sep = "")
N_CONSTRAINT_HITS <- 0L

## -------------------------------------------------------------------------
## 4. CONTAMINATION: how degraded is the conditioning variable?
## -------------------------------------------------------------------------
#  A number to cite in the paper to justify that the construction is MAR in
#  practice. We compare the true neighbour mean with the one an analyst could
#  compute from the observed data alone.
M_test <- simulate_mask_mechanism(CFG_MECH$missing_rate, 700099L,
                                  gap_length = 1, mechanism = MECH)
X_test <- REF; X_test[M_test] <- NA
NEIGH_obs <- t(sapply(seq_len(n_obs), function(i)
  sapply(seq_len(n_dis), function(j) mean(X_test[i, -j], na.rm = TRUE))))
cat("\nContamination of the conditioning variable:\n")
cat(sprintf("  correlation(true neighbours, observable neighbours) = %.4f\n",
            cor(as.numeric(NEIGH), as.numeric(NEIGH_obs), use = "complete.obs")))
cat(sprintf("  share of unavailable neighbours per cell           = %.1f %%\n",
            100 * mean(is.na(X_test[, -1]))))

## -------------------------------------------------------------------------
## 5. LOOP (identical to 05, a single mechanism)
## -------------------------------------------------------------------------
res_g <- list(); res_d <- list(); res_v <- list(); res_k <- list(); kk <- 0

for (gap in CFG_MECH$gap_lengths) {
for (rep in seq_len(CFG_MECH$n_rep)) {
  rate <- CFG_MECH$missing_rate
  seed <- 700000L + 1000L*rep + rate + 7L*gap   # disjoint from the other arms
  M <- simulate_mask_mechanism(rate, seed, gap, MECH)
  saveRDS(M, file.path(DOS_MECH, "masks",
                       sprintf("mask_%s_rep%02d_gap%d.rds", MECH, rep, gap)))
  X <- REF; X[M] <- NA

  res_v[[length(res_v) + 1]] <- diagnose_mechanism(M, MECH) |>
    mutate(rep = rep, gap_length = gap, .before = 1)

  if (rep <= CFG_MECH$tuning_reps) {
    k_opt <- tune_k(X, rate, B = CFG_MECH$tuning_B)
  } else {
    seen  <- bind_rows(res_k) |> filter(gap_length == !!gap)
    k_opt <- as.integer(names(sort(table(seen$value), decreasing = TRUE))[1])
  }
  res_k[[length(res_k) + 1]] <- tibble(
    mechanism = MECH, rep = rep, gap_length = gap,
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
        mutate(mechanism = MECH, rep = rep, missing_rate = rate,
               gap_length = gap, method = nm, variant = v,
               sec = as.numeric(tps$toc - tps$tic), .before = 1)
      res_d[[kk]] <- compute_metrics_by_district(Z, M) |>
        mutate(mechanism = MECH, rep = rep, gap_length = gap,
               method = nm, variant = v, .before = 1)
    }
  }
  cat(sprintf("  %s / gap %d / rep %02d done\n", MECH, gap, rep))
}}

metrics_strong <- bind_rows(res_g)
write_csv(metrics_strong, file.path(DOS_MECH, "metrics", "metrics_mar_strong.csv"))
write_csv(bind_rows(res_d), file.path(DOS_MECH, "metrics", "district_mar_strong.csv"))
write_csv(bind_rows(res_v) |> mutate(force_neighbour = FORCE_NEIGHBOUR),
          file.path(DOS_MECH, "metrics", "validation_mar_strong.csv"))
cat(sprintf("\nRow constraint active %d times across the %d masks of this arm.\n",
            N_CONSTRAINT_HITS, 2 * CFG_MECH$n_rep))

## -------------------------------------------------------------------------
## 6. ANALYSIS ACROSS THE FOUR ARMS
## -------------------------------------------------------------------------
mcar <- read_csv("output/metrics/metrics_long.csv", show_col_types = FALSE) |>
  filter(missing_rate == CFG_MECH$missing_rate,
         gap_length %in% CFG_MECH$gap_lengths, rep <= CFG_MECH$n_rep) |>
  mutate(mechanism = "MCAR")
others <- read_csv(file.path(DOS_MECH, "metrics", "metrics_mechanism.csv"),
                   show_col_types = FALSE)
ARM_ORDER <- c("MCAR", "MAR", "MAR_strong", "MNAR")
all_arms <- bind_rows(mcar |> select(any_of(names(metrics_strong))),
                      others, metrics_strong) |>
  mutate(mechanism = factor(mechanism, levels = ARM_ORDER))

best_variant <- all_arms |>
  group_by(mechanism, gap_length, method, variant) |>
  summarise(m = mean(MAE), .groups = "drop") |>
  group_by(mechanism, gap_length, method) |>
  slice_min(m, n = 1, with_ties = FALSE) |> ungroup() |>
  select(mechanism, gap_length, method, variant)
perf <- all_arms |> inner_join(best_variant,
                               by = c("mechanism", "gap_length", "method", "variant"))

ranks <- perf |> group_by(mechanism, gap_length, rep) |>
  mutate(rank = rank(MAE)) |> ungroup() |>
  group_by(mechanism, gap_length, method) |>
  summarise(mean_rank = mean(rank), MAE = mean(MAE), bias = mean(bias),
            sensitivity = mean(sensitivity), youden = mean(youden),
            .groups = "drop")

cat("\n=== MEAN RANK BY MECHANISM ===\n")
tab <- ranks |> select(mechanism, gap_length, method, mean_rank) |>
  pivot_wider(names_from = c(mechanism, gap_length), values_from = mean_rank)
ref_col <- if ("MCAR_1" %in% names(tab)) "MCAR_1" else names(tab)[2]
tab <- tab |> arrange(.data[[ref_col]])
print(as.data.frame(tab), row.names = FALSE, digits = 3)
write_csv(tab, file.path(DOS_MECH, "metrics", "ranks_four_mechanisms.csv"))

cat("\n=== AGREEMENT OF THE RANKINGS (Spearman) ===\n")
pairs <- expand.grid(a = ARM_ORDER, b = ARM_ORDER, gap = CFG_MECH$gap_lengths,
                     stringsAsFactors = FALSE) |>
  filter(match(a, ARM_ORDER) < match(b, ARM_ORDER))
for (i in seq_len(nrow(pairs))) {
  ra <- ranks |> filter(mechanism == pairs$a[i], gap_length == pairs$gap[i]) |>
    arrange(method) |> pull(mean_rank)
  rb <- ranks |> filter(mechanism == pairs$b[i], gap_length == pairs$gap[i]) |>
    arrange(method) |> pull(mean_rank)
  ct <- suppressWarnings(cor.test(ra, rb, method = "spearman"))
  mark <- if (pairs$a[i] == "MAR_strong" && pairs$b[i] == "MNAR")
    "   <-- THE KEY CONTRAST (equal intensity)" else ""
  cat(sprintf("  gap %d: %-11s vs %-11s rho = %+.3f  p = %.4f%s\n",
              pairs$gap[i], pairs$a[i], pairs$b[i],
              ct$estimate, ct$p.value, mark))
}

cat("\n=== REVERSAL WITH GAP LENGTH ===\n")
inv <- ranks |> select(mechanism, gap_length, method, mean_rank) |>
  pivot_wider(names_from = gap_length, values_from = mean_rank,
              names_prefix = "gap") |>
  mutate(delta = gap8 - gap1) |> arrange(mechanism, gap1)
print(as.data.frame(inv), row.names = FALSE, digits = 3)
write_csv(inv, file.path(DOS_MECH, "metrics", "reversal_four_mechanisms.csv"))
d_wide <- inv |> select(mechanism, method, delta) |>
  pivot_wider(names_from = mechanism, values_from = delta)
cat("\nConsistency of the reversal DIRECTION (Spearman on the deltas):\n")
for (m2 in setdiff(names(d_wide), c("method", "MCAR"))) {
  ct <- suppressWarnings(cor.test(d_wide$MCAR, d_wide[[m2]], method = "spearman"))
  cat(sprintf("  MCAR vs %-11s: rho = %+.3f  p = %.4f\n",
              m2, ct$estimate, ct$p.value))
}

cat("\n=== BIAS AND PEAK DETECTION (isolated gaps) ===\n")
tb <- ranks |> filter(gap_length == 1) |>
  select(mechanism, method, bias, sensitivity) |>
  pivot_wider(names_from = mechanism, values_from = c(bias, sensitivity))
print(as.data.frame(tb), row.names = FALSE, digits = 3)
write_csv(tb, file.path(DOS_MECH, "metrics", "bias_four_mechanisms.csv"))
cat("\nREADING: if bias(MAR_strong) stays small while bias(MNAR) is strongly\n",
    "negative AT EQUAL SELECTION STRENGTH, then the estimation bias comes\n",
    "from the MNAR NATURE of the missingness, not merely from its intensity.\n",
    "That is the argument the MAR arm alone could not support.\n", sep = "")

## -------------------------------------------------------------------------
## 7. FIGURE
## -------------------------------------------------------------------------
p <- ranks |>
  mutate(gap_lab = factor(paste0(gap_length, "-week gaps"),
                          levels = paste0(CFG_MECH$gap_lengths, "-week gaps")),
         method = reorder(method, mean_rank)) |>
  ggplot(aes(mean_rank, method, colour = mechanism)) +
  geom_line(aes(group = method), colour = "grey85", linewidth = 3, alpha = .6) +
  geom_point(size = 3) + facet_wrap(~ gap_lab) +
  scale_colour_manual(values = c(MCAR = "#1b9e77", MAR = "#d95f02",
                                 MAR_strong = "#e7298a", MNAR = "#7570b3"),
                      name = "Mechanism") +
  labs(x = "Mean rank (1 = best)", y = NULL,
       title = "Robustness of the ranking to the missingness mechanism",
       subtitle = sprintf("MAR_strong and MNAR are calibrated to equal selection strength (Cohen's d ~ %.2f)",
                          TARGET)) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(file.path(DOS_MECH, "figures", "fig9b_four_mechanisms.png"), p,
       width = 10, height = 6, dpi = 300)

cat("\n=== Done ===\n")
cat("Now re-run source(\"R/06_mechanism_diagnostics.R\"): it will also read\n",
    "the MAR_strong masks and complete the diagnostics table.\n", sep = "")
