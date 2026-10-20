###############################################################################
#  04 — FINAL ANALYSES, DECISION GUIDE AND FIGURES
#
#  Run AFTER 03 and 07: re-reads their result files, recomputes no imputation.
#  Runtime: a few minutes.
#
#  Produces:
#    - district-level analysis (spatial component, level 1)
#    - critical difference diagram (the summary figure)
#    - computational cost and the accuracy / time trade-off
#    - a quantified decision guide
#    - the paper's figures
###############################################################################

if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(tidyverse); library(lubridate)

CFG <- list(output_dir = "output")
MET <- file.path(CFG$output_dir, "metrics")
FIG <- file.path(CFG$output_dir, "figures")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

metrics         <- read_csv(file.path(MET, "metrics_long.csv"),        show_col_types = FALSE)
distr           <- read_csv(file.path(MET, "metrics_by_district.csv"), show_col_types = FALSE)
characteristics <- read_csv("data/derived/reference_characteristics.csv",
                            show_col_types = FALSE)
ranks           <- read_csv(file.path(MET, "friedman_ranks_main.csv"), show_col_types = FALSE)

# Number of blocks (missingness rate x repetition) expected in every
# per-cell summary below: 60 at n_rep = 20, 6 in a dry run at n_rep = 2.
# Read from the metrics rather than written as 60, so that the checks hold
# whatever CFG$n_rep was used in 02 (this script does not load functions.R).
N_BLOCKS <- n_distinct(metrics$missing_rate) * n_distinct(metrics$rep)

UNIVARIATE <- c("LinearInterp", "Kalman", "SeasonalDecomp")
FAMILY <- c(Mean = "Naive",
            MICE_inference = "Cross-sectional", MICE_precision = "Cross-sectional",
            MICE_deterministic = "Cross-sectional", KNN = "Cross-sectional",
            MissForest = "Cross-sectional",
            LinearInterp = "Temporal", Kalman = "Temporal",
            SeasonalDecomp = "Temporal",
            MissForest_ST = "Hybrid")


## =========================================================================
## 1. DISTRICT-LEVEL ANALYSIS  (spatial component, level 1)
## =========================================================================
#  Question: does the relative performance of the methods depend on district
#  characteristics? This is what turns a global ranking into an operational
#  recommendation.
dis <- distr |>
  filter(gap_length == 1) |>
  group_by(district, method, variant) |>
  summarise(MAE = mean(MAE), MAPE = mean(MAPE), .groups = "drop")

# MAPE table, district x method (best variant of each method).
# with_ties = FALSE is ESSENTIAL: Mean ignores the period covariates, so its
# three variants are strictly identical. Without it slice_min returns all
# three ties, the join triples the rows and pivot_wider produces list-columns.
best_variant <- dis |> group_by(method, variant) |>
  summarise(m = mean(MAE), .groups = "drop") |>
  group_by(method) |> slice_min(m, n = 1, with_ties = FALSE) |>
  ungroup() |> select(method, variant)

district_table <- dis |>
  inner_join(best_variant, by = c("method", "variant")) |>
  select(district, method, MAPE) |>
  pivot_wider(names_from = method, values_from = MAPE) |>
  left_join(characteristics |>
              select(district, mean_corr_other_districts, cv, seasonal_r2,
                     autocorr_lag1),
            by = "district") |>
  arrange(mean_corr_other_districts)
print(district_table, n = 20, width = Inf)
write_csv(district_table, file.path(MET, "district_analysis.csv"))

#  Two relationships to test (exploratory: only 11 districts)
#   (a) does the LEVEL of error track the district's own volatility?
#   (b) does the GAIN of multivariate methods track correlation with peers?
# Same precaution: keep ONE variant per method before pivoting.
gain <- dis |>
  inner_join(best_variant, by = c("method", "variant")) |>
  filter(method %in% c("Mean", "MissForest", "LinearInterp")) |>
  select(district, method, MAE) |>
  pivot_wider(names_from = method, values_from = MAE) |>
  mutate(gain_MissForest_pct   = 100 * (Mean - MissForest) / Mean,
         gain_LinearInterp_pct = 100 * (Mean - LinearInterp) / Mean) |>
  left_join(characteristics, by = "district")

cat("\n=== District relationships (Spearman, n = 11, EXPLORATORY) ===\n")
test_pair <- function(x, y, nx, ny) {
  ct <- suppressWarnings(cor.test(x, y, method = "spearman"))
  cat(sprintf("  %-24s vs %-26s rho = %+.2f   p = %.4f\n",
              nx, ny, ct$estimate, ct$p.value))
}
CHARS <- c("cv", "mean_corr_other_districts", "seasonal_r2", "autocorr_lag1")
for (v in CHARS) test_pair(gain$MissForest, gain[[v]], "MAE MissForest", v)
cat("\n")
for (v in CHARS) test_pair(gain$gain_MissForest_pct, gain[[v]],
                           "MissForest gain (%)", v)
# With 8 tests, a Bonferroni-corrected threshold is 0.05/8 = 0.00625.
cat("\n  Bonferroni threshold for 8 tests: p < 0.00625\n")
write_csv(gain, file.path(MET, "gain_by_district.csv"))


## =========================================================================
## 2. CRITICAL DIFFERENCE DIAGRAM
## =========================================================================
#  The figure that summarises Friedman + Nemenyi at a glance.
#  Two methods less than CD apart are NOT distinguishable.
n_blocks <- metrics |> filter(gap_length == 1, variant == "NoTime") |>
  distinct(rep, missing_rate) |> nrow()
k  <- nrow(ranks)
CD <- qtukey(0.95, k, Inf) / sqrt(2) * sqrt(k * (k + 1) / (6 * n_blocks))
cat(sprintf("\nCritical difference (Nemenyi, k = %d, N = %d): CD = %.2f\n",
            k, n_blocks, CD))

rg <- ranks |> arrange(mean_rank) |> mutate(y = row_number())
# maximal groups of indistinguishable methods
grp <- list()
for (i in seq_len(k)) {
  j <- max(which(rg$mean_rank - rg$mean_rank[i] < CD))
  if (j > i) grp[[length(grp) + 1]] <- c(i, j)
}
grp <- grp[!duplicated(map_chr(grp, ~ paste(.x, collapse = "-")))]
grp <- grp[!map_lgl(seq_along(grp), function(a)
  any(map_lgl(seq_along(grp), function(b)
    b != a && grp[[b]][1] <= grp[[a]][1] && grp[[b]][2] >= grp[[a]][2])))]
seg <- map_dfr(seq_along(grp), ~ tibble(
  x = rg$mean_rank[grp[[.x]][1]], xend = rg$mean_rank[grp[[.x]][2]],
  y = -0.6 - 0.35 * .x))

p_cd <- ggplot(rg, aes(mean_rank, 0)) +
  geom_segment(aes(xend = mean_rank, yend = -0.15), linewidth = .3) +
  geom_point(size = 2.4) +
  geom_text(aes(label = method, y = 0.18), angle = 35, hjust = 0, size = 3.2) +
  geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = y),
               linewidth = 1.6, colour = "grey25", inherit.aes = FALSE) +
  annotate("segment", x = min(rg$mean_rank), xend = min(rg$mean_rank) + CD,
           y = 0.75, yend = 0.75, linewidth = .8, colour = "#c0392b") +
  annotate("text", x = min(rg$mean_rank) + CD / 2, y = 0.92,
           label = sprintf("CD = %.2f", CD), size = 3.2, colour = "#c0392b") +
  scale_x_continuous("Mean rank (1 = best)", breaks = 1:10) +
  scale_y_continuous(NULL, limits = c(-0.6 - 0.35 * length(grp) - 0.2, 1.8)) +
  labs(title = "Comparison of imputation methods",
       subtitle = sprintf("Friedman + Nemenyi post-hoc, %d blocks. Connected methods are not distinguishable.",
                          n_blocks)) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())
ggsave(file.path(FIG, "fig1_critical_difference.png"), p_cd,
       width = 9, height = 4.5, dpi = 300)


## =========================================================================
## 3. COMPUTATIONAL COST AND TRADE-OFF
## =========================================================================
#  The central argument of the "limited resources" framing: how much accuracy
#  do we buy per additional second of computation?
ref_mae <- metrics |> filter(gap_length == 1, method == "Mean") |>
  pull(MAE) |> mean()
cost <- metrics |>
  filter(gap_length == 1) |>
  group_by(method, variant) |>
  summarise(MAE = mean(MAE), sec = mean(sec), .groups = "drop") |>
  group_by(method) |> slice_min(MAE, n = 1, with_ties = FALSE) |> ungroup() |>
  mutate(family        = FAMILY[method],
         gain_vs_mean  = ref_mae - MAE,
         gain_per_sec  = ifelse(sec > 0.005, gain_vs_mean / sec, NA)) |>
  arrange(MAE)
cat("\n=== Computational cost and trade-off ===\n")
print(cost |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = 20)
write_csv(cost, file.path(MET, "computational_cost.csv"))


## =========================================================================
## 4. QUANTIFIED DECISION GUIDE
## =========================================================================
peaks    <- read_csv(file.path(MET, "peak_detection.csv"), show_col_types = FALSE)
gap_sens <- read_csv(file.path(MET, "gap_structure_sensitivity.csv"),
                     show_col_types = FALSE)

# The "no peak detected" figure is read from the data rather than hard-coded,
# so the guide can never drift away from the results it describes.
mean_sensitivity <- peaks |> filter(method == "Mean") |> pull(sensitivity) |> mean()

#  VARIANCE PRESERVATION: read from the data, never asserted.
#  An earlier version of this guide claimed that MICE was "the only family that
#  preserves variance (ratio ~1)". The metrics say otherwise, and the reason is
#  in functions.R: the two pmm variants return an AVERAGED completion
#  (imp_mice_stochastic ends with Reduce(`+`, ...) / m) and MICE-norm.predict
#  imputes the regression prediction. Both converge to, or estimate, the
#  conditional expectation and shrink variance by construction. What the
#  var_ratio measures here is therefore that ESTIMAND, not the ability of the
#  MICE family to preserve a distribution. The row below is computed.
#  CONVENTION (decided 2026-10-10): each method is scored under its BEST
#  predictor block, the block selected on MAE by 03 (selected_variant.csv).
#  This is the convention of Table 1 and of ranks_by_metric.csv. An earlier
#  version averaged var_ratio over the three blocks, which gave a different
#  value for the best MICE variant (0.892 instead of 0.904) than the text.
selected <- read_csv(file.path(MET, "selected_variant.csv"), show_col_types = FALSE)
var_pres <- metrics |>
  filter(gap_length == 1) |>
  inner_join(selected, by = c("method", "variant")) |>
  group_by(method, variant) |>
  summarise(var_ratio = mean(var_ratio), .groups = "drop") |>
  arrange(desc(var_ratio))
stopifnot(nrow(var_pres) == nrow(selected))
#  Display names: the guide is a paper table, so it must not mix internal
#  code labels ("SeasonalDecomp") with prose labels ("Seasonal decomposition").
#  The MICE labels also state the ESTIMAND actually scored, not an inferential
#  use the implementation does not perform.
DISPLAY <- c(SeasonalDecomp = "Seasonal decomposition",
             LinearInterp = "Linear interpolation", Kalman = "Kalman smoothing",
             MissForest = "MissForest", MissForest_ST = "MissForest-ST",
             KNN = "KNN", Mean = "Mean imputation",
             MICE_inference = "MICE-pmm (m=5)", MICE_precision = "MICE-pmm (m=20)",
             MICE_deterministic = "MICE-norm.predict")
best_var      <- unname(DISPLAY[var_pres$method[1]])
best_var_val  <- var_pres$var_ratio[1]
best_mice     <- var_pres |> filter(grepl("^MICE", method)) |> slice(1)
cat("\n--- Variance preservation, isolated gaps (measured) ---\n")
print(var_pres, n = 20)

#  Every number and every method name below is READ FROM THE RESULTS.
#  An earlier version hard-coded "moves from rank 7 to rank 1" and a row
#  "District weakly correlated with its peers / multivariate methods gain
#  almost nothing"; both became false when the series was extended to
#  2018-2025 (MissForest 7 -> 2, MissForest-ST 3 -> 1; gain vs mean
#  correlation with the other districts no longer significant). The district
#  row is therefore removed: the district analysis (n = 11) is exploratory
#  and does not support a recommendation.
rank_of <- function(m, v, col) gap_sens[[col]][gap_sens$method == m & gap_sens$variant == v]
long_gap <- gap_sens |> slice_min(rank_gap8, n = 1, with_ties = FALSE)
mf  <- c(rank_of("MissForest", "SeasonalTrend", "rank_gap1"),
         rank_of("MissForest", "SeasonalTrend", "rank_gap8"))
mfs <- c(rank_of("MissForest_ST", "SeasonalTrend", "rank_gap1"),
         rank_of("MissForest_ST", "SeasonalTrend", "rank_gap8"))
fastest   <- cost |> filter(method == "LinearInterp")

#  PEAK DETECTION BY GAP LENGTH (added 2026-10-12, step 8 of the revision).
#  peak_detection.csv (03) covers 1-week gaps only, so the former "peaks" row
#  of this guide said nothing about prolonged gaps. Youden index of the ten
#  methods, each under its retained predictor block (selected_variant.csv,
#  the convention of Table 1), for every gap length; 60 blocks (repetition x
#  rate) per cell. Written here, not in R/11, because the guide below needs
#  it and run_all.R runs 04 before 11.
youden_by_gap <- metrics |>
  inner_join(selected, by = c("method", "variant")) |>
  group_by(gap_length, method, variant) |>
  summarise(youden = mean(youden, na.rm = TRUE),
            sensitivity = mean(sensitivity, na.rm = TRUE),
            specificity = mean(specificity, na.rm = TRUE),
            MAE = mean(MAE), n_blocks = n(), .groups = "drop") |>
  group_by(gap_length) |>
  mutate(rank_youden = rank(-youden), rank_MAE = rank(MAE)) |>
  ungroup() |> arrange(gap_length, rank_youden)
stopifnot(all(youden_by_gap$n_blocks == N_BLOCKS))
write_csv(youden_by_gap, file.path(MET, "youden_by_gap.csv"))

#  Paired comparisons of the Youden index over the 60 blocks (same model as
#  paired_gap8.csv in R/11): (i) 1-week gaps, the two best methods on Youden;
#  (ii) 8-week gaps, the best method on Youden against the best on MAE.
#  Both pairs are read from youden_by_gap, never hard-coded.
paired_youden <- function(gap, m_a, m_b) {
  va <- selected$variant[selected$method == m_a]
  vb <- selected$variant[selected$method == m_b]
  d <- metrics |>
    filter(gap_length == gap,
           (method == m_a & variant == va) | (method == m_b & variant == vb)) |>
    select(rep, missing_rate, method, youden, MAE) |>
    pivot_wider(names_from = method, values_from = c(youden, MAE))
  ya <- d[[paste0("youden_", m_a)]]; yb <- d[[paste0("youden_", m_b)]]
  ma <- d[[paste0("MAE_", m_a)]];    mb <- d[[paste0("MAE_", m_b)]]
  stopifnot(nrow(d) == N_BLOCKS, !anyNA(ya), !anyNA(yb))
  wt <- suppressWarnings(wilcox.test(ya, yb, paired = TRUE))
  nz <- sum(ya != yb)
  bt <- binom.test(sum(ya > yb), nz, p = 0.5)
  tibble(gap_length = gap, method_a = m_a, variant_a = va,
         method_b = m_b, variant_b = vb, n_blocks = nrow(d),
         youden_a = mean(ya), youden_b = mean(yb), mean_diff_youden = mean(ya - yb),
         blocks_a_higher = sum(ya > yb), blocks_tied = sum(ya == yb),
         wilcoxon_V = unname(wt$statistic), wilcoxon_p = wt$p.value,
         sign_test_p = bt$p.value,
         mae_a = mean(ma), mae_b = mean(mb))
}
yg1 <- youden_by_gap |> filter(gap_length == 1)
yg8 <- youden_by_gap |> filter(gap_length == 8)
pair1 <- paired_youden(1, yg1$method[1], yg1$method[2])
pair8 <- paired_youden(8, yg8$method[1],
                       yg8 |> slice_min(rank_MAE, n = 1, with_ties = FALSE) |> pull(method))
paired_youden_tbl <- bind_rows(pair1, pair8)
write_csv(paired_youden_tbl, file.path(MET, "paired_youden.csv"))
cat("\n--- Peak detection (Youden) by gap length ---\n")
print(youden_by_gap |> filter(rank_youden <= 4), n = 20, width = Inf)
print(paired_youden_tbl, width = Inf)

#  SITUATIONS (revised 2026-10-19, step 8b of the revision, batch 7).
#  The guide is written to match Table 5 of the paper, row by row. Every
#  number is read from a result file or computed here from one; the method
#  names that are not read from the results are fixed by the text of the
#  paper.
#  Row 2 quotes the reduced mechanism design (Table 4), so this script reads
#  the outputs of 05 and 07: run_all.R runs 04 after them (and before 11,
#  which reads youden_by_gap.csv written above).
MMET <- file.path("output_mechanism", "metrics")
fmt_p <- function(p) if (p < 0.001) "p < 0.001" else sprintf("p = %.2f", p)

#  (i) Row 1. Kalman smoothing against linear interpolation, 1-week gaps,
#  paired over the 60 blocks (repetition x rate). Same computation as
#  paired_kalman_linearinterp.csv in R/11: the MAE gap is expressed
#  relative to the MAE of the method that trails.
kl <- metrics |>
  filter(gap_length == 1, method %in% c("Kalman", "LinearInterp")) |>
  inner_join(selected, by = c("method", "variant")) |>
  select(rep, missing_rate, method, MAE) |>
  pivot_wider(names_from = method, values_from = MAE)
stopifnot(nrow(kl) == N_BLOCKS, !anyNA(kl$Kalman), !anyNA(kl$LinearInterp))
kal_leads <- mean(kl$Kalman) < mean(kl$LinearInterp)
kl_lead   <- if (kal_leads) kl$Kalman else kl$LinearInterp
kl_trail  <- if (kal_leads) kl$LinearInterp else kl$Kalman
kl_won    <- sum(kl_lead < kl_trail)
kl_pct    <- 100 * (mean(kl_trail) - mean(kl_lead)) / mean(kl_trail)

#  (ii) Row 2. MAE of the two forests at 8 weeks, seasonal and trend
#  predictors (the block quoted in the row, as for mf and mfs above).
forest_mae8 <- metrics |>
  filter(gap_length == 8, method %in% c("MissForest_ST", "MissForest"),
         variant == "SeasonalTrend") |>
  group_by(method) |> summarise(MAE = mean(MAE), n = n(), .groups = "drop")
stopifnot(nrow(forest_mae8) == 2, all(forest_mae8$n == N_BLOCKS))
mae8_of <- function(m) forest_mae8$MAE[forest_mae8$method == m]

#  (iii) Row 2. Reduced mechanism design, MNAR, 8-week gaps: mean ranks of
#  the two forests (ranks_four_mechanisms.csv, written by 07) and Nemenyi
#  critical difference with N = number of repetitions read from the metric
#  file (same formula as section 2 above and as R/11, section 3). The gap
#  compared with the CD is the distance from the lower-ranked forest to the
#  best method under MNAR at 8 weeks.
r4 <- read_csv(file.path(MMET, "ranks_four_mechanisms.csv"), show_col_types = FALSE)
n_mnar <- read_csv(file.path(MMET, "metrics_mechanism.csv"), show_col_types = FALSE) |>
  filter(mechanism == "MNAR", gap_length == 8) |>
  distinct(rep) |> nrow()
k_mech  <- nrow(r4)
cd_mnar <- qtukey(0.95, k_mech, Inf) / sqrt(2) *
           sqrt(k_mech * (k_mech + 1) / (6 * n_mnar))
mnar_of  <- function(m) r4$MNAR_8[r4$method == m]
mnar_gap <- max(mnar_of("MissForest_ST"), mnar_of("MissForest")) - min(r4$MNAR_8)
cat(sprintf("\nMNAR, 8-week gaps: N = %d, CD = %.2f, gap to the best = %.2f\n",
            n_mnar, cd_mnar, mnar_gap))

#  (iv) Row 5. Time per imputed dataset, 1-week gaps, best block of each
#  method (computational_cost.csv, section 3 above).
sec_of <- function(m) cost$sec[cost$method == m]

guide <- tibble(
  situation = c(
    "Isolated gaps (1 week)",
    "Prolonged gaps (8 weeks), accuracy of counts",
    "Priority on epidemic peak detection",
    "Distributional fidelity required (variance preserved), isolated gaps",
    "Modest hardware or limited technical skills",
    "To be avoided in all circumstances"),
  method = c(
    "Linear interpolation or Kalman smoothing",
    "MissForest-ST or MissForest, with seasonal and trend predictors",
    # tie at 1 week when the paired test does not separate the two best
    sprintf("%s (isolated gaps); %s (8-week gaps)",
            if (pair1$wilcoxon_p > 0.05)
              paste(unname(DISPLAY[pair1$method_a]), "or", unname(DISPLAY[pair1$method_b]))
            else unname(DISPLAY[pair1$method_a]),
            unname(DISPLAY[pair8$method_a])),
    best_var,
    "Linear interpolation (isolated gaps); Seasonal decomposition (8-week gaps)",
    "Mean imputation"),
  rationale = c(
    sprintf(paste("mean ranks %.2f and %.2f; %s beats %s in %d of %d blocks,",
                  "with an MAE gap of %.1f%%"),
            ranks$mean_rank[ranks$method == "LinearInterp"],
            ranks$mean_rank[ranks$method == "Kalman"],
            if (kal_leads) "Kalman smoothing" else "linear interpolation",
            if (kal_leads) "interpolation" else "Kalman smoothing",
            kl_won, nrow(kl), kl_pct),
    sprintf(paste("from 1 to 8 weeks, rank %d -> %d (MissForest-ST) and",
                  "%d -> %d (MissForest) of 24 combinations, with close errors",
                  "(%.1f and %.1f); under MNAR, descriptive ranks %.1f and %.1f",
                  "over %d repetitions, gap %s the CD; mechanisms: Table 4;",
                  "detection: next row"),
            mfs[1], mfs[2], mf[1], mf[2],
            mae8_of("MissForest_ST"), mae8_of("MissForest"),
            mnar_of("MissForest_ST"), mnar_of("MissForest"), n_mnar,
            if (mnar_gap < cd_mnar) "below" else "above"),
    sprintf(paste("isolated gaps: Youden index %.3f and %.3f (paired Wilcoxon, %s);",
                  "8-week gaps: %.3f against %.3f for %s (higher in %d of %d",
                  "blocks, %s)"),
            pair1$youden_a, pair1$youden_b, fmt_p(pair1$wilcoxon_p),
            pair8$youden_a, pair8$youden_b, unname(DISPLAY[pair8$method_b]),
            pair8$blocks_a_higher, pair8$n_blocks, fmt_p(pair8$wilcoxon_p)),
    sprintf(paste("variance ratio %.3f, closest to 1, computed over the eleven",
                  "districts pooled, so it includes the level differences",
                  "between districts; the best MICE variant reaches only %.3f,",
                  "because the two pmm variants average their completions and",
                  "MICE-norm.predict imputes the regression prediction, both of",
                  "which shrink the variance"),
            best_var_val, best_mice$var_ratio[1]),
    sprintf(paste("%.2f and %.2f s per imputed dataset, no hyperparameter to",
                  "tune, same R package (imputeTS)"),
            sec_of("LinearInterp"), sec_of("SeasonalDecomp")),
    sprintf("sensitivity = %.3f: no epidemic peak detected", mean_sensitivity))
)
cat("\n=== Decision guide ===\n"); print(guide, n = 20, width = Inf)
write_csv(guide, file.path(MET, "decision_guide.csv"))


## =========================================================================
## 5. PAPER FIGURES
## =========================================================================
pal <- c(Naive = "#d95f02", `Cross-sectional` = "#7570b3",
         Temporal = "#1b9e77", Hybrid = "#e7298a")

# Fig 2 — MAE by method and missingness rate
f2 <- metrics |> filter(gap_length == 1) |>
  group_by(method, variant, missing_rate) |>
  summarise(m = mean(MAE), s = sd(MAE), .groups = "drop") |>
  group_by(method, missing_rate) |> slice_min(m, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(family = FAMILY[method], method = fct_reorder(method, m))
p2 <- ggplot(f2, aes(m, method, colour = family)) +
  # geom_errorbarh is deprecated since ggplot2 4.0: use orientation = "y"
  geom_errorbar(aes(xmin = m - s, xmax = m + s), orientation = "y",
                width = .25, linewidth = .4) +
  geom_point(size = 2) +
  facet_wrap(~ missing_rate, labeller = label_both) +
  scale_colour_manual(values = pal) +
  labs(x = "MAE (mean ± SD over repetitions)", y = NULL, colour = NULL,
       title = "Imputation error by method and missingness rate",
       subtitle = "Isolated missingness (one-week gaps)") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig2_mae_by_rate.png"), p2, width = 10, height = 5, dpi = 300)

# Fig 3 — accuracy versus detection: the trade-off
p3 <- peaks |> mutate(family = FAMILY[method]) |>
  ggplot(aes(MAE, youden, colour = family)) +
  geom_point(size = 2.4) +
  ggrepel::geom_text_repel(aes(label = paste0(method, " (", variant, ")")),
                           size = 2.7, max.overlaps = 30, show.legend = FALSE) +
  scale_colour_manual(values = pal) +
  labs(x = "MAE (lower = better)", y = "Youden index (higher = better)",
       colour = NULL, title = "Average accuracy versus alert capability",
       subtitle = "Each point is a method x variant; the two criteria do not rank alike") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig3_accuracy_vs_detection.png"), p3,
       width = 8.5, height = 6, dpi = 300)

# Fig 4 — the reversal driven by the missingness structure
p4 <- gap_sens |>
  pivot_longer(starts_with("rank_gap"), names_to = "gap", values_to = "rank") |>
  mutate(gap = as.integer(str_remove(gap, "rank_gap")),
         family = FAMILY[method],
         key = paste(method, variant)) |>
  ggplot(aes(factor(gap), rank, group = key, colour = family)) +
  geom_line(linewidth = .5, alpha = .8) + geom_point(size = 1.6) +
  scale_y_reverse(breaks = c(1, 5, 10, 15, 20, 24)) +
  scale_colour_manual(values = pal) +
  labs(x = "Gap length (weeks)", y = "Rank (1 = best)", colour = NULL,
       title = "Gap length changes the leading methods",
       subtitle = "Same missingness rate, isolated versus clustered gaps") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig4_reversal.png"), p4, width = 7.5, height = 5.5, dpi = 300)

# Fig 5 — decomposition of the period effect
pe <- read_csv(file.path(MET, "period_effect.csv"), show_col_types = FALSE)
p5 <- pe |> filter(gap_length %in% c(1, 8), method != "Mean") |>
  group_by(method, gap_length) |>
  summarise(Seasonality = mean(gain_seasonal_pct),
            Trend       = mean(gain_trend_pct), .groups = "drop") |>
  pivot_longer(c(Seasonality, Trend), names_to = "component", values_to = "gain") |>
  ggplot(aes(gain, fct_reorder(method, gain, sum), fill = component)) +
  geom_col() +
  facet_wrap(~ gap_length,
             labeller = labeller(gap_length = \(x) paste(x, "week gaps"))) +
  scale_fill_manual(values = c(Seasonality = "#1b9e77", Trend = "#7570b3")) +
  labs(x = "MAE gain from the period covariates (%)", y = NULL, fill = NULL,
       title = "Where does the gain from the period covariates come from?",
       subtitle = "Gain from seasonal and trend predictors, by gap length") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig5_period_effect.png"), p5, width = 9, height = 5, dpi = 300)

# Fig 6 — heterogeneity between districts
p6 <- gain |>
  ggplot(aes(mean_corr_other_districts, gain_MissForest_pct)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70", linewidth = .5) +
  geom_point(aes(size = cv), colour = "#1b4965") +
  ggrepel::geom_text_repel(aes(label = district), size = 3) +
  labs(x = "Mean correlation of the district with its peers",
       y = "Gain of MissForest over mean imputation (%)",
       size = "Coefficient\nof variation",
       title = "Where is multivariate imputation useful?",
       subtitle = "Exploratory (n = 11 districts); no association passes the Bonferroni threshold") +
  theme_minimal(base_size = 10)
ggsave(file.path(FIG, "fig6_districts.png"), p6, width = 8, height = 6, dpi = 300)

cat("\n=== Final analyses complete. Figures in", FIG, "===\n")
cat("Files produced: district_analysis.csv, gain_by_district.csv,\n",
    "  computational_cost.csv, youden_by_gap.csv, paired_youden.csv,\n",
    "  decision_guide.csv + 6 figures.\n")
