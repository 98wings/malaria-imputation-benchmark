###############################################################################
#  06 — CORRECTED DIAGNOSTICS OF THE MISSINGNESS MECHANISMS
#
#  WHY THIS SCRIPT EXISTS
#  The positive control in 05 originally used the "p_time" test of
#  validate_mcar(). That test aggregates missingness BY YEAR:
#      obs <- tapply(rowSums(M), year, sum)
#  It therefore asks whether missingness is uniform from one YEAR to the next.
#  The programmed MAR mechanism is SEASONAL, i.e. INTRA-annual: it is entirely
#  invisible to a year-aggregated test. Hence an apparent "60% detection" for
#  MAR at gap 1 -- not a weak mechanism, simply the wrong test.
#
#  This script recomputes the diagnostics with the tests that actually match
#  the mechanisms, FROM THE MASKS ALREADY SAVED ON DISK. No imputation is
#  redone. Runtime: seconds.
#
#  It includes MCAR as a NEGATIVE CONTROL: the same tests must be
#  non-significant there. That is what makes the table publishable.
#
#  Run it again after 07 to pick up the MAR_strong masks.
###############################################################################

source("R/functions_mechanism.R")

GAPS  <- CFG_MECH$gap_lengths
REPS  <- seq_len(CFG_MECH$n_rep)
RATE  <- CFG_MECH$missing_rate

RAINY <- PERIOD_SEASONAL$rainy_season
cat(sprintf("Share of calendar weeks in the rainy season: %.1f %%\n",
            100 * mean(RAINY)))

## -------------------------------------------------------------------------
##  The tests, chosen to match what is actually programmed
## -------------------------------------------------------------------------
diagnose <- function(M, mechanism, rep, gap) {
  Mv        <- M & !EXCLUDE
  vec_rainy <- rep(RAINY, n_dis)
  vec_miss  <- as.vector(M)

  # (1) SEASON TEST — the one that was missing. 2x2 table: is missingness
  #     more frequent during the rainy season? That IS the definition of the
  #     programmed MAR, so this is the test that must reject.
  tb <- table(rainy = vec_rainy, missing = vec_miss)
  or_season <- if (all(dim(tb) == c(2, 2)))
    (tb[2, 2] * tb[1, 1]) / (tb[2, 1] * tb[1, 2]) else NA_real_
  p_season  <- if (all(dim(tb) == c(2, 2))) chisq.test(tb)$p.value else NA_real_

  # (2) DISTRICT TEST — the other component of the programmed MAR.
  p_district <- chisq.test(colSums(M))$p.value

  # (3) VALUE TESTS — these separate MAR from MNAR.
  #     Under strict MAR the masked values are NOT chosen for their value:
  #     any residual gap comes only from the correlation between season and
  #     incidence. Under MNAR the selection is direct and strong.
  Z    <- scale(REF)
  p_ks <- suppressWarnings(ks.test(Z[Mv], Z[!Mv])$p.value)
  p_mw <- suppressWarnings(wilcox.test(Z[Mv], Z[!Mv])$p.value)

  tibble(mechanism = mechanism, rep = rep, gap_length = gap,
         pct_rainy  = 100 * mean(vec_rainy[vec_miss] == 1),
         or_season  = or_season,
         p_season   = p_season,
         p_district = p_district,
         mean_ratio = mean(REF[Mv]) / mean(REF[!Mv]),
         cohen_d    = (mean(Z[Mv]) - mean(Z[!Mv])) / sd(Z[!Mv]),
         p_ks       = p_ks,
         p_mw       = p_mw)
}

## -------------------------------------------------------------------------
##  Read the masks already on disk
## -------------------------------------------------------------------------
res <- list()
for (gap in GAPS) for (rep in REPS) {
  # Negative control: the MCAR masks from the main run, same perimeter
  f <- file.path("output", "masks",
                 sprintf("mask_rep%02d_rate%02d_gap%d.rds", rep, RATE, gap))
  if (file.exists(f))
    res[[length(res) + 1]] <- diagnose(readRDS(f), "MCAR", rep, gap)

  # MAR_strong exists only after 07: missing files are ignored.
  for (mech in c("MAR", "MAR_strong", "MNAR")) {
    f <- file.path(DOS_MECH, "masks",
                   sprintf("mask_%s_rep%02d_gap%d.rds", mech, rep, gap))
    if (file.exists(f))
      res[[length(res) + 1]] <- diagnose(readRDS(f), mech, rep, gap)
  }
}
diag <- bind_rows(res)
stopifnot(nrow(diag) > 0)
write_csv(diag, file.path(DOS_MECH, "metrics", "mechanism_diagnostics.csv"))

## -------------------------------------------------------------------------
##  Summary table — intended for the supplementary material
## -------------------------------------------------------------------------
summary_tab <- diag |>
  group_by(mechanism, gap_length) |>
  summarise(n = n(),
            pct_rainy = mean(pct_rainy),
            or_season = mean(or_season),
            pct_reject_season   = 100 * mean(p_season   < 0.05),
            pct_reject_district = 100 * mean(p_district < 0.05),
            mean_ratio = mean(mean_ratio),
            cohen_d    = mean(cohen_d),
            pct_reject_ks = 100 * mean(p_ks < 0.05),
            pct_reject_mw = 100 * mean(p_mw < 0.05),
            .groups = "drop") |>
  arrange(gap_length,
          factor(mechanism, levels = c("MCAR", "MAR", "MAR_strong", "MNAR")))

## WARNING — THE p-VALUES ARE INVALID FOR gap_length > 1.
#  The tests assume independent cells. With 8-week gaps there are only about
#  36 independent blocks, not 286 cells: the effective sample size is divided
#  by 8 and the tests over-reject massively. The proof is in the table itself:
#  MCAR at gap 8 rejects the season test in ~70% of replications even though
#  or_season = 1.13 and cohen_d = -0.007, i.e. no real effect at all.
#  => For gap > 1, read ONLY the effect sizes (or_season, cohen_d,
#     mean_ratio), never the p-values. It was the MCAR negative control that
#     revealed this invalidity: without it, those rejections would have been
#     mistaken for a validation of the mechanism.

cat("\n=== CORRECTED MECHANISM DIAGNOSTICS ===\n")
cat("(pct_reject_* = percentage of the", length(REPS),
    "repetitions where the test rejects at 5%)\n\n")
print(as.data.frame(summary_tab), row.names = FALSE, digits = 3)

cat("\nEXPECTED READING\n")
cat("  MCAR       : all pct_reject close to 5%, or_season ~ 1, mean_ratio ~ 1\n")
cat("  MAR        : pct_reject_season = 100%, or_season ~ 3, BUT mean_ratio\n")
cat("               close to 1 -> missingness is DISPLACED in time without\n")
cat("               being selected for its value. That is exactly MAR.\n")
cat("  MAR_strong : like MAR for the season, but cohen_d matched to MNAR.\n")
cat("  MNAR       : pct_reject_season = 100% AND mean_ratio clearly > 1,\n")
cat("               cohen_d clearly > 0 -> direct selection on the missing value.\n")
cat("\n  For gap > 1: read the EFFECT SIZES only. See the comment above.\n")

## -------------------------------------------------------------------------
##  SELECTION STRENGTH — the number a reviewer will look for
## -------------------------------------------------------------------------
#  cohen_d measures how many standard deviations higher the masked cells are
#  than the observed ones. It is the only honest way to answer "was your MAR
#  strong enough to have been able to change the conclusions?".
cat("\n=== SELECTION STRENGTH ON THE VALUE (Cohen's d) ===\n")
print(as.data.frame(summary_tab |>
                      select(mechanism, gap_length, cohen_d, mean_ratio)),
      row.names = FALSE, digits = 3)
cat("\nIf cohen_d(MAR) is close to 0, the robustness observed under MAR proves\n",
    "little: the mechanism could not have broken the ranking. Either say so\n",
    "explicitly in the paper, or add a more selective MAR arm -- which is\n",
    "exactly what 07_mar_strong_arm.R does.\n", sep = "")

cat("\n=== Done ===\n")
cat("Output: ", file.path(DOS_MECH, "metrics", "mechanism_diagnostics.csv"),
    "\n", sep = "")
