###############################################################################
#  SHARED MACHINERY FOR THE MISSINGNESS-MECHANISM SENSITIVITY ANALYSIS
#
#  Sourced by 05, 06 and 07. It defines the propensity weights, the mask
#  sampler under those weights, and the mechanism diagnostics.
#
#  In the French version, scripts 06 and 07 recovered these definitions by
#  reading the pipeline file line by line and eval()-ing a prefix of it:
#      src <- readLines("01_pipeline_imputation.R")
#      eval(parse(text = src[1:(grep("BOUCLE PRINCIPALE", src) - 1)]))
#  That worked, but renaming a section heading silently broke it. Sourcing
#  two ordinary files is what a reader of the repository expects.
###############################################################################

source("R/functions.R")

DOS_MECH <- "output_mechanism"
for (d in c("masks", "metrics", "figures"))
  dir.create(file.path(DOS_MECH, d), recursive = TRUE, showWarnings = FALSE)

## -------------------------------------------------------------------------
## 1. THE MECHANISMS
## -------------------------------------------------------------------------
#  Rubin's definitions, applied to the letter:
#    MCAR : P(missing) depends on nothing.
#    MAR  : P(missing) depends on OBSERVED variables.
#    MNAR : P(missing) depends on the MISSING value itself.
#
#  The classic trap: building a "MAR" from the previous week's value. If that
#  week is itself masked, the propensity depends on an unobserved quantity --
#  and it is no longer MAR. We therefore rely on variables that are NEVER
#  punched with holes by construction:
#    (a) the rainy season -> a period covariate, never masked
#    (b) the district     -> known for every cell
#  Epidemiological reading: during the rainy season tracks are impassable and
#  health facilities report less well; and some districts structurally have
#  weaker reporting capacity. This is MAR in the strict sense, and plausible.
#
#  MAR_strong additionally uses the values of the OTHER districts in the same
#  week (observed, hence still MAR) -- see 07_mar_strong_arm.R for why.

FORCE_SEASON   <- log(3)   # a rainy-season cell is 3x more at risk
FORCE_DISTRICT <- log(2)   # SD of the district effect, in log-odds
FORCE_VALUE    <- log(2)   # MNAR: +1 SD of cases -> 2x more at risk
FORCE_NEIGHBOUR <- NA_real_  # MAR_strong: calibrated in 07

#  District effect DRAWN ONCE and frozen: it is a property of the design, not
#  replication noise. Without this, the district effect would average to zero
#  across the 10 repetitions and the mechanism would reduce to a season effect.
set.seed(20240614)
U_DISTRICT <- as.numeric(scale(rnorm(n_dis)))
names(U_DISTRICT) <- districts

Z_REF  <- scale(REF)                              # per-district z-scores, MNAR
NEIGH  <- (rowSums(REF) - REF) / (n_dis - 1)      # mean of the other districts
Z_NEIGH <- scale(NEIGH)                           # per-district z-scores, MAR_strong

mechanism_weights <- function(mechanism) {
  logw <- matrix(0, n_obs, n_dis, dimnames = list(NULL, districts))
  if (mechanism %in% c("MAR", "MNAR", "MAR_strong")) {
    logw <- logw + FORCE_SEASON * PERIOD_SEASONAL$rainy_season
    logw <- logw + rep(FORCE_DISTRICT * U_DISTRICT, each = n_obs)
  }
  if (mechanism == "MNAR")       logw <- logw + FORCE_VALUE     * as.numeric(Z_REF)
  if (mechanism == "MAR_strong") logw <- logw + FORCE_NEIGHBOUR * as.numeric(Z_NEIGH)
  exp(logw - max(logw))          # -max: avoids any numerical overflow
}

## -------------------------------------------------------------------------
## 2. THE MASK SAMPLER
## -------------------------------------------------------------------------
#  Sampling is WITHOUT REPLACEMENT with these weights, not independent
#  Bernoulli draws. Reason: this guarantees exactly the same NUMBER of masked
#  cells as under MCAR. Otherwise we would be changing both the mechanism AND
#  the rate, and could no longer attribute the observed difference.
#
#  THE ROW CONSTRAINT, AND WHY IT IS NOT A HACK.
#  A strong MAR built on spatial covariation concentrates missingness into
#  ENTIRE WEEKS: NEIGH[i,j] is the mean of the other ten districts in week i,
#  a quantity dominated by the week's overall level, so it is very nearly a
#  ROW variable. At high FORCE_NEIGHBOUR the sampler empties whole weeks.
#  A week in which all 11 districts are missing is not an imputation problem:
#  no cross-sectional method has a donor, and the question falls outside the
#  compared perimeter. We therefore cap the number of missing districts per
#  week. The counter below reports how often the cap actually bound.
#
#  DEFAULT = n_dis - 1L, which is exactly the guard the main pipeline already
#  enforces (rowSums(!M) > 0). With that default, and for a given seed, this
#  sampler produces EXACTLY the same masks as an unconstrained
#  sample.int(N, n_na, prob = w): drawing k items is the first k steps of the
#  same sequential procedure. Script 07 tightens it to n_dis - 3L.
MAX_MISSING_PER_WEEK <- n_dis - 1L
MAX_MISSING_PER_DISTRICT <- n_obs - 3L
N_CONSTRAINT_HITS <- 0L

simulate_mask_mechanism <- function(rate, seed, gap_length = 1,
                                    mechanism = "MAR") {
  set.seed(seed)
  W     <- mechanism_weights(mechanism)
  M     <- matrix(FALSE, n_obs, n_dis, dimnames = list(NULL, districts))
  n_na  <- round(n_obs * n_dis * rate / 100)
  hits  <- 0L

  if (gap_length == 1) {
    nl <- integer(n_obs); nc <- integer(n_dis); n <- 0L
    order_idx <- sample.int(n_obs * n_dis, n_obs * n_dis, prob = as.numeric(W))
    for (s in order_idx) {
      i <- ((s - 1L) %%  n_obs) + 1L
      j <- ((s - 1L) %/% n_obs) + 1L
      if (nl[i] >= MAX_MISSING_PER_WEEK || nc[j] >= MAX_MISSING_PER_DISTRICT) {
        hits <- hits + 1L; next
      }
      M[s] <- TRUE; nl[i] <- nl[i] + 1L; nc[j] <- nc[j] + 1L; n <- n + 1L
      if (n >= n_na) break
    }
    stopifnot(n == n_na)
  } else {
    starts <- expand.grid(i = seq_len(n_obs - gap_length), j = seq_len(n_dis))
    w_st   <- W[cbind(starts$i, starts$j)]
    guard  <- 0L
    while (sum(M) < n_na && guard < 50000L) {
      guard <- guard + 1L
      s  <- sample.int(nrow(starts), 1, prob = w_st)
      i0 <- starts$i[s]; j0 <- starts$j[s]
      ii <- i0:(i0 + gap_length - 1L)
      fresh <- ii[!M[ii, j0]]
      if (!length(fresh)) next
      if (any(rowSums(M[fresh, , drop = FALSE]) + 1L > MAX_MISSING_PER_WEEK) ||
          sum(M[, j0]) + length(fresh) > MAX_MISSING_PER_DISTRICT) {
        hits <- hits + 1L; next
      }
      M[ii, j0] <- TRUE
    }
    stopifnot(guard < 50000L)
  }
  N_CONSTRAINT_HITS <<- N_CONSTRAINT_HITS + hits
  stopifnot(all(colSums(!M) > 2), all(rowSums(!M) > 0))
  M
}

## -------------------------------------------------------------------------
## 3. MECHANISM DIAGNOSTICS — the positive control
## -------------------------------------------------------------------------
#  Here the tests of validate_mcar() are used IN REVERSE.
#  Under MCAR we wanted p > 0.05 (nothing distinguishes the missingness).
#  Under MAR/MNAR we want p < 0.05: if the tests do NOT reject, the programmed
#  mechanism is too weak to be detectable and the sensitivity analysis tests
#  nothing at all. This control is what separates a real sensitivity analysis
#  from a decorative figure.
#
#  IMPORTANT CORRECTION. The "p_time" test of validate_mcar() aggregates
#  missingness BY YEAR. The MAR mechanism programmed here is SEASONAL, hence
#  intra-annual: it is invisible to a year-aggregated test. Using p_time as
#  the positive control leads to wrongly concluding that the mechanism is
#  weak. The matching test is the 2x2 table missingness x rainy season, added
#  below.
diagnose_mechanism <- function(M, mechanism) {
  v  <- validate_mcar(M)
  Mv <- M & !EXCLUDE
  Z  <- scale(REF)
  tb <- table(rep(PERIOD_SEASONAL$rainy_season, n_dis), as.vector(M))
  or_s <- if (all(dim(tb) == c(2, 2)))
    (tb[2, 2] * tb[1, 1]) / (tb[2, 1] * tb[1, 2]) else NA_real_
  p_s  <- if (all(dim(tb) == c(2, 2))) chisq.test(tb)$p.value else NA_real_
  tibble(mechanism      = mechanism,
         p_season       = p_s,           # <- the real positive control
         or_season      = or_s,
         p_districts    = v$p_districts,
         p_time_annual  = v$p_time,      # kept, but NOT relevant here
         p_values_ks    = v$p_values_ks,
         # how much higher are the masked cells?
         mean_ratio     = mean(REF[Mv]) / mean(REF[!Mv]),
         # selection strength on the value: ~0 under strict MAR, >0 under MNAR
         cohen_d        = (mean(Z[Mv]) - mean(Z[!Mv])) / sd(Z[!Mv]),
         pct_missing_in_rainy =
           100 * mean(PERIOD_SEASONAL$rainy_season[row(M)[M]] == 1),
         detected       = p_s < 0.05 | v$p_values_ks < 0.05)
}

## -------------------------------------------------------------------------
## 4. REDUCED SCOPE FOR THE SENSITIVITY ANALYSIS
## -------------------------------------------------------------------------
#  We do NOT redo all 4320 imputations. A sensitivity analysis must target
#  what could break the conclusion:
#    rate = 10 %        -> the median rate, the one most often reported
#    gap = 1 AND 8      -> ESSENTIAL. The most original result of this study
#                          is the ranking reversal with gap length (Kalman
#                          1st -> 8th, MissForest 7th -> 1st). If that
#                          reversal does not survive a change of mechanism,
#                          the central contribution falls.
#    n_rep = 10         -> enough for a Friedman test on 10 methods
#                          (N = 10 blocks, k = 10: CD stays usable, though
#                          large -- see the caveat in the results document)
CFG_MECH <- list(
  missing_rate = 10,
  gap_lengths  = c(1, 8),
  n_rep        = 10,
  tuning_reps  = 3,
  tuning_B     = 5
)

cat("Mechanism machinery loaded.\n")
