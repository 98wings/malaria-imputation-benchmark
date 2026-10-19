###############################################################################
#  SHARED SETUP AND FUNCTIONS
#
#  Sourced by 02_run_benchmark.R and by every mechanism-sensitivity script.
#  It defines, in this order:
#     0. packages and configuration
#     1. the reference data and the period covariates
#     2. mask simulation and MCAR validation
#     3. the ten imputation methods
#     4. honest hyperparameter tuning
#     5. the evaluation metrics
#
#  WHY THIS FILE EXISTS
#  The sensitivity scripts need exactly the same imputation code as the main
#  benchmark. Copying the functions into each script guarantees that the
#  copies eventually diverge. Sourcing one file is what a reviewer means by
#  "the same code was used throughout".
#
#  Sourcing this file runs nothing expensive: it loads the reference,
#  builds the covariates and defines functions. Roughly one second.
###############################################################################

## ---------------------------------------------------------------- 0. SETUP
pkgs <- c("tidyverse", "lubridate", "mice", "VIM", "missRanger", "imputeTS",
          "forecast", "PMCMRplus", "tictoc")
missing_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(missing_pkgs)) install.packages(missing_pkgs)
# ORDER MATTERS: mice masks filter(), cbind() and rbind(). tidyverse is
# therefore loaded LAST so that dplyr::filter keeps precedence.
invisible(lapply(setdiff(pkgs, c("tidyverse", "lubridate")), library,
                 character.only = TRUE))
library(lubridate); library(tidyverse)

set.seed(2024)

CFG <- list(
  # ONE reference only: 417 ISO weeks, 2018-2025, assembled by
  # data/extend_reference_2018_2025.py and certified by 01_build_reference.R.
  # Never point at the raw extractions in data/raw/.
  # (The article's first version used reference_adamaoua_weekly_2018_2022.csv,
  #  260 rows, in which weeks 52 and 53 of 2020 were merged.)
  reference_file    = "data/derived/reference_adamaoua_weekly_2018_2025.csv",
  flagged_cells_file = "data/derived/reference_flagged_cells.csv",
  output_dir        = "output",
  missing_rates     = c(5, 10, 20),
  n_rep             = 20,          # 20 rather than 10: the cost is machine
                                   # time, not analyst time
  gap_lengths       = c(1, 4, 8),  # missingness structure: isolated / 4wk / 8wk
  tuning_reps       = 5,           # k is tuned on the first 5 reps only
  tuning_B          = 10,          # internal re-masking rounds per tuning
  forecast_reps     = 5,           # reps used for downstream forecast checks
  variants          = c("NoTime", "Seasonal", "SeasonalTrend"),
  date_col          = "periodid"
)
# >>> FOR A DRY RUN: set n_rep = 2 above, check the output, then restore 20.

for (d in c("masks", "imputed", "metrics", "figures"))
  dir.create(file.path(CFG$output_dir, d), recursive = TRUE, showWarnings = FALSE)

# --- Guard against STALE OUTPUTS ------------------------------------------
# Files left by a run on another reference (e.g. the 260-week series) have a
# different number of rows. Scripts 03 and 06 read whatever is on disk, so a
# stale file would be mixed silently with the new results, or crash later.
check_no_stale_outputs <- function(n_rows) {
  masks <- c(list.files(file.path(CFG$output_dir, "masks"), "\\.rds$", full.names = TRUE),
             list.files(file.path("output_mechanism", "masks"), "\\.rds$", full.names = TRUE))
  imputed <- list.files(file.path(CFG$output_dir, "imputed"), "\\.csv$", full.names = TRUE)
  bad <- c(masks[vapply(masks, function(f) nrow(readRDS(f)) != n_rows, logical(1))],
           imputed[vapply(imputed, function(f)
             length(readLines(f, warn = FALSE)) - 1L != n_rows, logical(1))])
  if (length(bad))
    stop(length(bad), " output file(s) come from a run on a different reference ",
         "(expected ", n_rows, " rows), e.g. ", bad[1],
         ".\nMove the folders 'output' and 'output_mechanism' away, then rerun ",
         "from 02_run_benchmark.R.", call. = FALSE)
}


## ------------------------------------------------------- 1. REFERENCE DATA
raw_ref <- read_csv(CFG$reference_file, show_col_types = FALSE)

dates     <- as.Date(raw_ref[[CFG$date_col]], format = "%Y-%m-%d")
districts <- setdiff(names(raw_ref), CFG$date_col)
REF       <- as.matrix(raw_ref[, districts])       # 417 x 11, no NA
n_obs     <- nrow(REF); n_dis <- ncol(REF)

stopifnot(!anyNA(REF), !any(is.na(dates)), all(diff(dates) == 7),
          n_dis == 11,
          # every ISO week of the period, none missing or merged
          n_obs == sum(isoweek(as.Date(paste0(isoyear(min(dates)):isoyear(max(dates)),
                                              "-12-28")))))
cat(sprintf("Reference: %d weeks x %d districts, from %s to %s\n",
            n_obs, n_dis, min(dates), max(dates)))
check_no_stale_outputs(n_obs)

# --- Cells EXCLUDED from the metrics (see 01_build_reference.R) ------------
# Belel 2019-09-01 was already a mean imputation in the source data: it is not
# ground truth. It is removed from every evaluation, while remaining available
# as a predictor and inside the reference series used for forecasting.
EXCLUDE <- matrix(FALSE, n_obs, n_dis, dimnames = list(NULL, districts))
# The file is REQUIRED: without it the reconstructed Belel cell would silently
# be scored as ground truth. Run 01_build_reference.R first.
stopifnot(file.exists(CFG$flagged_cells_file))
flagged  <- read_csv(CFG$flagged_cells_file, show_col_types = FALSE)
flag_row <- match(format(as.Date(flagged$periodid), "%Y-%m-%d"),
                  format(dates, "%Y-%m-%d"))
stopifnot(!anyNA(flag_row), all(flagged$district %in% districts))
for (i in seq_len(nrow(flagged))) EXCLUDE[flag_row[i], flagged$district[i]] <- TRUE
cat(sprintf("Cells excluded from the metrics: %d\n", sum(EXCLUDE)))

# --- Period covariates: encoded on the DAY OF YEAR ------------------------
# NOT on isoweek: the week ending 2021-01-03 is ISO week 53 of 2020, which
# would overflow a /52 encoding and create a discontinuity.
# yday/365.25 is bounded, cannot overflow, and absorbs calendar drift
# (the first observation of a year moves between 3 and 9 January).
# `dates` are week-ENDING Sundays.
ydn <- lubridate::yday(dates) / 365.25
#
# THREE VARIANTS, NOT TWO. This distinction matters.
#   NoTime        = no temporal predictor at all
#   Seasonal      = SEASONALITY only (Fourier terms + rainy-season indicator)
#   SeasonalTrend = seasonality + TREND (normalised time index)
#
# Why separate them: "trend" is a monotone index of time. Once standardised,
# it mechanically brings neighbouring weeks closer together in the KNN
# distance -- in other words it turns KNN into a disguised temporal method.
# Comparing only NoTime and SeasonalTrend would leave us unable to say whether
# the gain comes from SEASONALITY (interesting, epidemiological) or from a
# bare TIME INDEX (trivial, and already captured by interpolation).
# A reviewer would ask exactly that question.
PERIOD_SEASONAL <- data.frame(
  sin1 = sin(2*pi*ydn), cos1 = cos(2*pi*ydn),   # annual cycle
  sin2 = sin(4*pi*ydn), cos2 = cos(4*pi*ydn),   # semi-annual harmonic
  # Rainy season in Adamaoua: approx. April-October. TO BE CONFIRMED with
  # the supervisor or local climatology before publication. The Cameroon
  # ONACC climate profile gives May-October; the widest definition is used
  # here as the conservative choice.
  rainy_season = as.integer(month(dates) %in% 4:10)
)
PERIOD_SEASONAL_TREND <- cbind(PERIOD_SEASONAL, trend = seq_len(n_obs) / n_obs)

# Helper used everywhere: variant name -> predictor block
period_block <- function(variant)
  switch(variant,
         NoTime        = NULL,
         Seasonal      = PERIOD_SEASONAL,
         SeasonalTrend = PERIOD_SEASONAL_TREND,
         stop("unknown variant: ", variant))


## --------------------------------------------------------- 2. MASK SIMULATION
# The mask applies ONLY to the district columns. The period covariates are
# never punched with holes.
#
# gap_length = 1  -> isolated MCAR cells: the reference condition, standing
#                    for one-off omissions. Cells are drawn independently, so
#                    neighbouring cells can form short runs, more often at
#                    higher rates: see output/metrics/gap_length_distribution.csv
#                    (R/11, section 11).
# gap_length > 1  -> block MCAR: blocks of gap_length consecutive weeks in one
#                    district, standing for prolonged non-reporting. Blocks
#                    can touch or overlap, and start at i <= n_obs - gap_length,
#                    so the last week is never masked. How the ranking of the
#                    methods varies with gap length is analysed in 03 and R/11.
simulate_mask <- function(rate, seed, gap_length = 1) {
  set.seed(seed)
  M <- matrix(FALSE, n_obs, n_dis, dimnames = list(NULL, districts))
  n_na <- round(n_obs * n_dis * rate / 100)
  if (gap_length == 1) {
    M[sample.int(n_obs * n_dis, n_na)] <- TRUE
  } else {
    while (sum(M) < n_na) {
      j <- sample.int(n_dis, 1); i <- sample.int(n_obs - gap_length, 1)
      M[i:(i + gap_length - 1), j] <- TRUE
    }
  }
  stopifnot(all(colSums(!M) > 2), all(rowSums(!M) > 0))
  M
}

# --- MASK VALIDATION: four tests, run on every mask ------------------------
# Checks empirically that the simulation really is MCAR. Almost no published
# paper does this; the resulting table belongs in the supplementary material.
validate_mcar <- function(M) {
  yr <- lubridate::year(dates)
  # 1. uniformity across districts
  p1 <- chisq.test(colSums(M))$p.value
  # 2. uniformity over time (expected counts proportional to exposure)
  obs <- tapply(rowSums(M), yr, sum)
  exp <- as.numeric(table(yr)) * n_dis; exp <- exp / sum(exp) * sum(obs)
  p2 <- chisq.test(as.numeric(obs), p = exp / sum(exp))$p.value
  # 3. THE DECISIVE TEST: do the masked true values follow the same
  #    distribution as the observed ones? Under MCAR, yes.
  Z  <- scale(REF)
  p3 <- suppressWarnings(ks.test(Z[M], Z[!M])$p.value)
  # 4. independence of missingness between districts (~5% significant pairs)
  pr <- combn(n_dis, 2, function(ij) {
    tb <- table(M[, ij[1]], M[, ij[2]])
    if (!all(dim(tb) == c(2, 2))) return(NA_real_)
    fisher.test(tb)$p.value
  })
  tibble(p_districts = p1, p_time = p2, p_values_ks = p3,
         pct_pairs_significant = 100 * mean(pr < 0.05, na.rm = TRUE),
         valid = p3 > 0.05 & p2 > 0.05)
}


## ------------------------------------------------------------- 3. THE METHODS
# Every function takes X (the matrix with holes) plus P (a period data frame
# or NULL) and returns a complete matrix of the same dimensions as X.

imp_mean <- function(X, P = NULL) {
  for (j in seq_len(ncol(X))) X[is.na(X[, j]), j] <- mean(X[, j], na.rm = TRUE)
  X
}

# Stochastic MICE (pmm): m completions drawn by pmm, then AVERAGED; the average
# shrinks the variance (see the note on M_INFERENCE below).
imp_mice_stochastic <- function(X, P = NULL, m = 5) {
  D   <- as.data.frame(X); if (!is.null(P)) D <- cbind(D, P)
  fit <- mice(D, m = m, maxit = 20, method = "pmm", printFlag = FALSE, seed = 1)
  Reduce(`+`, lapply(1:m, function(i)
    as.matrix(complete(fit, i)[, seq_len(ncol(X))]))) / m
}

# Deterministic MICE: the conditional mean. Comparable on pointwise error.
imp_mice_deterministic <- function(X, P = NULL) {
  D   <- as.data.frame(X); if (!is.null(P)) D <- cbind(D, P)
  fit <- mice(D, m = 1, maxit = 20, method = "norm.predict",
              printFlag = FALSE, seed = 1)
  as.matrix(complete(fit, 1)[, seq_len(ncol(X))])
}

# KNN: k is CRITICAL (49% difference measured between k = 1 and k = 5).
# k is chosen by repeated internal re-masking -- see tune_k() below.
imp_knn <- function(X, P = NULL, k = 5) {
  D  <- as.data.frame(X); if (!is.null(P)) D <- cbind(D, P)
  mu <- colMeans(D, na.rm = TRUE)
  sd <- apply(D, 2, sd, na.rm = TRUE); sd[sd == 0] <- 1
  Z  <- sweep(sweep(D, 2, mu, "-"), 2, sd, "/")
  Zi <- VIM::kNN(as.data.frame(Z), k = k, imp_var = FALSE)
  as.matrix(sweep(sweep(Zi[, seq_len(ncol(X))], 2, sd[seq_len(ncol(X))], "*"),
                  2, mu[seq_len(ncol(X))], "+"))
}

# missRanger: equivalent to missForest, 10 to 50 times faster.
imp_missforest <- function(X, P = NULL) {
  D <- as.data.frame(X); if (!is.null(P)) D <- cbind(D, P)
  as.matrix(missRanger(D, num.trees = 100, pmm.k = 0, verbose = 0,
                       seed = 1)[, seq_len(ncol(X))])
}

# --- TEMPORAL BASELINES: the main omission of the original study plan ------
imp_linear_interp <- function(X, P = NULL)
  apply(X, 2, function(v)
    as.numeric(imputeTS::na_interpolation(v, option = "linear")))

# na_kalman(model = "StructTS") fails on some series (non-convergent
# optimisation, singular matrix). Fallback cascade: auto.arima, then linear
# interpolation. Fallbacks are counted and reported at the end: if there are
# many, it must be stated in the Methods section.
N_KALMAN_FALLBACKS <- 0L
imp_kalman_vector <- function(v) {
  r <- try(as.numeric(imputeTS::na_kalman(v, model = "StructTS")), silent = TRUE)
  if (!inherits(r, "try-error") && !anyNA(r)) return(r)
  N_KALMAN_FALLBACKS <<- N_KALMAN_FALLBACKS + 1L
  r <- try(as.numeric(imputeTS::na_kalman(v, model = "auto.arima")), silent = TRUE)
  if (!inherits(r, "try-error") && !anyNA(r)) return(r)
  as.numeric(imputeTS::na_interpolation(v))
}
imp_kalman <- function(X, P = NULL) apply(X, 2, imp_kalman_vector)

imp_seasonal_decomp <- function(X, P = NULL)
  apply(X, 2, function(v)
    as.numeric(imputeTS::na_seadec(ts(v, frequency = 52),
                                   algorithm = "interpolation")))

# --- FIFTH FAMILY: spatio-temporal, without any geostatistics --------------
# Idea: give missRanger the temporal lags of every district IN ADDITION to
# the other districts. The model learns on its own which neighbours (spatial
# AND temporal) matter, and where. No shapefile required.
imp_missforest_st <- function(X, P = NULL) {
  lagm <- function(v, k)
    if (k > 0) c(rep(NA, k), head(v, -k)) else c(tail(v, k), rep(NA, -k))
  L <- do.call(cbind, lapply(c(1, 2, -1, -2), function(k) apply(X, 2, lagm, k = k)))
  colnames(L) <- paste0("lag", rep(c(1, 2, -1, -2), each = ncol(X)), "_",
                        rep(districts, 4))
  D <- as.data.frame(cbind(X, L)); if (!is.null(P)) D <- cbind(D, P)
  as.matrix(missRanger(D, num.trees = 100, pmm.k = 0, verbose = 0,
                       seed = 1)[, seq_len(ncol(X))])
}

M_INFERENCE <- 5    # m = 5,  completions AVERAGED (label kept for compatibility;
                    # no Rubin's rules are applied in this benchmark)
M_PRECISION <- 20   # m = 20, completions averaged -> point precision
# Measured: m = 1 gives MAE 78, m = 50 gives 44.6. Averaging the completions
# converges to the conditional expectation, hence towards norm.predict.
# This is not merely a tuning choice, it is a choice of estimand.

METHODS <- list(
  Mean                = imp_mean,
  MICE_inference      = function(X, P) imp_mice_stochastic(X, P, m = M_INFERENCE),
  MICE_precision      = function(X, P) imp_mice_stochastic(X, P, m = M_PRECISION),
  MICE_deterministic  = imp_mice_deterministic,
  KNN                 = imp_knn,        # k injected dynamically, see the loop
  MissForest          = imp_missforest,
  LinearInterp        = imp_linear_interp,
  Kalman              = imp_kalman,
  SeasonalDecomp      = imp_seasonal_decomp,
  MissForest_ST       = imp_missforest_st
)
# Purely univariate methods: the period block cannot be passed to them, so
# they are evaluated under the NoTime variant only.
UNIVARIATE <- c("LinearInterp", "Kalman", "SeasonalDecomp")


## ------------------------------ 4. HONEST HYPERPARAMETER TUNING ------------
#
#  PRINCIPLE: never tune a parameter by looking at the ground truth.
#  A district-level practitioner has no ground truth -- if they cannot
#  recover your k, your recommendation is worthless. Tuning on the truth also
#  unfairly favours the most tunable methods (KNN) over those with nothing to
#  tune (interpolation): the comparison stops being fair.
#
#  PROCEDURE: re-mask a sample of the ALREADY OBSERVED cells (whose values we
#  therefore know without cheating), measure the error there, and repeat B
#  times to stabilise.
#
#  MEASURED ON THESE DATA: B = 1 leaves an optimism of 1.58 MAE (too noisy),
#  B = 10 brings it down to 0.73 MAE (~1.5%, negligible). Use B >= 10.
#  Worth reporting: k = 5 fixed in advance performs just as well.
tune_k <- function(X, rate, grid = c(3, 5, 7, 10, 15, 20), B = 10) {
  err <- matrix(NA_real_, B, length(grid), dimnames = list(NULL, grid))
  obs <- which(!is.na(X))
  for (b in seq_len(B)) {
    set.seed(4242 + b)
    sel   <- sample(obs, round(length(obs) * rate / 100))
    X2    <- X; truth <- X[sel]; X2[sel] <- NA
    for (i in seq_along(grid))
      err[b, i] <- mean(abs(imp_knn(X2, NULL, k = grid[i])[sel] - truth))
  }
  grid[which.min(colMeans(err))]
}


## ---------------------------------------------------------------- 5. METRICS
# NB: M & !EXCLUDE -> unreliable reference cells are removed from the
# evaluation even when they were masked and imputed normally.
#
# EPIDEMIC THRESHOLD: the 90th percentile of each district, computed on the
# complete reference. Used by the detection metric below.
THRESHOLD_P90 <- apply(REF, 2, quantile, probs = 0.90)

compute_metrics <- function(Z, M) {
  M <- M & !EXCLUDE
  truth <- REF[M]; est <- Z[M]; d <- est - truth

  # --- DETECTION OF HIGH-INCIDENCE WEEKS -------------------------------
  # This is THE surveillance question: does an imputed cell cross the
  # epidemic threshold when the true value crosses it?
  # It is indispensable because MAE, RMSE and MAPE all reward smoothing:
  # a method that pulls everything towards the mean scores a good average
  # error while erasing the peaks entirely. This metric penalises exactly
  # that, and it is the one that matters for an alert system.
  cell_threshold <- matrix(THRESHOLD_P90, n_obs, n_dis, byrow = TRUE)[M]
  truth_high <- truth > cell_threshold
  est_high   <- est   > cell_threshold
  sens <- if (any(truth_high))  mean(est_high[truth_high])   else NA_real_
  spec <- if (any(!truth_high)) mean(!est_high[!truth_high]) else NA_real_

  tibble(
    MAE         = mean(abs(d)),
    RMSE        = sqrt(mean(d^2)),
    MAPE        = mean(abs(d / truth)) * 100,  # safe: no zero in REF
    bias        = mean(d),                     # reveals underestimation
    var_ratio   = var(est) / var(truth),       # reveals smoothing
    R2          = cor(est, truth)^2,
    sensitivity = sens,                        # peaks recovered
    specificity = spec,                        # false alarms avoided
    youden      = sens + spec - 1,             # single criterion, 1 = perfect
    n_peaks     = sum(truth_high)
  )
}

compute_metrics_by_district <- function(Z, M) {
  M <- M & !EXCLUDE
  map_dfr(seq_len(n_dis), function(j) {
    if (!any(M[, j])) return(NULL)
    v <- REF[M[, j], j]; e <- Z[M[, j], j]; d <- e - v
    tibble(district = districts[j], n_imputed = length(v),
           MAE = mean(abs(d)), RMSE = sqrt(mean(d^2)),
           MAPE = mean(abs(d / v)) * 100, bias = mean(d))
  })
}

cat("Shared functions loaded.\n")
