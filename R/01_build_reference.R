###############################################################################
#  01 — CERTIFY THE REFERENCE DATASET
#
#  Run ONCE, before anything else.
#  Input   : data/derived/reference_adamaoua_weekly_2018_2025.csv (never modified)
#  Outputs : data/derived/reference_flagged_cells.csv
#            data/derived/reference_belel_sensitivity.csv
#            data/derived/reference_characteristics.csv
#
#  WHERE THE INPUT COMES FROM
#  The reference series is assembled by data/extend_reference_2018_2025.py
#  (Python) from the extractions kept in data/raw/. That script
#    - starts from the 2018-2022 series built from Palu_Adamaoua_2018-2022.csv,
#    - splits the row that held ISO weeks 52 AND 53 of 2020 added together,
#    - appends 2023-2025 from the supervisor's extraction.
#  Each row is one ISO week (Monday to Sunday) and periodid is the Sunday on
#  which it ends. See data/derived/reference_2018_2025_README.txt.
#
#  THIS SCRIPT NO LONGER BUILDS THE SERIES. It certifies it (dates, counts,
#  upstream imputations), records the cell excluded from the metrics, and
#  computes the characteristics reported in the paper.
#  The version that built the 2018-2022 series is kept in
#  _Archives/R/01_build_reference_2018_2022.R.
###############################################################################

library(tidyverse); library(lubridate)

SRC <- "data/derived/reference_adamaoua_weekly_2018_2025.csv"
OUT <- "data/derived"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

raw <- read_csv(SRC, show_col_types = FALSE)

## --------------------------------------------------------------------------
## 1. CHOICE OF SOURCE: the file carrying BOTH the date and the total
## --------------------------------------------------------------------------
#  We start from Palu_Adamaoua_2018-2022.csv because it is the SUPERSET: it
#  contains everything the other two contain. We then remove what we do not
#  want, rather than trying to add back what is missing.
#  Starting from the SP file would be a mistake: it has no date column, and a
#  lost time index cannot be recovered. That index carries the seasonality,
#  the lag structure, and the entire downstream forecasting analysis.

districts <- setdiff(names(raw), c("periodid", "Palu_Adamaoua"))
stopifnot(length(districts) == 11)

## --------------------------------------------------------------------------
## 2. THE TOTAL COLUMN: dropped from the file, not from the analysis
## --------------------------------------------------------------------------
#  Palu_Adamaoua is EXACTLY the sum of the 11 districts (verified: maximum
#  discrepancy = 0). It is therefore a derived column, not data.
#  We do not store it, for two reasons:
#    (a) once stored, it eventually finds its way into the predictor set ->
#        information leakage, since it contains the very value being imputed;
#    (b) after imputation the original total no longer matches the imputed
#        columns -> any plot then compares two different things without
#        anyone noticing.
#  Rule: the total is ALWAYS recomputed, on both sides.
#      total_ref <- rowSums(REF)      # for the reference
#      total_imp <- rowSums(Z)        # for each imputed dataset
#  (The 2018-2025 reference carries no total column: nothing to check or drop.
#  The test is kept for any future source file that would carry one.)
if ("Palu_Adamaoua" %in% names(raw)) {
  total_gap <- max(abs(raw$Palu_Adamaoua - rowSums(raw[, districts])))
  cat(sprintf("Check: total equals sum of districts, max discrepancy = %.6f\n",
              total_gap))
  stopifnot(total_gap < 1e-6)
}

## --------------------------------------------------------------------------
## 3. THE DATE: already normalised, but certified here
## --------------------------------------------------------------------------
#  periodid is already stored as ISO YYYY-MM-DD in the source file. The
#  "date problem" observed during exploratory work came from Excel, which
#  re-renders dates as DD/MM/YYYY on opening. Never re-save a CSV from Excel;
#  open it from R only. (Opening is harmless; saving is what corrupts.)
dates <- as.Date(raw$periodid, format = "%Y-%m-%d")

stopifnot(
  !any(is.na(dates)),                      # everything parsed
  !any(duplicated(dates)),                 # no duplicate week
  all(diff(dates) == 7),                   # strict 7-day step
  !is.unsorted(dates)                      # chronological order
)
cat(sprintf("Dates: %d weeks, from %s to %s, 7-day step, weekday = %s\n",
            length(dates), min(dates), max(dates),
            unique(weekdays(dates))))

#  Every periodid is the SUNDAY ending an ISO week, and the series holds every
#  ISO week of the period: 52 per year, 53 in 2020. This is the check that
#  would have caught the merged weeks 52 + 53 of 2020 in the first extraction.
iso_years  <- isoyear(min(dates)):isoyear(max(dates))
n_expected <- sum(isoweek(as.Date(paste0(iso_years, "-12-28"))))
stopifnot(all(wday(dates, week_start = 1) == 7),      # Sundays
          isoweek(min(dates)) == 1,                   # starts on ISO week 1
          length(dates) == n_expected)                # no week missing or merged
cat(sprintf("ISO weeks: %d expected for %d-%d, %d present\n",
            n_expected, min(iso_years), max(iso_years), length(dates)))

## --------------------------------------------------------------------------
## 4. THE ISO WEEK 53 TRAP
## --------------------------------------------------------------------------
#  isoweek() returns 53 for 2021-01-03 (the week 28 Dec 2020 - 3 Jan 2021). An encoding of the form
#  sin(2*pi*week/52) would then overflow the cycle and introduce a
#  discontinuity exactly at a year boundary.
#  Seasonality is therefore encoded on the DAY OF YEAR instead: it cannot
#  overflow, and it absorbs calendar drift (the first observation of a year
#  moves between 3 and 9 January).
n_iso53 <- sum(isoweek(dates) == 53)
cat(sprintf("ISO week 53 occurrences: %d  -> encode on yday, not on isoweek\n",
            n_iso53))

## --------------------------------------------------------------------------
## 5. THE ONLY UNOBSERVED VALUE IN THE SERIES: Belel, 2019-09-01
## --------------------------------------------------------------------------
#  CONTEXT (confirmed with the data providers): this is the only genuinely
#  missing value in the whole 2018-2022 dataset. It had been filled upstream
#  with Belel's 2019 annual mean, hence the stored value 156.098039.
#
#  WHY IT CANNOT BE KEPT AS IS
#  Its neighbouring weeks are 79 and 70. A value of 156 is an artificial
#  spike. The cell plays three roles in this study and corrupts two of them:
#    - ground truth for the error metrics  -> INVALID (never observed)
#    - predictor for other cells           -> distorted: a spike that
#                                             never happened
#    - reference series for forecasting    -> same distortion
#  Worse: if the cell falls inside a simulated mask, mean imputation
#  "recovers" it almost perfectly. We would then be measuring the internal
#  consistency of an artefact, not the performance of a method.
#
#  DECISION, IN THREE PARTS
#   (a) REPLACE with 74 = linear interpolation between adjacent weeks,
#       (79 + 70) / 2 = 74.5, rounded (these are counts). A simple,
#       conventional rule fixed a priori -- and therefore not circular with
#       respect to the benchmark that compares these same methods.
#   (b) EXCLUDE from all metrics: a reconstructed value cannot serve as
#       ground truth, however good it is (see EXCLUDE in 02_run_benchmark.R).
#   (c) DECLARE in the paper, together with the sensitivity table below.
#
#  SENSITIVITY (estimates obtained by setting the cell back to NA):
#       Linear interpolation      74.5   <- retained
#       Seasonal decomposition    69.7
#       MissForest + period       85.4
#       KNN (k = 7)               87.1
#       MissForest                87.9
#       MICE deterministic       121.2
#       Annual mean              156.1   <- original value, implausible
#  Every serious method lands between 70 and 88. Because the cell is excluded
#  from the metrics, its only residual effect is as a predictor: 1 cell out of
#  the whole matrix (0.035% for 2018-2022, 0.022% for 2018-2025). The exact choice has no measurable consequence -- but
#  documenting it is what separates careful work from the rest.
CELL_DATE <- "2019-09-01"; CELL_DISTRICT <- "Belel"; CELL_VALUE <- 74

#  In the 2018-2025 reference the cell already holds the retained value (74).
#  The value found in the original extraction is recorded here as a constant.
SOURCE_VALUE <- 156.098039215686     # Belel 2019 annual mean, Palu_Adamaoua_2018-2022.csv
idx        <- match(CELL_DATE, format(as.Date(raw$periodid), "%Y-%m-%d"))
src_value  <- SOURCE_VALUE
stopifnot(!is.na(idx), raw[[CELL_DISTRICT]][idx] == CELL_VALUE)
neighbours <- raw[[CELL_DISTRICT]][(idx - 2):(idx + 2)]
cat(sprintf("\nCell %s / %s: %.6f -> %d\n", CELL_DATE, CELL_DISTRICT,
            src_value, CELL_VALUE))
cat("  neighbouring weeks (W-2..W+2):", paste(neighbours, collapse = ", "), "\n")

flagged <- tibble(
  periodid          = CELL_DATE,
  district          = CELL_DISTRICT,
  row_index         = idx,
  source_value      = src_value,
  source_value_origin = "Belel 2019 annual mean (imputation applied upstream of this study)",
  retained_value    = CELL_VALUE,
  replacement_rule  = "linear interpolation between adjacent weeks: (79+70)/2 = 74.5, rounded",
  excluded_from_metrics = TRUE,
  exclusion_reason  = "never observed: cannot serve as ground truth",
  retained_role     = "predictor for other cells and reference series for forecasting"
)
write_csv(flagged, file.path(OUT, "reference_flagged_cells.csv"))

sensitivity <- tibble(
  method = c("Linear interpolation (RETAINED)", "Seasonal decomposition",
             "MissForest + period", "KNN (k=7)", "MissForest",
             "MICE deterministic", "Annual mean (source value)"),
  estimate = c(74.5, 69.7, 85.4, 87.1, 87.9, 121.2, 156.1)) |>
  mutate(difference_vs_retained = round(estimate - 74.5, 1))
write_csv(sensitivity, file.path(OUT, "reference_belel_sensitivity.csv"))
print(sensitivity)

#  Generic detection: are there OTHER upstream-imputed cells, invisible
#  because they are whole plausible numbers? Test: a value exactly equal to
#  its column's global or annual mean is almost certainly an imputation.
#  Result on the 2018-2022 extraction: exactly 1 cell (Belel, above), a
#  NON-INTEGER value far from its neighbours.
#  Re-run here on 2018-2025, where the Belel cell already holds 74.
detect_upstream_imputations <- function(df, dts) {
  map_dfr(districts, function(d) {
    v <- df[[d]]
    global_mean <- rep(mean(v), length(v))
    annual_mean <- ave(v, year(dts), FUN = mean)
    k <- which(abs(v - global_mean) < 1e-6 | abs(v - annual_mean) < 1e-6)
    if (!length(k)) return(NULL)
    tibble(district = d, periodid = format(dts[k], "%Y-%m-%d"), value = v[k])
  })
}
detected <- detect_upstream_imputations(raw, dates)
cat(sprintf("\nCells equal to a mean (global or annual): %d\n", nrow(detected)))
print(detected)

#  REVIEWED CANDIDATES. A cell listed here was examined by hand and judged to
#  be a genuine observation that happens to equal a mean.
#    Meiganga, 2020-05-24 (ISO week 21 of 2020): value 400, and the 52 weeks
#    dated 2020 average exactly 400. Kept as observed because (i) the original
#    2020 extraction, data/raw/data (3).csv, already gives the integer 400,
#    (ii) it lies close to its neighbours (425 before, 485 after), and (iii) an
#    upstream mean fill leaves a non-integer value, as it did for Belel.
#  TO EXCLUDE IT INSTEAD: delete its line below and add it to `flagged` above.
REVIEWED <- tribble(
  ~district,  ~periodid,
  "Meiganga", "2020-05-24"
)
unreviewed <- anti_join(detected, REVIEWED, by = c("district", "periodid"))
if (nrow(unreviewed)) {
  cat("\nUNREVIEWED cells equal to a mean. Inspect each one, then either add it\n",
      "to REVIEWED (genuine value) or to `flagged` (excluded from the metrics):\n")
  print(unreviewed)
}
stopifnot(nrow(unreviewed) == 0)   # halts if new candidates appear

## --------------------------------------------------------------------------
## 6. FINAL CHECKS ON THE REFERENCE (it is read, never rewritten, here)
## --------------------------------------------------------------------------
reference <- raw |>
  select(periodid, all_of(districts)) |>
  mutate(periodid = format(dates, "%Y-%m-%d"))

stopifnot(!anyNA(reference), nrow(reference) == n_expected, ncol(reference) == 12,
          all(map_lgl(reference[districts], \(x) all(x %% 1 == 0))),
          all(map_lgl(reference[districts], \(x) all(x > 0))))

## --------------------------------------------------------------------------
## 7. CHARACTERISATION: to be reported in the Data section of the paper
## --------------------------------------------------------------------------
#  These statistics are not decorative. They determine which method will win,
#  and they give the reader the criteria needed to judge whether our
#  conclusions transfer to their own data.
M  <- as.matrix(reference[, districts])
yd <- yday(dates)

#  NOTE (subtle, and it cost us a run): the statistics are computed BEFORE the
#  tibble, not inside it. tibble() builds its columns SEQUENTIALLY and makes
#  each finished column visible to the arguments that follow it. A column
#  named `mean` or `sd` therefore SHADOWS the functions mean() and sd() for
#  every later argument: apply(M, 2, sd) is then handed a numeric vector
#  instead of a function, and fails with
#        object 'FUN' of mode 'function' not found
#  The French version escaped this only because its columns were called
#  `moyenne` and `ecart_type`. We keep the English names and move the
#  computation out. base:: / stats:: prefixes make the intent explicit.
col_mean <- colMeans(M)
col_sd   <- apply(M, 2, stats::sd)
col_min  <- apply(M, 2, base::min)
col_max  <- apply(M, 2, base::max)
col_ac1  <- apply(M, 2, \(v) cor(v[-1], head(v, -1)))
col_ac52 <- apply(M, 2, \(v) cor(v[-(1:52)], head(v, -52)))
col_corr <- sapply(districts,
                   \(d) base::mean(cor(M)[d, setdiff(districts, d)]))
col_r2   <- apply(M, 2, \(v) {
  p <- ave(v, cut(yd, breaks = seq(0, 366, by = 7)), FUN = base::mean)
  1 - sum((v - p)^2) / sum((v - base::mean(v))^2)
})

characteristics <- tibble(
  district        = districts,
  mean            = col_mean,
  sd              = col_sd,
  cv              = col_sd / col_mean,
  min             = col_min,
  max             = col_max,
  n_zeros         = colSums(M == 0),
  autocorr_lag1   = col_ac1,
  autocorr_lag52  = col_ac52,
  mean_corr_other_districts = col_corr,
  seasonal_r2     = col_r2
)
write_csv(characteristics, file.path(OUT, "reference_characteristics.csv"))
print(characteristics |> mutate(across(where(is.numeric), \(x) round(x, 3))))

cat(sprintf("\nZeros in the whole matrix: %d  -> MAPE %s\n", sum(M == 0),
            ifelse(sum(M == 0) == 0, "usable without SMAPE",
                                     "MUST BE REPLACED by SMAPE")))
cat(sprintf("Mean between-district correlation: %.3f\n",
            mean(cor(M)[upper.tri(cor(M))])))
cat("\n=== Reference certified. Continue with 02_run_benchmark.R ===\n")
