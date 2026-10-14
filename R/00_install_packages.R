###############################################################################
#  00 — INSTALL PACKAGES
#
#  RUN ONCE, before anything else.
#  Installs the required packages. Allow 10 to 20 minutes the first time
#  (some packages compile from source). An internet connection is required.
#
#  In RStudio: open this file, then click "Source" (or Ctrl+Shift+S).
###############################################################################

pkgs <- c(
  # data manipulation and graphics
  "tidyverse", "lubridate",
  # imputation
  "mice", "VIM", "missRanger", "imputeTS",
  # time series and forecasting
  "forecast",
  # statistical comparison tests
  "PMCMRplus",
  # utilities
  "tictoc", "openxlsx", "ggrepel"
)

# Optional: only needed for the choropleth map in 08_district_map.R.
# Installing sf can be slow on Windows; skip it if you are not making the map.
optional_pkgs <- c("sf", "spdep", "httr")

missing <- pkgs[!pkgs %in% rownames(installed.packages())]

if (length(missing) == 0) {
  cat("All required packages are already installed.\n")
} else {
  cat("Installing:", paste(missing, collapse = ", "), "\n\n")
  install.packages(missing, dependencies = TRUE)
}

# Verification: does each package actually load?
cat("\n--- Verification ---\n")
ok <- vapply(pkgs, function(p) {
  r <- requireNamespace(p, quietly = TRUE)
  cat(sprintf("  %-14s %s\n", p, ifelse(r, "OK", "FAILED")))
  r
}, logical(1))

cat("\n--- Optional (map only) ---\n")
for (p in optional_pkgs)
  cat(sprintf("  %-14s %s\n", p,
              ifelse(requireNamespace(p, quietly = TRUE), "OK", "not installed")))

if (all(ok)) {
  cat("\nEverything is ready. You can now run R/01_build_reference.R\n")
} else {
  cat("\nFailed packages:", paste(pkgs[!ok], collapse = ", "), "\n")
  cat("Re-run install.packages() for those, and read the error message.\n")
  cat("On Windows, a compilation error is usually fixed by installing\n")
  cat("Rtools: https://cran.r-project.org/bin/windows/Rtools/\n")
}

cat("\nR version:", R.version.string, "\n")
if (getRversion() < "4.1.0")
  cat("WARNING: these scripts use the native pipe |>, which requires R >= 4.1.\n",
      "Update R, or replace |> with %>% throughout.\n")
