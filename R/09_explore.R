###############################################################################
#  09 — INSPECT THE DATA WITHOUT RISK  (helper functions, not part of the run)
#
#  Load it whenever you want to LOOK at the data:
#      source("R/09_explore.R")
#
#  THE RULE, stated correctly:
#    Opening a CSV in Excel damages nothing.
#    SAVING from Excel can break everything (';' separator, decimal comma,
#    reformatted dates).
#  The problem is that the Ctrl+S reflex and the "Do you want to save your
#  changes?" dialog both exist. These functions mean you never have to take
#  that risk.
###############################################################################

library(tidyverse); library(lubridate)

REF_PATH <- "data/derived/reference_adamaoua_weekly_2018_2025.csv"
dir.create("preview", showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a


## ==========================================================================
## OPTION 1 — Look at it inside RStudio (fastest, zero risk)
## ==========================================================================
#  Opens a sortable, filterable table in RStudio. The file on disk is never
#  touched: you are handling an in-memory copy.
inspect <- function(path = REF_PATH) {
  d <- read_csv(path, show_col_types = FALSE)
  View(d, title = basename(path))
  invisible(d)
}


## ==========================================================================
## OPTION 2 — A console summary, usually enough
## ==========================================================================
overview <- function(path = REF_PATH, n = 10) {
  d <- read_csv(path, show_col_types = FALSE)
  cat("File:", path, "\n")
  cat(sprintf("Dimensions: %d rows x %d columns\n", nrow(d), ncol(d)))
  cat(sprintf("Missing values: %d (%.2f %%)\n\n",
              sum(is.na(d)), 100 * mean(is.na(d))))
  print(head(d, n))
  cat("\n--- Column summary ---\n")
  num <- d |> select(where(is.numeric))
  print(tibble(
    column = names(num),
    n_na   = colSums(is.na(num)),
    min    = apply(num, 2, min, na.rm = TRUE),
    median = apply(num, 2, median, na.rm = TRUE),
    mean   = round(colMeans(num, na.rm = TRUE), 1),
    max    = apply(num, 2, max, na.rm = TRUE)
  ), n = 50)
  invisible(d)
}


## ==========================================================================
## OPTION 3 — Excel, but on a disposable .xlsx copy
## ==========================================================================
#  THE solution if you want Excel: export a copy to .xlsx in preview/. Open
#  it, sort it, colour it, delete it -- the working CSV is unaffected.
#  The .xlsx is also immune to the two CSV traps: it stores types (a date
#  stays a date) and depends on no separator.
export_for_excel <- function(path = REF_PATH, name = NULL) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")
  d <- read_csv(path, show_col_types = FALSE)
  if ("periodid" %in% names(d)) d$periodid <- as.Date(d$periodid)
  out <- file.path("preview",
                   name %||% paste0("PREVIEW_",
                                    tools::file_path_sans_ext(basename(path)),
                                    ".xlsx"))

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "data")
  openxlsx::writeDataTable(wb, "data", d)
  openxlsx::freezePane(wb, "data", firstActiveRow = 2, firstActiveCol = 2)
  openxlsx::setColWidths(wb, "data", 1:ncol(d), widths = "auto")
  openxlsx::conditionalFormatting(wb, "data", cols = 1:ncol(d),
    rows = 2:(nrow(d) + 1), type = "contains", rule = "NA",
    style = openxlsx::createStyle(bgFill = "#FFC7CE"))
  openxlsx::saveWorkbook(wb, out, overwrite = TRUE)

  cat("Disposable copy written:", out, "\n")
  cat("Open it, edit it, delete it: the working CSV is intact.\n")
  invisible(out)
}


## ==========================================================================
## OPTION 4 — See the STRUCTURE of missingness (more telling than a table)
## ==========================================================================
#  A 417-row table shows very little. This map does.
missingness_map <- function(path) {
  d  <- read_csv(path, show_col_types = FALSE)
  dt <- if ("periodid" %in% names(d)) as.Date(d$periodid) else seq_len(nrow(d))
  d |>
    select(where(is.numeric)) |>
    mutate(date = dt) |>
    pivot_longer(-date, names_to = "district") |>
    mutate(state = ifelse(is.na(value), "missing", "observed")) |>
    ggplot(aes(date, district, fill = state)) +
    geom_tile() +
    scale_fill_manual(values = c(observed = "grey88", missing = "#c0392b")) +
    labs(title = "Structure of the missing data", x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 11) + theme(legend.position = "bottom")
}


## ==========================================================================
## OPTION 5 — Plot the series
## ==========================================================================
plot_series <- function(path = REF_PATH, keep = NULL) {
  d  <- read_csv(path, show_col_types = FALSE)
  dt <- as.Date(d$periodid)
  lg <- d |> select(where(is.numeric)) |> mutate(date = dt) |>
    pivot_longer(-date, names_to = "district")
  if (!is.null(keep)) lg <- filter(lg, district %in% keep)
  ggplot(lg, aes(date, value)) +
    geom_line(linewidth = .35, colour = "#1b4965") +
    facet_wrap(~ district, scales = "free_y", ncol = 3) +
    labs(title = "Weekly malaria incidence by district", x = NULL, y = "Cases") +
    theme_minimal(base_size = 10)
}


## ==========================================================================
## OPTION 6 — Check that a CSV has not been damaged
## ==========================================================================
#  Run it if you have any doubt after handling a file. It detects the three
#  classic kinds of damage caused by saving from Excel.
check_integrity <- function(path = REF_PATH) {
  lines <- readLines(path, n = 3, warn = FALSE)
  ok <- TRUE

  if (grepl(";", lines[1])) {
    cat("X SEPARATOR: the file uses ';' instead of ','\n")
    cat("  -> Excel re-saved it in a French locale. Restore from data/raw/.\n")
    ok <- FALSE
  }
  d <- suppressWarnings(try(read_csv(path, show_col_types = FALSE), silent = TRUE))
  if (inherits(d, "try-error")) {
    cat("X The file can no longer be read.\n"); return(invisible(FALSE))
  }

  if ("periodid" %in% names(d)) {
    dt <- suppressWarnings(as.Date(d$periodid, format = "%Y-%m-%d"))
    if (any(is.na(dt))) {
      cat("X DATES: ISO format lost. Example found:", d$periodid[1], "\n")
      cat("  -> Excel reformatted to DD/MM/YYYY. Restore from data/raw/.\n")
      ok <- FALSE
    } else if (!all(diff(dt) == 7)) {
      cat("X DATES: the weekly step is broken (rows sorted or deleted?).\n")
      ok <- FALSE
    }
  }
  num <- d |> select(where(is.numeric))
  if (ncol(num) < ncol(d) - 1) {
    cat("X TYPES: numeric columns have become text.\n")
    cat("  -> decimal comma introduced by Excel. Restore from data/raw/.\n")
    ok <- FALSE
  }
  if (ok) cat("OK - separator, dates and types intact.\n")
  invisible(ok)
}


cat("Available functions:\n",
    "  overview()          console summary\n",
    "  inspect()           RStudio table viewer\n",
    "  export_for_excel()  disposable .xlsx copy in preview/\n",
    "  missingness_map(f)  map of the missing values\n",
    "  plot_series()       time series plots\n",
    "  check_integrity()   post-handling check\n")
