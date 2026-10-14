###############################################################################
#  08 — CHOROPLETH MAPS OF THE HEALTH DISTRICTS  (optional)
#
#  Prerequisites: 01, 02, 03, 04 must have run.
#  Run from the project root (open the .Rproj in RStudio).
#
#  THIS SCRIPT HAS TWO PARTS.
#  PART A downloads the health-district boundaries, inspects them, and builds
#         a name-correspondence table with the 11 districts.
#         ---> run PART A first and check its output.
#  PART B runs only if the correspondence is complete (11/11).
#
#  BOUNDARY SOURCE: eHealth Africa, via openAFRICA
#    https://open.africa/dataset/cameroon-health-district
#  Citation a reviewer will expect:
#    Rosencrans LC, et al. Mapping for Health in Cameroon: Polio Legacy and
#    Beyond. J Infect Dis. 2017;216(suppl_1):S105-S112.
#
#  DO NOT use geodata::gadm("CMR", level = 2): that returns the 5 ADMINISTRATIVE
#  departments of Adamaoua, not the 11 HEALTH districts. They are different
#  divisions and the join will silently fail.
###############################################################################

library(sf); library(dplyr); library(readr); library(tidyr); library(ggplot2)

GEO_DIR <- "data/geo"
FIG     <- "output/figures"
MET     <- "output/metrics"
dir.create(GEO_DIR, recursive = TRUE, showWarnings = FALSE)

GEO_FILE <- file.path(GEO_DIR, "cameroon_health_districts.geojson")
ADA_FILE <- file.path(GEO_DIR, "adamaoua_districts.gpkg")

## =========================================================================
## PART A.1 — DOWNLOAD
## =========================================================================
#  The eHealth Africa WFS server presents a self-signed SSL certificate: a
#  standard download fails with "certificate verify failed". We therefore
#  disable certificate verification FOR THIS URL ONLY. Acceptable here
#  (public, non-sensitive data) but not to be generalised.
WFS_URL <- paste0(
  "https://gis-a.ie.ehealthafrica.org/geoserver/eHA_db/ows",
  "?service=WFS&version=1.0.0&request=GetFeature",
  "&typeName=eHA_db:sv_health_districts_cameroon",
  "&outputFormat=application/json",
  "&authkey=fdfe9a37-d2d0-4210-9a15-25dab5d907fa")

#  NOTE: testing file.exists() alone is not enough. A failed download leaves a
#  ZERO-BYTE file behind, the script then skips the download and st_read()
#  fails with an opaque error. (This actually happened in the French version:
#  data/geo held a 0-byte geojson.) We therefore also require a plausible size.
if (!file.exists(GEO_FILE) || file.size(GEO_FILE) < 1e5) {
  if (file.exists(GEO_FILE))
    cat("Existing file is too small (failed download): downloading again.\n")
  cat("Downloading boundaries (this may take 1-2 minutes)...\n")
  ok <- FALSE
  # Attempt 1: command-line curl (-k ignores the certificate)
  try({
    download.file(WFS_URL, GEO_FILE, method = "curl",
                  extra = "-k -L --max-time 300", quiet = FALSE)
    ok <- file.exists(GEO_FILE) && file.size(GEO_FILE) > 1e5
  }, silent = TRUE)
  # Attempt 2: httr, in case curl is unavailable on the machine
  if (!ok) {
    cat("curl unavailable or failed -> trying httr...\n")
    try({
      httr::GET(WFS_URL, httr::config(ssl_verifypeer = 0L), httr::timeout(300),
                httr::write_disk(GEO_FILE, overwrite = TRUE))
      ok <- file.exists(GEO_FILE) && file.size(GEO_FILE) > 1e5
    }, silent = TRUE)
  }
  if (!ok)
    stop("Download failed.\nManual fallback: open this URL in a browser\n  ",
         WFS_URL, "\naccept the security warning, save the file as\n  ",
         GEO_FILE, "\nthen re-run this script.")
}
cat(sprintf("File present: %s (%.1f MB)\n", GEO_FILE, file.size(GEO_FILE) / 1e6))

## =========================================================================
## PART A.2 — INSPECTION
## =========================================================================
geo <- sf::st_read(GEO_FILE, quiet = TRUE)

cat("\n=== LAYER STRUCTURE ===\n")
cat("Number of features (189 health districts expected):", nrow(geo), "\n")
cat("Coordinate reference system:", sf::st_crs(geo)$input, "\n")
cat("Available columns:\n"); print(names(geo))

cat("\n=== ATTRIBUTE PREVIEW (first 5 rows, geometry dropped) ===\n")
print(as.data.frame(sf::st_drop_geometry(geo))[1:5, ], width = Inf)

#  We do not guess which column holds the district name: we test.
#  Heuristic: the region column has ~10 distinct values (Cameroon's regions),
#  the district column has ~189.
text_cols <- names(geo)[sapply(sf::st_drop_geometry(geo), is.character)]
cat("\n=== TEXT COLUMNS AND THEIR CARDINALITY ===\n")
for (cc in text_cols) {
  v <- sf::st_drop_geometry(geo)[[cc]]
  cat(sprintf("  %-28s %4d distinct values | e.g.: %s\n",
              cc, dplyr::n_distinct(v),
              paste(utils::head(unique(v), 3), collapse = " / ")))
}
card <- sapply(text_cols, \(cc) dplyr::n_distinct(sf::st_drop_geometry(geo)[[cc]]))
col_region   <- text_cols[which.min(abs(card - 10))]
col_district <- text_cols[which.min(abs(card - 189))]
cat(sprintf("\nDetected region column   : %s\n", col_region))
cat(sprintf("Detected district column : %s\n", col_district))
cat("\n=== VALUES OF THE REGION COLUMN ===\n")
print(sort(unique(sf::st_drop_geometry(geo)[[col_region]])))

## =========================================================================
## PART A.3 — NAME CORRESPONDENCE
## =========================================================================
#  A critical and underestimated step: a join on place names almost never
#  works first time (accents, hyphens, case, competing spellings). We
#  normalise aggressively BEFORE joining, then VERIFY that all 11 match. A
#  silently unmatched district would appear grey on the map, and a reviewer
#  would notice.
normalise <- function(x) {
  x |>
    iconv(from = "UTF-8", to = "ASCII//TRANSLIT") |>   # e-acute -> e
    toupper() |>
    gsub("[^A-Z]", "", x = _)                          # drop spaces, _, -
}

reference <- readr::read_csv(
  "data/derived/reference_adamaoua_weekly_2018_2025.csv", show_col_types = FALSE)
my_districts <- setdiff(names(reference), "periodid")
stopifnot(length(my_districts) == 11)

geo <- geo |> mutate(.name = .data[[col_district]],
                     .region = .data[[col_region]],
                     .key = normalise(.name))

correspondence <- tibble::tibble(district = my_districts,
                                 .key = normalise(my_districts)) |>
  left_join(sf::st_drop_geometry(geo) |>
              select(.key, geo_name = .name, geo_region = .region) |>
              distinct(.key, .keep_all = TRUE),
            by = ".key")

cat("\n=== CORRESPONDENCE TABLE (check by eye) ===\n")
print(as.data.frame(correspondence |> select(district, geo_name, geo_region)),
      row.names = FALSE)

unmatched <- correspondence |> filter(is.na(geo_name)) |> pull(district)
if (length(unmatched) > 0) {
  cat("\n!!! ", length(unmatched), " unmatched district(s): ",
      paste(unmatched, collapse = ", "), "\n", sep = "")
  cat("Candidates available in the most likely region:\n")
  print(sort(sf::st_drop_geometry(geo) |>
               filter(grepl("ADAMA", normalise(.region))) |> pull(.name)))
  cat("\nSTOPPING HERE. Complete the correspondence table by hand with the\n",
      "real names rather than guessing them.\n", sep = "")
  stop("Incomplete correspondence -- PART B not executed.")
}
cat("\nCorrespondence complete: 11/11. We can map.\n")

ada <- geo |>
  inner_join(correspondence |> select(.key, district), by = ".key") |>
  select(district, geometry) |> sf::st_make_valid()
stopifnot(nrow(ada) == 11)

#  Save the filtered layer: the script stays reproducible even if the
#  eHealth Africa server goes down (which happens).
sf::st_write(ada, ADA_FILE, delete_dsn = TRUE, quiet = TRUE)
cat("Adamaoua layer saved:", ADA_FILE, "\n")

## =========================================================================
## PART B — THE MAPS
## =========================================================================

## Map 1 — mean weekly incidence (goes in the Methods section, to describe
##         the study area). A descriptive map is expected by any reviewer of
##         a paper on spatially indexed data.
incidence <- tibble::tibble(
  district = my_districts,
  mean_cases = sapply(my_districts, \(d) mean(reference[[d]])))
m1 <- ada |> left_join(incidence, by = "district")

p1 <- ggplot(m1) +
  geom_sf(aes(fill = mean_cases), colour = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "magma", direction = -1,
                       name = "Mean weekly\ncases") +
  geom_sf_text(aes(label = district), size = 2.6, colour = "grey20") +
  labs(title = "Mean weekly malaria incidence, Adamaoua, 2018-2025",
       caption = "Boundaries: eHealth Africa (Rosencrans et al., 2017)") +
  theme_void(base_size = 11)
ggsave(file.path(FIG, "fig7_incidence_map.png"), p1,
       width = 8, height = 6.5, dpi = 300)

## Map 2 — gain of MissForest over mean imputation, by district (goes in the
##         Results: the spatial translation of the district-level finding).
f_gain <- file.path(MET, "gain_by_district.csv")
if (file.exists(f_gain)) {
  gain <- readr::read_csv(f_gain, show_col_types = FALSE)
  m2 <- ada |> left_join(gain |> select(district, gain_MissForest_pct,
                                        mean_corr_other_districts),
                         by = "district")
  stopifnot(!anyNA(m2$gain_MissForest_pct))

  p2 <- ggplot(m2) +
    geom_sf(aes(fill = gain_MissForest_pct), colour = "white", linewidth = 0.3) +
    scale_fill_viridis_c(option = "viridis", name = "Gain vs\nmean (%)") +
    geom_sf_text(aes(label = sprintf("%s\n%.0f%%", district,
                                     gain_MissForest_pct)),
                 size = 2.4, colour = "white", fontface = "bold") +
    labs(title = "Benefit of MissForest over mean imputation",
         subtitle = "The gain is largest where a district correlates with its peers",
         caption = "Boundaries: eHealth Africa (Rosencrans et al., 2017)") +
    theme_void(base_size = 11)
  ggsave(file.path(FIG, "fig8_gain_map.png"), p2, width = 8, height = 6.5, dpi = 300)

  ## Spatial coherence check: is the "peer correlation" computed in script 04
  ## consistent with REAL contiguity? Until now "neighbour" was a statistical
  ## notion, never a geographical one. This is the only genuine scientific
  ## added value of the map, beyond illustration.
  nb <- sf::st_touches(ada)
  n_nb <- lengths(nb)
  cat("\n=== GEOGRAPHICAL CONTIGUITY ===\n")
  for (i in seq_len(nrow(ada)))
    cat(sprintf("  %-20s %d neighbour(s): %s\n", ada$district[i], n_nb[i],
                paste(ada$district[nb[[i]]], collapse = ", ")))
  if (any(n_nb == 0))
    cat("\nWarning: district(s) with no neighbour -> check st_make_valid.\n")

  ## Moran's I on mean incidence: is there spatial structure?
  ## Without it, a reviewer will ask "did you test spatial autocorrelation?".
  if (requireNamespace("spdep", quietly = TRUE) && all(n_nb > 0)) {
    lw <- spdep::nb2listw(spdep::poly2nb(as(ada, "Spatial")), style = "W",
                          zero.policy = TRUE)
    mt <- spdep::moran.test(m1$mean_cases, lw, zero.policy = TRUE)
    cat(sprintf("\nMoran's I on mean incidence: I = %+.3f, p = %.4f\n",
                mt$estimate[1], mt$p.value))
    cat("  (n = 11 districts: EXPLORATORY, present it as such)\n")
  } else {
    cat("\nFor Moran's test: install.packages('spdep')\n")
  }
} else {
  cat("\n", f_gain, " missing: map 2 skipped. Run R/04_final_analysis.R first.\n")
}

cat("\n=== Done ===\n")
cat("Figures written to ", FIG, "\n", sep = "")
