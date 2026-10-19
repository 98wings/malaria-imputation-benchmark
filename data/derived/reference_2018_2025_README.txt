Corrected reference series
==========================
Built by data/extend_reference_2018_2025.py. Do not edit by hand.

reference_adamaoua_weekly_2018_2022_corrected.csv : 261 weeks, 2018-01-07 to 2023-01-01
reference_adamaoua_weekly_2018_2025.csv           : 417 weeks, 2018-01-07 to 2025-12-28
11 districts. Every row is one ISO week (Monday to Sunday); periodid is the
Sunday on which the week ends.

Defect corrected
  In reference_adamaoua_weekly_2018_2022.csv (used by the article) the row dated
  2020-12-27 is the sum of ISO weeks 52 and 53 of 2020, and every date from 2021
  onwards is one week early. Verified against raw/data (3).csv: weeks 1 to 51 of
  2020 are identical to the reference, and weeks 52 + 53 equal the row 2020-12-27
  in all 11 districts.
  Correction: the row is split into week 52 (2020-12-27) and week 53 (2021-01-03),
  and all later rows are moved forward by seven days. The total number of cases
  over 2018-2022 is unchanged (907588).

2021-2022, comparison with raw/malaria_2021_2025_fusionne.csv (103 weeks, 1133 cells)
  identical cells : 1128
  different cells : 5, listed in reference_2018_2025_overlap_differences.csv;
                    the value of the article's reference was kept.
  weeks absent from that file : 2023-01-01 (taken from the article's reference).

2023-2025 : appended from raw/malaria_2021_2025_fusionne.csv (156 weeks).
Dropped   : district Mayo-Baleo (no value in 2021-2022, zeros in 2023-2025).

Status: the article and all its results use reference_adamaoua_weekly_2018_2025.csv (417 weeks).
