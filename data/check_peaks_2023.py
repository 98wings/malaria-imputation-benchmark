"""Checks on the autumn 2023 peaks and on a possible level break in early 2023.

Reads  data/derived/reference_adamaoua_weekly_2018_2025.csv (never modified).
Writes data/derived/peaks_2023_checks.csv        (one row per district x year)
       data/derived/level_break_2023_check.csv   (one row per district)

Peak checks, applied to every week from September to December of every year,
for the three districts whose 2023 autumn values are the series maxima
(Ngaoundere Rural, Bankim, Djohong). A week is flagged when its value:
  equal_mean_median : equals the rounded mean or median of its district,
                      computed over its calendar year or over the whole series;
  repeat_previous   : equals the previous week (last value carried forward);
  interpolated      : equals the rounded mean of the previous and next weeks;
  double_count      : equals the sum of the previous and next weeks, or twice
                      one of them.
The same checks on the other years give the rate expected by chance.

Level break: for each district, mean weekly count in January-June of 2022,
2023 and 2024, and the ratios 2023/2022 and 2024/2023. A district carved out
of another in early 2023 would show a drop in the mother district.
"""
import pandas as pd, numpy as np, os

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "derived", "reference_adamaoua_weekly_2018_2025.csv")
r = pd.read_csv(REF, parse_dates=["periodid"])
districts = [c for c in r.columns if c != "periodid"]
r["year"] = r.periodid.dt.year
r["month"] = r.periodid.dt.month

rows = []
for d in ["Ngaoundere_Rural", "Bankim", "Djohong"]:
    v = r[d].to_numpy()
    prev = np.r_[np.nan, v[:-1]]; nxt = np.r_[v[1:], np.nan]
    glob = {round(np.mean(v)), round(np.median(v))}
    for y, g in r.groupby("year"):
        ann = {round(g[d].mean()), round(g[d].median())}
        idx = g.index[g.month >= 9]
        f_mean = [v[i] in glob | ann for i in idx]
        f_rep = [v[i] == prev[i] for i in idx]
        f_int = [(not np.isnan(prev[i]) and not np.isnan(nxt[i])
                  and v[i] == round((prev[i] + nxt[i]) / 2)) for i in idx]
        f_dbl = [(not np.isnan(prev[i]) and not np.isnan(nxt[i])
                  and (v[i] == prev[i] + nxt[i] or v[i] == 2 * prev[i] or v[i] == 2 * nxt[i]))
                 for i in idx]
        anyf = [a or b or c or e for a, b, c, e in zip(f_mean, f_rep, f_int, f_dbl)]
        rows.append(dict(district=d, year=y, n_weeks_sep_dec=len(idx),
                         max_sep_dec=int(v[idx].max()),
                         equal_mean_median=sum(f_mean), repeat_previous=sum(f_rep),
                         interpolated=sum(f_int), double_count=sum(f_dbl),
                         any_flag=sum(anyf)))
pk = pd.DataFrame(rows)
pk.to_csv(os.path.join(HERE, "derived", "peaks_2023_checks.csv"), index=False)

h1 = r[r.month <= 6]
lb = []
for d in districts:
    m = {y: h1.loc[h1.year == y, d].mean() for y in (2022, 2023, 2024)}
    lb.append(dict(district=d, mean_h1_2022=m[2022], mean_h1_2023=m[2023], mean_h1_2024=m[2024],
                   ratio_2023_2022=m[2023] / m[2022], ratio_2024_2023=m[2024] / m[2023]))
lb = pd.DataFrame(lb)
lb.to_csv(os.path.join(HERE, "derived", "level_break_2023_check.csv"), index=False)
print(pk.to_string(index=False)); print(); print(lb.round(3).to_string(index=False))
