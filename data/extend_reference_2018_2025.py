#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Builds the corrected weekly reference series, 2018-2022 and 2018-2025.

    python extend_reference_2018_2025.py        (run from anywhere)

Reads (never modified)
    data/derived/reference_adamaoua_weekly_2018_2022.csv   reference used by the article
    data/raw/data (3).csv                                  2020, ISO weeks 1 to 53
    data/raw/malaria_2021_2025_fusionne.csv                2021-2025, ISO weeks
Writes
    data/derived/reference_adamaoua_weekly_2018_2022_corrected.csv   261 weeks
    data/derived/reference_adamaoua_weekly_2018_2025.csv             417 weeks
    data/derived/reference_2018_2025_overlap_differences.csv
    data/derived/reference_2018_2025_README.txt

What is corrected
  In the reference used by the article, the row dated 2020-12-27 is the SUM of
  ISO weeks 52 and 53 of 2020, and every date from 2021 onwards is one week
  early as a consequence. Here:
    1. the row 2020-12-27 is split into week 52 (2020-12-27) and week 53
       (2021-01-03), taken from the 2020 extraction;
    2. every later row of the reference is moved forward by seven days;
    3. 2023-2025 is appended from the new file.
  After correction every row is one ISO week (Monday to Sunday) and `periodid`
  is the Sunday on which it ends.
  All other values of the article's reference are kept as they are.
The script stops with an error if any check fails; it never writes a partial result.
"""
import datetime as dt
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OLD = os.path.join(HERE, 'derived', 'reference_adamaoua_weekly_2018_2022.csv')
RAW2020 = os.path.join(HERE, 'raw', 'data (3).csv')
RAWNEW = os.path.join(HERE, 'raw', 'malaria_2021_2025_fusionne.csv')
OUT22 = os.path.join(HERE, 'derived', 'reference_adamaoua_weekly_2018_2022_corrected.csv')
OUT25 = os.path.join(HERE, 'derived', 'reference_adamaoua_weekly_2018_2025.csv')
DIFF = os.path.join(HERE, 'derived', 'reference_2018_2025_overlap_differences.csv')
NOTE = os.path.join(HERE, 'derived', 'reference_2018_2025_README.txt')

MERGED = pd.Timestamp('2020-12-27')      # row of the reference that holds weeks 52 + 53
WEEK = pd.Timedelta(days=7)
DROP = ['Mayo_Baleo']
HEADER = re.compile(r'^Week (\d+) (\d{4}-\d{2}-\d{2}) - (\d{4}-\d{2}-\d{2})$')


def stop(msg):
    sys.exit('STOP: ' + msg)


def clean(name):
    return re.sub(r'[ \-]+', '_', re.sub(r'^.*District\s+', '', name).strip())


# ---------------------------------------------------------------- article reference
old = pd.read_csv(OLD, parse_dates=['periodid']).set_index('periodid')
districts = list(old.columns)

# ---------------------------------------------------------------- 2020, weeks 1 to 53
r20 = pd.read_csv(RAW2020)
r20 = r20.rename(columns={c: clean(c) for c in r20.columns if 'District' in c})
weeks = r20['periodid'].str.extract(r'^2020W(\d+)$')[0]
if weeks.isna().any() or list(weeks.astype(int)) != list(range(1, 54)):
    stop('the 2020 file does not hold ISO weeks 1 to 53 in order')
if sorted(c for c in r20.columns if c in districts) != sorted(districts):
    stop('districts of the 2020 file differ from the reference')
r20.index = pd.DatetimeIndex([pd.Timestamp(dt.date.fromisocalendar(2020, w, 7)) for w in range(1, 54)])
r20 = r20[districts].astype(int)

before = old[old.index < MERGED]
if not r20.iloc[:51].equals(before.loc[r20.index[:51]]):
    stop('weeks 1 to 51 of 2020 do not match the reference')
if list(old.loc[MERGED]) != list(r20.iloc[51] + r20.iloc[52]):
    stop('row 2020-12-27 of the reference is not the sum of weeks 52 and 53')

# ---------------------------------------------------------------- corrected 2018-2022
after = old[old.index > MERGED].copy()
after.index = after.index + WEEK
ref22 = pd.concat([before, r20.iloc[51:53], after])
ref22.index.name = 'periodid'

# ---------------------------------------------------------------- 2021-2025 file
wide = pd.read_csv(RAWNEW)
wide['District'] = wide['District'].map(clean)
wide = wide[~wide['District'].isin(DROP)]
if sorted(wide['District']) != sorted(districts):
    stop('districts of the 2021-2025 file differ from the reference')
ends = {}
for col in wide.columns[1:]:
    m = HEADER.match(col)
    if not m:
        stop('unexpected column name: %r' % col)
    start, end = pd.Timestamp(m.group(2)), pd.Timestamp(m.group(3))
    if start.dayofweek != 0 or (end - start).days != 6 or start.isocalendar()[1] != int(m.group(1)):
        stop('column is not a Monday-to-Sunday ISO week: %r' % col)
    ends[col] = end
new = wide.set_index('District').T
new.index = pd.DatetimeIndex([ends[c] for c in new.index], name='periodid')
new = new.sort_index()[districts]
if new.index.duplicated().any() or new.isna().any().any() or ((new % 1) != 0).any().any():
    stop('duplicated weeks, missing values or non-integer counts in the 2021-2025 file')
new = new.astype(int)

# ---------------------------------------------------------------- overlap 2021-2022
common = new.index.intersection(ref22.index)
delta = (new.loc[common] - ref22.loc[common]).stack()
delta = delta[delta != 0]
n_cells = len(common) * len(districts)
if len(delta) > 0.02 * n_cells:
    stop('%d of %d overlapping cells differ: the files are not aligned' % (len(delta), n_cells))
diff = pd.DataFrame({
    'periodid': [d.strftime('%Y-%m-%d') for d, _ in delta.index],
    'district': [k for _, k in delta.index],
    'value_new_file': [int(new.loc[d, k]) for d, k in delta.index],
    'value_reference_kept': [int(ref22.loc[d, k]) for d, k in delta.index],
})
diff['difference'] = diff['value_new_file'] - diff['value_reference_kept']
absent = [d.strftime('%Y-%m-%d') for d in ref22.index[ref22.index >= new.index.min()] if d not in new.index]

# ---------------------------------------------------------------- 2018-2025
ref25 = pd.concat([ref22, new[new.index > ref22.index.max()]])
ref25.index.name = 'periodid'


def check(df, label):
    if set(pd.Series(df.index).diff().dropna().dt.days) != {7}:
        stop(label + ': not one row every seven days')
    if set(df.index.dayofweek) != {6}:
        stop(label + ': a date is not a Sunday')
    first = df.index[0] - pd.Timedelta(days=6)
    expected = sum(dt.date(y, 12, 28).isocalendar()[1] for y in range(first.isocalendar()[0], df.index[-1].isocalendar()[0] + 1))
    if first.isocalendar()[1] != 1 or len(df) != expected:
        stop('%s: %d rows for %d ISO weeks' % (label, len(df), expected))
    if df.isna().any().any():
        stop(label + ': missing values')


check(ref22, '2018-2022')
check(ref25, '2018-2025')
if int(ref22.values.sum()) != int(old.values.sum()):
    stop('the total number of cases 2018-2022 changed')        # splitting a row must not create or lose cases

# ---------------------------------------------------------------- write
for df, path in ((ref22, OUT22), (ref25, OUT25)):
    out = df.copy()
    out.index = out.index.strftime('%Y-%m-%d')
    out.to_csv(path)
diff.to_csv(DIFF, index=False)
fmt = lambda d: d.strftime('%Y-%m-%d')
with open(NOTE, 'w', encoding='utf-8') as f:
    f.write(
        "Corrected reference series\n"
        "==========================\n"
        "Built by data/extend_reference_2018_2025.py. Do not edit by hand.\n\n"
        "reference_adamaoua_weekly_2018_2022_corrected.csv : %d weeks, %s to %s\n"
        "reference_adamaoua_weekly_2018_2025.csv           : %d weeks, %s to %s\n"
        "%d districts. Every row is one ISO week (Monday to Sunday); periodid is the\n"
        "Sunday on which the week ends.\n\n"
        "Defect corrected\n"
        "  In reference_adamaoua_weekly_2018_2022.csv (used by the article) the row dated\n"
        "  2020-12-27 is the sum of ISO weeks 52 and 53 of 2020, and every date from 2021\n"
        "  onwards is one week early. Verified against raw/data (3).csv: weeks 1 to 51 of\n"
        "  2020 are identical to the reference, and weeks 52 + 53 equal the row 2020-12-27\n"
        "  in all %d districts.\n"
        "  Correction: the row is split into week 52 (2020-12-27) and week 53 (2021-01-03),\n"
        "  and all later rows are moved forward by seven days. The total number of cases\n"
        "  over 2018-2022 is unchanged (%d).\n\n"
        "2021-2022, comparison with raw/malaria_2021_2025_fusionne.csv (%d weeks, %d cells)\n"
        "  identical cells : %d\n"
        "  different cells : %d, listed in reference_2018_2025_overlap_differences.csv;\n"
        "                    the value of the article's reference was kept.\n"
        "  weeks absent from that file : %s (taken from the article's reference).\n\n"
        "2023-2025 : appended from raw/malaria_2021_2025_fusionne.csv (%d weeks).\n"
        "Dropped   : district Mayo-Baleo (no value in 2021-2022, zeros in 2023-2025).\n\n"
        "Status: the article and its results still use reference_adamaoua_weekly_2018_2022.csv.\n"
        "The benchmark has NOT been rerun on the corrected series.\n"
        % (len(ref22), fmt(ref22.index[0]), fmt(ref22.index[-1]),
           len(ref25), fmt(ref25.index[0]), fmt(ref25.index[-1]),
           len(districts), len(districts), int(ref22.values.sum()),
           len(common), n_cells, n_cells - len(delta), len(delta), ', '.join(absent) or 'none',
           len(ref25) - len(ref22)))

print('written :')
print('  %s : %d weeks, %s to %s' % (os.path.basename(OUT22), len(ref22), fmt(ref22.index[0]), fmt(ref22.index[-1])))
print('  %s : %d weeks, %s to %s' % (os.path.basename(OUT25), len(ref25), fmt(ref25.index[0]), fmt(ref25.index[-1])))
print('  row 2020-12-27 split into weeks 52 and 53; dates from 2021 moved forward by 7 days')
print('  overlap 2021-2022 : %d of %d cells identical, %d different (article values kept)' % (n_cells - len(delta), n_cells, len(delta)))
print('  weeks absent from the 2021-2025 file :', ', '.join(absent) or 'none')
