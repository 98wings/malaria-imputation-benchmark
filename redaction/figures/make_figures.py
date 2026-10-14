#!/usr/bin/env python3
"""
Manuscript figures, generated ONLY from the outputs of the R pipeline.

    python3 make_figures.py          # EN + FR, main figures + supplement

Sources (no value is typed by hand):
  data/derived/reference_adamaoua_weekly_2018_2025.csv   reference series (417 weeks)
  output/masks/mask_rep01_rate10_gap{1,4,8}.rds          actual masks (read with pyreadr)
  output/imputed/, output/imputed_examples/              imputed series of the example figure
  output/metrics/*.csv, output_mechanism/metrics/*.csv   metrics
The figure -> file -> column mapping is in redaction/data_provenance.csv.

Outputs:
  redaction/tex/figN_*.pdf and figN_*_fr.pdf  vector, embedded fonts, called
                                              without a path by the .tex
                                              (the sn-jnl template forbids subfolders)
  redaction/figures/png/*.png                 300 dpi, for proofreading or upload
  redaction/supplements/S5_fig*.pdf|.png      figures of Additional file 5

Conventions (decided by the panel, see redaction/plan/rapport_figures_panel.md):
  - width 170 mm (BMC full page); sans-serif fonts of 7 to 8.5 pt at final
    size; panel letters in bold lower case; no title inside the image (BMC:
    title and legend belong in the manuscript).
  - colour = method FAMILY, identical on every figure; Okabe-Ito palette,
    colour-blind safe; the marker shape doubles the colour for greyscale
    reading.
  - Fig. 4: colour encodes the MECHANISM (orange -> vermilion -> black, from
    the weakest to the strongest selection), with one shape per mechanism.
Dependencies: pandas, numpy, matplotlib, pyreadr (pip install pyreadr).
"""
import math
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap

HERE = os.path.dirname(os.path.abspath(__file__))
RED = os.path.dirname(HERE)
ROOT = os.path.dirname(RED)
MET = os.path.join(ROOT, 'output', 'metrics')
MECH = os.path.join(ROOT, 'output_mechanism', 'metrics')
MASKS = os.path.join(ROOT, 'output', 'masks')
REF = os.path.join(ROOT, 'data', 'derived', 'reference_adamaoua_weekly_2018_2025.csv')
OUT_TEX = os.path.join(RED, 'tex')
OUT_PNG = os.path.join(RED, 'figures', 'png')
OUT_SUP = os.path.join(RED, 'supplements')
os.makedirs(OUT_PNG, exist_ok=True)

MM = 1 / 25.4
FULL = 170 * MM

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
    'font.size': 8, 'axes.labelsize': 8, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'legend.fontsize': 7, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
    'ytick.major.width': 0.6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
    'lines.linewidth': 1.0, 'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 300,
    'legend.frameon': False, 'axes.unicode_minus': True,
})

# ------------------------------------------------------------- conventions
FAMILY = {'Kalman': 'Temporal', 'LinearInterp': 'Temporal', 'SeasonalDecomp': 'Temporal',
          'MissForest_ST': 'Hybrid', 'MissForest': 'Cross-sectional', 'KNN': 'Cross-sectional',
          'MICE_inference': 'Cross-sectional', 'MICE_precision': 'Cross-sectional',
          'MICE_deterministic': 'Cross-sectional', 'Mean': 'Naive'}
FAM_COL = {'Temporal': '#009E73', 'Hybrid': '#CC79A7', 'Cross-sectional': '#0072B2',
           'Naive': '#7F7F7F'}
FAM_MRK = {'Temporal': 'o', 'Hybrid': 'D', 'Cross-sectional': 's', 'Naive': 'X'}
FAM_ORDER = ['Temporal', 'Hybrid', 'Cross-sectional', 'Naive']

MECHS = ['MCAR', 'MAR', 'MAR_strong', 'MNAR']
MECH_STYLE = {'MCAR': dict(color='#999999', marker='o', mfc='white'),
              'MAR': dict(color='#E69F00', marker='^', mfc='#E69F00'),
              'MAR_strong': dict(color='#D55E00', marker='D', mfc='#D55E00'),
              'MNAR': dict(color='#000000', marker='s', mfc='#000000')}

DISTRICT = {'Ngaoundere_Rural': 'Ngaoundéré Rural', 'Ngaoundere_Urbain': 'Ngaoundéré Urbain',
            'Tignere': 'Tignère'}

T = {
 'en': dict(
   method={'SeasonalDecomp': 'Seasonal decomposition', 'LinearInterp': 'Linear interpolation',
           'Kalman': 'Kalman smoothing', 'MissForest': 'MissForest', 'MissForest_ST': 'MissForest-ST',
           'KNN': 'KNN', 'Mean': 'Mean imputation', 'MICE_inference': 'MICE-pmm (m=5)',
           'MICE_precision': 'MICE-pmm (m=20)', 'MICE_deterministic': 'MICE-norm.predict'},
   family={'Temporal': 'Temporal', 'Hybrid': 'Hybrid', 'Cross-sectional': 'Cross-sectional',
           'Naive': 'Naive'},
   mech={'MCAR': 'MCAR', 'MAR': 'MAR', 'MAR_strong': 'Strong MAR', 'MNAR': 'MNAR'},
   ratio='Weekly cases relative to\nthe district median',
   gaps='{g}-week gaps: {n} masked cells', gap1='1-week gaps: {n} masked cells',
   week='Epidemiological week', rank_axis='Rank among 24 combinations (1 = best)',
   gap_axis='Gap length (weeks)', mae_axis='MAE (cases per week)',
   mean_rank='Mean rank (1 = best)', cd='CD = {cd:.2f}',
   drank='Change in mean rank, 8-week minus 1-week gaps',
   worse='worse with long gaps', better='better with long gaps',
   bias='Mean bias at 1-week gaps (cases per week)',
   mae_gap1='MAE at 1-week gaps (cases per week)', youden="Youden's index",
   h4='Forecast, 4-week horizon', h52='Forecast, 52-week horizon',
   deg='RMSE degradation\n(cases per week)',
   rho='ρ = {r:+.2f}',
   seas='Seasonality', trend='Trend', gain='MAE gain over the NoTime block (%)',
   gapt='{g}-week gaps', corr='Mean correlation with the other districts',
   gainmf='Gain of MissForest\nover mean imputation (%)',
   ref='Reference', masked='masked weeks', cases='Weekly cases, {d}',
   ex1='1-week gaps', ex8='8-week gap', todo='run R/10_example_series.R',
 ),
 'fr': dict(
   method={'SeasonalDecomp': 'Décomposition saisonnière', 'LinearInterp': 'Interpolation linéaire',
           'Kalman': 'Lissage de Kalman', 'MissForest': 'MissForest', 'MissForest_ST': 'MissForest-ST',
           'KNN': 'KNN', 'Mean': 'Imputation par la moyenne', 'MICE_inference': 'MICE-pmm (m=5)',
           'MICE_precision': 'MICE-pmm (m=20)', 'MICE_deterministic': 'MICE-norm.predict'},
   family={'Temporal': 'Temporelle', 'Hybrid': 'Hybride', 'Cross-sectional': 'Transversale',
           'Naive': 'Naïve'},
   mech={'MCAR': 'MCAR', 'MAR': 'MAR', 'MAR_strong': 'MAR renforcé', 'MNAR': 'MNAR'},
   ratio='Cas hebdomadaires rapportés\nà la médiane du district',
   gaps='Lacunes de {g} semaines : {n} cellules masquées', gap1='Lacunes de 1 semaine : {n} cellules masquées',
   week='Semaine épidémiologique', rank_axis='Rang parmi 24 combinaisons (1 = meilleur)',
   gap_axis='Longueur des lacunes (semaines)', mae_axis='MAE (cas par semaine)',
   mean_rank='Rang moyen (1 = meilleur)', cd='DC = {cd:.2f}',
   drank='Variation du rang moyen\n(lacunes de 8 semaines moins lacunes de 1 semaine)',
   worse='moins bon sur lacunes longues', better='meilleur sur lacunes longues',
   bias='Biais moyen, lacunes de 1 semaine\n(cas par semaine)',
   mae_gap1='MAE, lacunes de 1 semaine (cas par semaine)', youden='Indice de Youden',
   h4='Prévision, horizon de 4 semaines', h52='Prévision, horizon de 52 semaines',
   deg='Dégradation de la RMSE\n(cas par semaine)',
   rho='ρ = {r:+.2f}',
   seas='Saisonnalité', trend='Tendance', gain='Gain de MAE par rapport au bloc NoTime (%)',
   gapt='Lacunes de {g} semaines', corr='Corrélation moyenne avec les autres districts',
   gainmf="Gain de MissForest sur\nl'imputation par la moyenne (%)",
   ref='Référence', masked='semaines masquées', cases='Cas hebdomadaires, {d}',
   ex1='Lacunes de 1 semaine', ex8='Lacune de 8 semaines', todo='lancer R/10_example_series.R',
 ),
}


def fmt(x, lang):
    """Number in the language's format (decimal comma in French)."""
    return x.replace('.', ',') if lang == 'fr' else x


def panel(ax, letter, x=-0.02, y=1.02):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10, fontweight='bold',
            ha='right', va='bottom')


def save(fig, name, lang, where='tex'):
    suffix = '' if lang == 'en' else '_fr'
    if where == 'tex':
        pdf = os.path.join(OUT_TEX, name + suffix + '.pdf')
        png = os.path.join(OUT_PNG, name + suffix + '.png')
    else:
        pdf = os.path.join(OUT_SUP, name + suffix + '.pdf')
        png = os.path.join(OUT_SUP, name + suffix + '.png')
    fig.savefig(pdf, bbox_inches='tight', pad_inches=0.02)
    fig.savefig(png, bbox_inches='tight', pad_inches=0.02, dpi=300)
    plt.close(fig)
    print('  ', os.path.relpath(pdf, ROOT))


def load(*p):
    return pd.read_csv(os.path.join(*p))


def family_handles(lang, fams=FAM_ORDER):
    return [Line2D([], [], ls='', marker=FAM_MRK[f], color=FAM_COL[f], ms=4.5,
                   label=T[lang]['family'][f]) for f in fams]


def district_characteristics():
    """Series characteristics of the 11 districts, as written by the R pipeline
    (R/04_final_analysis.R, gain_by_district.csv) from the current reference."""
    return load(MET, 'gain_by_district.csv')


# ================================================================ FIGURE 1
def fig_design(lang):
    import pyreadr
    t = T[lang]
    ref = pd.read_csv(REF, parse_dates=['periodid']).set_index('periodid')
    districts = list(ref.columns)
    rel = np.log2(ref / ref.median())
    masks = {g: pyreadr.read_r(os.path.join(MASKS, 'mask_rep01_rate10_gap%d.rds' % g))[None]
             for g in (1, 4, 8)}
    for g, m in masks.items():
        assert m.shape == ref.shape and list(m.columns) == districts, g
    dates = mdates.date2num(ref.index.to_pydatetime())
    step = dates[1] - dates[0]
    ext = [dates[0] - step / 2, dates[-1] + step / 2, len(districts) - 0.5, -0.5]
    names = [DISTRICT.get(d, d) for d in districts]

    fig = plt.figure(figsize=(FULL, 118 * MM))
    gs = fig.add_gridspec(4, 2, height_ratios=[1.6, 1, 1, 1], width_ratios=[1, 0.018],
                          hspace=0.38, wspace=0.02)
    ax = fig.add_subplot(gs[0, 0])
    lim = 1.5
    im = ax.imshow(rel.T.values, aspect='auto', cmap='viridis', vmin=-lim, vmax=lim,
                   extent=ext, interpolation='nearest')
    cax = fig.add_subplot(gs[0, 1])
    cb = fig.colorbar(im, cax=cax, ticks=[-1, 0, 1])
    cb.ax.set_yticklabels(['×0.5' if lang == 'en' else '×0,5', '×1', '×2'])
    cb.set_label(t['ratio'], fontsize=7)
    cb.outline.set_linewidth(0.4)
    panel(ax, 'a')
    axes = [ax]
    for i, g in enumerate((1, 4, 8)):
        a = fig.add_subplot(gs[i + 1, 0], sharex=ax)
        m = masks[g].values.astype(float).T
        a.imshow(m, aspect='auto', cmap=ListedColormap(['#F2F2F2', '#000000']),
                 extent=ext, interpolation='nearest', vmin=0, vmax=1)
        head = t['gap1'] if g == 1 else t['gaps']
        a.text(0.0, 1.03, head.format(g=g, n=int(m.sum())), transform=a.transAxes,
               fontsize=7, va='bottom')
        panel(a, 'bcd'[i])
        axes.append(a)
    for i, a in enumerate(axes):
        a.set_yticks(range(len(districts)))
        # District names on panel a only: masks b-d follow the same row
        # order (stated in the legend).
        a.set_yticklabels(names if i == 0 else [''] * len(names), fontsize=6)
        a.tick_params(axis='y', length=0, pad=1.5)
        a.xaxis_date()
        # ticks on 1 January, year labels centred on mid-year
        a.xaxis.set_major_locator(mdates.YearLocator())
        a.xaxis.set_major_formatter(plt.NullFormatter())
        a.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=7, bymonthday=1))
        a.xaxis.set_minor_formatter(mdates.DateFormatter('%Y'))
        a.tick_params(axis='x', which='minor', length=0)
        for s in a.spines.values():
            s.set_visible(False)
    for a in axes[:-1]:
        plt.setp(a.get_xticklabels(which='both'), visible=False)
    axes[-1].set_xlabel(t['week'])
    save(fig, 'fig1_design', lang)
    return len(ref), str(ref.index[0].date()), str(ref.index[-1].date())


# ================================================================ FIGURE 2
# Nemenyi critical value q_0.05 = qtukey(0.95, k, Inf) / sqrt(2), the formula
# used by R/04_final_analysis.R. It is COMPUTED here (studentized range of k
# standard normals, infinite degrees of freedom), not copied from a table.
def _phi(z):
    return np.exp(-0.5 * z * z) / math.sqrt(2 * math.pi)


_erf = np.frompyfunc(math.erf, 1, 1)


def _Phi(z):
    return 0.5 * (1.0 + _erf(z / math.sqrt(2)).astype(float))


def ptukey_inf(q, k):
    """P(range of k iid N(0,1) <= q)."""
    z = np.linspace(-9.0, 9.0, 18001)
    f = k * _phi(z) * (_Phi(z + q) - _Phi(z)) ** (k - 1)
    return float(np.trapezoid(f, z) if hasattr(np, 'trapezoid') else np.trapz(f, z))


def qtukey_inf(p, k):
    lo, hi = 0.0, 10.0
    for _ in range(60):
        mid = (lo + hi) / 2
        if ptukey_inf(mid, k) < p:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def nemenyi_cd(k, n):
    return qtukey_inf(0.95, k) / math.sqrt(2) * math.sqrt(k * (k + 1) / (6.0 * n))


def fig_cd(lang):
    t = T[lang]
    r = load(MET, 'friedman_ranks_main.csv').sort_values('mean_rank').reset_index(drop=True)
    m = load(MET, 'metrics_long.csv')
    n_blocks = len(m[(m.gap_length == 1) & (m.variant == 'NoTime')][['rep', 'missing_rate']]
                   .drop_duplicates())
    k = len(r)
    cd = nemenyi_cd(k, n_blocks)
    ranks = r.mean_rank.values
    # maximal cliques: methods whose rank difference is below the CD
    cliques = []
    for i in range(k):
        j = max(jj for jj in range(k) if ranks[jj] - ranks[i] < cd)
        if j > i and not any(a <= i and j <= b for a, b in cliques):
            cliques.append((i, j))

    fig, ax = plt.subplots(figsize=(FULL, 52 * MM))
    ax.set_xlim(0.2, 10.8)
    ax.set_ylim(-4.25, 1.3)
    ax.axis('off')
    y0 = 0
    ax.plot([1, 10], [y0, y0], color='k', lw=0.8)
    for v in range(1, 11):
        ax.plot([v, v], [y0, y0 + 0.12], color='k', lw=0.8)
        ax.text(v, y0 + 0.2, str(v), ha='center', va='bottom', fontsize=7)
    ax.text(5.5, y0 + 0.62, t['mean_rank'], ha='center', va='bottom', fontsize=8)
    # critical difference bar
    ax.plot([1, 1 + cd], [y0 + 1.0, y0 + 1.0], color='k', lw=1.2)
    for v in (1, 1 + cd):
        ax.plot([v, v], [y0 + 0.92, y0 + 1.08], color='k', lw=1.2)
    ax.text(1 + cd / 2, y0 + 1.12, fmt(t['cd'].format(cd=cd), lang), ha='center',
            va='bottom', fontsize=7)
    half = (k + 1) // 2
    rows = np.linspace(-0.95, -3.35, half)
    for i, (meth, rk) in enumerate(zip(r.method, ranks)):
        fam = FAMILY[meth]
        left = i < half
        yy = rows[i] if left else rows[k - 1 - i]
        xe = 0.9 if left else 10.1
        ax.plot([rk, rk, xe], [y0, yy, yy], color='k', lw=0.6)
        ax.plot(rk, y0, marker=FAM_MRK[fam], color=FAM_COL[fam], ms=4.5, zorder=3,
                mec='k', mew=0.3)
        lab = '%s (%s)' % (t['method'][meth], fmt('%.2f' % rk, lang))
        ax.text(xe - 0.05 if left else xe + 0.05, yy, lab, ha='right' if left else 'left',
                va='center', fontsize=7)
    for c, (a, b) in enumerate(cliques):
        yy = -0.28 - 0.17 * c
        ax.plot([ranks[a] - 0.05, ranks[b] + 0.05], [yy, yy], color='k', lw=2.2,
                solid_capstyle='butt')
    ax.legend(handles=family_handles(lang), loc='lower center', ncol=4,
              bbox_to_anchor=(0.5, 0.0), handletextpad=0.3, columnspacing=1.2)
    save(fig, 'fig2_critical_difference', lang)
    return cd, n_blocks, cliques


# ================================================================ FIGURE 3
# The six combinations discussed in the text (Results, reversal table).
HIGHLIGHT = [('Kalman', 'NoTime'), ('LinearInterp', 'NoTime'), ('SeasonalDecomp', 'NoTime'),
             ('MissForest_ST', 'SeasonalTrend'), ('MissForest', 'SeasonalTrend'),
             ('MICE_precision', 'SeasonalTrend')]
HL_LS = {('Kalman', 'NoTime'): ('-', 'o'), ('LinearInterp', 'NoTime'): ('--', '^'),
         ('SeasonalDecomp', 'NoTime'): (':', 's'), ('MissForest_ST', 'SeasonalTrend'): ('-', 'D'),
         ('MissForest', 'SeasonalTrend'): ('-', 'o'), ('MICE_precision', 'SeasonalTrend'): ('--', '^')}


def fig_reversal(lang):
    t = T[lang]
    s = load(MET, 'gap_structure_sensitivity.csv')
    m = load(MET, 'metrics_long.csv')
    gaps = [1, 4, 8]
    fig = plt.figure(figsize=(FULL, 84 * MM))
    # panel a, area for direct labels, panel b
    a = fig.add_axes([0.075, 0.25, 0.31, 0.70])
    b = fig.add_axes([0.715, 0.25, 0.275, 0.70])
    hl = set(HIGHLIGHT)
    seen_mean = False
    for _, x in s.iterrows():
        key = (x.method, x.variant)
        if key in hl:
            continue
        if x.method == 'Mean':
            if seen_mean:
                continue
            seen_mean = True
        a.plot(gaps, [x.rank_gap1, x.rank_gap4, x.rank_gap8], color='#C8C8C8', lw=0.6,
               zorder=1)
    for key in HIGHLIGHT:
        x = s[(s.method == key[0]) & (s.variant == key[1])].iloc[0]
        y = [x.rank_gap1, x.rank_gap4, x.rank_gap8]
        col = FAM_COL[FAMILY[key[0]]]
        ls, mk = HL_LS[key]
        a.plot(gaps, y, color=col, ls=ls, marker=mk, ms=4, lw=1.3, zorder=3,
               mec='white', mew=0.4)
        lab = '%s, %s' % (t['method'][key[0]], key[1])
        a.text(8.35, y[-1], lab, va='center', ha='left', fontsize=6.5)
        a.text(0.65, y[0], '%d' % y[0], va='center', ha='right', fontsize=6.5)
    a.set_ylim(len(s) + 0.5, 0.3)
    a.set_yticks([1, 5, 10, 15, 20, len(s)])
    a.set_xlim(0.2, 8.3)
    a.set_xticks(gaps)
    a.set_xlabel(t['gap_axis'])
    a.set_ylabel(t['rank_axis'])
    panel(a, 'a', x=-0.12)

    for key in HIGHLIGHT:
        d = m[(m.method == key[0]) & (m.variant == key[1])]
        g = d.groupby('gap_length').MAE.agg(['mean', 'std', 'count'])
        ci = 1.96 * g['std'] / np.sqrt(g['count'])
        col = FAM_COL[FAMILY[key[0]]]
        ls, mk = HL_LS[key]
        b.errorbar(g.index, g['mean'], yerr=ci, color=col, ls=ls, marker=mk, ms=4, lw=1.3,
                   capsize=2, elinewidth=0.8, mec='white', mew=0.4,
                   label='%s, %s' % (t['method'][key[0]], key[1]))
    b.set_xticks(gaps)
    b.set_xlim(0.3, 8.7)
    b.set_xlabel(t['gap_axis'])
    b.set_ylabel(t['mae_axis'])
    h, l = b.get_legend_handles_labels()
    fig.legend(h, l, loc='lower center', ncol=3, bbox_to_anchor=(0.53, -0.01),
               fontsize=6.5, handlelength=3.0, columnspacing=1.5)
    panel(b, 'b', x=-0.2)
    save(fig, 'fig3_reversal', lang)


# ================================================================ FIGURE 4
def fig_mechanism(lang):
    t = T[lang]
    rev = load(MECH, 'reversal_four_mechanisms.csv')
    bias = load(MECH, 'bias_four_mechanisms.csv').set_index('method')
    order = (rev[rev.mechanism == 'MCAR'].sort_values('delta').method.tolist())
    ypos = {mth: i for i, mth in enumerate(order)}
    dodge = dict(zip(MECHS, [-0.27, -0.09, 0.09, 0.27]))
    fig, (a, b) = plt.subplots(1, 2, figsize=(FULL, 78 * MM), sharey=True,
                               gridspec_kw=dict(wspace=0.08))
    for mech in MECHS:
        st = MECH_STYLE[mech]
        d = rev[rev.mechanism == mech]
        a.plot(d.delta, [ypos[x] + dodge[mech] for x in d.method], ls='', marker=st['marker'],
               color=st['color'], mfc=st['mfc'], ms=4, mew=0.8, label=t['mech'][mech])
        b.plot([bias.loc[x, 'bias_' + mech] for x in order],
               [ypos[x] + dodge[mech] for x in order], ls='', marker=st['marker'],
               color=st['color'], mfc=st['mfc'], ms=4, mew=0.8)
    for ax in (a, b):
        ax.axvline(0, color='#555555', lw=0.6, ls='--', zorder=0)
        for i in range(len(order)):
            ax.axhline(i, color='#EBEBEB', lw=0.5, zorder=-1)
    a.set_yticks(range(len(order)))
    a.set_yticklabels([t['method'][x] for x in order])
    a.set_ylim(-0.6, len(order) - 0.4)
    a.set_xlabel(t['drank'])
    b.set_xlabel(t['bias'])
    b.tick_params(axis='y', length=0)
    panel(a, 'a', x=-0.03)
    panel(b, 'b', x=-0.03)
    a.legend(loc='lower center', bbox_to_anchor=(1.04, 1.03), ncol=4, handletextpad=0.2,
             columnspacing=1.4)
    save(fig, 'fig4_mechanism', lang)


# ================================================================ FIGURE 5
# labels: the combinations cited in the text (Results, peak detection);
# values are offsets in points (layout only).
LABEL5 = {('SeasonalDecomp', 'NoTime'): (10, 3), ('LinearInterp', 'NoTime'): (14, -9),
          ('KNN', 'SeasonalTrend'): (10, -4), ('Mean', 'NoTime'): (-4, 12)}


def spearman(x, y):
    return pd.Series(x).rank().corr(pd.Series(y).rank())


def fig_criteria(lang):
    t = T[lang]
    c = load(MET, 'criteria_comparison.csv')
    c = c.drop_duplicates(subset=['MAE', 'youden', 'degradation_h4', 'degradation_h52'])
    full = load(MET, 'criteria_comparison.csv')
    fig, axs = plt.subplots(1, 3, figsize=(FULL, 70 * MM), gridspec_kw=dict(wspace=0.5))
    fig.subplots_adjust(bottom=0.27, top=0.88, left=0.07, right=0.99)
    cols = [('youden', t['youden'], t['youden']),
            ('degradation_h4', t['h4'], t['deg']),
            ('degradation_h52', t['h52'], t['deg'])]
    rhos = {}
    for i, (ax, (col, head, lab)) in enumerate(zip(axs, cols)):
        for fam in FAM_ORDER:
            d = c[c.method.map(FAMILY) == fam]
            ax.plot(d.MAE, d[col], ls='', marker=FAM_MRK[fam], color=FAM_COL[fam], ms=4,
                    mec='white', mew=0.3, alpha=0.95)
        ax.set_ylabel(lab)
        ax.set_title(head, fontsize=7.5, loc='left', pad=4)
        if col != 'youden':
            ax.axhline(0, color='#555555', lw=0.6, ls='--', zorder=0)
            rho = spearman(full.MAE, full[col])
            rhos[col] = rho
            txt = t['rho'].format(r=rho).replace('-', '−')
            ax.text(0.97, 0.05 if col == 'degradation_h4' else 0.86, fmt(txt, lang),
                    transform=ax.transAxes, ha='right',
                    va='bottom' if col == 'degradation_h4' else 'top', fontsize=7)
        panel(ax, 'abc'[i], x=-0.2 if i == 0 else -0.3, y=1.03)
    fig.text(0.53, 0.13, t['mae_gap1'], ha='center', va='center', fontsize=8)
    for key, off in LABEL5.items():
        x = c[(c.method == key[0]) & (c.variant == key[1])].iloc[0]
        lab = t['method'][key[0]] + ('' if key[0] in ('Mean', 'Kalman', 'LinearInterp',
                                                       'SeasonalDecomp') else ', ' + key[1])
        axs[0].annotate(lab, (x.MAE, x.youden), xytext=off, textcoords='offset points',
                        fontsize=6, ha='left' if off[0] > 0 else 'right', va='center',
                        arrowprops=dict(arrowstyle='-', lw=0.4, color='#555555',
                                        shrinkA=0, shrinkB=2))
    fig.legend(handles=family_handles(lang), loc='lower center', ncol=4,
               bbox_to_anchor=(0.53, 0.0), handletextpad=0.3, columnspacing=1.2)
    save(fig, 'fig5_beyond_mae', lang)
    return rhos


# ======================================================= SUPPLEMENT S5
def figS_period(lang):
    t = T[lang]
    p = load(MET, 'period_effect.csv')
    p = p[(p.method != 'Mean') & (p.gap_length.isin([1, 8]))]
    g = p.groupby(['gap_length', 'method'])[['gain_seasonal_pct', 'gain_trend_pct']].mean()
    order = (g.loc[1].sum(axis=1).sort_values().index.tolist())
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 62 * MM), sharey=True,
                            gridspec_kw=dict(wspace=0.08))
    # Two bars per method (seasonal increment, then trend increment) rather
    # than a stack: an increment can be negative (KNN, 8-week gaps), and a
    # stacked bar would hide it.
    hgt = 0.36
    for ax, gap, L in zip(axs, (1, 8), 'ab'):
        d = g.loc[gap].reindex(order)
        y = np.arange(len(order))
        ax.barh(y + hgt / 2, d.gain_seasonal_pct, color='#4D4D4D', height=hgt, label=t['seas'])
        ax.barh(y - hgt / 2, d.gain_trend_pct, color='#BDBDBD', height=hgt, label=t['trend'])
        ax.axvline(0, color='k', lw=0.6)
        ax.set_yticks(y)
        ax.set_yticklabels([t['method'][x] for x in order])
        ax.set_xlabel(t['gain'])
        ax.set_title(t['ex1'] if gap == 1 else t['gapt'].format(g=gap), fontsize=7.5, loc='left')
        ax.tick_params(axis='y', length=0)
        panel(ax, L, x=-0.02, y=1.06)
    axs[1].legend(loc='lower right')
    xmin = min(0.0, min(a.get_xlim()[0] for a in axs))
    xmax = max(a.get_xlim()[1] for a in axs)
    for a in axs:
        a.set_xlim(xmin, xmax)
    save(fig, 'S5_figS1_period_effect', lang, where='sup')


# label offsets in points (layout only), for districts whose labels would collide
DISTRICT_LABEL_OFFSET = {'Ngaoundere_Urbain': (-4, 5, 'right'), 'Banyo': (0, 6, 'center'),
                         'Tignere': (5, -2, 'left'), 'Dang': (-4, 3, 'right'),
                         'Bankim': (4, -4, 'left')}


def figS_districts(lang):
    t = T[lang]
    g = load(MET, 'gain_by_district.csv')
    fig, ax = plt.subplots(figsize=(110 * MM, 78 * MM))
    ax.plot(g.mean_corr_other_districts, g.gain_MissForest_pct, ls='', marker='o', ms=4.5,
            color=FAM_COL['Cross-sectional'], mec='white', mew=0.4)
    for _, x in g.iterrows():
        dx, dy, ha = DISTRICT_LABEL_OFFSET.get(x.district, (4, 2, 'left'))
        ax.annotate(DISTRICT.get(x.district, x.district),
                    (x.mean_corr_other_districts, x.gain_MissForest_pct), xytext=(dx, dy),
                    textcoords='offset points', fontsize=6.3, ha=ha)
    rho = spearman(g.mean_corr_other_districts, g.gain_MissForest_pct)
    ax.text(0.03, 0.95, fmt(t['rho'].format(r=rho), lang) + ' (n = %d)' % len(g),
            transform=ax.transAxes, va='top', fontsize=7)
    ax.set_xlabel(t['corr'])
    ax.set_ylabel(t['gainmf'])
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0, x1 + 0.04)
    save(fig, 'S5_figS2_districts', lang, where='sup')
    return rho


# Illustrative example: district with the median mean correlation with the
# other districts, first 8-week gap of the mask (repetition 1, rate 10%),
# +/- 12 weeks. Rule fixed BEFORE looking at the curves, so as not to pick a
# convenient case. MissForest-ST is drawn (best combination at 8 weeks).
EX_COMBOS = [('Kalman', 'NoTime', '-', 'o'), ('LinearInterp', 'NoTime', '--', '^'),
             ('MissForest_ST', 'SeasonalTrend', '-', 'D'), ('Mean', 'NoTime', ':', 'X')]
EX_DIR = os.path.join(ROOT, 'output', 'imputed_examples')


def figS_example(lang):
    import pyreadr
    t = T[lang]
    ref = pd.read_csv(REF, parse_dates=['periodid']).set_index('periodid')
    ch = district_characteristics()
    ch = ch.sort_values('mean_corr_other_districts').reset_index(drop=True)
    d = ch.district.iloc[len(ch) // 2]                      # median district
    m8 = pyreadr.read_r(os.path.join(MASKS, 'mask_rep01_rate10_gap8.rds'))[None][d].values
    m1 = pyreadr.read_r(os.path.join(MASKS, 'mask_rep01_rate10_gap1.rds'))[None][d].values
    start = int(np.argmax(m8)); end = start
    while end + 1 < len(m8) and m8[end + 1]:
        end += 1
    lo, hi = max(0, start - 12), min(len(m8), end + 13)
    x = ref.index[lo:hi]
    step = (ref.index[1] - ref.index[0]) / 2
    fig, axs = plt.subplots(1, 2, figsize=(FULL, 66 * MM), sharey=True,
                            gridspec_kw=dict(wspace=0.14))
    fig.subplots_adjust(bottom=0.24, left=0.08, right=0.98)
    complete = True
    for ax, mask, gap, head, L in zip(axs, (m1, m8), (1, 8), (t['ex1'], t['ex8']), 'ab'):
        for i in np.where(mask[lo:hi])[0]:
            ax.axvspan(x[i] - step, x[i] + step, color='#E6E6E6', lw=0, zorder=0)
        ax.plot(x, ref[d].values[lo:hi], color='k', lw=1.1, marker='o', ms=2.2, zorder=4,
                label=t['ref'])
        for nm, v, ls, mk in EX_COMBOS:
            f = (os.path.join(ROOT, 'output', 'imputed', 'imp_rep01_rate10_%s_%s.csv' % (nm, v))
                 if gap == 1 else
                 os.path.join(EX_DIR, 'imp_rep01_rate10_gap8_%s_%s.csv' % (nm, v)))
            if not os.path.exists(f):
                complete = False
                continue
            z = pd.read_csv(f, parse_dates=['periodid']).set_index('periodid')[d].values
            # check: outside the mask, the imputed series is the reference
            assert np.allclose(z[~mask], ref[d].values[~mask]), f
            y = np.where(mask, z, np.nan)[lo:hi]
            col = FAM_COL[FAMILY[nm]]
            if gap == 1:
                ax.plot(x, y, ls='', marker=mk, color=col, ms=4.5, mec='white', mew=0.4,
                        zorder=5, label=t['method'][nm])
            else:
                yy = y.copy()                 # join to the neighbouring observed weeks
                for j in (start - lo - 1, end - lo + 1):
                    if 0 <= j < len(yy):
                        yy[j] = ref[d].values[lo + j]
                ax.plot(x, yy, ls=ls, lw=1.2, color=col, zorder=5)
                ax.plot(x, y, ls='', marker=mk, color=col, ms=4, mec='white', mew=0.4,
                        zorder=6, label=t['method'][nm])
        if gap == 8 and not complete:
            ax.text(0.5, 0.9, t['todo'], transform=ax.transAxes, ha='center', fontsize=7,
                    color='#B00000')
        ax.set_title(head, fontsize=7.5, loc='left', pad=4)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
        panel(ax, L, x=-0.02, y=1.03)
    axs[1].tick_params(axis='y', length=0)
    axs[0].set_ylabel(t['cases'].format(d=DISTRICT.get(d, d)))
    h, l = axs[0].get_legend_handles_labels()
    from matplotlib.patches import Patch
    h.append(Patch(color='#E6E6E6')); l.append(t['masked'])
    fig.legend(h, l, loc='lower center', ncol=6, bbox_to_anchor=(0.52, 0.02),
               handletextpad=0.4, columnspacing=1.2, fontsize=6.5)
    save(fig, 'fig_example_series', lang)
    if not complete:
        print('   WARNING: panel b incomplete, run R/10_example_series.R')
    return d, str(ref.index[lo].date()), str(ref.index[hi - 1].date()), \
        str(ref.index[start].date()), str(ref.index[end].date()), complete


if __name__ == '__main__':
    for lang in ('en', 'fr'):
        print('[%s]' % lang)
        info1 = fig_design(lang)
        cd, n, cl = fig_cd(lang)
        fig_reversal(lang)
        fig_mechanism(lang)
        rhos = fig_criteria(lang)
        figS_period(lang)
        rho_s2 = figS_districts(lang)
        ex = figS_example(lang)
    print('Reference: %d weeks, %s to %s' % info1)
    print('Nemenyi CD recomputed: %.3f (k = 10, N = %d), cliques %s' % (cd, n, cl))
    print('Fig. 5 rho (24 combinations): h4 %+.3f, h52 %+.3f'
          % (rhos['degradation_h4'], rhos['degradation_h52']))
    print('Fig. S2 rho: %+.3f' % rho_s2)
    print('Example: district %s, window %s to %s, 8-week gap %s to %s, panel b complete: %s'
          % ex)
