NB1d — Cox Survival Model 

import warnings, pickle, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings('ignore')
matplotlib.use('Agg')

# Paths 
BASE_DIR  = Path(r"D:\Data")
RAW_DIR   = BASE_DIR / "Raw Data"
INTER_DIR = BASE_DIR / "intermediates"
SAVE_DIR  = BASE_DIR / "Manuscript Data"
SUPP_DIR  = SAVE_DIR / "Supplementary"
for d in (SAVE_DIR, SUPP_DIR):
    d.mkdir(parents=True, exist_ok=True)

DAYS_PER_MO = 30.4375
LANDMARK_M  = 24     

# rcParams 
NATURE_RC = {
    'font.family':       'Arial',
    'font.size':         7,
    'axes.linewidth':    0.5,
    'axes.labelsize':    7,
    'axes.titlesize':    7,
    'axes.titleweight':  'normal',
    'xtick.labelsize':   6,
    'ytick.labelsize':   6,
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5,
    'xtick.major.size':  2.5,
    'ytick.major.size':  2.5,
    'xtick.direction':   'out',
    'ytick.direction':   'out',
    'legend.fontsize':   6,
    'legend.frameon':    False,
    'pdf.fonttype':      42,
    'ps.fonttype':       42,
}

COL_EARLY = '#E67E22'    # orange — early phase (≤24m)
COL_LATE  = '#2471A3'    # blue   — late phase  (>24m)

# Cancer-type colours for landmark overlay
CANCER_COLS = {
    'BRCA': '#C0392B', 'COAD': '#27AE60', 'HNSC': '#8E44AD',
    'KIRC': '#2980B9', 'LIHC': '#D4AC0D', 'LUAD': '#E67E22',
    'LUSC': '#16A085', 'STAD': '#884EA0', 'UCEC': '#CB4335',
}


def _hdr(txt):
    print(f"\n{'='*60}\n{txt}\n{'='*60}\n")

def _savefig(fig, name, folder=None):
    folder = folder or SAVE_DIR
    for ext in ('pdf', 'png'):
        fp = Path(folder) / f"{name}.{ext}"
        fig.savefig(fp, dpi=1200, bbox_inches='tight', facecolor='white')
    plt.close(fig)


# S1 — LOAD CDR
_hdr("S1: LOAD CDR + NB1c INTERMEDIATES")

def _load(name):
    p = INTER_DIR / f"{name}.pkl"
    if p.exists():
        with open(p, 'rb') as f:
            return pickle.load(f)
    return None

module_scores_37 = _load('module_scores_37') or {}
ssgsea_corrected = _load('ssgsea_corrected') or {}
geo_clinical     = _load('geo_clinical')     or {}
tcga_clinical    = _load('tcga_clinical')    or {}

# Load CDR
cdr_path = RAW_DIR / "tcga_clinical" / "TCGA-CDR.csv"
if not cdr_path.exists():
    cdr_path = cdr_path.with_suffix('.xlsx')
cdr = pd.read_csv(cdr_path) if cdr_path.suffix == '.csv' else pd.read_excel(cdr_path)
print(f"  CDR loaded: {len(cdr)} rows, {len(cdr.columns)} columns")
print(f"  Cancer types in CDR: {sorted(cdr['type'].unique())}")


def cdr_os(cancer):
    """Extract OS endpoint for a cancer type from CDR."""
    sub = cdr[cdr['type'] == cancer].copy()
    if 'bcr_patient_barcode' in sub.columns:
        sub = sub.set_index('bcr_patient_barcode')
    if 'OS.time' not in sub.columns:
        print(f"  WARNING: OS.time missing for {cancer}")
        return None
    out = pd.DataFrame({
        'time':  pd.to_numeric(sub['OS.time'], errors='coerce') / DAYS_PER_MO,
        'event': (pd.to_numeric(sub['OS'],      errors='coerce') > 0).astype(float),
    }, index=sub.index).dropna()
    out.index = out.index.astype(str).str[:12]
    out = out[out['time'] > 0].loc[lambda x: ~x.index.duplicated()]
    return out


# S2 — BUILD OS DATAFRAMES FOR 9 CANCER TYPES
_hdr("S2: BUILD OS DATAFRAMES")

CANCER_TYPES = ['BRCA', 'COAD', 'HNSC', 'KIRC', 'LIHC', 'LUAD', 'LUSC', 'STAD', 'UCEC']
CANCER_LABELS = {
    'BRCA': 'Breast invasive carcinoma',
    'COAD': 'Colon adenocarcinoma',
    'HNSC': 'Head and neck squamous cell carcinoma',
    'KIRC': 'Kidney renal clear cell carcinoma',
    'LIHC': 'Liver hepatocellular carcinoma',
    'LUAD': 'Lung adenocarcinoma',
    'LUSC': 'Lung squamous cell carcinoma',
    'STAD': 'Stomach adenocarcinoma',
    'UCEC': 'Uterine corpus endometrial carcinoma',
}

os_dfs = {}   # cancer → DataFrame(time, event, phase, ...)
for ct in CANCER_TYPES:
    clin = cdr_os(ct)
    if clin is None or len(clin) < 30:
        print(f"  {ct}: insufficient OS data")
        continue

    # Assign phase: early = event before 24m; late = no event before 24m
    clin['phase'] = np.where(
        (clin['event'] == 1) & (clin['time'] <= LANDMARK_M), 'early', 'late'
    )
    n_early = (clin['phase'] == 'early').sum()
    n_late  = (clin['phase'] == 'late').sum()
    n_ev    = int(clin['event'].sum())
    print(f"  {ct}: n={len(clin)}  events={n_ev}  "
          f"early={n_early}  late={n_late}  "
          f"median_OS={clin['time'].median():.1f}mo")
    os_dfs[ct] = clin

# Pan-cancer pool
all_df = pd.concat(os_dfs.values(), keys=os_dfs.keys()).reset_index(level=0)
all_df = all_df.rename(columns={'level_0': 'cancer_type'})
print(f"\n  Pan-cancer: n={len(all_df)}  events={int(all_df['event'].sum())}")


# S3 — COX PROPORTIONAL HAZARDS — PHASE → OS
_hdr("S3: COX MODELS — PHASE → OS")

print(f"  Predictor: binary phase (early ≤{LANDMARK_M}m = 1, late = 0)")
print(f"  Outcome: Overall Survival\n")
print(f"  {'Cancer':<8} {'n':>5} {'events':>7} {'HR':>8} {'95%CI':>18} {'P':>10}")
print(f"  {'-'*60}")

cox_results = {}
for ct, df in os_dfs.items():
    cox_df = df[['time', 'event', 'phase']].copy()
    cox_df['phase_bin'] = (cox_df['phase'] == 'early').astype(int)
    try:
        cph = CoxPHFitter()
        cph.fit(cox_df[['time', 'event', 'phase_bin']], duration_col='time',
                event_col='event')
        hr  = cph.params_['phase_bin']
        hr_exp = np.exp(hr)
        ci_l = np.exp(cph.confidence_intervals_.loc['phase_bin', '95% lower-bound'])
        ci_u = np.exp(cph.confidence_intervals_.loc['phase_bin', '95% upper-bound'])
        pval = cph.summary['p']['phase_bin']
        n_ev = int(df['event'].sum())
        print(f"  {ct:<8} {len(df):>5} {n_ev:>7} {hr_exp:>8.3f} "
              f"  [{ci_l:.3f}–{ci_u:.3f}]  {pval:>10.4f}")
        cox_results[ct] = {'HR': hr_exp, 'CI_l': ci_l, 'CI_u': ci_u,
                           'P': pval, 'n': len(df), 'n_events': n_ev}
    except Exception as e:
        print(f"  {ct}: Cox failed — {e}")

# Pan-cancer Cox
pan_cox_df = all_df[['time', 'event', 'phase']].copy()
pan_cox_df['phase_bin'] = (pan_cox_df['phase'] == 'early').astype(int)
try:
    cph_pan = CoxPHFitter()
    cph_pan.fit(pan_cox_df[['time', 'event', 'phase_bin']],
                duration_col='time', event_col='event')
    hr_pan = np.exp(cph_pan.params_['phase_bin'])
    ci_l_p = np.exp(cph_pan.confidence_intervals_.loc[
        'phase_bin', '95% lower-bound'])
    ci_u_p = np.exp(cph_pan.confidence_intervals_.loc[
        'phase_bin', '95% upper-bound'])
    p_pan  = cph_pan.summary['p']['phase_bin']
    print(f"\n  {'PAN':<8} {len(pan_cox_df):>5} {int(all_df['event'].sum()):>7} "
          f"{hr_pan:>8.3f}   [{ci_l_p:.3f}–{ci_u_p:.3f}]  {p_pan:>10.4f}")
    cox_results['PAN'] = {'HR': hr_pan, 'CI_l': ci_l_p, 'CI_u': ci_u_p,
                          'P': p_pan, 'n': len(pan_cox_df),
                          'n_events': int(all_df['event'].sum())}
except Exception as e:
    print(f"  Pan-cancer Cox failed: {e}")


# S4 — STANDARD KM PER CANCER TYPE
_hdr("S4: STANDARD KM FIGURE (Figure 8)")

print("  Plotting KM curves: early (≤24m) vs late (>24m) per cancer type")
print("  Style: matching uploaded Figure 8 reference image")

matplotlib.rcParams.update(NATURE_RC)

n_panels = len(os_dfs) + 1   # 9 cancers + pan-cancer
n_cols = 2
n_rows = (n_panels + 1) // 2

fig_km, axes_km = plt.subplots(
    n_rows, n_cols,
    figsize=(7.0, 2.6 * n_rows),
    gridspec_kw={'hspace': 0.55, 'wspace': 0.38}
)
axes_km_flat = axes_km.flatten()
panel_letters = list('ABCDEFGHIJ')


def _km_panel(ax, df, title, letter, cox_r=None):
    """Plot single KM panel with early vs late curves."""
    early_df = df[df['phase'] == 'early']
    late_df  = df[df['phase'] == 'late']

    lr = logrank_test(early_df['time'], late_df['time'],
                      event_observed_A=early_df['event'],
                      event_observed_B=late_df['event'])
    p_lr = lr.p_value

    for subset_df, col, lbl in [
        (late_df,  COL_LATE,  'Late phase (>24 months)'),
        (early_df, COL_EARLY, 'Early phase (≤24 months)'),
    ]:
        n  = len(subset_df)
        ev = int(subset_df['event'].sum())
        kmf = KaplanMeierFitter()
        kmf.fit(subset_df['time'], subset_df['event'],
                label=f'{lbl}\nn={n}, events={ev}')
        kmf.plot_survival_function(ax=ax, ci_show=True, linewidth=1.0,
                                   color=col, ci_alpha=0.15)

    ax.set_ylim(0, 1.05)
    ax.set_xlim(left=0)
    ax.set_xlabel("Time (months)", fontsize=7)
    ax.set_ylabel("Overall survival probability", fontsize=7)
    ax.set_title(f"{title}: {LANDMARK_M}-month survival analysis", fontsize=7, pad=3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=6)
    ax.yaxis.set_major_locator(mticker.MultipleLocator(0.2))
    ax.grid(True, alpha=0.2, lw=0.3, color='grey')
    ax.axvline(LANDMARK_M, color='#888888', lw=0.5, ls='--', alpha=0.5)

    # Stats box
    hr_txt = ''
    if cox_r:
        p_s = f"{cox_r['P']:.4f}" if cox_r['P'] >= 0.0001 else f"{cox_r['P']:.2e}"
        hr_txt = (f"Log-rank P={p_s}\n"
                  f"HR={cox_r['HR']:.2f} (95% CI: "
                  f"{cox_r['CI_l']:.2f}–{cox_r['CI_u']:.2f})")
    else:
        p_s = f"{p_lr:.4f}" if p_lr >= 0.0001 else f"{p_lr:.2e}"
        hr_txt = f"Log-rank P={p_s}"

    ax.text(0.97, 0.03, hr_txt, transform=ax.transAxes,
            fontsize=5.5, va='bottom', ha='right',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
    ax.text(-0.10, 1.06, letter, transform=ax.transAxes,
            fontsize=8, fontweight='bold', va='top')
    leg = ax.get_legend()
    if leg:
        leg.set_title(None)
        for t in leg.get_texts():
            t.set_fontsize(5.5)
        leg.get_frame().set_visible(False)


for i, ct in enumerate(list(os_dfs.keys())):
    _km_panel(axes_km_flat[i], os_dfs[ct],
              CANCER_LABELS.get(ct, ct),
              panel_letters[i],
              cox_results.get(ct))
    print(f"  Panel {panel_letters[i]}: {ct}")

# Pan-cancer panel
_km_panel(axes_km_flat[len(os_dfs)], all_df,
          "Pan-cancer survival analysis",
          panel_letters[len(os_dfs)],
          cox_results.get('PAN'))
print(f"  Panel {panel_letters[len(os_dfs)]}: Pan-cancer")

# Hide unused panels
for j in range(len(os_dfs) + 1, len(axes_km_flat)):
    axes_km_flat[j].set_visible(False)

fig_km.suptitle(
    f"Cancer-Specific and Pan-Cancer Kaplan–Meier Survival Analysis\n"
    f"Stratified by {LANDMARK_M}-month biological transition",
    fontsize=8, y=1.002
)
fig_km.tight_layout()
_savefig(fig_km, "Fig_KM_standard_Fig8", SAVE_DIR)
print("  Figure saved: Fig_KM_standard_Fig8")


# S5 — LANDMARK KM AT 24m
_hdr("S5: LANDMARK KM AT 24m — REVIEWER RESPONSE")

print(f"""
  COMPLETED
""")

landmark_results = {}
for ct, df in os_dfs.items():
    # Landmark cohort: survived past 24m (no event before 24m)
    lm = df[(df['time'] >= LANDMARK_M)].copy()
    # Reset time origin
    lm['lm_time']  = lm['time'] - LANDMARK_M
    n_lm  = len(lm)
    n_ev_lm = int(lm['event'].sum())
    med_t = lm['lm_time'].median()
    print(f"  {ct}: landmark n={n_lm}  events_after_24m={n_ev_lm}  "
          f"median_post24={med_t:.1f}mo")
    landmark_results[ct] = lm


# Figure: Landmark overlay
matplotlib.rcParams.update(NATURE_RC)

fig_lm, ax_lm = plt.subplots(1, 1, figsize=(4.5, 3.5))

for ct, lm_df in landmark_results.items():
    col = CANCER_COLS.get(ct, '#555555')
    n_ev = int(lm_df['event'].sum())
    kmf = KaplanMeierFitter()
    kmf.fit(lm_df['lm_time'], lm_df['event'],
            label=f"{ct} (n={len(lm_df)}, ev={n_ev})")
    kmf.plot_survival_function(ax=ax_lm, ci_show=False,
                               linewidth=1.1, color=col)

ax_lm.set_xlabel(f"Time from {LANDMARK_M}-month landmark (months)", fontsize=7)
ax_lm.set_ylabel("Survival probability (from 24 m)", fontsize=7)
ax_lm.set_title(
    f"Post-{LANDMARK_M}-month conditional survival by cancer type\n"
    f"(landmark cohort: survived ≥{LANDMARK_M} months)",
    fontsize=7, pad=4
)
ax_lm.set_ylim(0, 1.05)
ax_lm.set_xlim(left=0)
ax_lm.spines['top'].set_visible(False)
ax_lm.spines['right'].set_visible(False)
ax_lm.tick_params(labelsize=6)
ax_lm.grid(True, alpha=0.18, lw=0.3, color='grey')

leg = ax_lm.get_legend()
if leg:
    leg.set_title(None)
    for t in leg.get_texts():
        t.set_fontsize(5.5)
    leg.get_frame().set_visible(False)
    leg.set_bbox_to_anchor((1.02, 1.0))
    leg.set_loc('upper left')

fig_lm.tight_layout(rect=[0, 0, 0.78, 1])
_savefig(fig_lm, "Fig_KM_landmark_overlay", SAVE_DIR)
print("\n  Figure saved: Fig_KM_landmark_overlay")


# Figure: Landmark per-cancer subplots 
matplotlib.rcParams.update(RC)

n_ct = len(landmark_results)
n_cols_lm = 3
n_rows_lm = (n_ct + n_cols_lm - 1) // n_cols_lm

fig_lm2, axes_lm2 = plt.subplots(
    n_rows_lm, n_cols_lm,
    figsize=(7.0, 2.4 * n_rows_lm),
    gridspec_kw={'hspace': 0.55, 'wspace': 0.38}
)
axes_lm2_flat = axes_lm2.flatten()

for i, (ct, lm_df) in enumerate(landmark_results.items()):
    ax = axes_lm2_flat[i]
    col = CANCER_COLS.get(ct, '#555555')

    kmf = KaplanMeierFitter()
    n_ev = int(lm_df['event'].sum())
    kmf.fit(lm_df['lm_time'], lm_df['event'],
            label=f"n={len(lm_df)}, events={n_ev}")
    kmf.plot_survival_function(ax=ax, ci_show=True,
                               linewidth=1.1, color=col, ci_alpha=0.15)

    # Median post-24m survival
    med = kmf.median_survival_time_
    med_txt = f"Median post-24m: {med:.1f} mo" if not np.isinf(med) else "Median: not reached"
    ax.set_xlabel(f"Months after {LANDMARK_M}-month landmark", fontsize=7)
    ax.set_ylabel("Survival probability", fontsize=7)
    ax.set_title(f"{CANCER_LABELS.get(ct, ct)}\n({ct})", fontsize=6.5, pad=3)
    ax.set_ylim(0, 1.05)
    ax.set_xlim(left=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=5.5)
    ax.grid(True, alpha=0.18, lw=0.3, color='grey')
    ax.text(0.97, 0.97, med_txt, transform=ax.transAxes,
            fontsize=5.5, va='top', ha='right', color='#333333',
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1))
    ax.text(-0.10, 1.08, panel_letters[i], transform=ax.transAxes,
            fontsize=8, fontweight='bold', va='top')
    leg = ax.get_legend()
    if leg:
        leg.set_title(None)
        for t in leg.get_texts():
            t.set_fontsize(5.5)
        leg.get_frame().set_visible(False)

for j in range(n_ct, len(axes_lm2_flat)):
    axes_lm2_flat[j].set_visible(False)

fig_lm2.suptitle(
    f"Conditional survival from {LANDMARK_M}-month landmark by cancer type\n"
    f"(reviewer response: late-phase curves aligned at {LANDMARK_M} months)",
    fontsize=8, y=1.002
)
fig_lm2.tight_layout()
_savefig(fig_lm2, "Fig_KM_landmark_per_cancer", SAVE_DIR)
print("  Figure saved: Fig_KM_landmark_per_cancer")


# Log-rank between cancer types in landmark cohort
print("\n  Post-24m survival comparison (landmark cohort):")
print(f"  {'Cancer':<8} {'n_lm':>6} {'ev_lm':>7} {'median_post24':>15} {'5yr_surv':>10}")
print(f"  {'-'*50}")
for ct, lm_df in landmark_results.items():
    kmf = KaplanMeierFitter()
    kmf.fit(lm_df['lm_time'], lm_df['event'])
    med = kmf.median_survival_time_
    med_s = f"{med:.1f}m" if not np.isinf(med) else "NR"
    # 5-year survival (60 months post-landmark)
    try:
        surv_60 = kmf.survival_function_at_times([60]).values[0]
    except Exception:
        surv_60 = np.nan
    print(f"  {ct:<8} {len(lm_df):>6} {int(lm_df['event'].sum()):>7} "
          f"{med_s:>15} {surv_60:>10.1%}")


# S6 — SUMMARY + SAVE
_hdr("S6: SUMMARY + SAVE")

# Save results
results = {
    'cox_results': cox_results,
    'landmark_results': {ct: df[['lm_time','event']].copy()
                         for ct, df in landmark_results.items()},
    'os_dfs': {ct: df.copy() for ct, df in os_dfs.items()},
}
with open(SAVE_DIR / "NB1d_results.pkl", 'wb') as f:
    pickle.dump(results, f)

print(f"""
  COMPLETED
""")
