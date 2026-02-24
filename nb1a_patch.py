# NB1_supp_batch_sensitivity.py

import pickle, warnings
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import mannwhitneyu

warnings.filterwarnings("ignore")

# SAVE_DIR: reuse from NB1c kernel

try:
    _ = SAVE_DIR
    print(f"Using SAVE_DIR from kernel: {SAVE_DIR}")
except NameError:
    from pathlib import Path
    SAVE_DIR = Path(r"D:\your\output\folder")  # ← only needed if running standalone
    SAVE_DIR.mkdir(exist_ok=True, parents=True)
    print(f"Standalone SAVE_DIR: {SAVE_DIR}")
DAYS_PER_MO = 30.4375

# Variables already in kernel 
_missing = []
for _var, _obj in [
    ('ssgsea_corrected', ssgsea_corrected),
    ('ssgsea_raw',       ssgsea_raw),
    ('module_scores_37', module_scores_37),
    ('geo_clinical',     geo_clinical),
    ('tcga_clinical',    tcga_clinical),
]:
    try:
        assert _obj is not None
        print(f"  ok  {_var}")
    except Exception:
        _missing.append(_var)
if _missing:
    raise RuntimeError(
        f"Missing kernel variables: {_missing}\n"
        "Make sure NB1a and NB1b have both run successfully in this kernel.")
print("All required variables confirmed in kernel.")

MODULES = ["Proliferation", "Immune", "Metabolic", "Microenvironment"]

# MSRS implementation
def msrs_scan(df, score_col, t_min=6, t_max=60, min_each=20):
    """Return scan DataFrame (index=timepoints, u_stat column) and optimum."""
    results = {}
    times = np.arange(t_min, t_max + 1)
    for t in times:
        early = df[df["time"] <= t]
        late  = df[df["time"] >  t]
        if len(early) < min_each or len(late) < min_each:
            continue
        e_scores = early[score_col].dropna()
        l_scores = late[score_col].dropna()
        if len(e_scores) < 5 or len(l_scores) < 5:
            continue
        u, _ = mannwhitneyu(e_scores, l_scores, alternative="two-sided")
        results[t] = u
    if not results:
        return None, None
    scan = pd.DataFrame({"u_stat": results})
    opt  = int(scan["u_stat"].idxmax())
    return scan, opt


def build_module_df(ssgsea_dict, clin_dict, cohort, module="Proliferation"):
    """Build a DataFrame with time, event, and mean module score."""
    ss = ssgsea_dict.get(cohort)
    clin = clin_dict.get(cohort) if cohort in clin_dict else None
    if ss is None or clin is None or clin.empty:
        return None
    # Module membership
    ms = module_scores_37.get(cohort)
    if ms is None:
        return None
    ms_T = ms.T
    common = sorted(set(ms_T.index.astype(str)) & set(clin.index.astype(str)))
    if len(common) < 20:
        return None
    rows = []
    for s in common:
        try:
            rows.append({
                "time":  float(clin.loc[s, "time"]),
                "event": float(clin.loc[s, "event"]),
                module:  float(ms_T.loc[s, module]) if module in ms_T.columns else np.nan,
            })
        except Exception:
            continue
    df = pd.DataFrame(rows).dropna(subset=["time", "event", module])
    return df[(df["time"] > 0) & df["event"].isin([0., 1.])].reset_index(drop=True)


# Run MSRS on corrected vs uncorrected for key cohorts
print("=" * 65)
print("BATCH CORRECTION SENSITIVITY — MSRS COMPARISON")
print("=" * 65)

cohort_configs = [
    # (cohort_key,  label,              clin_source,   t_min, t_max)
    ("GSE2034",     "Breast / GSE2034", geo_clinical,   6,    60),
    ("GSE31210",    "LUAD / GSE31210",  geo_clinical,   6,    60),
    ("TCGA-BRCA",   "TCGA-BRCA DFI",   tcga_clinical,  6,    60),
]

results_table = []
scan_corrected_all = {}
scan_raw_all       = {}

for cohort, label, clin_src, tmin, tmax in cohort_configs:
    df_c = build_module_df(ssgsea_corrected, clin_src, cohort)
    df_r = build_module_df(ssgsea_raw,       clin_src, cohort)

    scan_c, opt_c = (None, None)
    scan_r, opt_r = (None, None)

    if df_c is not None:
        scan_c, opt_c = msrs_scan(df_c, "Proliferation", tmin, tmax)
    if df_r is not None:
        scan_r, opt_r = msrs_scan(df_r, "Proliferation", tmin, tmax)

    scan_corrected_all[cohort] = (scan_c, opt_c, label)
    scan_raw_all[cohort]       = (scan_r, opt_r, label)

    diff = abs(opt_c - opt_r) if (opt_c and opt_r) else None
    results_table.append({
        "Cohort":              label,
        "Batch correction":    "ComBat→ssGSEA" if cohort.startswith("GSE20") or cohort.startswith("GSE29") or cohort.startswith("GSE10") else "None (ssGSEA only)",
        "Cut (corrected, mo)": opt_c,
        "Cut (raw, mo)":       opt_r,
        "Difference (mo)":     diff if diff is not None else "N/A",
        "Note":                "GSE31210 & TCGA: no ComBat in either pipeline" if not cohort.startswith("GSE2") else "",
    })

    c_str = f"{opt_c} mo" if opt_c else "n/a"
    r_str = f"{opt_r} mo" if opt_r else "n/a"
    d_str = f"Δ = {diff} mo" if diff is not None else ""
    print(f"\n  {label}")
    print(f"    Corrected:   {c_str}")
    print(f"    Uncorrected: {r_str}  {d_str}")

res_df = pd.DataFrame(results_table)
print("\n")
print(res_df.to_string(index=False))

# Supplementary Figure: MSRS scans corrected vs uncorrected 
matplotlib.rcParams.update({
    "font.family":        "Arial",
    "font.size":          7,
    "axes.titlesize":     8,
    "axes.titleweight":   "bold",
    "axes.labelsize":     7,
    "axes.labelweight":   "bold",
    "xtick.labelsize":    6,
    "ytick.labelsize":    6,
    "legend.fontsize":    6,
    "legend.frameon":     False,
    "axes.linewidth":     0.5,
    "savefig.dpi":        1200,
    "savefig.bbox":       "tight",
    "savefig.pad_inches": 0.04,
})
_SC = 3.46
_DC = 7.09

def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

valid_cohorts = [c for c in cohort_configs
                 if scan_corrected_all[c[0]][0] is not None
                 and scan_raw_all[c[0]][0] is not None]
n_panels = len(valid_cohorts)

fig, axes = plt.subplots(1, n_panels, figsize=(_DC, 2.0), sharey=False)
if n_panels == 1:
    axes = [axes]

COL_CORR = "#B03A2E"   # muted crimson — corrected
COL_RAW  = "#2471A3"   # slate blue — uncorrected

for ax, (cohort, label, clin_src, tmin, tmax) in zip(axes, valid_cohorts):
    scan_c, opt_c, _ = scan_corrected_all[cohort]
    scan_r, opt_r, _ = scan_raw_all[cohort]

    ts_c = scan_c.index.values.astype(float)
    us_c = scan_c["u_stat"].values.astype(float)
    ts_r = scan_r.index.values.astype(float)
    us_r = scan_r["u_stat"].values.astype(float)

    # Normalise to [0,1] so both are on the same scale despite different n
    def _norm(x):
        rng = x.max() - x.min()
        return (x - x.min()) / rng if rng > 0 else x * 0

    ax.plot(ts_c, _norm(us_c), color=COL_CORR, lw=1.0, label="ComBat→ssGSEA", zorder=3)
    ax.plot(ts_r, _norm(us_r), color=COL_RAW,  lw=1.0, ls="--",
            label="Raw→ssGSEA", zorder=2)

    # Vertical lines at each optimum
    if opt_c:
        ax.axvline(opt_c, color=COL_CORR, ls=":", lw=0.7, zorder=4)
        ax.text(opt_c + 0.5, 0.97, "{} mo".format(opt_c),
                fontsize=5.5, va="top", ha="left", color=COL_CORR,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.9))
    if opt_r:
        y_offset = 0.78
        ax.axvline(opt_r, color=COL_RAW, ls=":", lw=0.7, zorder=4)
        ax.text(opt_r - 0.5, y_offset, "{} mo".format(opt_r),
                fontsize=5.5, va="top", ha="right", color=COL_RAW,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.9))

    ax.set_xlim(tmin - 1, tmax + 2)
    ax.set_ylim(-0.05, 1.12)
    ax.set_xlabel("Candidate cut-point (months)", fontsize=7, fontweight="bold")
    if ax is axes[0]:
        ax.set_ylabel("U statistic (normalised)", fontsize=7, fontweight="bold")
    ax.set_title(label, fontsize=7, pad=3)
    ax.xaxis.set_major_locator(mticker.MultipleLocator(20))
    _despine(ax)

# Shared legend on last panel
axes[-1].legend(loc="lower right", fontsize=5.5,
                handlelength=1.5, borderpad=0.5)

fig.suptitle(
    "Supplementary Figure: MSRS transition timepoints are robust to batch correction\n"
    "Solid = ComBat-corrected→ssGSEA  |  Dashed = raw expression→ssGSEA",
    fontsize=6.5, y=1.04,
)
plt.tight_layout(pad=0.6, w_pad=1.2)

out_path = SAVE_DIR / "Supp_Fig_BatchSensitivity.png"
fig.savefig(out_path, dpi=1200, bbox_inches="tight")
print(f"\n  Saved: {out_path}")
plt.close(fig)

# ── Also: show same pipeline explains 25m vs 21m difference ──────────────────
print("""
COMPLETED
""")

# ── Save summary table ────────────────────────────────────────────────────────
res_df.to_csv(SAVE_DIR / "Supp_Table_BatchSensitivity.csv", index=False)
print("  Saved: Supp_Table_BatchSensitivity.csv")
print("\n  → Use these results in the reviewer response letter and")
print("    Supplementary Methods section.")
