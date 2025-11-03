
"""
================================================================================
24-MONTH BIOLOGICAL TRANSITION: VISUALIZATION CODE
Generate Publication-Quality Figures
================================================================================

Description:
    Generates Kaplan-Meier curves, forest plots, and summary tables for
    multi-cancer validation results.

Authors: Shanghai Jiao Tong University - Fuxilab
Date: November 2025
License: MIT
Repository: https://github.com/Sjtu-Fuxilab/COSMOS
================================================================================
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import warnings
warnings.filterwarnings('ignore')

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

# Nature journal style
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['pdf.fonttype'] = 42

COLORS = {
    'early_phase': '#D55E00',
    'late_phase': '#0072B2',
    'significant': '#009E73',
    'nonsignificant': '#999999',
}


def parse_arguments():
    parser = argparse.ArgumentParser(description='Generate validation figures')
    parser.add_argument('--results_dir', type=str, required=True)
    parser.add_argument('--output_dir', type=str, default='./figures')
    parser.add_argument('--dpi', type=int, default=1200)
    return parser.parse_args()


def load_results(results_dir):
    results_path = Path(results_dir)
    patient_data = pd.read_csv(results_path / "results" / "patient_level_data.csv")
    cancer_results = pd.read_csv(results_path / "results" / "cancer_specific_results.csv")
    print(f"✅ Loaded: {len(patient_data)} patients, {len(cancer_results)} cancer types")
    return patient_data, cancer_results


def create_kaplan_meier_figure(patient_data, cancer_results, output_dir, dpi):
    print(f"\n📊 Creating Kaplan-Meier curves...")

    cancer_types = sorted(patient_data['cancer_type'].unique())
    fig = plt.figure(figsize=(18, 16))

    cancer_names = {
        'BRCA': 'Breast invasive carcinoma',
        'LUAD': 'Lung adenocarcinoma',
        'LUSC': 'Lung squamous cell carcinoma',
        'KIRC': 'Kidney renal clear cell carcinoma',
        'LIHC': 'Liver hepatocellular carcinoma',
        'STAD': 'Stomach adenocarcinoma',
        'COAD': 'Colon adenocarcinoma',
        'HNSC': 'Head and neck squamous cell carcinoma',
        'UCEC': 'Uterine corpus endometrial carcinoma'
    }

    # Individual cancer types
    for idx, ct in enumerate(cancer_types):
        ax = plt.subplot(4, 3, idx + 1)

        ct_data = patient_data[patient_data['cancer_type'] == ct].copy()
        ct_data['phase_early'] = (ct_data['time'] <= 24).astype(int)

        kmf = KaplanMeierFitter()

        # Late phase
        mask_late = ct_data['phase_early'] == 0
        if mask_late.sum() > 0:
            kmf.fit(ct_data.loc[mask_late, 'time'], ct_data.loc[mask_late, 'status'],
                   label=f"Late (>24mo) n={mask_late.sum()}")
            kmf.plot_survival_function(ax=ax, ci_show=True, color=COLORS['late_phase'], linewidth=2)

        # Early phase
        mask_early = ct_data['phase_early'] == 1
        if mask_early.sum() > 0:
            kmf.fit(ct_data.loc[mask_early, 'time'], ct_data.loc[mask_early, 'status'],
                   label=f"Early (≤24mo) n={mask_early.sum()}")
            kmf.plot_survival_function(ax=ax, ci_show=True, color=COLORS['early_phase'], linewidth=2)

        # Add statistics
        ct_result = cancer_results[cancer_results['Cancer_Type'] == ct]
        if len(ct_result) > 0:
            hr = ct_result['HR'].iloc[0]
            p = ct_result['P'].iloc[0]
            stats_text = f"HR={hr:.2f}\nP={p:.4f}"
            ax.text(0.98, 0.02, stats_text, transform=ax.transAxes, fontsize=7,
                   verticalalignment='bottom', horizontalalignment='right',
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

        ax.set_xlabel('Time (months)', fontsize=9, fontweight='bold')
        ax.set_ylabel('Survival probability', fontsize=9, fontweight='bold')
        ax.set_title(cancer_names[ct], fontsize=10, fontweight='bold')
        ax.legend(loc='upper right', fontsize=7)
        ax.grid(alpha=0.25)

    # Pan-cancer analysis
    ax_pan = plt.subplot(4, 3, 10)
    ax_pan.set_position([0.1, 0.05, 0.8, 0.22])

    patient_data['phase_early'] = (patient_data['time'] <= 24).astype(int)
    kmf = KaplanMeierFitter()

    mask_late = patient_data['phase_early'] == 0
    kmf.fit(patient_data.loc[mask_late, 'time'], patient_data.loc[mask_late, 'status'],
           label=f"Late phase n={mask_late.sum()}")
    kmf.plot_survival_function(ax=ax_pan, ci_show=True, color=COLORS['late_phase'], linewidth=3)

    mask_early = patient_data['phase_early'] == 1
    kmf.fit(patient_data.loc[mask_early, 'time'], patient_data.loc[mask_early, 'status'],
           label=f"Early phase n={mask_early.sum()}")
    kmf.plot_survival_function(ax=ax_pan, ci_show=True, color=COLORS['early_phase'], linewidth=3)

    # Cox regression for pan-cancer
    cox_data = patient_data[['time', 'status', 'phase_early']].dropna()
    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(cox_data, duration_col='time', event_col='status')

    hr_pan = np.exp(cph.params_['phase_early'])
    ci_l = np.exp(cph.confidence_intervals_.iloc[0, 0])
    ci_u = np.exp(cph.confidence_intervals_.iloc[0, 1])
    c_pan = cph.concordance_index_

    early = patient_data[patient_data['phase_early'] == 1]
    late = patient_data[patient_data['phase_early'] == 0]
    lr = logrank_test(early['time'], late['time'], early['status'], late['status'])

    stats_text = f"Log-rank P={lr.p_value:.2e}\nHR={hr_pan:.2f} [{ci_l:.2f}-{ci_u:.2f}]\nC-index={c_pan:.3f}"
    ax_pan.text(0.98, 0.02, stats_text, transform=ax_pan.transAxes, fontsize=9,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.9),
               verticalalignment='bottom', horizontalalignment='right', fontweight='bold')

    ax_pan.set_xlabel('Time (months)', fontsize=11, fontweight='bold')
    ax_pan.set_ylabel('Survival probability', fontsize=11, fontweight='bold')
    ax_pan.set_title('Pan-cancer: 24-month biological transition', fontsize=12, fontweight='bold')
    ax_pan.legend(fontsize=9)
    ax_pan.grid(alpha=0.25)

    plt.tight_layout()

    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)

    plt.savefig(output_path / "Figure_KaplanMeier.png", dpi=dpi, bbox_inches='tight')
    plt.savefig(output_path / "Figure_KaplanMeier.pdf", bbox_inches='tight')
    plt.close()

    print(f"   ✅ Saved Kaplan-Meier figure")


def create_forest_plot(cancer_results, output_dir, dpi):
    print(f"\n📊 Creating forest plot...")

    fig, ax = plt.subplots(figsize=(10, 7))

    cancer_results = cancer_results.sort_values('HR', ascending=True)
    y_pos = np.arange(len(cancer_results))

    for i, (_, row) in enumerate(cancer_results.iterrows()):
        is_sig = row['P'] < 0.05
        color = COLORS['significant'] if is_sig else COLORS['nonsignificant']
        marker = 'D' if is_sig else 'o'

        ax.plot([row['CI_Lower'], row['CI_Upper']], [i, i], color=color, linewidth=2.5)
        ax.scatter(row['HR'], i, s=100, color=color, marker=marker,
                  edgecolor='black', linewidth=2, zorder=3)

    ax.axvline(x=1, color='red', linestyle='--', linewidth=2, alpha=0.7)

    ax.set_yticks(y_pos)
    ax.set_yticklabels([f"{row['Cancer_Type']} (n={row['N']})" 
                        for _, row in cancer_results.iterrows()], fontsize=9)
    ax.set_xlabel('Hazard Ratio (95% CI)', fontsize=11, fontweight='bold')
    ax.set_title('24-Month Phase Effect Across Cancer Types', fontsize=12, fontweight='bold')
    ax.grid(axis='x', alpha=0.25)
    ax.set_xscale('log')

    plt.tight_layout()

    output_path = Path(output_dir)
    plt.savefig(output_path / "Figure_ForestPlot.png", dpi=dpi, bbox_inches='tight')
    plt.savefig(output_path / "Figure_ForestPlot.pdf", bbox_inches='tight')
    plt.close()

    print(f"   ✅ Saved forest plot")


def create_summary_table(patient_data, cancer_results, output_dir):
    print(f"\n📊 Creating summary table...")

    # Pan-cancer
    patient_data['phase_early'] = (patient_data['time'] <= 24).astype(int)
    cox_data = patient_data[['time', 'status', 'phase_early']].dropna()

    cph = CoxPHFitter(penalizer=0.01)
    cph.fit(cox_data, duration_col='time', event_col='status')

    hr_pan = np.exp(cph.params_['phase_early'])
    ci_l_pan = np.exp(cph.confidence_intervals_.iloc[0, 0])
    ci_u_pan = np.exp(cph.confidence_intervals_.iloc[0, 1])
    p_pan = cph.summary.loc['phase_early', 'p']
    c_pan = cph.concordance_index_

    table = cancer_results.copy()
    table['HR_CI'] = table.apply(
        lambda row: f"{row['HR']:.2f} ({row['CI_Lower']:.2f}-{row['CI_Upper']:.2f})", axis=1
    )

    pan_row = pd.DataFrame({
        'Cancer_Type': ['Pan-Cancer'],
        'N': [len(cox_data)],
        'Events': [int(cox_data['status'].sum())],
        'HR_CI': [f"{hr_pan:.2f} ({ci_l_pan:.2f}-{ci_u_pan:.2f})"],
        'P': [f"{p_pan:.2e}"],
        'C_index': [f"{c_pan:.3f}"]
    })

    table['P'] = table['P'].apply(lambda x: f"{x:.4f}" if x >= 0.0001 else f"{x:.2e}")
    table['C_index'] = table['C_index'].apply(lambda x: f"{x:.3f}")

    table_display = table[['Cancer_Type', 'N', 'Events', 'HR_CI', 'P', 'C_index']]
    table_final = pd.concat([table_display, pan_row], ignore_index=True)

    table_final.columns = ['Cancer Type', 'N', 'Events', 'HR (95% CI)', 'P-value', 'C-index']

    output_path = Path(output_dir)
    table_final.to_csv(output_path / "Table_ValidationResults.csv", index=False)

    print(f"   ✅ Saved summary table")
    print(f"\n{table_final.to_string(index=False)}")


def main():
    args = parse_arguments()

    print("="*80)
    print("GENERATING FIGURES")
    print("="*80)

    patient_data, cancer_results = load_results(args.results_dir)

    create_kaplan_meier_figure(patient_data, cancer_results, args.output_dir, args.dpi)
    create_forest_plot(cancer_results, args.output_dir, args.dpi)
    create_summary_table(patient_data, cancer_results, args.output_dir)

    print(f"\n{'='*80}")
    print(f"✅ All figures saved to: {args.output_dir}")
    print("="*80)


if __name__ == "__main__":
    main()
