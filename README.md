# COSMOS

**Co**mputational **S**emi-supervised **MO**dular **S**election for temporal recurrence biology in solid tumors.

This repository contains the reproducibility pipeline for the COSMOS framework, which identifies a data-driven 21 to 25 month transition window separating early proliferation-driven recurrence from late immune- and microenvironment-driven recurrence across multiple solid tumor types.

## Overview

COSMOS combines unsupervised modular feature reduction with a Maximally Selected Rank Statistics (MSRS) scan to discover temporal transition cut-points in recurrence biology, then validates the resulting phase definition against independent overall-survival endpoints across nine TCGA cancer types and against whole-slide histopathology in 9,816 diagnostic images.

Key findings reproduced by this pipeline:

- A 21 to 25 month convergence window across four cancer types (breast, lung adenocarcinoma, glioblastoma, renal clear cell)
- Proliferation module Cohen's d = 0.76 (95% CI: +0.48 to +1.06) in pooled discovery cohorts
- Pan-cancer Cox HR = 4.25 (95% CI: 3.61 to 5.01) for phase-stratified overall survival, n = 2,871
- 24-month landmark HR = 2.33 (95% CI: 1.82 to 2.98), n = 2,440
- Whole-slide image validation HR = 1.80 (95% CI: 1.55 to 2.08), n = 4,881 patients, 9,816 slides
- 36 of 37 pre-specified pathways significant at FDR < 0.05 after Benjamini-Hochberg correction
- Spearman correlation rho = +0.60 between PAM50 subtype Ki-67 rank and median recurrence time

## Repository contents

```
.
├── 01_preprocessing.ipynb              GEO + TCGA loading, ComBat, ssGSEA
├── 02_cosmos_feature_selection.ipynb   Hopkins, K-means, IQR + IsoForest, MSLR
├── 03_msrs_phase_analysis.ipynb        MSRS scan, phase comparison, JSD, BH, bootstrap
├── 04_pancancer_cox.ipynb              Per-cancer + pan-cancer Cox, landmark
├── 05_wsi_validation.ipynb             ResNet50 features, patient-level Cox
├── 06_subtype_analysis.ipynb           PAM50 subtype-stratified Spearman
├── README.md
└── requirements.txt
```

Each notebook is self-contained and writes intermediate pickled outputs that downstream notebooks consume.

## Installation

Python 3.10 or higher is required. Create a fresh environment:

```bash
python -m venv cosmos-env
source cosmos-env/bin/activate
pip install -r requirements.txt
```

Core dependencies:

| Package | Tested version |
|---------|---------------|
| numpy | 1.26 |
| pandas | 2.1 |
| scipy | 1.11 |
| scikit-learn | 1.3 |
| statsmodels | 0.14 |
| lifelines | 0.27 |
| gseapy | 1.1 |
| inmoose (ComBat) | 0.7 |
| torch | 2.1 |
| torchvision | 0.16 |
| openslide-python | 1.3 |
| tqdm, joblib | latest |

A CUDA-enabled GPU is recommended for notebook 05 (whole-slide feature extraction). Notebooks 01 through 04 and 06 run on CPU.

## Data layout

Set the `COSMOS_BASE_DIR` environment variable to the root of your data directory:

```bash
export COSMOS_BASE_DIR=/path/to/cosmos_data
```

Expected directory structure:

```
$COSMOS_BASE_DIR/
├── Raw Data/
│   ├── geo/
│   │   ├── GSE2034/
│   │   │   ├── GSE2034_series_matrix.txt
│   │   │   ├── GSE2034_family.soft.gz
│   │   │   └── GSE2034_clinical.csv         # manual curation of DMFS
│   │   ├── GSE2990/
│   │   │   ├── GSE2990_series_matrix.txt
│   │   │   ├── GSE2990_family.soft.gz
│   │   │   └── GSE2990_suppl_info.txt
│   │   ├── GSE103746/
│   │   └── GSE31210/
│   ├── tcga_rnaseq/
│   │   ├── TCGA-BRCA/
│   │   │   ├── counts/                      # STAR-aligned gene counts
│   │   │   └── patient_file_map.csv
│   │   ├── TCGA-LUAD/
│   │   ├── TCGA-KIRC/
│   │   ├── TCGA-GBM/
│   │   └── gencode.v36.annotation.gtf.gz
│   ├── tcga_clinical/
│   │   ├── TCGA-CDR.csv                     # https://gdc.cancer.gov
│   │   └── brca_tcga_pan_can_atlas_2018_clinical_data.tsv  # cBioPortal
│   └── pathway_databases/
│       ├── h.all.gmt                        # MSigDB v2023.2 Hallmark
│       ├── c2.cp.kegg_legacy.gmt
│       ├── c2.cp.reactome.gmt
│       └── c5.go.bp.gmt
├── intermediates/                            # created by notebooks
└── wsi_results/                              # created by notebook 05
```

Whole-slide images for notebook 05 are located under a separate directory specified by `COSMOS_WSI_DIR`. Files should be named with TCGA patient barcodes (e.g., `TCGA-XX-YYYY-...svs`). Both flat and per-cancer-type subdirectory layouts are supported.

## Running the pipeline

Notebooks are designed to run sequentially. Each notebook writes its outputs to `$COSMOS_BASE_DIR/intermediates/` for downstream notebooks to consume.

```bash
jupyter nbconvert --to notebook --execute 01_preprocessing.ipynb
jupyter nbconvert --to notebook --execute 02_cosmos_feature_selection.ipynb
jupyter nbconvert --to notebook --execute 03_msrs_phase_analysis.ipynb
jupyter nbconvert --to notebook --execute 04_pancancer_cox.ipynb
jupyter nbconvert --to notebook --execute 05_wsi_validation.ipynb
jupyter nbconvert --to notebook --execute 06_subtype_analysis.ipynb
```

Approximate runtimes on a workstation with 384 GB RAM and an RTX 4090:

| Notebook | CPU only | With GPU |
|----------|----------|----------|
| 01 preprocessing | 35 min | n/a |
| 02 COSMOS | 12 min | n/a |
| 03 MSRS + BH + bootstrap | 18 min | n/a |
| 04 pan-cancer Cox | 2 min | n/a |
| 05 WSI validation | impractical | 9 to 14 hr |
| 06 subtype analysis | < 1 min | n/a |

## Reproducibility

All randomized components use a fixed seed of 42:

- `numpy`, `random`, and `PYTHONHASHSEED` are seeded at the top of each notebook
- `sklearn` estimators receive `random_state=42`
- ComBat and ssGSEA permutation tests use `np.random.RandomState(seed)` for reproducible streams
- PyTorch (notebook 05) seeds both CPU and CUDA streams

The pipeline produces deterministic outputs given identical input data and software versions. Small numerical differences (within reported confidence intervals) may occur across `scikit-learn` or `gseapy` versions.

## Key technical choices

- **Endpoint priority:** DFI as primary recurrence endpoint, PFI as fallback. GBM uses PFI as primary because DFI is not meaningful for an infiltrative primary brain tumor.
- **Phase definition:** Binary indicator of recurrence event observed within 24 months. Censored observations and late events both belong to the "late" group.
- **Cox covariates:** All pan-cancer models are stratified by cancer type. The whole-slide validation model includes `phase_early + PC1` with cancer-type stratification.
- **Pathway panel:** 37 pathways across four biological modules (Proliferation, Immune, Metabolic, Microenvironment), pre-specified before any survival analysis to prevent post-hoc selection.
- **WSI feature extractor:** ResNet50 pretrained on ImageNet. Pathology foundation models (UNI, CONCH, Virchow) are not used here because they were pretrained on TCGA, which would create circular dependence with the validation cohort.

## Outputs

Each notebook writes one or more pickle files to `$COSMOS_BASE_DIR/intermediates/`:

| Notebook | Outputs |
|----------|---------|
| 01 | `geo_clinical.pkl`, `tcga_clinical.pkl`, `geo_expr_corrected.pkl`, `tcga_expr_corrected.pkl`, `ssgsea_corrected.pkl`, `ssgsea_raw.pkl`, `combined_genesets.pkl`, `tcga_endpoint_used.pkl` |
| 02 | `stage1_results.pkl`, `stage2_results.pkl`, `stage3_results.pkl`, `core_pathways_37.pkl`, `module_scores_37.pkl`, `brca_specific_modules.pkl`, `luad_specific_modules.pkl` |
| 03 | `transition_results.pkl` (cut-points, phase comparisons, JSD, BH-corrected pathway table, bootstrap CIs) |
| 04 | `cox_summary.pkl` (per-cancer, pan-cancer, leave-one-out, landmark) |
| 05 | `patient_level_data.csv`, `cancer_specific_results.csv`, `summary.json`, `wsi_pca_model.pkl`, `wsi_scaler.pkl` |
| 06 | `subtype_results.pkl` (per-subtype metrics, Spearman rho and p-value) |

## Manuscript

If you use this code, please cite:

> Zafar SA, Khan AAK, Nazir A, Khalid F, Faisal MS, Qin W. *COmputational Semi-supervised MOdular Selection (COSMOS): A pan-cancer framework for discovering temporal transition windows in tumor recurrence biology.* Manuscript under review.

## Funding

This work was supported by the Shanghai Jiao Tong University Yang Generation Award (SJTU YG2025QNA31).

## License

MIT License. See `LICENSE` for details. Use of TCGA, GEO, and cBioPortal data is subject to their respective data use agreements.

## Contact

For questions about the code, open an issue on GitHub. For research collaboration, contact the corresponding author through the affiliated institution.
