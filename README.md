# COSMOS: 24-Month Biological Transition Validation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

**Official validation code for:** "Temporal Stratification Reveals Discrete Proliferation-to-Immune Transition at 24 Months in Cancer Recurrence"

---

## 📋 Overview

This repository provides complete, reproducible code for validating the **24-month biological transition** in cancer recurrence using whole slide imaging (WSI) data from 2,871 TCGA patients across 9 cancer types.

### What This Code Does

✅ **Cox Proportional Hazards Regression** - Proper survival analysis with censored data  
✅ **Patient-Level Aggregation** - Ensures statistical independence  
✅ **Multi-Cancer Validation** - Analyzes 9 TCGA cancer types (1,362 events)  
✅ **Publication Figures** - Generates Kaplan-Meier curves and forest plots  
✅ **Full Reproducibility** - Fixed random seeds, documented methods  

### Key Results

| Metric | Value |
|--------|-------|
| **Patients** | 2,871 |
| **Events** | 1,362 |
| **Hazard Ratio** | 56.48 (95% CI: 39.57-80.62) |
| **P-value** | <0.001 |
| **C-index** | 0.76 |

---

## 🚀 Quick Start

### Installation

```bash
# Clone repository
git clone https://github.com/Sjtu-Fuxilab/COSMOS.git
cd COSMOS

# Install dependencies
pip install -r requirements.txt
```

### Run Validation

```bash
python full_validation.py \
    --data_dir /path/to/tcga \
    --cdr_file /path/to/TCGA-CDR.xlsx \
    --output_dir ./results
```

---

## 📦 Repository Contents

```
COSMOS/
├── README.md              # This file
├── LICENSE                # MIT License
├── requirements.txt       # Dependencies
├── full_validation.py     # Main validation script
└── docs/
    ├── METHODS.md        # Detailed methodology
    └── INSTALLATION.md   # Installation guide
```

---

## 🔬 Methodology

### Cox Proportional Hazards Regression

This analysis uses **gold-standard survival analysis methods**:

1. **Patient-Level Aggregation**: Multiple slides averaged per patient
2. **Cox Regression**: Time-to-event analysis with censoring
3. **24-Month Phase**: Binary predictor (early ≤24 vs late >24 months)
4. **Metrics**: Hazard ratios, 95% CI, C-index, log-rank tests

---

## 📊 Data Access

### TCGA Data (Public)

**Whole Slide Images:**
- Source: [NCI Genomic Data Commons](https://portal.gdc.cancer.gov)
- Cancer types: BRCA, COAD, HNSC, KIRC, LIHC, LUAD, LUSC, STAD, UCEC

**Clinical Data:**
- Source: TCGA-CDR (Liu et al., Cell 2018)
- DOI: [10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052)

---

## 💻 System Requirements

**Minimum:**
- Python 3.9+
- 16GB RAM
- NVIDIA GPU with 6GB VRAM
- 500GB storage

**Recommended:**
- 32GB RAM
- NVIDIA RTX 3090 (24GB)
- 1TB SSD

---

## 📚 Citation

```bibtex
@article{cosmos2025,
  title={Temporal Stratification Reveals Discrete Proliferation-to-Immune 
         Transition at 24 Months in Cancer Recurrence},
  author={Fuxilab, Shanghai Jiao Tong University},
  year={2025}
}
```

---

## 📝 License

MIT License - see [LICENSE](LICENSE) file.

---

## 📧 Contact

**Fuxilab**  
Shanghai Jiao Tong University

**Issues:** [GitHub Issues](https://github.com/Sjtu-Fuxilab/COSMOS/issues)
