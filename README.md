# COSMOS: 24-Month Biological Transition Validation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

**Official validation code for:** "Temporal Stratification Reveals Discrete Proliferation-to-Immune Transition at 24 Months in Cancer Recurrence"

---

## 📋 Overview

Complete, reproducible code for validating the **24-month biological transition** in cancer recurrence using WSI data from 2,871 TCGA patients across 9 cancer types.

### Features

✅ **Cox Proportional Hazards Regression** - Proper survival analysis  
✅ **Patient-Level Aggregation** - Statistical independence  
✅ **Multi-Cancer Validation** - 9 TCGA cancer types  
✅ **Publication Figures** - Kaplan-Meier curves & forest plots  
✅ **Full Reproducibility** - Fixed seeds, documented methods  

### Key Results

| Metric | Value |
|--------|-------|
| **Patients** | 2,871 |
| **Events** | 1,362 |
| **HR** | 56.48 (95% CI: 39.57-80.62) |
| **P-value** | <0.001 |
| **C-index** | 0.76 |

---

## 🚀 Quick Start

```bash
# Clone & install
git clone https://github.com/Sjtu-Fuxilab/COSMOS.git
cd COSMOS
pip install -r requirements.txt

# Run validation
python full_validation.py \
    --data_dir /path/to/tcga \
    --cdr_file /path/to/TCGA-CDR.xlsx \
    --output_dir ./results

# Generate figures
python visualization.py \
    --results_dir ./results \
    --output_dir ./figures
```

---

## 📦 Repository Contents

```
COSMOS/
├── README.md                # This file
├── LICENSE                  # MIT License
├── requirements.txt         # Dependencies
├── full_validation.py       # Main validation
├── visualization.py         # Figure generation
└── docs/
    ├── METHODS.md          # Methodology
    └── INSTALLATION.md     # Installation
```

---

## 🔬 Methodology

Uses **gold-standard Cox proportional hazards regression**:

1. **Patient-Level Aggregation** - Multiple slides averaged
2. **Cox Regression** - Time-to-event with censoring
3. **24-Month Phase** - Binary predictor (≤24 vs >24 months)
4. **Metrics** - HR, 95% CI, C-index, log-rank tests

---

## 📊 Data

**TCGA Data (Public):**
- **WSI:** [GDC Portal](https://portal.gdc.cancer.gov)
- **Clinical:** TCGA-CDR (Liu et al., Cell 2018)
- **Cancer types:** BRCA, COAD, HNSC, KIRC, LIHC, LUAD, LUSC, STAD, UCEC

---

## 💻 Requirements

**Minimum:** Python 3.9+, 16GB RAM, NVIDIA GPU (6GB)

**Recommended:** 32GB RAM, RTX 3090 (24GB)

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

MIT License - see [LICENSE](LICENSE)

---

## 📧 Contact

**Fuxilab**  
Shanghai Jiao Tong University

**Issues:** [GitHub Issues](https://github.com/Sjtu-Fuxilab/COSMOS/issues)
