# 24-Month Biological Transition: Multi-Cancer WSI Validation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

**Official validation code for:** "Temporal Stratification Reveals Discrete Proliferation-to-Immune Transition at 24 Months in Cancer Recurrence"

## 📋 Overview

This repository provides complete, reproducible code for validating the **24-month biological transition** in cancer recurrence using whole slide imaging (WSI) data from 2,871 TCGA patients across 9 cancer types.

### What This Code Does

✅ **Cox Proportional Hazards Regression** - Proper survival analysis with censored data  
✅ **Patient-Level Aggregation** - Ensures statistical independence  
✅ **Multi-Cancer Validation** - Analyzes 9 TCGA cancer types (1,362 events)  
✅ **Publication Figures** - Generates Figure 8 (Kaplan-Meier curves) and forest plots  
✅ **Full Reproducibility** - Fixed random seeds, documented methods  

### Key Results

| Metric | Value |
|--------|-------|
| **Patients** | 2,871 |
| **Events** | 1,362 |
| **Hazard Ratio** | 56.48 (95% CI: 39.57-80.62) |
| **P-value** | <0.001 |
| **C-index** | 0.76 |
| **Power** | >99% |

---

## 🚀 Quick Start

### Installation
```bash
