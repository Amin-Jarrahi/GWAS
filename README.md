# Genomic Causal Inference Pipeline

[![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A modular pipeline that combines **machine learning prediction** with **causal inference** to estimate the effect of genetic variants (SNPs) on human height using GTEx genotype data.

---

## Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Installation](#installation)
- [Data Preparation](#data-preparation)
- [Usage](#usage)
- [Methods](#methods)
- [Configuration](#configuration)
- [Output](#output)
- [Citation](#citation)
- [License](#license)

---

## Overview

This pipeline takes **VCF genotype files** and a **phenotype spreadsheet** (height, sex, age) as input and runs:

| Category | Methods |
|---|---|
| **Prediction** | OLS, Random Forest (RF), K-Nearest Neighbors (KNN) |
| **Instrumental Variable** | Two-Stage Least Squares (2SLS), Two-Stage RF (2SRF), Two-Stage KNN (2SKNN) |
| **Causal Inference** | Mendelian Randomization (Wald, IVW, MR-Egger, Weighted Median), Double Machine Learning (DML), DoWhy Causal Model |

All methods are evaluated with cross-validation and results are exported as Excel spreadsheets and plots.

---

## Repository Structure

```
genomic-causal-inference/
├── README.md                  # This file
├── LICENSE                    # MIT License
├── pyproject.toml             # Build system & project metadata
├── requirements.txt           # Pinned dependencies
├── setup.cfg                  # Package configuration
├── Makefile                   # Convenience commands
├── .gitignore                 # Git ignore rules
├── CITATION.cff               # Citation metadata
│
├── src/                       # Source package
│   └── genomic_causal/
│       ├── __init__.py        # Package init & version
│       ├── config.py          # Central Config dataclass + CLI parser
│       ├── data_loader.py     # VCF & phenotype loading, encoding
│       ├── feature_selection.py  # Correlation-based feature ranking
│       ├── models.py          # OLS, RF, KNN & two-stage variants
│       ├── causal.py          # MR, DML, DoWhy wrappers
│       ├── evaluation.py      # Metrics, cross-validation harness
│       ├── visualization.py   # Plotting utilities
│       └── pipeline.py        # End-to-end orchestration
│
├── scripts/
│   └── run_pipeline.py        # CLI entry point
│
├── tests/
│   ├── __init__.py
│   ├── conftest.py            # Shared pytest fixtures
│   ├── test_data_loader.py    # Data loading tests
│   ├── test_models.py         # Model tests
│   └── test_causal.py         # Causal inference tests
│
├── notebooks/
│   └── example_analysis.ipynb # Interactive walkthrough
│
├── docs/
│   └── METHODS.md             # Detailed method descriptions
│
├── data/                      # ⬇ User places input data here
│   └── .gitkeep
│
└── results/                   # ⬇ Pipeline writes output here
    └── .gitkeep
```

---

## Installation

### Option A — pip install (editable)

```bash
git clone https://github.com/Amin-Jarrahi/genomic-causal-inference.git
cd genomic-causal-inference
pip install -e ".[dev]"
```

### Option B — requirements only

```bash
git clone https://github.com/Amin-Jarrahi/genomic-causal-inference.git
cd genomic-causal-inference
pip install -r requirements.txt
```

### Option C — Make

```bash
make install      # editable install
make install-dev  # editable install + dev/test dependencies
```

---

## Data Preparation

1. **VCF files** — Place gzipped VCF files (one per gene) in `data/vcf/`:
   ```
   data/vcf/ENSG00000143476.17.vcf.gz
   data/vcf/ENSG00000146197.8.vcf.gz
   ...
   ```

2. **Phenotype file** — Place the GTEx height demographics spreadsheet in `data/`:
   ```
   data/GTEx-Height-Demographics.xlsx
   ```
   Expected columns: `SUBJID`, `HGHT`, `SEX`, `AGE`

---

## Usage

### Command Line

```bash
# Run with defaults
python scripts/run_pipeline.py

# Specify paths explicitly
python scripts/run_pipeline.py \
    --vcf_dir ./data/vcf \
    --height_file ./data/GTEx-Height-Demographics.xlsx \
    --output_dir ./results

# Skip optional causal methods
python scripts/run_pipeline.py --no_dowhy --no_mr

# Change model hyperparameters
python scripts/run_pipeline.py --rf_n_estimators 200 --knn_n_neighbors 7 --seed 42
```

### Python API

```python
from genomic_causal import Config, run_pipeline

cfg = Config(
    vcf_dir="./data/vcf",
    height_file="./data/GTEx-Height-Demographics.xlsx",
    output_dir="./results",
    rf_n_estimators=200,
)
run_pipeline(cfg)
```

### Makefile shortcuts

```bash
make run            # Run pipeline with defaults
make test           # Run test suite
make lint           # Lint with ruff
make clean          # Remove results & caches
```

---

## Methods

See [`docs/METHODS.md`](docs/METHODS.md) for full mathematical descriptions. Summary:

- **OLS / RF / KNN** — Standard supervised regression from SNP features → height.
- **2SLS / 2SRF / 2SKNN** — Two-stage instrumental variable methods. Stage 1 predicts the treatment SNP from instruments; Stage 2 regresses height on the predicted treatment.
- **Mendelian Randomization** — Uses genetic variants as instruments to estimate causal SNP → height effects. Four estimators: Wald Ratio, Inverse-Variance Weighted (IVW), MR-Egger, Weighted Median.
- **Double Machine Learning** — Chernozhukov et al. (2018) framework with cross-fitting and heteroskedasticity-robust standard errors.
- **DoWhy** — Microsoft's causal inference library for graphical causal model specification, identification, estimation, and refutation.

---

## Configuration

All parameters are controlled via the `Config` dataclass. Override via CLI flags or by passing values directly in Python.

| Parameter | Default | Description |
|---|---|---|
| `vcf_dir` | `./data/vcf` | Directory containing `.vcf.gz` files |
| `height_file` | `./data/GTEx-Height-Demographics.xlsx` | Phenotype spreadsheet |
| `output_dir` | `./results` | Where to write results |
| `vcf_header_row` | `3385` | 0-indexed row of the VCF header line |
| `gene_list` | 9 ENSG IDs | Which genes to load |
| `rf_n_estimators` | `100` | Random Forest trees |
| `knn_n_neighbors` | `5` | KNN neighbors |
| `dml_n_splits` | `5` | DML cross-fitting folds |
| `random_seed` | `2022` | Global random seed |
| `run_dowhy` | `True` | Enable DoWhy analysis |
| `run_mendelian_randomization` | `True` | Enable MR |
| `run_double_ml` | `True` | Enable DML |

---

## Output

The pipeline writes to `output_dir/`:

| File | Contents |
|---|---|
| `OLS_Results.xlsx` | OLS metrics per feature count |
| `2SLS_Results.xlsx` | Two-Stage Least Squares metrics |
| `RF_Results.xlsx` | Random Forest metrics |
| `2SRF_Results.xlsx` | Two-Stage Random Forest metrics |
| `KNN_Results.xlsx` | KNN metrics |
| `2SKNN_Results.xlsx` | Two-Stage KNN metrics |
| `MR_Results.xlsx` | Mendelian Randomization estimates |
| `DML_Results.xlsx` | Double ML estimates |
| `DoWhy_Results.xlsx` | DoWhy estimates & refutations |
| `*_plot.png` | Metrics vs. feature-count plots |

---

## Citation

If you use this pipeline in your research, please cite:

```bibtex
@software{jarrahi2026genomic_causal,
  author    = {Jarrahi, Alex(Amin)},
  title     = {Genomic Causal Inference Pipeline},
  year      = {2026},
  url       = {https://github.com/Amin-Jarrahi/genomic-causal-inference}
}
```

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
