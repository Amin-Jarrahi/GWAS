"""
Genomic Causal Inference Pipeline
==================================

A modular pipeline combining machine learning prediction with causal inference
methods to estimate the effect of genetic variants (SNPs) on human height.

Modules:
    config            – Central Config dataclass and CLI argument parser
    data_loader       – VCF and phenotype data loading / encoding
    feature_selection – Correlation-based feature ranking
    models            – OLS, RF, KNN and two-stage IV variants
    causal            – Mendelian Randomization, DML, DoWhy
    evaluation        – Metrics computation and cross-validation harness
    visualization     – Plotting utilities
    pipeline          – End-to-end orchestration
"""

__version__ = "1.0.0"
__author__ = "Alex Jarrahi"
__email__ = "ajarrahi@vols.utk.edu"

from genomic_causal.config import Config
from genomic_causal.pipeline import run_pipeline

__all__ = ["Config", "run_pipeline", "__version__"]
