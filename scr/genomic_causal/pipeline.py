"""
End-to-end pipeline orchestration.

Ties together data loading, feature selection, model training,
causal inference, and result export.
"""

from __future__ import annotations

import logging
import os
import random
import warnings
from typing import Optional

import numpy as np
import pandas as pd

from genomic_causal.config import Config, parse_args
from genomic_causal.data_loader import build_dataset
from genomic_causal.feature_selection import rank_features_by_correlation
from genomic_causal.evaluation import generate_cv_splits, run_cv_sweep, save_results
from genomic_causal.visualization import plot_all_metrics, plot_mr_forest
from genomic_causal.causal import (
    run_mendelian_randomization,
    run_double_ml,
    run_dowhy_analysis,
)

logger = logging.getLogger(__name__)


def _set_seed(seed: int) -> None:
    """Set global random seeds for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)


def run_pipeline(cfg: Config) -> None:
    """Execute the full genomic causal inference pipeline.

    Steps:
        1. Load VCF genotype data and phenotype data.
        2. Rank features by correlation with height.
        3. Run cross-validated prediction models (OLS, RF, KNN + two-stage).
        4. Run causal inference (MR, DML, DoWhy).
        5. Save results and plots.

    Parameters
    ----------
    cfg : Config
        Pipeline configuration.
    """
    _set_seed(cfg.random_seed)
    os.makedirs(cfg.output_dir, exist_ok=True)

    # ── Step 1: Load data ──────────────────────────────────────────────────
    logger.info("=" * 70)
    logger.info("STEP 1 — Loading data")
    logger.info("=" * 70)

    X, Y, height_df, genotype_df, all_variant_ids = build_dataset(cfg)

    # Confounders: Sex & Age
    confounders = height_df[["Sex", "Age"]].values.astype(np.float64)

    # ── Step 2: Feature ranking ────────────────────────────────────────────
    logger.info("=" * 70)
    logger.info("STEP 2 — Feature ranking (correlation with height)")
    logger.info("=" * 70)

    column_names = list(genotype_df.columns) if genotype_df is not None else None
    sorted_indices, sorted_names = rank_features_by_correlation(X, Y, column_names)

    # Reorder X so column 0 = most correlated (treatment for causal methods)
    X_ranked = X[:, sorted_indices]

    # ── Step 3: Cross-validated prediction models ──────────────────────────
    logger.info("=" * 70)
    logger.info("STEP 3 — Cross-validated prediction models")
    logger.info("=" * 70)

    splits = generate_cv_splits(len(Y), cfg.cv_fold_size, cfg.random_seed)
    prediction_results = run_cv_sweep(
        X_ranked, Y, np.arange(X_ranked.shape[1]),  # already sorted
        splits, cfg.feature_counts, cfg,
    )

    save_results(prediction_results, cfg.output_dir)
    plot_all_metrics(prediction_results, cfg.output_dir)

    # ── Step 4: Causal Inference ───────────────────────────────────────────
    logger.info("=" * 70)
    logger.info("STEP 4 — Causal inference")
    logger.info("=" * 70)

    # Use top-100 features for causal methods (or all if fewer)
    n_causal_features = min(100, X_ranked.shape[1])
    X_causal = X_ranked[:, :n_causal_features]
    causal_feature_names = sorted_names[:n_causal_features]

    # 4a) Mendelian Randomization
    mr_results = None
    if cfg.run_mendelian_randomization:
        try:
            mr_results = run_mendelian_randomization(
                X_causal, Y, confounders, causal_feature_names,
                n_instruments=min(20, n_causal_features - 1),
                seed=cfg.random_seed,
            )
            # Save MR results
            if mr_results:
                mr_df = pd.DataFrame([
                    {
                        "Method": r.method,
                        "Estimate": r.estimate,
                        "SE": r.se,
                        "P-value": r.pvalue,
                        "CI_Lower": r.ci_lower,
                        "CI_Upper": r.ci_upper,
                        "N_Instruments": r.n_instruments,
                    }
                    for r in mr_results
                ])
                mr_path = os.path.join(cfg.output_dir, "MR_Results.xlsx")
                mr_df.to_excel(mr_path, index=False)
                logger.info("MR results saved → %s", mr_path)

                plot_mr_forest(mr_results, cfg.output_dir)
        except Exception as e:
            logger.error("Mendelian Randomization failed: %s", e)

    # 4b) Double Machine Learning
    dml_result = None
    if cfg.run_double_ml:
        try:
            dml_result = run_double_ml(X_causal, Y, confounders, cfg)
            dml_df = pd.DataFrame([{
                "Estimate": dml_result.estimate,
                "SE": dml_result.se,
                "P-value": dml_result.pvalue,
                "CI_Lower": dml_result.ci_lower,
                "CI_Upper": dml_result.ci_upper,
                "N_Splits": dml_result.n_splits,
            }])
            dml_path = os.path.join(cfg.output_dir, "DML_Results.xlsx")
            dml_df.to_excel(dml_path, index=False)
            logger.info("DML results saved → %s", dml_path)
        except Exception as e:
            logger.error("Double ML failed: %s", e)

    # 4c) DoWhy
    dowhy_result = None
    if cfg.run_dowhy:
        try:
            dowhy_result = run_dowhy_analysis(
                X_causal, Y, confounders, causal_feature_names,
                n_instruments=min(20, n_causal_features - 1),
            )
            if dowhy_result is not None:
                dowhy_df = pd.DataFrame([{
                    "Estimate": dowhy_result.estimate,
                    "Refutation_Placebo": dowhy_result.refutation_placebo,
                    "Refutation_Random_Common_Cause": dowhy_result.refutation_random_common_cause,
                    "Refutation_Subset": dowhy_result.refutation_subset,
                }])
                dowhy_path = os.path.join(cfg.output_dir, "DoWhy_Results.xlsx")
                dowhy_df.to_excel(dowhy_path, index=False)
                logger.info("DoWhy results saved → %s", dowhy_path)
        except Exception as e:
            logger.error("DoWhy analysis failed: %s", e)

    # ── Summary ────────────────────────────────────────────────────────────
    logger.info("=" * 70)
    logger.info("PIPELINE COMPLETE — results written to: %s", cfg.output_dir)
    logger.info("=" * 70)


def main() -> None:
    """CLI entry point."""
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s  %(levelname)-8s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    cfg = parse_args()
    run_pipeline(cfg)


if __name__ == "__main__":
    main()
