"""
Cross-validation harness and results aggregation.
"""

from __future__ import annotations

import logging
import random
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from genomic_causal.config import Config
from genomic_causal.feature_selection import select_top_k
from genomic_causal.models import (
    FoldResult,
    fit_2sknn,
    fit_2sls,
    fit_2srf,
    fit_knn,
    fit_ols,
    fit_rf,
)

logger = logging.getLogger(__name__)


def generate_cv_splits(
    n: int,
    fold_size: int,
    seed: int,
) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Generate non-overlapping train/test splits.

    Each fold uses ``fold_size`` samples as the test set and the remainder
    as training.  Folds are drawn without replacement until the pool is
    exhausted.

    Parameters
    ----------
    n : int
        Total number of samples.
    fold_size : int
        Number of test samples per fold.
    seed : int
        Random seed.

    Returns
    -------
    list of (train_indices, test_indices) tuples.
    """
    rng = random.Random(seed)
    indices = list(range(n))
    rng.shuffle(indices)

    splits = []
    for start in range(0, n - fold_size + 1, fold_size):
        test_idx = np.array(indices[start : start + fold_size])
        train_idx = np.array(indices[:start] + indices[start + fold_size :])
        splits.append((train_idx, test_idx))

    logger.info("Generated %d CV folds (fold_size=%d, n=%d)", len(splits), fold_size, n)
    return splits


def run_cv_sweep(
    X: np.ndarray,
    Y: np.ndarray,
    sorted_indices: np.ndarray,
    splits: List[Tuple[np.ndarray, np.ndarray]],
    feature_counts: List[int],
    cfg: Config,
) -> Dict[str, pd.DataFrame]:
    """Run all six models across feature-count sweep with cross-validation.

    Returns
    -------
    dict mapping model name → DataFrame with columns:
        [n_features, MAE, MSE, RMSE, R2]
    """
    model_names = ["OLS", "2SLS", "RF", "2SRF", "KNN", "2SKNN"]
    all_results: Dict[str, List[Dict]] = {name: [] for name in model_names}

    max_features = X.shape[1]

    for k in feature_counts:
        if k > max_features:
            logger.debug("Skipping k=%d (only %d features available)", k, max_features)
            continue

        X_k = select_top_k(X, sorted_indices, k)
        logger.info("Feature count k=%d  —  running CV", k)

        # Accumulators per model
        fold_metrics: Dict[str, List[FoldResult]] = {name: [] for name in model_names}

        for train_idx, test_idx in splits:
            X_tr, X_te = X_k[train_idx], X_k[test_idx]
            y_tr, y_te = Y[train_idx], Y[test_idx]

            # Standard models
            fold_metrics["OLS"].append(fit_ols(X_tr, y_tr, X_te, y_te))
            fold_metrics["RF"].append(fit_rf(X_tr, y_tr, X_te, y_te, cfg))
            fold_metrics["KNN"].append(fit_knn(X_tr, y_tr, X_te, y_te, cfg))

            # Two-stage models (need ≥2 features)
            if k >= 2:
                fold_metrics["2SLS"].append(fit_2sls(X_tr, y_tr, X_te, y_te))
                fold_metrics["2SRF"].append(fit_2srf(X_tr, y_tr, X_te, y_te, cfg))
                fold_metrics["2SKNN"].append(fit_2sknn(X_tr, y_tr, X_te, y_te, cfg))

        # Average across folds
        for name in model_names:
            folds = fold_metrics[name]
            if not folds:
                continue
            avg = {
                "n_features": k,
                "MAE": np.mean([f.mae for f in folds]),
                "MSE": np.mean([f.mse for f in folds]),
                "RMSE": np.mean([f.rmse for f in folds]),
                "R2": np.mean([f.r2 for f in folds]),
            }
            all_results[name].append(avg)

    return {name: pd.DataFrame(rows) for name, rows in all_results.items() if rows}


def save_results(
    results: Dict[str, pd.DataFrame],
    output_dir: str,
) -> None:
    """Save each model's results to an Excel file."""
    import os
    os.makedirs(output_dir, exist_ok=True)

    for name, df in results.items():
        path = os.path.join(output_dir, f"{name}_Results.xlsx")
        df.to_excel(path, index=False)
        logger.info("Saved %s → %s", name, path)
