"""
Correlation-based feature ranking for SNP selection.
"""

from __future__ import annotations

import logging
from typing import List, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def rank_features_by_correlation(
    X: np.ndarray,
    Y: np.ndarray,
    column_names: List[str] | None = None,
) -> Tuple[np.ndarray, List[str]]:
    """Rank features by absolute Pearson correlation with the target.

    Parameters
    ----------
    X : np.ndarray, shape (n_samples, n_features)
        Feature matrix.
    Y : np.ndarray, shape (n_samples,)
        Target vector.
    column_names : list[str], optional
        Feature names (for logging). If None, indices are used.

    Returns
    -------
    sorted_indices : np.ndarray
        Feature indices sorted by descending |correlation|.
    sorted_names : list[str]
        Corresponding feature names.
    """
    n_features = X.shape[1]
    correlations = np.zeros(n_features)

    for j in range(n_features):
        col = X[:, j]
        # Skip constant columns (zero variance)
        if np.std(col) < 1e-12:
            correlations[j] = 0.0
        else:
            correlations[j] = abs(np.corrcoef(col, Y)[0, 1])

    # Replace NaN correlations with 0
    correlations = np.nan_to_num(correlations, nan=0.0)

    sorted_indices = np.argsort(correlations)[::-1]

    if column_names is None:
        column_names = [str(i) for i in range(n_features)]

    sorted_names = [column_names[i] for i in sorted_indices]

    logger.info(
        "Feature ranking complete — top-5 correlations: %s",
        [(sorted_names[i], f"{correlations[sorted_indices[i]]:.4f}") for i in range(min(5, len(sorted_names)))],
    )
    return sorted_indices, sorted_names


def select_top_k(
    X: np.ndarray,
    sorted_indices: np.ndarray,
    k: int,
) -> np.ndarray:
    """Select the top-k features from X based on pre-computed ranking.

    Parameters
    ----------
    X : np.ndarray
        Full feature matrix.
    sorted_indices : np.ndarray
        Feature indices sorted by importance (descending).
    k : int
        Number of features to keep. Clamped to the available count.

    Returns
    -------
    np.ndarray
        Reduced feature matrix with k columns.
    """
    k = min(k, X.shape[1])
    return X[:, sorted_indices[:k]]
