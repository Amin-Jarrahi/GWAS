"""
Prediction models: OLS, Random Forest, KNN, and their two-stage IV variants.

The two-stage methods implement instrumental-variable estimation:
  Stage 1 — Predict the endogenous treatment from instruments (using train data).
  Stage 2 — Regress outcome on the predicted treatment (using train data).
  Predict — Apply both stages to held-out test data.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Dict, Optional, Tuple

import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.neighbors import KNeighborsRegressor

from genomic_causal.config import Config

logger = logging.getLogger(__name__)


# ═══════════════════════════════════════════════════════════════════════════════
#  Result Container
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class FoldResult:
    """Metrics from a single cross-validation fold."""
    y_true: np.ndarray
    y_pred: np.ndarray
    mae: float
    mse: float
    rmse: float
    r2: float


# ═══════════════════════════════════════════════════════════════════════════════
#  Standard Models
# ═══════════════════════════════════════════════════════════════════════════════

def fit_ols(
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray, y_test: np.ndarray,
) -> FoldResult:
    """Ordinary Least Squares regression."""
    model = LinearRegression()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    return _make_result(y_test, y_pred)


def fit_rf(
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray, y_test: np.ndarray,
    cfg: Config,
) -> FoldResult:
    """Random Forest regression."""
    model = RandomForestRegressor(
        n_estimators=cfg.rf_n_estimators,
        random_state=cfg.random_seed,
    )
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    return _make_result(y_test, y_pred)


def fit_knn(
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray, y_test: np.ndarray,
    cfg: Config,
) -> FoldResult:
    """K-Nearest Neighbors regression."""
    model = KNeighborsRegressor(n_neighbors=cfg.knn_n_neighbors)
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    return _make_result(y_test, y_pred)


# ═══════════════════════════════════════════════════════════════════════════════
#  Two-Stage Instrumental Variable Models
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Causal structure assumed:
#     Instruments (Z) → Treatment (T) → Outcome (Y)
#
#  We split X into:
#     T = X[:, 0]    (first / most-correlated SNP = treatment)
#     Z = X[:, 1:]   (remaining SNPs = instruments)
#
#  Stage 1:  T̂ = f(Z)        — fit on TRAIN only
#  Stage 2:  Y = g(T̂)        — fit on TRAIN only
#  Test:     T̂_test = f(Z_test),  Ŷ_test = g(T̂_test)
# ═══════════════════════════════════════════════════════════════════════════════


def _split_treatment_instruments(
    X: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Split feature matrix into treatment (col 0) and instruments (cols 1+)."""
    if X.shape[1] < 2:
        raise ValueError(
            "Two-stage models require at least 2 features "
            "(1 treatment + 1 instrument)."
        )
    return X[:, 0], X[:, 1:]


def fit_2sls(
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray, y_test: np.ndarray,
) -> FoldResult:
    """Two-Stage Least Squares (2SLS).

    Stage 1: OLS  Z_train → T_train
    Stage 2: OLS  T̂_train → Y_train
    Test:    T̂_test = Stage1(Z_test),  Ŷ_test = Stage2(T̂_test)
    """
    T_train, Z_train = _split_treatment_instruments(X_train)
    T_test, Z_test = _split_treatment_instruments(X_test)

    # Stage 1 — predict treatment from instruments
    stage1 = LinearRegression()
    stage1.fit(Z_train, T_train)
    T_hat_train = stage1.predict(Z_train)
    T_hat_test = stage1.predict(Z_test)

    # Stage 2 — regress outcome on predicted treatment
    stage2 = LinearRegression()
    stage2.fit(T_hat_train.reshape(-1, 1), y_train)
    y_pred = stage2.predict(T_hat_test.reshape(-1, 1))

    return _make_result(y_test, y_pred)


def fit_2srf(
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray, y_test: np.ndarray,
    cfg: Config,
) -> FoldResult:
    """Two-Stage Random Forest (2SRF).

    Stage 1: RF   Z_train → T_train
    Stage 2: RF   T̂_train → Y_train
    Test:    T̂_test = Stage1(Z_test),  Ŷ_test = Stage2(T̂_test)
    """
    T_train, Z_train = _split_treatment_instruments(X_train)
    T_test, Z_test = _split_treatment_instruments(X_test)

    # Stage 1
    stage1 = RandomForestRegressor(
        n_estimators=cfg.rf_n_estimators,
        random_state=cfg.random_seed,
    )
    stage1.fit(Z_train, T_train)
    T_hat_train = stage1.predict(Z_train)
    T_hat_test = stage1.predict(Z_test)

    # Stage 2
    stage2 = RandomForestRegressor(
        n_estimators=cfg.rf_n_estimators,
        random_state=cfg.random_seed,
    )
    stage2.fit(T_hat_train.reshape(-1, 1), y_train)
    y_pred = stage2.predict(T_hat_test.reshape(-1, 1))

    return _make_result(y_test, y_pred)


def fit_2sknn(
    X_train: np.ndarray, y_train: np.ndarray,
    X_test: np.ndarray, y_test: np.ndarray,
    cfg: Config,
) -> FoldResult:
    """Two-Stage K-Nearest Neighbors (2SKNN).

    Stage 1: KNN  Z_train → T_train
    Stage 2: KNN  T̂_train → Y_train
    Test:    T̂_test = Stage1(Z_test),  Ŷ_test = Stage2(T̂_test)
    """
    T_train, Z_train = _split_treatment_instruments(X_train)
    T_test, Z_test = _split_treatment_instruments(X_test)

    # Stage 1
    stage1 = KNeighborsRegressor(n_neighbors=cfg.knn_n_neighbors)
    stage1.fit(Z_train, T_train)
    T_hat_train = stage1.predict(Z_train)
    T_hat_test = stage1.predict(Z_test)

    # Stage 2
    stage2 = KNeighborsRegressor(n_neighbors=cfg.knn_n_neighbors)
    stage2.fit(T_hat_train.reshape(-1, 1), y_train)
    y_pred = stage2.predict(T_hat_test.reshape(-1, 1))

    return _make_result(y_test, y_pred)


# ═══════════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _make_result(y_true: np.ndarray, y_pred: np.ndarray) -> FoldResult:
    """Compute regression metrics and wrap in a FoldResult."""
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

    mae = mean_absolute_error(y_true, y_pred)
    mse = mean_squared_error(y_true, y_pred)
    rmse = np.sqrt(mse)
    r2 = r2_score(y_true, y_pred)

    return FoldResult(
        y_true=y_true,
        y_pred=y_pred,
        mae=mae,
        mse=mse,
        rmse=rmse,
        r2=r2,
    )
