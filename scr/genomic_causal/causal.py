"""
Causal inference methods:
  - Mendelian Randomization (Wald Ratio, IVW, MR-Egger, Weighted Median)
  - Double Machine Learning (DML) with cross-fitting
  - DoWhy causal model with robustness refutations

Confounders used: Sex, Age (from phenotype data).
Treatment:        First (most-correlated) SNP.
Instruments:      Remaining top SNPs.
Outcome:          Height.
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold

from genomic_causal.config import Config

logger = logging.getLogger(__name__)

# DoWhy — optional
try:
    from dowhy import CausalModel
    DOWHY_AVAILABLE = True
except ImportError:
    DOWHY_AVAILABLE = False


# ═══════════════════════════════════════════════════════════════════════════════
#  Mendelian Randomization
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class MRResult:
    """Container for Mendelian Randomization estimates."""
    method: str
    estimate: float
    se: float
    pvalue: float
    ci_lower: float
    ci_upper: float
    n_instruments: int


def _wald_ratio(beta_zy: float, se_zy: float,
                beta_zx: float, se_zx: float) -> Tuple[float, float]:
    """Single-instrument Wald ratio: β_causal = β_ZY / β_ZX."""
    ratio = beta_zy / beta_zx
    # Delta-method SE
    se = abs(ratio) * np.sqrt((se_zy / beta_zy) ** 2 + (se_zx / beta_zx) ** 2)
    return ratio, se


def _ivw_estimate(
    beta_zy: np.ndarray, se_zy: np.ndarray,
    beta_zx: np.ndarray, se_zx: np.ndarray,
) -> Tuple[float, float]:
    """Inverse-Variance Weighted (IVW) MR — fixed-effect meta-analysis."""
    weights = beta_zx ** 2 / se_zy ** 2
    beta_ivw = np.sum(weights * (beta_zy / beta_zx)) / np.sum(weights)
    se_ivw = 1.0 / np.sqrt(np.sum(weights))
    return beta_ivw, se_ivw


def _mr_egger(
    beta_zy: np.ndarray, se_zy: np.ndarray,
    beta_zx: np.ndarray, se_zx: np.ndarray,
) -> Tuple[float, float, float, float]:
    """MR-Egger regression — allows for directional pleiotropy.

    Returns: (slope, slope_se, intercept, intercept_se)
    """
    from sklearn.linear_model import LinearRegression

    weights = 1.0 / se_zy ** 2
    W_sqrt = np.sqrt(weights)

    # Weighted regression of β_ZY on β_ZX (with intercept)
    X_w = (beta_zx * W_sqrt).reshape(-1, 1)
    y_w = beta_zy * W_sqrt

    reg = LinearRegression(fit_intercept=True)
    reg.fit(X_w, y_w)

    slope = reg.coef_[0]
    intercept = reg.intercept_

    # Residual-based SE
    y_pred = reg.predict(X_w)
    n = len(beta_zy)
    residuals = y_w - y_pred
    s2 = np.sum(residuals ** 2) / max(n - 2, 1)

    X_w_mat = np.column_stack([np.ones(n), X_w.ravel()])
    cov_matrix = s2 * np.linalg.pinv(X_w_mat.T @ X_w_mat)
    intercept_se = np.sqrt(max(cov_matrix[0, 0], 0))
    slope_se = np.sqrt(max(cov_matrix[1, 1], 0))

    return slope, slope_se, intercept, intercept_se


def _weighted_median(
    beta_zy: np.ndarray, se_zy: np.ndarray,
    beta_zx: np.ndarray, se_zx: np.ndarray,
    n_boot: int = 1000, seed: int = 42,
) -> Tuple[float, float]:
    """Weighted median MR estimator — robust if ≤50% instruments are invalid."""
    ratios = beta_zy / beta_zx
    weights = beta_zx ** 2 / se_zy ** 2
    weights = weights / weights.sum()

    # Weighted median
    order = np.argsort(ratios)
    cumw = np.cumsum(weights[order])
    idx = np.searchsorted(cumw, 0.5)
    idx = min(idx, len(ratios) - 1)
    beta_wm = ratios[order[idx]]

    # Bootstrap SE
    rng = np.random.RandomState(seed)
    boot_estimates = []
    for _ in range(n_boot):
        b_zy = beta_zy + rng.normal(0, se_zy)
        b_zx = beta_zx + rng.normal(0, se_zx)
        # Avoid division by zero
        mask = np.abs(b_zx) > 1e-12
        if mask.sum() == 0:
            continue
        r = b_zy[mask] / b_zx[mask]
        w = b_zx[mask] ** 2 / se_zy[mask] ** 2
        w = w / w.sum()
        o = np.argsort(r)
        cw = np.cumsum(w[o])
        i = min(np.searchsorted(cw, 0.5), len(r) - 1)
        boot_estimates.append(r[o[i]])

    se_wm = np.std(boot_estimates) if boot_estimates else np.nan
    return beta_wm, se_wm


def run_mendelian_randomization(
    X: np.ndarray,
    Y: np.ndarray,
    confounders: np.ndarray,
    feature_names: List[str],
    n_instruments: int = 10,
    seed: int = 42,
) -> List[MRResult]:
    """Run all four MR estimators.

    Parameters
    ----------
    X : array, shape (n, p)
        Genotype matrix (columns ranked by correlation).
    Y : array, shape (n,)
        Outcome (height).
    confounders : array, shape (n, q)
        Confounder matrix (Sex, Age).
    feature_names : list[str]
        Column names for X.
    n_instruments : int
        Number of top SNPs to use as instruments (after the treatment SNP).
    seed : int
        Random seed for bootstrap.

    Returns
    -------
    list[MRResult]
    """
    logger.info("Running Mendelian Randomization with %d instruments", n_instruments)

    treatment = X[:, 0]  # Most-correlated SNP
    n_inst = min(n_instruments, X.shape[1] - 1)
    instruments = X[:, 1 : 1 + n_inst]

    # Residualise treatment and outcome on confounders
    from sklearn.linear_model import LinearRegression

    if confounders.shape[1] > 0:
        reg_t = LinearRegression().fit(confounders, treatment)
        treatment_resid = treatment - reg_t.predict(confounders)
        reg_y = LinearRegression().fit(confounders, Y)
        Y_resid = Y - reg_y.predict(confounders)
    else:
        treatment_resid = treatment
        Y_resid = Y

    # Per-instrument regression coefficients
    n = len(Y)
    beta_zx = np.zeros(n_inst)
    se_zx = np.zeros(n_inst)
    beta_zy = np.zeros(n_inst)
    se_zy = np.zeros(n_inst)

    for j in range(n_inst):
        z = instruments[:, j]
        if confounders.shape[1] > 0:
            reg_z = LinearRegression().fit(confounders, z)
            z_resid = z - reg_z.predict(confounders)
        else:
            z_resid = z

        # Z → X
        slope_zx, intercept_zx, r_zx, p_zx, se_slope_zx = stats.linregress(z_resid, treatment_resid)
        beta_zx[j] = slope_zx
        se_zx[j] = se_slope_zx

        # Z → Y
        slope_zy, intercept_zy, r_zy, p_zy, se_slope_zy = stats.linregress(z_resid, Y_resid)
        beta_zy[j] = slope_zy
        se_zy[j] = se_slope_zy

    results: List[MRResult] = []

    # 1) Wald Ratio (strongest instrument)
    best_j = np.argmax(np.abs(beta_zx))
    if abs(beta_zx[best_j]) > 1e-12:
        est, se = _wald_ratio(beta_zy[best_j], se_zy[best_j],
                              beta_zx[best_j], se_zx[best_j])
        z_stat = est / se if se > 0 else 0.0
        pval = 2 * (1 - stats.norm.cdf(abs(z_stat)))
        results.append(MRResult("Wald Ratio", est, se, pval,
                                est - 1.96 * se, est + 1.96 * se, 1))

    # 2) IVW
    valid = np.abs(beta_zx) > 1e-12
    if valid.sum() >= 2:
        est, se = _ivw_estimate(beta_zy[valid], se_zy[valid],
                                beta_zx[valid], se_zx[valid])
        z_stat = est / se if se > 0 else 0.0
        pval = 2 * (1 - stats.norm.cdf(abs(z_stat)))
        results.append(MRResult("IVW", est, se, pval,
                                est - 1.96 * se, est + 1.96 * se, int(valid.sum())))

    # 3) MR-Egger
    if valid.sum() >= 3:
        slope, slope_se, intercept, intercept_se = _mr_egger(
            beta_zy[valid], se_zy[valid], beta_zx[valid], se_zx[valid],
        )
        z_stat = slope / slope_se if slope_se > 0 else 0.0
        pval = 2 * (1 - stats.norm.cdf(abs(z_stat)))
        results.append(MRResult("MR-Egger", slope, slope_se, pval,
                                slope - 1.96 * slope_se, slope + 1.96 * slope_se,
                                int(valid.sum())))

    # 4) Weighted Median
    if valid.sum() >= 3:
        est, se = _weighted_median(beta_zy[valid], se_zy[valid],
                                   beta_zx[valid], se_zx[valid],
                                   seed=seed)
        z_stat = est / se if se > 0 else 0.0
        pval = 2 * (1 - stats.norm.cdf(abs(z_stat)))
        results.append(MRResult("Weighted Median", est, se, pval,
                                est - 1.96 * se, est + 1.96 * se, int(valid.sum())))

    for r in results:
        logger.info(
            "  MR %-16s  β=%.4f  SE=%.4f  p=%.4g  95%%CI [%.4f, %.4f]",
            r.method, r.estimate, r.se, r.pvalue, r.ci_lower, r.ci_upper,
        )
    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  Double Machine Learning
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class DMLResult:
    """Container for Double Machine Learning estimate."""
    estimate: float
    se: float
    pvalue: float
    ci_lower: float
    ci_upper: float
    n_splits: int


def run_double_ml(
    X: np.ndarray,
    Y: np.ndarray,
    confounders: np.ndarray,
    cfg: Config,
) -> DMLResult:
    """Partially-linear Double Machine Learning (Chernozhukov et al., 2018).

    Model:
        Y = θ·T + g(W) + ε
        T = m(W) + η

    where T = X[:, 0] (treatment SNP), W = [X[:, 1:], confounders].

    Uses K-fold cross-fitting with Random Forest nuisance models.

    Parameters
    ----------
    X : array, shape (n, p)
        Genotype features (columns ranked by correlation).
    Y : array, shape (n,)
        Outcome (height).
    confounders : array, shape (n, q)
        Confounder matrix (Sex, Age).
    cfg : Config
        Pipeline configuration (for n_splits, rf params, seed).

    Returns
    -------
    DMLResult
    """
    logger.info("Running Double Machine Learning (K=%d)", cfg.dml_n_splits)

    treatment = X[:, 0]
    # Covariates W = remaining SNPs + confounders
    W = np.hstack([X[:, 1:], confounders]) if confounders.shape[1] > 0 else X[:, 1:]

    n = len(Y)
    residuals_y = np.zeros(n)
    residuals_t = np.zeros(n)

    kf = KFold(n_splits=cfg.dml_n_splits, shuffle=True, random_state=cfg.random_seed)

    for fold_idx, (train_idx, test_idx) in enumerate(kf.split(W)):
        W_train, W_test = W[train_idx], W[test_idx]
        Y_train, Y_test = Y[train_idx], Y[test_idx]
        T_train, T_test = treatment[train_idx], treatment[test_idx]

        # Nuisance model for Y ~ W
        model_y = RandomForestRegressor(
            n_estimators=cfg.rf_n_estimators, random_state=cfg.random_seed,
        )
        model_y.fit(W_train, Y_train)
        residuals_y[test_idx] = Y_test - model_y.predict(W_test)

        # Nuisance model for T ~ W
        model_t = RandomForestRegressor(
            n_estimators=cfg.rf_n_estimators, random_state=cfg.random_seed,
        )
        model_t.fit(W_train, T_train)
        residuals_t[test_idx] = T_test - model_t.predict(W_test)

    # Final stage: regress residualised Y on residualised T
    denom = np.sum(residuals_t ** 2)
    if denom < 1e-12:
        logger.warning("DML: treatment residuals near zero — estimate unreliable.")
        return DMLResult(0.0, np.nan, np.nan, np.nan, np.nan, cfg.dml_n_splits)

    theta = np.sum(residuals_t * residuals_y) / denom

    # Heteroskedasticity-robust (HC0) standard error
    eps = residuals_y - theta * residuals_t
    var_theta = np.sum((residuals_t ** 2) * (eps ** 2)) / (denom ** 2)
    se_theta = np.sqrt(max(var_theta, 0))

    z_stat = theta / se_theta if se_theta > 0 else 0.0
    pval = 2 * (1 - stats.norm.cdf(abs(z_stat)))

    result = DMLResult(
        estimate=theta,
        se=se_theta,
        pvalue=pval,
        ci_lower=theta - 1.96 * se_theta,
        ci_upper=theta + 1.96 * se_theta,
        n_splits=cfg.dml_n_splits,
    )

    logger.info(
        "  DML  θ=%.4f  SE=%.4f  p=%.4g  95%%CI [%.4f, %.4f]",
        result.estimate, result.se, result.pvalue, result.ci_lower, result.ci_upper,
    )
    return result


# ═══════════════════════════════════════════════════════════════════════════════
#  DoWhy Causal Model
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class DoWhyResult:
    """Container for DoWhy causal analysis."""
    estimate: float
    refutation_placebo: Optional[float] = None
    refutation_random_common_cause: Optional[float] = None
    refutation_subset: Optional[float] = None


def run_dowhy_analysis(
    X: np.ndarray,
    Y: np.ndarray,
    confounders: np.ndarray,
    feature_names: List[str],
    n_instruments: int = 10,
) -> Optional[DoWhyResult]:
    """Run DoWhy causal model: estimate + three refutations.

    Parameters
    ----------
    X : array, shape (n, p)
        Genotype matrix (ranked by correlation).
    Y : array, shape (n,)
        Height.
    confounders : array, shape (n, q)
        Sex, Age.
    feature_names : list[str]
        Column names for X.
    n_instruments : int
        How many top SNPs to use as instruments.

    Returns
    -------
    DoWhyResult or None if DoWhy is unavailable.
    """
    if not DOWHY_AVAILABLE:
        logger.warning("DoWhy not installed — skipping.")
        return None

    logger.info("Running DoWhy causal model")

    n_inst = min(n_instruments, X.shape[1] - 1)
    treatment_name = feature_names[0] if feature_names else "treatment"
    instrument_names = feature_names[1 : 1 + n_inst] if len(feature_names) > 1 else []
    confounder_names = ["Sex", "Age"]

    # Build DataFrame
    data = pd.DataFrame(X[:, : 1 + n_inst], columns=[treatment_name] + instrument_names)
    for i, cn in enumerate(confounder_names):
        if i < confounders.shape[1]:
            data[cn] = confounders[:, i]
    data["Height"] = Y

    available_confounders = [cn for i, cn in enumerate(confounder_names)
                            if i < confounders.shape[1]]

    try:
        model = CausalModel(
            data=data,
            treatment=treatment_name,
            outcome="Height",
            common_causes=available_confounders if available_confounders else None,
            instruments=instrument_names if instrument_names else None,
        )

        identified = model.identify_effect(proceed_when_unidentifiable=True)

        estimate = model.estimate_effect(
            identified,
            method_name="iv.instrumental_variable" if instrument_names
                        else "backdoor.linear_regression",
        )
        est_value = float(estimate.value)
        logger.info("  DoWhy estimate: %.4f", est_value)

        result = DoWhyResult(estimate=est_value)

        # Refutations
        try:
            ref1 = model.refute_estimate(
                identified, estimate, method_name="placebo_treatment_refuter",
            )
            result.refutation_placebo = float(ref1.new_effect)
            logger.info("  Placebo refutation: %.4f", result.refutation_placebo)
        except Exception as e:
            logger.warning("  Placebo refutation failed: %s", e)

        try:
            ref2 = model.refute_estimate(
                identified, estimate, method_name="random_common_cause",
            )
            result.refutation_random_common_cause = float(ref2.new_effect)
            logger.info("  Random common cause: %.4f", result.refutation_random_common_cause)
        except Exception as e:
            logger.warning("  Random common cause refutation failed: %s", e)

        try:
            ref3 = model.refute_estimate(
                identified, estimate, method_name="data_subset_refuter",
            )
            result.refutation_subset = float(ref3.new_effect)
            logger.info("  Subset refutation: %.4f", result.refutation_subset)
        except Exception as e:
            logger.warning("  Subset refutation failed: %s", e)

        return result

    except Exception as e:
        logger.error("DoWhy analysis failed: %s", e)
        return None
