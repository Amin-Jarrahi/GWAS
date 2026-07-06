"""Tests for prediction models (OLS, RF, KNN, two-stage variants)."""

import numpy as np
import pytest

from genomic_causal.models import (
    FoldResult,
    fit_ols,
    fit_rf,
    fit_knn,
    fit_2sls,
    fit_2srf,
    fit_2sknn,
)


class TestStandardModels:
    """Basic smoke tests for standard prediction models."""

    def test_ols_returns_fold_result(self, synthetic_data):
        d = synthetic_data
        result = fit_ols(d["X"][:80], d["Y"][:80], d["X"][80:], d["Y"][80:])
        assert isinstance(result, FoldResult)
        assert result.y_pred.shape == (20,)
        assert np.isfinite(result.mae)
        assert np.isfinite(result.r2)

    def test_rf_returns_fold_result(self, synthetic_data, config):
        d = synthetic_data
        result = fit_rf(d["X"][:80], d["Y"][:80], d["X"][80:], d["Y"][80:], config)
        assert isinstance(result, FoldResult)
        assert result.mae >= 0

    def test_knn_returns_fold_result(self, synthetic_data, config):
        d = synthetic_data
        result = fit_knn(d["X"][:80], d["Y"][:80], d["X"][80:], d["Y"][80:], config)
        assert isinstance(result, FoldResult)
        assert result.mae >= 0


class TestTwoStageModels:
    """Smoke tests for two-stage IV models."""

    def test_2sls_runs(self, synthetic_data):
        d = synthetic_data
        result = fit_2sls(d["X"][:80], d["Y"][:80], d["X"][80:], d["Y"][80:])
        assert isinstance(result, FoldResult)
        assert np.isfinite(result.mse)

    def test_2srf_runs(self, synthetic_data, config):
        d = synthetic_data
        result = fit_2srf(d["X"][:80], d["Y"][:80], d["X"][80:], d["Y"][80:], config)
        assert isinstance(result, FoldResult)

    def test_2sknn_runs(self, synthetic_data, config):
        d = synthetic_data
        result = fit_2sknn(d["X"][:80], d["Y"][:80], d["X"][80:], d["Y"][80:], config)
        assert isinstance(result, FoldResult)

    def test_2sls_requires_two_features(self, synthetic_data):
        d = synthetic_data
        X_single = d["X"][:, :1]  # Only 1 feature
        with pytest.raises(ValueError, match="at least 2 features"):
            fit_2sls(X_single[:80], d["Y"][:80], X_single[80:], d["Y"][80:])


class TestMetricsSanity:
    """Ensure metric calculations are correct."""

    def test_mae_is_nonnegative(self, synthetic_data):
        d = synthetic_data
        result = fit_ols(d["X"][:80], d["Y"][:80], d["X"][80:], d["Y"][80:])
        assert result.mae >= 0
        assert result.mse >= 0
        assert result.rmse >= 0

    def test_rmse_equals_sqrt_mse(self, synthetic_data):
        d = synthetic_data
        result = fit_ols(d["X"][:80], d["Y"][:80], d["X"][80:], d["Y"][80:])
        assert abs(result.rmse - np.sqrt(result.mse)) < 1e-10
