"""Tests for causal inference methods (MR, DML, DoWhy)."""

import numpy as np
import pytest

from genomic_causal.causal import (
    run_mendelian_randomization,
    run_double_ml,
    MRResult,
    DMLResult,
)


class TestMendelianRandomization:
    """Tests for MR estimators."""

    def test_mr_returns_results(self, synthetic_data):
        d = synthetic_data
        results = run_mendelian_randomization(
            d["X"], d["Y"], d["confounders"], d["feature_names"],
            n_instruments=10, seed=42,
        )
        assert isinstance(results, list)
        assert len(results) >= 1  # At least Wald ratio should work
        assert all(isinstance(r, MRResult) for r in results)

    def test_mr_has_all_four_methods(self, synthetic_data):
        d = synthetic_data
        results = run_mendelian_randomization(
            d["X"], d["Y"], d["confounders"], d["feature_names"],
            n_instruments=10, seed=42,
        )
        methods = {r.method for r in results}
        # With 10 instruments, all four methods should run
        assert "Wald Ratio" in methods
        assert "IVW" in methods
        assert "MR-Egger" in methods
        assert "Weighted Median" in methods

    def test_mr_confidence_intervals(self, synthetic_data):
        d = synthetic_data
        results = run_mendelian_randomization(
            d["X"], d["Y"], d["confounders"], d["feature_names"],
            n_instruments=10, seed=42,
        )
        for r in results:
            assert r.ci_lower <= r.estimate <= r.ci_upper

    def test_mr_pvalues_valid(self, synthetic_data):
        d = synthetic_data
        results = run_mendelian_randomization(
            d["X"], d["Y"], d["confounders"], d["feature_names"],
            n_instruments=10, seed=42,
        )
        for r in results:
            assert 0 <= r.pvalue <= 1


class TestDoubleML:
    """Tests for Double Machine Learning."""

    def test_dml_returns_result(self, synthetic_data, config):
        d = synthetic_data
        result = run_double_ml(d["X"], d["Y"], d["confounders"], config)
        assert isinstance(result, DMLResult)
        assert np.isfinite(result.estimate)

    def test_dml_confidence_interval(self, synthetic_data, config):
        d = synthetic_data
        result = run_double_ml(d["X"], d["Y"], d["confounders"], config)
        assert result.ci_lower <= result.estimate <= result.ci_upper

    def test_dml_pvalue_valid(self, synthetic_data, config):
        d = synthetic_data
        result = run_double_ml(d["X"], d["Y"], d["confounders"], config)
        assert 0 <= result.pvalue <= 1

    def test_dml_detects_true_effect(self, synthetic_data, config):
        """With a strong true effect, DML should estimate something non-zero."""
        d = synthetic_data
        result = run_double_ml(d["X"], d["Y"], d["confounders"], config)
        # true_effect = 2.5, estimate should be in the right ballpark
        # (not testing exact value — just that it's meaningfully non-zero)
        assert abs(result.estimate) > 0.1
