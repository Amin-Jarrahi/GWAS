"""
Shared pytest fixtures for the test suite.

These fixtures provide synthetic data so tests run without real VCF files.
"""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def synthetic_data():
    """Generate a small synthetic dataset mimicking the pipeline's data shapes.

    Returns a dict with keys: X, Y, confounders, feature_names
    """
    rng = np.random.RandomState(42)
    n_samples = 100
    n_features = 50

    # Genotype matrix (0, 1, 2 encoding)
    X = rng.choice([0.0, 1.0, 2.0], size=(n_samples, n_features), p=[0.6, 0.3, 0.1])

    # True causal effect from first SNP + noise
    true_effect = 2.5
    Y = true_effect * X[:, 0] + 0.5 * X[:, 1] + rng.normal(65, 4, n_samples)

    # Confounders: Sex (1 or 2), Age (0–4)
    sex = rng.choice([1.0, 2.0], size=n_samples)
    age = rng.choice([0.0, 1.0, 2.0, 3.0, 4.0], size=n_samples)
    confounders = np.column_stack([sex, age])

    feature_names = [f"SNP_{i}" for i in range(n_features)]

    return {
        "X": X.astype(np.float32),
        "Y": Y.astype(np.float64),
        "confounders": confounders,
        "feature_names": feature_names,
        "true_effect": true_effect,
    }


@pytest.fixture
def config():
    """Minimal Config for testing."""
    from genomic_causal.config import Config
    return Config(
        rf_n_estimators=10,  # Small for speed
        knn_n_neighbors=3,
        dml_n_splits=3,
        random_seed=42,
    )
