"""
Plotting utilities for pipeline results.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def plot_metric_vs_features(
    results: Dict[str, pd.DataFrame],
    metric: str,
    output_dir: str,
    title: Optional[str] = None,
    filename: Optional[str] = None,
) -> None:
    """Plot a single metric vs. number of features for all models.

    Parameters
    ----------
    results : dict
        Model name → DataFrame with columns [n_features, MAE, MSE, RMSE, R2].
    metric : str
        Which metric column to plot (e.g. "MAE", "R2").
    output_dir : str
        Directory to save the plot.
    title : str, optional
        Plot title. Defaults to ``metric + " vs. Number of Features"``.
    filename : str, optional
        Output filename. Defaults to ``metric + "_plot.png"``.
    """
    fig, ax = plt.subplots(figsize=(12, 6))

    for model_name, df in results.items():
        if metric not in df.columns:
            continue
        ax.plot(df["n_features"], df[metric], marker="o", markersize=3,
                label=model_name, linewidth=1.5)

    ax.set_xlabel("Number of Features (SNPs)", fontsize=12)
    ax.set_ylabel(metric, fontsize=12)
    ax.set_title(title or f"{metric} vs. Number of Features", fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)

    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, filename or f"{metric}_plot.png")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Plot saved → %s", path)


def plot_all_metrics(
    results: Dict[str, pd.DataFrame],
    output_dir: str,
) -> None:
    """Generate plots for all standard metrics."""
    for metric in ["MAE", "MSE", "RMSE", "R2"]:
        plot_metric_vs_features(results, metric, output_dir)


def plot_mr_forest(
    mr_results: list,
    output_dir: str,
    filename: str = "MR_forest_plot.png",
) -> None:
    """Forest plot for Mendelian Randomization estimates.

    Parameters
    ----------
    mr_results : list of MRResult
        Mendelian Randomization results.
    output_dir : str
        Output directory.
    filename : str
        Output file name.
    """
    if not mr_results:
        return

    methods = [r.method for r in mr_results]
    estimates = [r.estimate for r in mr_results]
    ci_lower = [r.ci_lower for r in mr_results]
    ci_upper = [r.ci_upper for r in mr_results]

    fig, ax = plt.subplots(figsize=(8, max(3, len(methods) * 0.8)))

    y_pos = range(len(methods))
    xerr_lower = [e - lo for e, lo in zip(estimates, ci_lower)]
    xerr_upper = [hi - e for e, hi in zip(estimates, ci_upper)]

    ax.errorbar(
        estimates, y_pos,
        xerr=[xerr_lower, xerr_upper],
        fmt="o", color="steelblue", capsize=5, markersize=8,
    )
    ax.axvline(0, color="grey", linestyle="--", linewidth=0.8)
    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(methods, fontsize=11)
    ax.set_xlabel("Causal Estimate (β)", fontsize=12)
    ax.set_title("Mendelian Randomization — Forest Plot", fontsize=13)
    ax.grid(True, axis="x", alpha=0.3)

    os.makedirs(output_dir, exist_ok=True)
    path = os.path.join(output_dir, filename)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("MR forest plot saved → %s", path)
