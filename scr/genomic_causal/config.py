"""
Central configuration and CLI argument parsing.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from typing import List


@dataclass
class Config:
    """All pipeline settings in one place.

    Override via CLI flags (see ``parse_args``) or by instantiating directly::

        cfg = Config(vcf_dir="/my/vcfs", rf_n_estimators=200)
    """

    # ── Paths ──────────────────────────────────────────────────────────────
    vcf_dir: str = "./data/vcf"
    height_file: str = "./data/GTEx-Height-Demographics.xlsx"
    output_dir: str = "./results"

    # ── VCF parsing ────────────────────────────────────────────────────────
    vcf_header_row: int = 3385  # 0-indexed row where the #CHROM header starts

    # ── Gene list ──────────────────────────────────────────────────────────
    gene_list: List[str] = field(default_factory=lambda: [
        "ENSG00000143476.17",
        "ENSG00000146197.8",
        "ENSG00000119681.11",
        "ENSG00000149257.13",
        "ENSG00000157766.15",
        "ENSG00000159899.14",
        "ENSG00000169047.5",
        "ENSG00000182752.9",
        "ENSG00000185960.13",
    ])

    # ── Feature-count sweeps ───────────────────────────────────────────────
    feature_counts: List[int] = field(default_factory=lambda: [
        2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100,
        200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500,
        3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 7000, 7500, 8000,
        8500, 9000, 9500, 10000, 12000, 14000, 16000, 18000, 20000,
    ])

    # ── Cross-validation ───────────────────────────────────────────────────
    n_subjects: int = 838
    cv_fold_size: int = 140
    random_seed: int = 2022

    # ── Model hyper-parameters ─────────────────────────────────────────────
    rf_n_estimators: int = 100
    knn_n_neighbors: int = 5
    dml_n_splits: int = 5

    # ── Flags ──────────────────────────────────────────────────────────────
    run_dowhy: bool = True
    run_mendelian_randomization: bool = True
    run_double_ml: bool = True


def parse_args(argv: list[str] | None = None) -> Config:
    """Build a ``Config`` from command-line arguments."""
    p = argparse.ArgumentParser(
        description="Genomic Causal Inference Pipeline — VCF → Height Prediction",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Paths
    p.add_argument("--vcf_dir", default=Config.vcf_dir, help="Directory with .vcf.gz files")
    p.add_argument("--height_file", default=Config.height_file, help="Phenotype spreadsheet")
    p.add_argument("--output_dir", default=Config.output_dir, help="Output directory")
    p.add_argument("--vcf_header_row", type=int, default=Config.vcf_header_row,
                    help="0-indexed row of the VCF header line")

    # Hyperparameters
    p.add_argument("--rf_n_estimators", type=int, default=Config.rf_n_estimators)
    p.add_argument("--knn_n_neighbors", type=int, default=Config.knn_n_neighbors)
    p.add_argument("--dml_n_splits", type=int, default=Config.dml_n_splits)
    p.add_argument("--seed", type=int, default=Config.random_seed, dest="random_seed")
    p.add_argument("--cv_fold_size", type=int, default=Config.cv_fold_size)

    # Flags
    p.add_argument("--no_dowhy", action="store_true", help="Skip DoWhy analysis")
    p.add_argument("--no_mr", action="store_true", help="Skip Mendelian Randomization")
    p.add_argument("--no_dml", action="store_true", help="Skip Double ML")

    args = p.parse_args(argv)

    return Config(
        vcf_dir=args.vcf_dir,
        height_file=args.height_file,
        output_dir=args.output_dir,
        vcf_header_row=args.vcf_header_row,
        rf_n_estimators=args.rf_n_estimators,
        knn_n_neighbors=args.knn_n_neighbors,
        dml_n_splits=args.dml_n_splits,
        random_seed=args.random_seed,
        cv_fold_size=args.cv_fold_size,
        run_dowhy=not args.no_dowhy,
        run_mendelian_randomization=not args.no_mr,
        run_double_ml=not args.no_dml,
    )
