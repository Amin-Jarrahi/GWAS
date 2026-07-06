"""
Data loading and preprocessing for VCF genotype files and phenotype data.
"""

from __future__ import annotations

import logging
import os
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from genomic_causal.config import Config

logger = logging.getLogger(__name__)

# ── Constants ──────────────────────────────────────────────────────────────────

VCF_META_COLS = [
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
]

GENOTYPE_MAP = {
    "0|0": 0.0, "0|1": 1.0, "1|0": 1.0, "1|1": 2.0, ".|.": np.nan,
    "0/0": 0.0, "0/1": 1.0, "1/0": 1.0, "1/1": 2.0, "./.": np.nan,
}

AGE_MAP = {
    "20-29": 0.0, "30-39": 0.0, "40-49": 1.0,
    "50-59": 2.0, "60-69": 3.0, "70-79": 4.0,
}

SEX_MAP = {"1": 2.0, "2": 1.0, 1: 2.0, 2: 1.0}


# ── VCF Loading ───────────────────────────────────────────────────────────────

def load_vcf(gene_name: str, cfg: Config) -> Tuple[pd.DataFrame, List[str]]:
    """Load a compressed VCF file and encode genotypes as 0/1/2.

    Parameters
    ----------
    gene_name : str
        Ensembl gene ID (used to find the file and name columns).
    cfg : Config
        Pipeline configuration.

    Returns
    -------
    vcf_encoded : pd.DataFrame
        Samples × variants matrix with additive genotype encoding.
    variant_ids : list[str]
        Raw variant IDs from the VCF ID column.
    """
    filepath = os.path.join(cfg.vcf_dir, f"{gene_name}.vcf.gz")
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"VCF file not found: {filepath}")

    vcf = pd.read_csv(
        filepath, compression="gzip", header=cfg.vcf_header_row, sep="\t",
    )
    variant_ids = vcf["ID"].tolist()

    # Keep only sample columns (drop VCF metadata columns)
    sample_cols = [c for c in vcf.columns if c not in VCF_META_COLS]
    vcf_samples = vcf[sample_cols].replace(GENOTYPE_MAP)

    # Transpose so rows = samples, columns = variants
    vcf_encoded = vcf_samples.T
    vcf_encoded.columns = [f"{gene_name}_{vid}" for vid in variant_ids]

    logger.info(
        "Loaded %s — %d variants × %d samples",
        gene_name, len(variant_ids), len(sample_cols),
    )
    return vcf_encoded, variant_ids


# ── Phenotype Loading ─────────────────────────────────────────────────────────

def load_height_data(cfg: Config) -> pd.DataFrame:
    """Load and encode the height / demographics spreadsheet.

    Expected columns in the Excel file: ``SUBJID``, ``HGHT``, ``SEX``, ``AGE``.

    Returns
    -------
    pd.DataFrame
        Columns: ``SUBJID``, ``Height``, ``Sex``, ``Age``.
    """
    if not os.path.isfile(cfg.height_file):
        raise FileNotFoundError(f"Height file not found: {cfg.height_file}")

    df = pd.read_excel(cfg.height_file)

    # Encode categorical variables
    df["AGE"] = df["AGE"].map(AGE_MAP).fillna(df["AGE"]).astype(float)
    df["SEX"] = df["SEX"].map(SEX_MAP).fillna(df["SEX"]).astype(float)
    df = df.rename(columns={"HGHT": "Height", "SEX": "Sex", "AGE": "Age"})

    logger.info("Height data loaded — %d subjects", len(df))
    return df


# ── Full Dataset Builder ──────────────────────────────────────────────────────

def build_dataset(
    cfg: Config,
) -> Tuple[np.ndarray, np.ndarray, pd.DataFrame, pd.DataFrame, List[str]]:
    """Combine all gene VCFs with phenotype data into a single dataset.

    Returns
    -------
    X : np.ndarray, shape (n_samples, n_variants)
        Genotype matrix (NaN-imputed with column medians).
    Y : np.ndarray, shape (n_samples,)
        Height values.
    height_df : pd.DataFrame
        Phenotype DataFrame with ``Height``, ``Sex``, ``Age`` columns.
    genotype_df : pd.DataFrame
        Full samples × variants DataFrame (pre-imputation).
    all_variant_ids : list[str]
        Flat list of all variant ID strings.
    """
    genotype_df: Optional[pd.DataFrame] = None
    height_df: Optional[pd.DataFrame] = None
    all_variant_ids: List[str] = []

    for i, gene in enumerate(cfg.gene_list):
        vcf_encoded, variant_ids = load_vcf(gene, cfg)
        all_variant_ids.extend(variant_ids)

        if genotype_df is None:
            genotype_df = vcf_encoded
        else:
            genotype_df = pd.concat([genotype_df, vcf_encoded], axis=1)

        # Build phenotype vectors by matching sample IDs (first gene only)
        if i == 0:
            height_raw = load_height_data(cfg)
            subj_lookup: Dict[str, Dict[str, float]] = {}
            for _, row in height_raw.iterrows():
                subj_lookup[str(row["SUBJID"])] = {
                    "Height": float(row["Height"]),
                    "Sex": float(row["Sex"]),
                    "Age": float(row["Age"]),
                }

            records = []
            for sid in vcf_encoded.index:
                sid_clean = str(sid).strip("(),' \"")
                if sid_clean in subj_lookup:
                    records.append(subj_lookup[sid_clean])

            height_df = pd.DataFrame(records)

    assert genotype_df is not None and height_df is not None

    # Convert to numpy and impute NaNs with column medians
    X = genotype_df.values.astype(np.float32)
    col_medians = np.nanmedian(X, axis=0)
    nan_mask = np.isnan(X)
    X[nan_mask] = np.take(col_medians, np.where(nan_mask)[1])

    Y = height_df["Height"].values.astype(np.float64)

    logger.info(
        "Dataset built — X shape: %s, Y shape: %s", X.shape, Y.shape,
    )
    return X, Y, height_df, genotype_df, all_variant_ids
