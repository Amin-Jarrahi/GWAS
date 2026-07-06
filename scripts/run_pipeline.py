#!/usr/bin/env python3
"""
CLI entry point for the Genomic Causal Inference Pipeline.

Usage:
    python scripts/run_pipeline.py --vcf_dir ./data/vcf --height_file ./data/heights.xlsx
    python scripts/run_pipeline.py --help
"""

import sys
from pathlib import Path

# Allow running from the repo root without installing
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from genomic_causal.pipeline import main

if __name__ == "__main__":
    main()
