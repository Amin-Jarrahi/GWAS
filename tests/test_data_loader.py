"""Tests for data loading and preprocessing utilities."""

import numpy as np
import pytest

from genomic_causal.data_loader import GENOTYPE_MAP, AGE_MAP, SEX_MAP


class TestGenotypeMaps:
    """Verify genotype encoding constants."""

    def test_homozygous_ref(self):
        assert GENOTYPE_MAP["0|0"] == 0.0
        assert GENOTYPE_MAP["0/0"] == 0.0

    def test_heterozygous(self):
        assert GENOTYPE_MAP["0|1"] == 1.0
        assert GENOTYPE_MAP["1|0"] == 1.0
        assert GENOTYPE_MAP["0/1"] == 1.0

    def test_homozygous_alt(self):
        assert GENOTYPE_MAP["1|1"] == 2.0
        assert GENOTYPE_MAP["1/1"] == 2.0

    def test_missing(self):
        assert np.isnan(GENOTYPE_MAP[".|."])
        assert np.isnan(GENOTYPE_MAP["./."])

    def test_age_encoding(self):
        assert AGE_MAP["20-29"] == 0.0
        assert AGE_MAP["50-59"] == 2.0
        assert AGE_MAP["70-79"] == 4.0

    def test_sex_encoding(self):
        # GTEx: 1 = male → 2.0,  2 = female → 1.0
        assert SEX_MAP[1] == 2.0
        assert SEX_MAP[2] == 1.0
