"""Tests pour l'algorithme de région partagée basé sur l'allèle mineur."""

import numpy as np
import pandas as pd
import pytest

from lodlink.config import GENO_AA, GENO_AB, GENO_BB, GENO_MISSING
from lodlink.html_viz import _has_minor_allele, analyze_shared_region


class FakePedigree:
    """Pedigree minimal pour les tests."""
    def __init__(self, affected, unaffected):
        self.affected = set(affected)
        self.unaffected = set(unaffected)


def _make_chr_results(marker_names, bp_positions, cm_positions):
    """Crée un chr_results minimal avec un markers_df."""
    mdf = pd.DataFrame({
        'name': marker_names,
        'bp': bp_positions,
        'cm': cm_positions,
    })
    return {'markers_df': mdf}


class TestHasMinorAllele:
    def test_missing_returns_none(self):
        assert _has_minor_allele(GENO_MISSING, True) is None
        assert _has_minor_allele(GENO_MISSING, False) is None

    def test_minor_is_a(self):
        # Minor = A → présent dans AA et AB
        assert _has_minor_allele(GENO_AA, minor_is_a=True) is True
        assert _has_minor_allele(GENO_AB, minor_is_a=True) is True
        assert _has_minor_allele(GENO_BB, minor_is_a=True) is False

    def test_minor_is_b(self):
        # Minor = B → présent dans BB et AB
        assert _has_minor_allele(GENO_BB, minor_is_a=False) is True
        assert _has_minor_allele(GENO_AB, minor_is_a=False) is True
        assert _has_minor_allele(GENO_AA, minor_is_a=False) is False


class TestAnalyzeSharedRegion:
    def test_clear_discrimination(self):
        """Les affectés portent l'allèle mineur, pas les non-affectés."""
        markers = ['m1', 'm2', 'm3', 'm4', 'm5']
        # freq_A < 0.5 → A est mineur
        freq_dict = {m: 0.1 for m in markers}

        # Affectés : tous AA (portent le mineur A)
        # Non-affectés : tous BB (ne portent pas le mineur A)
        genotypes = {}
        affected = [1, 2, 3]
        unaffected = [4, 5, 6]
        for m in markers:
            genotypes[m] = {}
            for a in affected:
                genotypes[m][a] = GENO_AA
            for u in unaffected:
                genotypes[m][u] = GENO_BB

        pedigree = FakePedigree(affected, unaffected)
        region = {
            'marker_names': markers,
            'start_bp': 1000,
            'end_bp': 5000,
        }
        chr_results = _make_chr_results(
            markers,
            [1000, 2000, 3000, 4000, 5000],
            [1.0, 2.0, 3.0, 4.0, 5.0],
        )

        result = analyze_shared_region(region, chr_results, pedigree, genotypes, freq_dict)

        assert result['has_shared'] is True
        assert result['n_affected_sharing'] == 3
        assert result['n_unaffected_sharing'] == 0
        assert set(result['affected_sharing_ids']) == {1, 2, 3}
        assert result['unaffected_sharing_ids'] == []

    def test_no_discrimination(self):
        """Tout le monde a le même génotype → pas de discrimination."""
        markers = ['m1', 'm2', 'm3']
        freq_dict = {m: 0.1 for m in markers}

        # Tout le monde est AB → tout le monde porte le mineur
        genotypes = {}
        affected = [1, 2]
        unaffected = [3, 4, 5]
        for m in markers:
            genotypes[m] = {}
            for ind in affected + unaffected:
                genotypes[m][ind] = GENO_AB

        pedigree = FakePedigree(affected, unaffected)
        region = {'marker_names': markers}
        chr_results = _make_chr_results(markers, [1000, 2000, 3000], [1.0, 2.0, 3.0])

        result = analyze_shared_region(region, chr_results, pedigree, genotypes, freq_dict)

        # Score devrait être négatif ou nul (2/2 - 3/3 = 0 → pas de segment positif robuste)
        # Ici score = (2/2)-(3/3) = 0 par marqueur, sum = 0 → has_shared = False
        assert result['has_shared'] is False

    def test_partial_sharing(self):
        """Certains affectés partagent, certains non-affectés aussi."""
        markers = ['m1', 'm2', 'm3', 'm4']
        # A est mineur (freq < 0.5)
        freq_dict = {m: 0.2 for m in markers}

        affected = [1, 2, 3]
        unaffected = [4, 5]

        # m1-m4: 2/3 affectés ont AA, 1 a BB. 0/2 non-affectés ont A.
        genotypes = {}
        for m in markers:
            genotypes[m] = {
                1: GENO_AA,  # affected, has minor
                2: GENO_AA,  # affected, has minor
                3: GENO_BB,  # affected, no minor
                4: GENO_BB,  # unaffected, no minor
                5: GENO_BB,  # unaffected, no minor
            }

        pedigree = FakePedigree(affected, unaffected)
        region = {'marker_names': markers}
        chr_results = _make_chr_results(markers, [1000, 2000, 3000, 4000], [1.0, 2.0, 3.0, 4.0])

        result = analyze_shared_region(region, chr_results, pedigree, genotypes, freq_dict)

        assert result['has_shared'] is True
        # 2/3 affectés partagent, 0/2 non-affectés
        assert result['n_affected_sharing'] == 2
        assert result['n_unaffected_sharing'] == 0

    def test_empty_markers(self):
        """Région sans marqueurs → None."""
        pedigree = FakePedigree([1, 2], [3])
        result = analyze_shared_region(
            {'marker_names': []}, {}, pedigree, {}, {}
        )
        assert result is None

    def test_single_affected(self):
        """Un seul affecté → None."""
        pedigree = FakePedigree([1], [2, 3])
        result = analyze_shared_region(
            {'marker_names': ['m1']}, {}, pedigree, {}, {}
        )
        assert result is None

    def test_kadane_finds_best_subsegment(self):
        """Vérifie que Kadane trouve le meilleur sous-segment, pas toute la région."""
        # 10 marqueurs:
        # - 3 premiers: non-affectés ont mineur, pas affectés (score négatif)
        # - 4 suivants: affectés ont mineur, pas non-affectés (score positif = signal)
        # - 3 derniers: non-affectés ont mineur, pas affectés (score négatif)
        markers = [f'm{i}' for i in range(10)]
        freq_dict = {m: 0.2 for m in markers}  # A est mineur

        affected = [1, 2, 3, 4]
        unaffected = [5, 6, 7, 8]

        genotypes = {}
        for i, m in enumerate(markers):
            genotypes[m] = {}
            if i < 3 or i >= 7:
                # Score négatif : non-affectés portent le mineur, pas les affectés
                for a in affected:
                    genotypes[m][a] = GENO_BB
                for u in unaffected:
                    genotypes[m][u] = GENO_AA
            else:
                # Signal : affectés ont mineur, non-affectés non
                for a in affected:
                    genotypes[m][a] = GENO_AA
                for u in unaffected:
                    genotypes[m][u] = GENO_BB

        pedigree = FakePedigree(affected, unaffected)
        region = {'marker_names': markers}
        bp = list(range(1000, 11000, 1000))
        cm = [float(x) for x in range(10)]
        chr_results = _make_chr_results(markers, bp, cm)

        result = analyze_shared_region(region, chr_results, pedigree, genotypes, freq_dict)

        assert result['has_shared'] is True
        # Le segment optimal devrait être m3-m6 (indices 3-6)
        assert result['start_marker'] == 'm3'
        assert result['end_marker'] == 'm6'
        assert result['n_affected_sharing'] == 4
        assert result['n_unaffected_sharing'] == 0
