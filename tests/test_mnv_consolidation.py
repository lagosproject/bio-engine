# SPDX-License-Identifier: MIT
"""
Fixture battery for conservative MNV consolidation.

Synthetic codon-aligned reference (treated as a coding transcript, CDS at c.1):

    c.:  1  2  3 | 4  5  6 | 7  8  9 |10 11 12|13 14 15|16 17 18|19 20 21|22 23 24
    nt:  A  T  G | C  A  C | A  A  A | G  G  T | T  T  T | C  G  A | C  C  C | T  A  A
    aa:    M     |   H     |   K     |   G     |   F     |   R     |   P     |   *

Cases exercise the merge policy and the biological correctness win.
"""
import os
import sys
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utilities.mnv_consolidation import consolidate_snv_runs

REFERENCE = "ATGCACAAAGGTTTTCGACCCTAA"  # 1-based c. positions
AC = "TEST_MNV.1"
COLUMNS = ["pos", "ref", "alt", "type", "genotype"]

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def _codon(c_pos):
    """Reference codon containing 1-based coding position c_pos."""
    start = ((c_pos - 1) // 3) * 3
    return REFERENCE[start:start + 3]


def _apply(c_pos, alt):
    cod = list(_codon(c_pos))
    cod[(c_pos - 1) % 3] = alt
    return CODON_TABLE["".join(cod)]


def _hgvs_set(columns, rows):
    i = columns.index("hgvs")
    return {r[i] for r in rows if r[i]}


def _run(rows):
    return consolidate_snv_runs(COLUMNS, rows, AC, "c", cds_start=0)


class TestMNVConsolidation(unittest.TestCase):
    # 1) Contiguous hom run -> single delins (the core behavior).
    def test_contiguous_hom_run_merges_to_delins(self):
        rows = [
            [4, "C", "T", "SNV", "hom. ALT"],
            [5, "A", "G", "SNV", "hom. ALT"],
            [6, "C", "A", "SNV", "hom. ALT"],
        ]
        cols, out = _run(rows)
        self.assertEqual(len(out), 1)
        self.assertEqual(_hgvs_set(cols, out), {"TEST_MNV.1:c.4_6delinsTGA"})
        self.assertEqual(out[0][cols.index("type")], "MNV")

    # 2) Two-base hom run merges.
    def test_two_base_hom_run_merges(self):
        rows = [[10, "G", "A", "SNV", "hom. ALT"], [11, "G", "T", "SNV", "hom. ALT"]]
        cols, out = _run(rows)
        self.assertEqual(_hgvs_set(cols, out), {"TEST_MNV.1:c.10_11delinsAT"})

    # 3) Heterozygous contiguous run is NOT merged (unphased).
    def test_het_run_not_merged(self):
        rows = [[13, "T", "A", "SNV", "het."], [14, "T", "C", "SNV", "het."]]
        cols, out = _run(rows)
        self.assertEqual(len(out), 2)

    # 4) Long het noise run is NOT merged (the CACNA1A-style artifact guard).
    def test_het_noise_not_merged(self):
        rows = [
            [16, "C", "A", "SNV", "het."],
            [17, "G", "T", "SNV", "het."],
            [18, "A", "C", "SNV", "het."],
        ]
        cols, out = _run(rows)
        self.assertEqual(len(out), 3)

    # 5) Isolated SNV passes through untouched.
    def test_isolated_snv_intact(self):
        rows = [[19, "C", "A", "SNV", "hom. ALT"]]
        cols, out = _run(rows)
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0][:5], [19, "C", "A", "SNV", "hom. ALT"])

    # 6) A gap breaks the run (no non-contiguous merge in this conservative pass).
    def test_gap_breaks_run(self):
        rows = [[22, "T", "A", "SNV", "hom. ALT"], [24, "A", "T", "SNV", "hom. ALT"]]
        cols, out = _run(rows)
        self.assertEqual(len(out), 2)

    # 7) Multi-allelic position (two SNVs at same pos) is never merged.
    def test_multiallelic_position_not_merged(self):
        rows = [
            [4, "C", "T", "SNV", "hom. ALT"],
            [5, "A", "G", "SNV", "hom. ALT"],
            [5, "A", "C", "SNV", "hom. ALT"],  # second allele at pos 5
        ]
        cols, out = _run(rows)
        self.assertEqual(len(out), 3)  # pos 5 ineligible -> nothing merges

    # 8) Indels pass through; an adjacent hom SNV pair around them still merges.
    def test_indels_passthrough(self):
        rows = [
            [10, "G", "A", "SNV", "hom. ALT"],
            [11, "G", "T", "SNV", "hom. ALT"],
            [13, "TT", "T", "Deletion", "hom. ALT"],
        ]
        cols, out = _run(rows)
        hgvs = _hgvs_set(cols, out)
        self.assertIn("TEST_MNV.1:c.10_11delinsAT", hgvs)
        self.assertEqual(len(out), 2)  # merged pair + untouched deletion

    # 9) THE CORRECTNESS WIN: per-SNV annotation vs merged consequence.
    def test_biological_correctness_stop_codon(self):
        # Codon c.4-6 = CAC (His). Three contiguous hom changes -> TGA (stop).
        per_snv = {_apply(4, "T"), _apply(5, "G"), _apply(6, "A")}
        self.assertEqual(per_snv, {"Y", "R", "Q"})          # 3 sense missense (wrong)
        merged_codon = "TGA"
        self.assertEqual(CODON_TABLE[merged_codon], "*")    # true effect: STOP
        # And the consolidator produces exactly that delins:
        rows = [[4, "C", "T", "SNV", "hom. ALT"],
                [5, "A", "G", "SNV", "hom. ALT"],
                [6, "C", "A", "SNV", "hom. ALT"]]
        cols, out = _run(rows)
        self.assertEqual(_hgvs_set(cols, out), {"TEST_MNV.1:c.4_6delinsTGA"})


if __name__ == "__main__":
    unittest.main(verbosity=2)
