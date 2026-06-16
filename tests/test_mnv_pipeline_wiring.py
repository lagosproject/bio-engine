# SPDX-License-Identifier: MIT
"""
Wiring test: TracyPipeline._ensure_hgvs_column applies MNV consolidation only
when opt-in (consolidate_mnv=True). Uses an NC_ genomic accession so no
reference-feature lookup (network/disk) is needed.
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utilities.tracy_pipeline import TracyPipeline

AC = "NC_000019.10"  # genomic -> h_type 'g', positions map directly


def _data():
    return {
        "variants": {
            "columns": ["pos", "ref", "alt", "type", "genotype"],
            "rows": [
                [100, "C", "T", "SNV", "hom. ALT"],
                [101, "A", "G", "SNV", "hom. ALT"],
                [102, "C", "A", "SNV", "hom. ALT"],
                [200, "G", "A", "SNV", "het."],   # isolated het -> untouched
            ],
        }
    }


class TestMNVPipelineWiring(unittest.TestCase):
    def setUp(self):
        self.tp = TracyPipeline(output_dir=tempfile.mkdtemp())

    def _hgvs(self, data):
        cols = data["variants"]["columns"]
        i = cols.index("hgvs")
        return [r[i] for r in data["variants"]["rows"]]

    def test_off_by_default_no_merge(self):
        data = _data()
        self.tp._ensure_hgvs_column(data, source_ac=AC, consolidate_mnv=False)
        self.assertEqual(len(data["variants"]["rows"]), 4)
        self.assertIn("NC_000019.10:g.100C>T", self._hgvs(data))

    def test_opt_in_merges_contiguous_hom_run(self):
        data = _data()
        self.tp._ensure_hgvs_column(data, source_ac=AC, consolidate_mnv=True)
        rows = data["variants"]["rows"]
        self.assertEqual(len(rows), 2)  # merged run + the het
        hgvs = self._hgvs(data)
        self.assertIn("NC_000019.10:g.100_102delinsTGA", hgvs)
        self.assertIn("NC_000019.10:g.200G>A", hgvs)


if __name__ == "__main__":
    unittest.main(verbosity=2)
