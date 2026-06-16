# SPDX-License-Identifier: MIT
"""
Battery for LocalCoordinateMapper (fully-offline NG_ -> chromosome via NCBI
RefSeqGene alignments). We inject a pre-built SQLite DB so the math/strand
logic is verified deterministically without the large FTP download.

RYR1 ground truth (from real jobs / Ensembl):
    NG_008866.1:g.65631C>T -> 19:38494330:C:T   (plus strand)
"""
import os
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utilities.local_coordinate_mapper import LocalCoordinateMapper


def _build_db(path, rows):
    conn = sqlite3.connect(str(path))
    conn.execute("""CREATE TABLE mappings (
        ng_accession TEXT PRIMARY KEY, chromosome TEXT, chrom_start INTEGER,
        chrom_end INTEGER, strand TEXT, target_strand TEXT)""")
    conn.executemany(
        "INSERT INTO mappings VALUES (?,?,?,?,?,?)", rows)
    conn.commit()
    conn.close()


class TestLocalCoordinateMapper(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def _mapper_with_db(self, rows):
        m = LocalCoordinateMapper(assembly="GRCh38")
        m.db_path = Path(self.tmp.name) / "db.sqlite"
        _build_db(m.db_path, rows)
        return m

    def test_plus_strand_matches_ground_truth(self):
        # chrom_start + local_pos - 1 = 38428700 + 65631 - 1 = 38494330
        m = self._mapper_with_db(
            [("NG_008866", "19", 38428700, 38587565, "+", "+")])
        self.assertEqual(
            m.resolve_ng_variant("NG_008866.1:g.65631C>T"), "19:38494330:C:T")
        m.close()

    def test_minus_strand_reverses_and_complements(self):
        # reverse: genomic = chrom_end - local_pos + 1 ; bases complemented
        # chrom_end 5000, local 100 -> 4901 ; C>T -> G>A
        m = self._mapper_with_db(
            [("NG_999999", "7", 4000, 5000, "-", "+")])
        self.assertEqual(
            m.resolve_ng_variant("NG_999999.1:g.100C>T"), "7:4901:G:A")
        m.close()

    def test_unknown_accession_returns_none(self):
        m = self._mapper_with_db(
            [("NG_008866", "19", 38428700, 38587565, "+", "+")])
        self.assertIsNone(m.resolve_ng_variant("NG_111111.1:g.10A>G"))
        m.close()

    def test_non_ng_returns_none(self):
        m = self._mapper_with_db(
            [("NG_008866", "19", 38428700, 38587565, "+", "+")])
        self.assertIsNone(m.resolve_ng_variant("NM_023035.2:c.100A>G"))
        m.close()


if __name__ == "__main__":
    unittest.main(verbosity=2)
