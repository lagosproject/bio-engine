# SPDX-License-Identifier: MIT
"""
Battery for TranscriptCoordinateMapper (offline NM_/NR_ -> chromosome, spliced).

Uses small synthetic transcripts with an injected SQLite DB so the exon-walk,
CDS-offset and strand/complement math are verified deterministically without the
UCSC download. Two exons of 10 bp each let us exercise the splice boundary.
"""
import os
import sqlite3
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utilities.transcript_coordinate_mapper import TranscriptCoordinateMapper

# exon arrays are 0-based half-open (UCSC convention)
ROWS = [
    # name, chrom, strand, cds_start, cds_end, exon_starts, exon_ends
    ("NM_000001", "1", "+", 1003, 2007, "1000,2000", "1010,2010"),
    ("NM_000002", "2", "-", 1003, 2007, "1000,2000", "1010,2010"),
]


def _build_db(path, rows):
    conn = sqlite3.connect(str(path))
    conn.execute("""CREATE TABLE transcripts (
        name TEXT PRIMARY KEY, chrom TEXT, strand TEXT,
        cds_start INTEGER, cds_end INTEGER, exon_starts TEXT, exon_ends TEXT)""")
    conn.executemany("INSERT INTO transcripts VALUES (?,?,?,?,?,?,?)", rows)
    conn.commit()
    conn.close()


class TestTranscriptCoordinateMapper(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.m = TranscriptCoordinateMapper("GRCh38")
        self.m.db_path = Path(self.tmp.name) / "tx.sqlite"
        _build_db(self.m.db_path, ROWS)

    def tearDown(self):
        self.m.close()
        self.tmp.cleanup()

    def r(self, v):
        return self.m.resolve_transcript_variant(v)

    # --- plus strand ---
    def test_plus_c1_is_cds_start(self):
        self.assertEqual(self.r("NM_000001.3:c.1A>G"), "1:1004:A:G")

    def test_plus_last_base_of_exon1(self):
        # c.7 -> genomic 1009(0-based) -> 1010
        self.assertEqual(self.r("NM_000001.3:c.7C>T"), "1:1010:C:T")

    def test_plus_splice_into_exon2(self):
        # c.8 jumps the intron to exon2 start 2000(0-based) -> 2001
        self.assertEqual(self.r("NM_000001.3:c.8G>A"), "1:2001:G:A")

    # --- minus strand (positions mirror, bases complement) ---
    def test_minus_c1_is_cds_end(self):
        # first coding base = highest coding coord = cds_end-1 = 2006 -> 2007 ; A>G -> T>C
        self.assertEqual(self.r("NM_000002.3:c.1A>G"), "2:2007:T:C")

    def test_minus_splice_into_exon1(self):
        # c.8 crosses to exon1 high base 1009(0-based) -> 1010 ; G>A -> C>T
        self.assertEqual(self.r("NM_000002.3:c.8G>A"), "2:1010:C:T")

    # --- non-coding transcript uses transcript position directly ---
    def test_noncoding_n_position(self):
        # n.1 = first transcript base (plus) = exon1 start 1000 -> 1001
        self.assertEqual(self.r("NM_000001.3:n.1A>G"), "1:1001:A:G")

    # --- version-agnostic lookup ---
    def test_version_stripped(self):
        self.assertEqual(self.r("NM_000001:c.1A>G"), "1:1004:A:G")

    # --- 5'UTR (c.-N): N bases 5' of c.1 ---
    def test_plus_5utr(self):
        self.assertEqual(self.r("NM_000001.3:c.-1A>G"), "1:1003:A:G")  # just before c.1
        self.assertEqual(self.r("NM_000001.3:c.-3A>G"), "1:1001:A:G")  # exon1 start

    def test_plus_5utr_out_of_bounds(self):
        self.assertIsNone(self.r("NM_000001.3:c.-4A>G"))  # before transcript start

    # --- 3'UTR (c.*N): N bases 3' of the stop ---
    def test_plus_3utr(self):
        self.assertEqual(self.r("NM_000001.3:c.*1A>G"), "1:2008:A:G")  # after last coding base
        self.assertEqual(self.r("NM_000001.3:c.*3A>G"), "1:2010:A:G")  # last transcript base

    def test_plus_3utr_out_of_bounds(self):
        self.assertIsNone(self.r("NM_000001.3:c.*4A>G"))  # past transcript end

    # --- UTR on minus strand: positions mirror and bases complement ---
    def test_minus_5utr(self):
        self.assertEqual(self.r("NM_000002.3:c.-1G>A"), "2:2008:C:T")

    def test_minus_3utr(self):
        self.assertEqual(self.r("NM_000002.3:c.*1G>A"), "2:1003:C:T")

    # --- conservative: intronic offsets and indels -> None (network fallback) ---
    def test_intronic_offset_none(self):
        self.assertIsNone(self.r("NM_000001.3:c.7+3A>G"))
        self.assertIsNone(self.r("NM_000001.3:c.7-2A>G"))
        self.assertIsNone(self.r("NM_000001.3:c.*5+1A>G"))

    def test_noncoding_rejects_utr_markers(self):
        self.assertIsNone(self.r("NM_000001.3:n.-1A>G"))
        self.assertIsNone(self.r("NM_000001.3:n.*1A>G"))

    def test_indel_none(self):
        self.assertIsNone(self.r("NM_000001.3:c.7_9del"))

    def test_unknown_transcript_none(self):
        self.assertIsNone(self.r("NM_999999.3:c.1A>G"))

    def test_non_transcript_none(self):
        self.assertIsNone(self.r("NG_008866.1:g.65631C>T"))


if __name__ == "__main__":
    unittest.main(verbosity=2)
