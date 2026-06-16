# SPDX-License-Identifier: MIT
"""
Integration battery for the combined NG_ recoding path inside
OpenCRAVATAnnotator._recode_to_vcf.

Proves the key property: with the offline mapper available, NG_ variants are
resolved correctly and Ensembl is NEVER called. If the offline mapper is down,
exactly one Ensembl recode (the anchor) is used, then the rest stay local.
"""
import os
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utilities.vep_utils import OpenCRAVATAnnotator
from utilities.ng_anchor import NGAnchorResolver

GROUND_TRUTH = {
    "NG_008866.1:g.65631C>T": "19:38494330:C:T",
    "NG_008866.1:g.65880G>A": "19:38494579:G:A",
}
BATCH = list(GROUND_TRUTH.keys()) + ["NG_008866.1:g.5131A>G"]  # extra, local-only


class TestNGRecodeIntegration(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.ann = OpenCRAVATAnnotator(assembly="GRCh38", oc_path="oc")
        # Fresh disk anchor cache per test.
        self.ann._ng_anchor = NGAnchorResolver(
            os.path.join(self.tmp.name, "ng_anchors.json"))

    def tearDown(self):
        self.tmp.cleanup()

    @patch.object(OpenCRAVATAnnotator, "_recode_to_vcf_via_ensembl")
    @patch.object(OpenCRAVATAnnotator, "_recode_to_vcf_via_opencravat")
    def test_offline_mapper_avoids_ensembl(self, mock_oc, mock_ens):
        mock_oc.return_value = {}
        mock_ens.return_value = {}
        # Offline mapper resolves everything.
        with patch.object(OpenCRAVATAnnotator, "_ng_resolve_offline",
                          side_effect=lambda v: GROUND_TRUTH.get(v) or self._linear(v)):
            res = self.ann._recode_to_vcf(BATCH)
        for v in BATCH:
            self.assertIn(v, res)
        self.assertEqual(res["NG_008866.1:g.65631C>T"], "19:38494330:C:T")
        mock_ens.assert_not_called()       # <-- the whole point
        mock_oc.assert_not_called()

    @patch.object(OpenCRAVATAnnotator, "_recode_to_vcf_via_ensembl")
    @patch.object(OpenCRAVATAnnotator, "_recode_to_vcf_via_opencravat")
    def test_ensembl_used_once_when_mapper_down(self, mock_oc, mock_ens):
        mock_oc.return_value = {}
        # One Ensembl recode answers the anchor seed; the rest go local.
        mock_ens.side_effect = lambda seeds: {s: GROUND_TRUTH[s] for s in seeds if s in GROUND_TRUTH}
        with patch.object(OpenCRAVATAnnotator, "_ng_resolve_offline", return_value=None):
            res = self.ann._recode_to_vcf(BATCH)
        # Anchor learned via exactly one Ensembl seed call, then arithmetic.
        self.assertEqual(mock_ens.call_count, 1)
        self.assertEqual(res["NG_008866.1:g.5131A>G"], "19:38433830:A:G")

    @patch.object(OpenCRAVATAnnotator, "_recode_to_vcf_via_ensembl")
    @patch.object(OpenCRAVATAnnotator, "_recode_to_vcf_via_opencravat")
    def test_transcript_resolved_locally_no_network(self, mock_oc, mock_ens):
        mock_oc.return_value = {}
        mock_ens.return_value = {}
        tx_gt = {"NM_000540.3:c.1A>G": "19:38433830:A:G",
                 "NM_023035.3:c.1A>G": "19:13506224:T:C"}
        with patch.object(OpenCRAVATAnnotator, "_tx_resolve_offline",
                          side_effect=lambda v: tx_gt.get(v)):
            res = self.ann._recode_to_vcf(list(tx_gt.keys()))
        self.assertEqual(res, tx_gt)
        mock_ens.assert_not_called()
        mock_oc.assert_not_called()

    @staticmethod
    def _linear(v):  # offset 38428699 for any RYR1 substitution
        import re
        m = re.match(r'NG_008866(?:\.\d+)?:g\.(\d+)([ACGT])>([ACGT])', v)
        if not m:
            return None
        return f"19:{int(m.group(1)) + 38428699}:{m.group(2)}:{m.group(3)}"


if __name__ == "__main__":
    unittest.main(verbosity=2)
