# SPDX-License-Identifier: MIT
import os
import sys
import unittest
from unittest.mock import patch

# Add project root to python path to resolve imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utilities.vep_utils import OpenCRAVATAnnotator, VEPAnnotator


class TestVEPAnnotatorModes(unittest.TestCase):
    @patch("utilities.vep_utils.OpenCRAVATAnnotator.get_annotations")
    @patch("utilities.vep_utils.OnlineVEPAnnotator.get_annotations")
    @patch("core.cache.cache.get_many")
    @patch("core.cache.cache.set_many")
    def test_mode_opencravat_bypasses_online(self, mock_set_many, mock_get_many, mock_online_annotate, mock_oc_annotate):
        mock_get_many.return_value = {}
        mock_oc_annotate.return_value = {
            "var1": {"gene_symbol": "GENE_OC", "clin_sig": ["Pathogenic"]}
        }

        annotator = VEPAnnotator(mode="opencravat")
        results = annotator.get_annotations(["var1"])

        # Verify OpenCRAVAT was called
        mock_oc_annotate.assert_called_once_with(["var1"])
        # Verify Online VEP was bypassed
        mock_online_annotate.assert_not_called()

        self.assertEqual(results["var1"]["gene_symbol"], "GENE_OC")
        self.assertEqual(results["var1"]["clin_sig"], ["Pathogenic"])

    @patch("utilities.vep_utils.OpenCRAVATAnnotator.get_annotations")
    @patch("utilities.vep_utils.OnlineVEPAnnotator.get_annotations")
    @patch("core.cache.cache.get_many")
    @patch("core.cache.cache.set_many")
    def test_mode_online_runs_both_and_merges(self, mock_set_many, mock_get_many, mock_online_annotate, mock_oc_annotate):
        mock_get_many.return_value = {}
        mock_oc_annotate.return_value = {
            "var1": {"gene_symbol": "GENE_OC", "clin_sig": ["Pathogenic"], "consequence": "missense"}
        }
        mock_online_annotate.return_value = {
            "var1": {"gene_symbol": "GENE_VEP", "clin_sig": ["Benign"], "polyphen": "probably damaging", "vep_raw": {"key": "val"}}
        }

        annotator = VEPAnnotator(mode="online")
        results = annotator.get_annotations(["var1"])

        # Verify both were called
        mock_oc_annotate.assert_called_once_with(["var1"])
        mock_online_annotate.assert_called_once_with(["var1"])

        # Verify merge logic
        self.assertEqual(results["var1"]["gene_symbol"], "GENE_VEP")  # VEP preferred/overwritten or merged
        self.assertEqual(results["var1"]["clin_sig"], ["Benign", "Pathogenic"])  # Merged and sorted list
        self.assertEqual(results["var1"]["polyphen"], "probably damaging")
        self.assertEqual(results["var1"]["consequence"], "missense")
        self.assertEqual(results["var1"]["vep_raw"], {"key": "val"})

    @patch("utilities.vep_utils.OpenCRAVATAnnotator._recode_to_vcf_via_opencravat")
    @patch("utilities.vep_utils.OpenCRAVATAnnotator._recode_to_vcf_via_ensembl")
    @patch("core.cache.cache.get_many")
    @patch("core.cache.cache.set")
    def test_recode_to_vcf_fallback(self, mock_cache_set, mock_get_many, mock_ensembl_recode, mock_oc_recode):
        # Cache misses all
        mock_get_many.return_value = {}

        # OC fails to resolve var2, resolves var1
        mock_oc_recode.return_value = {"var1": "1:100:A:T"}
        # Ensembl fallback resolves var2
        mock_ensembl_recode.return_value = {"var2": "2:200:G:C"}

        annotator = OpenCRAVATAnnotator()
        recode_map = annotator._recode_to_vcf(["var1", "var2"])

        mock_oc_recode.assert_called_once_with(["var1", "var2"])
        # Fallback should only query missing variant "var2"
        mock_ensembl_recode.assert_called_once_with(["var2"])

        self.assertEqual(recode_map, {
            "var1": "1:100:A:T",
            "var2": "2:200:G:C"
        })

if __name__ == "__main__":
    unittest.main()
