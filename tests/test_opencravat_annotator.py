import os
import sys
import unittest
from unittest.mock import patch, MagicMock

# Add project root to python path to resolve imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utilities.vep_utils import OpenCRAVATAnnotator

class TestOpenCRAVATAnnotator(unittest.TestCase):
    def setUp(self):
        # We can use the default "oc" binary installed in the conda environment
        self.annotator = OpenCRAVATAnnotator(assembly="GRCh38", oc_path="oc")

    @patch("utilities.vep_utils.OpenCRAVATAnnotator._recode_to_vcf")
    def test_clinvar_annotation_success(self, mock_recode):
        # Mock the HGVS to VCF recoding step to bypass the Ensembl API network calls
        # 17:43045712:T:C corresponds to the BRCA1 Tyr1853Cys variant (Pathogenic in ClinVar)
        mock_recode.return_value = {
            "NC_000017.11:g.43045712T>C": "17:43045712:T:C"
        }

        variants = ["NC_000017.11:g.43045712T>C"]
        results = self.annotator.get_annotations(variants)

        # Verify that we got results back
        self.assertIn("NC_000017.11:g.43045712T>C", results)
        anno = results["NC_000017.11:g.43045712T>C"]

        # Check mapped core features
        self.assertEqual(anno["gene_symbol"], "BRCA1")
        self.assertEqual(anno["consequence"], "MIS")
        self.assertEqual(anno["hgvs_p"], "p.Tyr1853Cys")

        # Verify ClinVar specific data
        self.assertIn("Pathogenic/Likely pathogenic", anno["clin_sig"])
        self.assertTrue(any("Hereditary breast ovarian cancer syndrome" in p for p in anno["phenotype"]))

        # Check raw oc_data presence
        self.assertIn("oc_data", anno)
        self.assertEqual(anno["oc_data"].get("clinvar__id"), "55627")
        self.assertEqual(anno["oc_data"].get("clinvar__germline_or_somatic"), "germline")

    def test_empty_variants(self):
        results = self.annotator.get_annotations([])
        self.assertEqual(results, {})

    @patch("utilities.vep_utils.OpenCRAVATAnnotator._recode_to_vcf")
    def test_unresolvable_variants(self, mock_recode):
        mock_recode.return_value = {}
        results = self.annotator.get_annotations(["invalid_variant"])
        self.assertEqual(results, {})

    @patch("utilities.vep_utils.OpenCRAVATAnnotator._recode_to_vcf")
    @patch("subprocess.run")
    def test_binary_not_found(self, mock_run, mock_recode):
        mock_recode.return_value = {
            "NC_000017.11:g.43045712T>C": "17:43045712:T:C"
        }
        # Force a FileNotFoundError to simulate missing binary
        mock_run.side_effect = FileNotFoundError()
        
        results = self.annotator.get_annotations(["NC_000017.11:g.43045712T>C"])
        self.assertEqual(results, {})

if __name__ == "__main__":
    unittest.main()
