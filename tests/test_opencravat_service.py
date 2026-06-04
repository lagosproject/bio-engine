import os
import sys
import unittest
from unittest.mock import patch, MagicMock

# Add project root to python path to resolve imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from services import opencravat

class TestOpenCRAVATService(unittest.TestCase):
    def test_parse_module_table_ls(self):
        # Mock stdout of "oc module ls"
        mock_stdout = """
Name                        Title                              Type             Version     Data source ver  Size      
casecontrol                 Case-Control                       postaggregator   1.2.0                        31.0 kB   
clingen-converter           Clingen Allele Registry Converter  converter        1.0.1                        14.4 kB   
excelreporter               Excel Reporter                     reporter         2.1.1                        1.6 MB    
go                          Gene Ontology                      annotator        2025.11.01  2025.11.01       47.3 MB   
hg38                        UCSC hg38 Gene Mapper              mapper           1.47.2                       732.1 MB  
"""
        modules = opencravat._parse_module_table(mock_stdout, default_installed=True)
        self.assertEqual(len(modules), 5)
        
        # Verify casecontrol
        m1 = next(m for m in modules if m["name"] == "casecontrol")
        self.assertEqual(m1["title"], "Case-Control")
        self.assertEqual(m1["type"], "postaggregator")
        self.assertEqual(m1["version"], "1.2.0")
        self.assertEqual(m1["size_bytes"], 31.0 * 1024)
        self.assertTrue(m1["installed"])
        
        # Verify clingen-converter
        m2 = next(m for m in modules if m["name"] == "clingen-converter")
        self.assertEqual(m2["title"], "Clingen Allele Registry Converter")
        self.assertEqual(m2["type"], "converter")
        self.assertEqual(m2["version"], "1.0.1")
        self.assertEqual(m2["size_bytes"], int(14.4 * 1024))
        
        # Verify hg38
        m3 = next(m for m in modules if m["name"] == "hg38")
        self.assertEqual(m3["title"], "UCSC hg38 Gene Mapper")
        self.assertEqual(m3["type"], "mapper")
        self.assertEqual(m3["version"], "1.47.2")
        self.assertEqual(m3["size_bytes"], int(732.1 * 1024 * 1024))

    def test_parse_module_table_ls_a(self):
        # Mock stdout of "oc module ls -a"
        mock_stdout = """
Name                               Title                                                 Type             Installed  Store ver   Store data ver                 Local ver   Local data ver  Size      
23andme-converter                  23andMe Converter                                     converter                   1.5.1                                                                  17.6 kB   
abraom                             ABRaOM                                                annotator                   1.0.0                                                                  113.6 MB  
excelreporter                      Excel Reporter                                        reporter         yes        2.1.1                                      2.1.1                       1.6 MB    
go                                 Gene Ontology                                         annotator        yes        2025.11.01  2025.11.01                 2025.11.01  2025.11.01      47.3 MB   
"""
        modules = opencravat._parse_module_table(mock_stdout, default_installed=None)
        self.assertEqual(len(modules), 4)

        # 23andme-converter (Not installed)
        m1 = next(m for m in modules if m["name"] == "23andme-converter")
        self.assertEqual(m1["title"], "23andMe Converter")
        self.assertEqual(m1["type"], "converter")
        self.assertEqual(m1["version"], "1.5.1")
        self.assertFalse(m1["installed"])
        
        # excelreporter (Installed)
        m2 = next(m for m in modules if m["name"] == "excelreporter")
        self.assertEqual(m2["title"], "Excel Reporter")
        self.assertEqual(m2["type"], "reporter")
        self.assertEqual(m2["version"], "2.1.1")
        self.assertTrue(m2["installed"])

    @patch("services.opencravat._run")
    @patch("services.opencravat._get_data_dir")
    def test_get_status_installed(self, mock_get_data_dir, mock_run):
        mock_run.return_value = (0, "3.1.1", "")
        mock_get_data_dir.return_value = None
        
        status = opencravat.get_status("oc")
        self.assertTrue(status["installed"])
        self.assertEqual(status["version"], "3.1.1")
        self.assertIsNone(status["error"])

    @patch("services.opencravat._run")
    def test_get_status_not_installed(self, mock_run):
        mock_run.return_value = (-1, "", "Command not found")
        
        status = opencravat.get_status("/invalid/path/oc")
        self.assertFalse(status["installed"])
        self.assertIsNotNone(status["error"])

if __name__ == "__main__":
    unittest.main()
