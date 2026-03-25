"""
Verification test for Windows binaries (bgzip and samtools).
This script exercises the reference loading and indexing logic.
"""

import os
import sys
import logging
import shutil
import subprocess

# Add the project root to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from services.reference import load_reference, ensure_indexed
from core.config import settings

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_ng008866_indexing():
    accession = "NG_008866"
    logger.info(f"Testing with accession: {accession}")
    
    # Ensure cache dir is clean for this test
    if os.path.exists(settings.cache_dir):
        # We don't want to delete the whole cache if it's production, 
        # but in CI it's fine.
        logger.info(f"Using cache dir: {settings.cache_dir}")

    try:
        # 1. Fetch and Index
        # This will call ensure_indexed internally if sequence length >= 50kb
        ref_path = load_reference(accession)
        logger.info(f"Reference loaded at: {ref_path}")
        
        # 2. Verify files exist
        if not ref_path.endswith(".gz"):
            raise RuntimeError(f"Expected compressed reference for {accession}, but got {ref_path}")
            
        fai_index = ref_path + ".fai"
        if not os.path.exists(fai_index):
            raise RuntimeError(f"FAI index missing: {fai_index}")
            
        logger.info("Successfully fetched, compressed, and indexed NG_008866")
        
    except Exception as e:
        logger.error(f"Test failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # In CI, we might want to override the paths to binaries if they are not in PATH
    logger.info(f"OS: {sys.platform}")
    logger.info(f"bgzip_path: {settings.bgzip_path}")
    logger.info(f"samtools_path: {settings.samtools_path}")
    
    # Check if binaries are executable
    for tool in ["bgzip", "samtools"]:
        path = getattr(settings, f"{tool}_path")
        if shutil.which(path):
            logger.info(f"Found {tool} at: {shutil.which(path)}")
        else:
            logger.warning(f"{tool} ('{path}') NOT found in PATH")

    test_ng008866_indexing()
