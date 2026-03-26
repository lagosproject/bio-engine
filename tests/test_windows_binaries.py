import os
import sys
import logging
import shutil
import tempfile
import subprocess

# Add the project root to sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from services.reference import ensure_indexed
from core.config import settings

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_windows_binary_functionality():
    logger.info(f"OS: {sys.platform}")
    logger.info(f"bgzip_path: {settings.bgzip_path}")
    logger.info(f"samtools_path: {settings.samtools_path}")
    
    # 1. Check if binaries are executable
    for tool in ["bgzip", "samtools"]:
        path = getattr(settings, f"{tool}_path")
        resolved_path = shutil.which(path)
        if resolved_path:
            logger.info(f"Found {tool} at: {resolved_path}")
            # Try to run it with --version or similar
            try:
                cmd = [resolved_path, "--version" if tool == "samtools" else "-h"]
                result = subprocess.run(cmd, capture_output=True, text=True)
                logger.info(f"{tool} version/help check: {result.returncode == 0}")
            except Exception as e:
                logger.warning(f"Could not execute {tool}: {e}")
        else:
            logger.error(f"{tool} ('{path}') NOT found in PATH")
            sys.exit(1)

    # 2. Test actual compression and indexing
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "test.fasta")
        gz_path = fasta_path + ".gz"
        
        logger.info("Creating dummy FASTA file...")
        with open(fasta_path, "w") as f:
            f.write(">test_seq\n")
            f.write("ATGC" * 12500 + "\n") # 50kb sequence to trigger ensure_indexed
            
        logger.info(f"Testing ensure_indexed on {fasta_path}...")
        try:
            # ensure_indexed will try pysam first, then fall back to binaries
            indexed_path = ensure_indexed(fasta_path)
            logger.info(f"Successfully processed reference at: {indexed_path}")
            
            if not indexed_path.endswith(".gz"):
                raise RuntimeError("Expected compressed output")
                
            fai_path = indexed_path + ".fai"
            if not os.path.exists(fai_path):
                raise RuntimeError("Index file (.fai) was not created")
            
            logger.info("Binary evaluation successful!")
        except Exception as e:
            logger.error(f"Functionality test failed: {e}")
            sys.exit(1)

    logger.info("WINDOWS BINARY VERIFICATION PASSED")

if __name__ == "__main__":
    test_windows_binary_functionality()
