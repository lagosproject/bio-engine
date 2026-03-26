import os
import sys
import logging
import tempfile
import shutil

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_pysam_capabilities():
    try:
        import pysam
        logger.info(f"pysam version: {pysam.__version__}")
    except ImportError:
        logger.error("pysam not installed")
        sys.exit(1)

    # Create a temporary directory for tests
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, "test.fasta")
        gz_path = fasta_path + ".gz"
        
        # 1. Create a dummy FASTA file
        logger.info("Creating dummy FASTA file...")
        with open(fasta_path, "w") as f:
            f.write(">test_seq\n")
            f.write("ATGC" * 1000 + "\n") # 4kb sequence
            
        # 2. Test bgzip compression via pysam
        logger.info(f"Testing pysam.tabix_compress on {fasta_path}...")
        try:
            pysam.tabix_compress(fasta_path, gz_path, force=True)
            if not os.path.exists(gz_path):
                raise RuntimeError("pysam.tabix_compress failed to create .gz file")
            logger.info(f"Successfully created compressed file: {gz_path}")
        except Exception as e:
            logger.error(f"pysam.tabix_compress failed: {e}")
            sys.exit(1)
            
        # 3. Test faidx indexing via pysam
        logger.info(f"Testing pysam.faidx on {gz_path}...")
        try:
            pysam.faidx(gz_path)
            fai_path = gz_path + ".fai"
            if not os.path.exists(fai_path):
                raise RuntimeError("pysam.faidx failed to create .fai file")
            logger.info(f"Successfully created index file: {fai_path}")
        except Exception as e:
            logger.error(f"pysam.faidx failed: {e}")
            sys.exit(1)
            
        # 4. Verify we can read from the indexed file
        logger.info("Verifying sequence retrieval...")
        with pysam.FastaFile(gz_path) as fa:
            seq = fa.fetch("test_seq")
            if seq != "ATGC" * 1000:
                raise RuntimeError("Sequence fetched from indexed file does not match original")
            logger.info("Sequence verification passed!")

    logger.info("ALL PYSAM TESTS PASSED ON LINUX")

if __name__ == "__main__":
    test_pysam_capabilities()
