import os
import sys
import tempfile
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_bgzip():
    bgzip_path = os.environ.get("BIO_BGZIP_PATH")
    if not bgzip_path:
        logger.error("BIO_BGZIP_PATH must be set for this test")
        sys.exit(1)
        
    logger.info(f"Testing bgzip at: {bgzip_path}")
    
    if not os.path.exists(bgzip_path):
        logger.warning(f"File not found exactly at: {bgzip_path}. It might rely on PATH resolution.")
        
    try:
        # Create a tiny dummy fasta
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta = os.path.join(tmpdir, "test.fasta")
            with open(fasta, "w") as f:
                f.write(">1\nATGC\n")
                
            cmd = [bgzip_path, "-c", fasta]
            logger.info(f"Running command: {cmd}")
            
            # Remove any msys/mingw/git paths from PATH to simulate the raw production environment
            env = os.environ.copy()
            clean_path = []
            for p in env.get("PATH", "").split(os.pathsep):
                p_lower = p.lower()
                if "msys" not in p_lower and "mingw" not in p_lower and "git" not in p_lower:
                    clean_path.append(p)
            env["PATH"] = os.pathsep.join(clean_path)
            
            logger.info(f"Cleaned PATH: {env['PATH']}")
            logger.info("Executing...")
            res = subprocess.run(cmd, capture_output=True, env=env)
            
            logger.info(f"STDOUT: {res.stdout}")
            logger.info(f"STDERR: {res.stderr}")
            
            if res.returncode != 0:
                logger.error(f"Command failed with exit code: {res.returncode}")
                # Print hex for missing DLL debugging
                hex_code = hex(res.returncode & 0xFFFFFFFF)
                logger.error(f"Hex exit code: {hex_code}")
                if hex_code == "0xc0000135" or str(res.returncode) == "3221225781":
                    logger.error("0xc0000135 indicates STATUS_DLL_NOT_FOUND. bgzip is missing a required DLL.")
                sys.exit(1)
            else:
                logger.info("bgzip executed successfully.")
    except Exception as e:
        logger.error(f"Exception during execution: {e}")
        sys.exit(1)

if __name__ == "__main__":
    test_bgzip()
