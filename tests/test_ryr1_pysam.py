import os
import sys
import logging
import pysam

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from services import reference as ref_service

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_ryr1_pysam_indexing():
    logger.info("Starting RYR1 Gene Search and pysam Indexing Test...")
    
    # 1. Search for RYR1
    query = "RYR1"
    logger.info(f"Searching for gene: {query}")
    try:
        results = ref_service.search_reference(query)
    except Exception as e:
        logger.warning(f"Search failed with error: {e}")
        results = []
    
    found_ng = False
    target_accession = "NG_008866"
    
    if results:
        for res in results:
            acc = res['accession']
            title = res['title']
            length = res['length']
            logger.info(f"Result: {acc} - {title} ({length} bp)")
            if target_accession in acc:
                found_ng = True
                # Use the full accession (with version if present)
                target_accession = acc
                break
    else:
        logger.warning("No search results returned. This might be a transient NCBI issue.")
            
    if not found_ng:
        if results:
            logger.warning(f"{target_accession} not found in search results. Using it directly anyway.")
        else:
            logger.info(f"Using known accession directly: {target_accession}")
    else:
        logger.info(f"Successfully found target accession in search results: {target_accession}")

    # 2. Load reference (this triggers download + ensure_indexed)
    logger.info(f"Loading and indexing reference: {target_accession}")
    try:
        ref_path = ref_service.load_reference(target_accession)
        logger.info(f"Reference loaded at: {ref_path}")
    except Exception as e:
        logger.error(f"Failed to load reference {target_accession}: {e}")
        sys.exit(1)

    # 3. Verify files exist and are compressed
    # Note: ensure_indexed will only compress if sequence length >= 50,000
    # NG_008866 is ~160kb, so it should be compressed.
    
    if not ref_path.endswith(".gz"):
        # If the downloaded file is small, it won't be compressed by default in bio-engine
        # But NG_008866 is large.
        logger.warning(f"Reference path {ref_path} does not end with .gz. Checking file size...")
        if os.path.exists(ref_path):
            size = os.path.getsize(ref_path)
            if size > 50000:
                logger.error(f"Error: Large file ({size} bytes) was not compressed.")
                sys.exit(1)
            else:
                logger.info(f"File is small ({size} bytes), compression not expected.")
        else:
            logger.error(f"File {ref_path} does not exist.")
            sys.exit(1)
    else:
        logger.info("Verified: Reference is compressed (.gz)")

    index_path = ref_path + ".fai"
    if not os.path.exists(index_path):
        # We try to index even if small if tracy/pysam is available.
        logger.error(f"Error: Index file not found: {index_path}")
        sys.exit(1)
    
    logger.info(f"Verified: Index file exists at {index_path}")

    # 4. Verify we can read from the indexed file using pysam
    logger.info("Verifying sequence retrieval via pysam...")
    try:
        # FastaFile handles both .fasta and .fasta.gz (if indexed)
        with pysam.FastaFile(ref_path) as fa:
            # Get the first sequence name
            contig = fa.references[0]
            logger.info(f"Fetching sequence from contig: {contig}")
            seq = fa.fetch(contig, 0, 100)
            if len(seq) == 100:
                logger.info(f"Successfully fetched 100bp from start: {seq}")
                print("\nRYR1 PYSAM VERIFICATION TEST PASSED")
            else:
                logger.error(f"Error: Fetched sequence length is {len(seq)}, expected 100")
                sys.exit(1)
    except Exception as e:
        logger.error(f"pysam retrieval failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    test_ryr1_pysam_indexing()
