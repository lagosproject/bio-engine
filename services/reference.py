"""
Reference Biology Service
=========================

This module manages DNA reference sequences. It handles downloading
sequences from NCBI (Entrez), local caching, feature extraction, 
format conversion (GenBank to FASTA), and indexing via samtools/tracy.
"""

import gzip
import io
import logging
import os
import shutil
import subprocess
import tempfile

from Bio import Entrez, SeqIO

from core.config import settings
from core.exceptions import ReferenceError

# Configure Entrez
Entrez.email = settings.entrez_email
logger = logging.getLogger(__name__)

def check_ncbi_reference_exists(accession: str) -> bool:
    """Verifies if an accession code exists in the NCBI nucleotide database without downloading."""
    try:
        with Entrez.esearch(db="nucleotide", term=accession, retmax=1) as handle:
            record = Entrez.read(handle)
            return int(record.get("Count", 0)) > 0
    except Exception as e:
        logger.error(f"NCBI existence check failed for {accession}: {e}")
        return False

def search_reference(query: str, retmax: int = 20, assembly: str | None = None) -> list[dict]:
    """
    Searches for references matching the query.
    Prioritizes NCBI for HG38 and Ensembl for HG19/GRCh37.
    """
    results = []
    is_hg19 = assembly and assembly.upper() in ["GRCH37", "HG19"]
    is_accession = any(query.upper().startswith(prefix) for prefix in ["NM_", "NC_", "NG_", "NR_", "XM_", "XR_"])

    def search_ensembl():
        ensembl_results = []
        try:
            from core.proxy_manager import proxy_manager
            base_url = "https://grch37.rest.ensembl.org" if is_hg19 else "https://rest.ensembl.org"
            client = proxy_manager.get_client("ensembl")
            
            # Lookup symbol to get transcripts
            url = f"{base_url}/lookup/symbol/human/{query}?expand=1"
            response = client.get(url, headers={"Content-Type": "application/json"})
            if response.status_code == 200:
                data = response.json()
                for transcript in data.get("Transcript", []):
                    # Try to map to RefSeq if possible (check multiple RefSeq DB names)
                    xref_ac = ""
                    refseq_dbs = ["RefSeq_mRNA", "RefSeq_genomic", "RefSeq_mRNA_predicted", "RefSeq_ncRNA"]
                    for xref in transcript.get("ExternalReference", []):
                        if xref.get("dbname") in refseq_dbs:
                            xref_ac = xref.get("primary_id")
                            break
                    
                    ensembl_results.append({
                        "accession": xref_ac or transcript.get("id"),
                        "title": f"{transcript.get('display_name', '')} (Ensembl/RefSeq mapping for {assembly})",
                        "length": transcript.get("end", 0) - transcript.get("start", 0),
                        "source": "Ensembl"
                    })
        except Exception as e:
            logger.warning(f"Ensembl search failed for {query}: {e}")
        return ensembl_results

    def search_ncbi():
        ncbi_results = []
        try:
            # We use a cleaner term without assembly over-filtering for transcripts
            # because RefSeq transcripts (NM_) are often not tagged with build version in their metadata
            term = query
            if not is_accession:
                term = f"{query} AND human[Organism] AND refseq[Filter]"

            with Entrez.esearch(db="nucleotide", term=term, retmax=retmax) as handle:
                record = Entrez.read(handle)
                id_list = record.get("IdList", [])
                
            if id_list:
                with Entrez.esummary(db="nucleotide", id=",".join(id_list)) as sum_handle:
                    summaries = Entrez.read(sum_handle)
                    for s in summaries:
                        ncbi_results.append({
                            "accession": s.get("Caption", ""),
                            "title": s.get("Title", ""),
                            "length": s.get("Length", 0),
                            "source": "NCBI"
                        })
        except Exception as e:
            logger.error(f"Failed to search NCBI reference for {query}: {e}")
        return ncbi_results

    # Execution order based on user preference
    if is_hg19 and not is_accession:
        # For HG19, Ensembl is much better for finding mapped transcripts
        results.extend(search_ensembl())
        ncbi_res = search_ncbi()
        for r in ncbi_res:
            if not any(res["accession"] == r["accession"] for res in results):
                results.append(r)
    else:
        # Default or HG38: NCBI first
        results.extend(search_ncbi())
        if assembly and not is_accession:
            ensembl_res = search_ensembl()
            for r in ensembl_res:
                if not any(res["accession"] == r["accession"] for res in results):
                    results.append(r)
    
    return results

def load_reference(ref_input: str, assembly: str | None = None) -> str:
    """Handles both local files and NCBI IDs, returning a valid FASTA path."""
    if os.path.exists(ref_input):
        return _process_local_fasta(ref_input)
    return _fetch_ncbi_reference(ref_input, assembly=assembly)

def _process_local_fasta(file_path: str) -> str:
    """Read and sanitize a local FASTA file."""
    try:
        with open(file_path, encoding="utf-8-sig") as f:
            content = f.read().strip()

        if not content.startswith(">"):
            if all(c.upper() in "ATGCN" for c in content.replace("\n", "").replace("\r", "")):
                content = ">Reference_Sequence\n" + content
                temp_fasta = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w")
                temp_fasta.write(content)
                temp_fasta.close()
                return temp_fasta.name
            raise ReferenceError("File does not start with '>' and is not raw DNA.")

        handle = io.StringIO(content)
        if not list(SeqIO.parse(handle, "fasta")):
            raise ReferenceError("No valid FASTA sequence found.")

        return file_path
    except Exception as e:
        raise ReferenceError(f"Failed to process local FASTA: {e}")


def _fetch_ncbi_reference(accession: str, assembly: str | None = None) -> str:
    """Fetch reference from NCBI and cache it. Use .gz for large files (>=50kb)."""
    cache_dir = settings.get_cache_dir(assembly=assembly)
    os.makedirs(cache_dir, exist_ok=True)
    
    fasta_cache_file = os.path.join(cache_dir, f"{accession}.fasta")
    gb_cache_file = os.path.join(cache_dir, f"{accession}.gb")

    if os.path.exists(fasta_cache_file) and os.path.exists(gb_cache_file):
        # Check size of existing fasta to decide if we should gzip it now
        try:
            record = SeqIO.read(fasta_cache_file, "fasta")
            if len(record.seq) >= 50000:
                # If large, return the compressed version via ensure_indexed
                # ensure_indexed will create .fasta.gz (bgzip) if needed
                return ensure_indexed(fasta_cache_file)
        except Exception:
            pass # Fasta might be corrupt or huge, proceed to re-download
        return fasta_cache_file

    try:
        # Fetch GenBank format (contains features + sequence)
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            record = SeqIO.read(handle, "genbank")

            # Save GenBank file
            with open(gb_cache_file, "w") as f:
                SeqIO.write(record, f, "genbank")

            # Save as regular .fasta
            with open(fasta_cache_file, "w") as f:
                SeqIO.write(record, f, "fasta")

            # Decide on format based on length
            if len(record.seq) >= 50000:
                 return ensure_indexed(fasta_cache_file)
            else:
                return fasta_cache_file

    except Exception as e:
        raise ReferenceError(f"NCBI fetch failed for {accession}: {e}")

def get_reference_features(ref_input: str, assembly: str | None = None) -> list:
    """
    Extracts features (exons, CDS, introns, etc.) from the reference.
    Ref_input can be an accession code or a file path.
    """
    # Ensure reference is loaded/fetched
    try:
        loaded_path = load_reference(ref_input, assembly=assembly)
    except ReferenceError:
        logger.warning(f"Could not load reference {ref_input} for feature extraction.")
        return []

    # Determine potential GenBank path
    gb_path = None

    # Case 1: ref_input is an existing local file
    if os.path.exists(ref_input):
        base, _ = os.path.splitext(ref_input)
        potential_gb = base + ".gb"
        if os.path.exists(potential_gb):
            gb_path = potential_gb

    # Case 2: check if it's in the cache (accession)
    if not gb_path:
        # Check based on loaded path from cache
        if "ncbi_cache" in loaded_path: # Simple heuristic check if it's from our cache
             # loaded_path could be .fasta or .fasta.gz
             # Extract accession from filename
             filename = os.path.basename(loaded_path)
             # remove extensions
             accession = filename.replace(".fasta.gz", "").replace(".fasta", "")
             cache_dir = settings.get_cache_dir(assembly=assembly)
             potential_gb = os.path.join(cache_dir, f"{accession}.gb")
             if os.path.exists(potential_gb):
                 gb_path = potential_gb

    if not gb_path or not os.path.exists(gb_path):
        return []

    try:
        record = SeqIO.read(gb_path, "genbank")
        features_data = []

        for feature in record.features:
            if feature.type in ["source", "gene", "CDS", "exon", "mRNA", "tRNA", "rRNA"]:
                feature_info = {
                    "type": feature.type,
                    "start": int(feature.location.start), # 0-based
                    "end": int(feature.location.end),     # 0-based
                    "strand": feature.location.strand,
                    "qualifiers": {k: v[0] if len(v) == 1 else v for k, v in feature.qualifiers.items()}
                }
                features_data.append(feature_info)

        return features_data
    except Exception as e:
        logger.error(f"Error parsing features from {gb_path}: {e}")
        return []

def get_fasta_sequence(ref_input: str, assembly: str | None = None) -> str:
    """Returns the sequence string from a FASTA file or NCBI accession."""
    if os.path.exists(ref_input) and ref_input.endswith(".gz"):
        path = ref_input
    elif os.path.exists(ref_input) and ref_input.endswith((".fasta", ".fa", ".fna")):
        path = ref_input
    else:
        path = load_reference(ref_input, assembly=assembly)

    open_func = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"

    with open_func(path, mode) as f:
        record = SeqIO.read(f, "fasta")
        return str(record.seq)


def get_lrg_mapping(accession: str, assembly: str | None = None) -> str | None:
    """
    Attempts to resolve an LRG (Locus Reference Genomic) equivalent for a RefSeq accession.
    Many NG_ records are mirrors of LRG_ records which Ensembl understands better.
    """
    if not accession.startswith("NG_"):
        return None

    # Determine potential GenBank path in cache
    cache_dir = settings.get_cache_dir(assembly=assembly)
    gb_path = os.path.join(cache_dir, f"{accession}.gb")
    if not os.path.exists(gb_path):
        # Try to load it if possible
        try:
            load_reference(accession, assembly=assembly)
        except:
            return None
    
    if not os.path.exists(gb_path):
        return None

    try:
        import re
        record = SeqIO.read(gb_path, "genbank")
        desc = record.description
        # Pattern: ... (LRG_766) ...
        match = re.search(r"\(LRG_(\d+)\)", desc)
        if match:
            return f"LRG_{match.group(1)}"
        
        # Also check db_xref in source features
        for feature in record.features:
            if feature.type == "source":
                for xref in feature.qualifiers.get("db_xref", []):
                    if xref.startswith("LRG:"):
                        return xref.replace("LRG:", "LRG_")
    except Exception as e:
        logger.error(f"Failed to parse LRG mapping from {gb_path}: {e}")
    
    return None

def ensure_indexed(ref_path: str) -> str:
    """Compresses and indexes a large reference file. Handles recovery from bad indices."""
    
    # Check for binaries, but don't fail immediately if pysam is available
    bgzip_available = shutil.which(settings.bgzip_path) if settings.bgzip_path else None
    samtools_available = shutil.which(settings.samtools_path) if settings.samtools_path else None

    # Determine paths
    if ref_path.endswith(".gz"):
        compressed_ref = ref_path
        source_ref = ref_path[:-3] # remove .gz
    else:
        compressed_ref = ref_path + ".gz"
        source_ref = ref_path

    # Step 1: Ensure compressed file exists
    if not os.path.exists(compressed_ref):
        if os.path.exists(source_ref):
            success = False
            # Try pysam first
            try:
                import pysam
                logger.info(f"Compressing {source_ref} to {compressed_ref} using pysam.tabix_compress")
                pysam.tabix_compress(source_ref, compressed_ref, force=True)
                success = True
            except (ImportError, Exception) as e:
                logger.warning(f"pysam.tabix_compress failed or unavailable: {e}")
                
            # Fallback to binary
            if not success:
                if bgzip_available:
                    logger.info(f"Falling back to {settings.bgzip_path} for compression")
                    try:
                        with open(compressed_ref, "wb") as f:
                            subprocess.run([settings.bgzip_path, "-c", source_ref], stdout=f, check=True, capture_output=False)
                        success = True
                    except Exception as sub_e:
                        logger.error(f"bgzip fallback failed: {sub_e}")
                else:
                    logger.error("bgzip binary not found and pysam failed.")

            if not success:
                if os.path.exists(compressed_ref):
                     os.remove(compressed_ref)
                raise RuntimeError(f"Failed to compress reference {source_ref}")
        else:
             raise FileNotFoundError(f"Cannot create {compressed_ref}, source {source_ref} not found.")

    # Step 2: Indexing with recovery
    index_fm9 = compressed_ref.replace(".gz", ".fm9") if compressed_ref.endswith(".gz") else compressed_ref + ".fm9"
    fai_index = compressed_ref + ".fai"

    try:
        # Tracy indexing (required for tracy functionality)
        if not os.path.exists(index_fm9):
            if shutil.which(settings.tracy_path):
                subprocess.run([settings.tracy_path, "index", "-o", index_fm9, compressed_ref], check=True, capture_output=True)
            else:
                logger.warning(f"Tracy binary not found at {settings.tracy_path}, skipping .fm9 index")

        # FAI indexing
        if not os.path.exists(fai_index):
            success = False
            # Try pysam
            try:
                import pysam
                logger.info(f"Indexing {compressed_ref} using pysam...")
                pysam.faidx(compressed_ref)
                success = True
            except (ImportError, Exception) as e:
                logger.warning(f"pysam indexing failed or unavailable: {e}")

            # Fallback to samtools
            if not success:
                if samtools_available:
                    logger.info(f"Falling back to {settings.samtools_path} for indexing")
                    subprocess.run([settings.samtools_path, "faidx", compressed_ref], check=True, capture_output=True)
                    success = True
                else:
                    logger.error("samtools binary not found and pysam failed.")

            if not success:
                logger.error(f"Failed to index reference {compressed_ref}. Continuing without index.")
                return compressed_ref

    except subprocess.CalledProcessError as e:
        logger.warning(f"Indexing failed for {compressed_ref}. Attempting to rebuild... Error: {e.stderr.decode() if e.stderr else str(e)}")

        # Cleanup potentially bad files
        for f in [compressed_ref, index_fm9, fai_index, compressed_ref + ".gzi"]:
            if os.path.exists(f):
                os.remove(f)

        # Retry logic (recursive-ish)
        if os.path.exists(source_ref):
             try:
                 return ensure_indexed(source_ref)
             except Exception as retry_e:
                 logger.error(f"Rebuild failed: {retry_e}. Returning original path.")
                 return source_ref
        else:
            logger.error(f"Indexing failed and cannot rebuild {compressed_ref} because source {source_ref} is missing. Returning original path.")
            return ref_path

    return compressed_ref
