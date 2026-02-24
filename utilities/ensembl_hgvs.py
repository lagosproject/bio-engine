import logging
import time
import httpx
from typing import List, Dict, Any
from services.reference import get_lrg_mapping

logger = logging.getLogger(__name__)

class EnsemblHGVS:
    """
    HGVS Annotator using Ensembl's REST API variant_recoder.
    Optimized for batch operations.
    """
    @staticmethod
    def format_hgvs(ac: str, hgvs_type: str, pos: int, ref: str, alt: str) -> str:
        """
        Formats a primary HGVS string.
        Example: NC_000007.14:g.55181378G>A
        """
        # Basic SNP/Indel formatting for primary ID
        # For simplicity, we use the standard > for SNPs and delins for others if needed
        # but Ensembl recoder is very flexible with input.
        if len(ref) == 1 and len(alt) == 1:
            return f"{ac}:{hgvs_type}.{pos}{ref}>{alt}"
        elif not ref: # Insertion
            return f"{ac}:{hgvs_type}.{pos}_{pos+1}ins{alt}"
        elif not alt: # Deletion
            return f"{ac}:{hgvs_type}.{pos}_{pos+len(ref)-1}del"
        else: # delins
            return f"{ac}:{hgvs_type}.{pos}_{pos+len(ref)-1}delins{alt}"

    def __init__(self, assembly: str = "GRCh38", timeout: float = 120.0):
        self.assembly = assembly.upper()
        if self.assembly == "GRCh37":
            self.base_url = "https://grch37.rest.ensembl.org"
        else:
            self.base_url = "https://rest.ensembl.org"
            
        # We increase the timeout because batch recoding can be heavy on Ensembl's side
        self.client = httpx.Client(timeout=timeout)
        self.headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }

    def __del__(self):
        try:
            self.client.close()
        except:
            pass

    def get_equivalents_batch(self, hgvs_variants: List[str], chunk_size: int = 5) -> Dict[str, List[str]]:
        """
        Main entry point for batch HGVS lookup. handles NG_ -> LRG_ mapping.
        """
        if not hgvs_variants:
            return {}

        # Pre-mapping: Detect NG_ and swap for LRG_ if available
        # This is a surgical fix because Ensembl rejects NG_ with 400 errors.
        id_map = {} # original -> used
        variants_to_lookup = []

        for v in hgvs_variants:
            ac_match = v.split(':')[0] if ':' in v else v
            if ac_match.startswith("NG_"):
                lrg = get_lrg_mapping(ac_match)
                if lrg:
                    new_v = v.replace(ac_match, lrg)
                    id_map[new_v] = v
                    variants_to_lookup.append(new_v)
                else:
                    variants_to_lookup.append(v)
            else:
                variants_to_lookup.append(v)

        final_results = {}
        for i in range(0, len(variants_to_lookup), chunk_size):
            chunk = variants_to_lookup[i:i + chunk_size]
            chunk_results = self._get_chunk_results(chunk)
            
            # Post-mapping: Restore original NG_ identifiers
            remapped_chunk_results = {}
            for lookup_v, alternatives in chunk_results.items():
                original_v = id_map.get(lookup_v, lookup_v)
                remapped_chunk_results[original_v] = alternatives

            # Ensure every input variant from the ORIGINAL list has at least itself
            current_chunk_originals = hgvs_variants[i:i + chunk_size]
            for v in current_chunk_originals:
                if v not in remapped_chunk_results or not remapped_chunk_results[v]:
                    remapped_chunk_results[v] = [v]

            final_results.update(remapped_chunk_results)
        
        return final_results

    def _get_chunk_results(self, chunk: List[str]) -> Dict[str, List[str]]:
        endpoint = f"{self.base_url}/variant_recoder/human"
        payload = {
            "ids": chunk,
            "fields": "hgvsg,hgvsc,hgvsp"
        }

        results_map = {}
        try:
            response = self.client.post(endpoint, json=payload, headers=self.headers)
            
            if response.status_code == 429:
                retry_after = int(response.headers.get("Retry-After", 1))
                time.sleep(retry_after)
                response = self.client.post(endpoint, json=payload, headers=self.headers)

            response.raise_for_status()
            data = response.json()

            for item in data:
                # Each item in the response list corresponds to one of our requested IDs
                for allele, info in item.items():
                    if not isinstance(info, dict):
                        continue
                    
                    input_id = info.get("input")
                    if not input_id:
                        continue
                    
                    equivalents = set()
                    # Always include the input itself if it's a valid notation
                    equivalents.add(input_id)

                    for field in ["hgvsg", "hgvsc", "hgvsp"]:
                        vals = info.get(field, [])
                        if isinstance(vals, list):
                            equivalents.update(vals)
                    
                    if input_id not in results_map:
                        results_map[input_id] = set()
                    results_map[input_id].update(equivalents)

            return {k: sorted(list(v)) for k, v in results_map.items()}
                    
        except Exception as e:
            logger.error(f"Ensembl chunk lookup failed for {chunk}: {e}")
            return {v: [v] for v in chunk}
