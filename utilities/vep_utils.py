import json
import logging
import re
import subprocess
import tempfile
import time
from abc import ABC, abstractmethod
from datetime import datetime
from pathlib import Path
from typing import Any

import httpx

# Configure traceability and logging levels
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

from services.reference import get_lrg_mapping
from core.proxy_manager import proxy_manager

class BaseVEPAnnotator(ABC):
    """
    Abstract base class for VEP annotation engines.
    """
    def __init__(self, assembly: str = "GRCh38"):
        self.assembly = assembly.upper()

    @abstractmethod
    def get_annotations(self, hgvs_variants: list[str]) -> dict[str, Any]:
        """
        Retrieves annotations for a list of HGVS variants.
        """
        pass

    def _rank_impact(self, impact: str) -> int:
        """
        Hierarchical heuristic classifier to prioritize transcriptional damage.
        """
        ranking = {
            "HIGH": 4,       # Truncations, stop codon loss, severe splicing alteration
            "MODERATE": 3,   # Missense mutations, in-frame deletions
            "LOW": 2,        # Synonymous mutations marginally altering translation
            "MODIFIER": 1    # Intronic, UTR variants or intergenic regions
        }
        return ranking.get(impact.upper(), 0)

    def _structure_results(self, raw_results: list[dict[str, Any]]) -> dict[str, Any]:
        """
        Common logic to structure raw VEP results into the application format.
        """
        structured_annotations: dict[str, Any] = {}

        for record in raw_results:
            # Unbreakable link via user's original identifier
            identifier = record.get("original_input", record.get("input"))

            # Taxonomic template for final variant report
            annotation = {
                "gene_symbol": "",
                "consequence": "",
                "impact": "",
                "hgvs_c": "",
                "hgvs_p": "",
                "sift": "",
                "polyphen": "",
                "clin_sig": [],
                "phenotype": [],
                "vep_raw": record,
                "retrieved_at": datetime.now().isoformat()
            }

            transcripts = record.get("transcript_consequences", [])
            most_severe_transcript = None
            highest_impact_score = -1

            # Transcript prioritization algorithm: scan all affected isoforms
            for transcript in transcripts:
                impact_str = transcript.get("impact", "")
                current_score = self._rank_impact(impact_str)

                # Clinical Modifier: Promote MANE Select universal transcripts
                is_mane = 'mane_select' in transcript

                # Elite selection of most severe consequence
                if current_score > highest_impact_score or (current_score == highest_impact_score and is_mane):
                    highest_impact_score = current_score
                    most_severe_transcript = transcript

            if most_severe_transcript:
                # Extraction and aesthetic mapping for human analyst
                annotation["gene_symbol"] = most_severe_transcript.get("gene_symbol", "")
                annotation["impact"] = most_severe_transcript.get("impact", "")

                consequence_terms = most_severe_transcript.get("consequence_terms", [])
                annotation["consequence"] = ", ".join(consequence_terms)

                annotation["hgvs_c"] = most_severe_transcript.get("hgvsc", "")
                annotation["hgvs_p"] = most_severe_transcript.get("hgvsp", "")

                # Consolidation of in-silico pathogenicity matrices
                if "sift_prediction" in most_severe_transcript and "sift_score" in most_severe_transcript:
                    annotation["sift"] = f"{most_severe_transcript['sift_prediction']} ({most_severe_transcript['sift_score']})"

                if "polyphen_prediction" in most_severe_transcript and "polyphen_score" in most_severe_transcript:
                    annotation["polyphen"] = f"{most_severe_transcript['polyphen_prediction']} ({most_severe_transcript['polyphen_score']})"

            # Extract clinical significance and phenotypes from colocated variants
            colocated_variants = record.get("colocated_variants", [])
            clin_sigs = set()
            phenotypes = []

            for cv in colocated_variants:
                if "clin_sig" in cv:
                    for sig in cv["clin_sig"]:
                        clin_sigs.add(sig.replace("_", " "))

                if "phenotype_or_disease" in cv and cv["phenotype_or_disease"] == 1:
                    if "id" in cv and not cv["id"].startswith("COS"):
                        phenotypes.append(cv["id"])
                    
                    if "pubmed" in cv:
                         for pmid in cv["pubmed"]:
                              if len(phenotypes) < 5:
                                  phenotypes.append(f"PMID:{pmid}")

            annotation["clin_sig"] = sorted(list(clin_sigs))
            annotation["phenotype"] = sorted(list(set(phenotypes)))[:10]

            if not identifier or not isinstance(identifier, str):
                continue

            structured_annotations[identifier] = annotation

        return structured_annotations

class OnlineVEPAnnotator(BaseVEPAnnotator):
    """
    Advanced genomic variant annotation engine interacting with Ensembl's REST API.
    """
    def __init__(self, assembly: str = "GRCh38", timeout: float = 30.0):
        super().__init__(assembly)
        if self.assembly == "GRCh37":
            self.base_url = "https://grch37.rest.ensembl.org"
        else:
            self.base_url = "https://rest.ensembl.org"

        self.client = proxy_manager.get_client("ensembl", timeout=timeout)
        self.headers = {
            "Content-Type": "application/json",
            "Accept": "application/json"
        }
        self.refseq_pattern = re.compile(r'^(NM_|NP_|NC_|NG_)', re.IGNORECASE)

    # Removed __del__ as client is managed by ProxyManager

    def _handle_rate_limit(self, response: httpx.Response) -> bool:
        if response.status_code == 429:
            retry_after = int(response.headers.get("Retry-After", 1))
            logging.warning(f"Restricción de tráfico detectada (HTTP 429). Suspendiendo ejecución por {retry_after} segundos...")
            time.sleep(retry_after)
            return True
        return False

    def recode_variants_batch(self, variants: list[str]) -> dict[str, str]:
        if not variants:
            return {}
        endpoint = f"{self.base_url}/variant_recoder/human"
        payload = {"ids": variants, "fields": "hgvsg,id"}
        translation_map = {}
        retries = 3

        for attempt in range(retries):
            try:
                response = self.client.post(endpoint, json=payload, headers=self.headers)
                if self._handle_rate_limit(response): continue
                response.raise_for_status()
                data = response.json()
                for item in data:
                    for allele_info in item.values():
                        if not isinstance(allele_info, dict): continue
                        original_input = allele_info.get("input")
                        if not original_input: continue
                        hgvsg_list = allele_info.get("hgvsg", [])
                        if hgvsg_list and isinstance(hgvsg_list, list) and len(hgvsg_list) > 0:
                            selected_hgvsg = hgvsg_list[0]
                            for h in hgvsg_list:
                                if h.startswith("NC_"):
                                    selected_hgvsg = h
                                    break
                            translation_map[original_input] = selected_hgvsg
                        else:
                            translation_map[original_input] = original_input
                        break
                break
            except httpx.HTTPError as e:
                logging.error(f"Network failure during genomic recoding (Attempt {attempt + 1}): {e}")
                if attempt == retries - 1:
                    translation_map = {v: v for v in variants}
                time.sleep(2)
        return translation_map

    def annotate_hgvs_batch(self, hgvs_variants: list[str]) -> list[dict[str, Any]]:
        if not hgvs_variants: return []
        variants_to_recode = []
        clean_variants = []
        id_map = {}
        for variant in hgvs_variants:
            ac_match = variant.split(':')[0] if ':' in variant else variant
            if ac_match.startswith("NG_"):
                lrg = get_lrg_mapping(ac_match)
                if lrg:
                    new_v = variant.replace(ac_match, lrg)
                    id_map[new_v] = variant
                    variants_to_recode.append(new_v)
                else:
                    variants_to_recode.append(variant)
            elif self.refseq_pattern.match(variant):
                variants_to_recode.append(variant)
            else:
                clean_variants.append(variant)

        recoded_map = {}
        if variants_to_recode:
            logging.info(f"Routing {len(variants_to_recode)} ambiguous/RefSeq sequences via Variant Recoder...")
            recoded_map = self.recode_variants_batch(variants_to_recode)

        final_vep_payload = list(set(clean_variants + list(recoded_map.values())))
        endpoint = f"{self.base_url}/vep/human/hgvs"
        payload = {
            "hgvs_notations": final_vep_payload,
            "ignore_invalid": 1, "refseq": 1, "mane": 1, "sift": "b", "polyphen": "b", "symbol": 1, "transcript_version": 1
        }
        results = []
        retries = 3
        for attempt in range(retries):
            try:
                response = self.client.post(endpoint, json=payload, headers=self.headers)
                if self._handle_rate_limit(response): continue
                if response.status_code == 400 and len(final_vep_payload) > 1:
                    return self._fallback_individual_requests(final_vep_payload, payload)
                response.raise_for_status()
                results = response.json()
                break
            except httpx.HTTPError as e:
                logging.error(f"Direct connection failure with VEP Engine (Attempt {attempt + 1}): {e}")
                if attempt == retries - 1: return []
                time.sleep(2)

        full_reverse_map = {}
        for used_v, hgvsg in recoded_map.items():
            original_v = id_map.get(used_v, used_v)
            full_reverse_map[hgvsg] = original_v

        for record in results:
            input_variant = record.get("input")
            record["original_input"] = full_reverse_map.get(input_variant, input_variant)
        return results

    def _fallback_individual_requests(self, variants: list[str], base_payload: dict) -> list[dict[str, Any]]:
        endpoint = f"{self.base_url}/vep/human/hgvs"
        results = []
        for variant in variants:
            payload = base_payload.copy()
            payload["hgvs_notations"] = [variant]
            try:
                response = self.client.post(endpoint, json=payload, headers=self.headers)
                if self._handle_rate_limit(response):
                    response = self.client.post(endpoint, json=payload, headers=self.headers)
                if response.status_code == 200:
                    results.extend(response.json())
            except Exception as e:
                logging.debug(f"Irresolvable variant isolated ({variant}): {e}")
                continue
        return results

    def get_annotations(self, hgvs_variants: list[str]) -> dict[str, Any]:
        raw_results = self.annotate_hgvs_batch(hgvs_variants)
        return self._structure_results(raw_results)

class LocalVEPAnnotator(BaseVEPAnnotator):
    """
    Annotator that uses a local VEP installation.
    """
    def __init__(self, assembly: str = "GRCh38", vep_path: str = "vep", vep_data: str | None = None):
        super().__init__(assembly)
        self.vep_path = vep_path
        self.vep_data = vep_data

    def get_annotations(self, hgvs_variants: list[str]) -> dict[str, Any]:
        if not hgvs_variants: return {}

        results = []
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as tmp_input:
            for variant in hgvs_variants:
                tmp_input.write(f"{variant}\n")
            tmp_input_path = tmp_input.name

        try:
            cmd = [
                self.vep_path,
                "--input_file", tmp_input_path,
                "--format", "hgvs",
                "--output_file", "STDOUT",
                "--json",
                "--database",
                "--assembly", self.assembly,
                "--everything"
            ]
            if self.vep_data:
                cmd.extend(["--dir", self.vep_data])
            
            logging.info(f"Running local VEP: {' '.join(cmd)}")
            process = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            for line in process.stdout.splitlines():
                if line.strip():
                    try:
                        results.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue
        except subprocess.CalledProcessError as e:
            logging.error(f"Local VEP failed: {e.stderr}")
        finally:
            Path(tmp_input_path).unlink(missing_ok=True)

        return self._structure_results(results)

class DockerVEPAnnotator(BaseVEPAnnotator):
    """
    Annotator that uses VEP via Docker.
    """
    def __init__(self, assembly: str = "GRCh38", image: str = "ensemblorg/ensembl-vep", vep_data: str | None = None):
        super().__init__(assembly)
        self.image = image
        self.vep_data = vep_data

    def get_annotations(self, hgvs_variants: list[str]) -> dict[str, Any]:
        if not hgvs_variants: return {}

        results = []
        # Create a temporary directory to share with Docker
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_dir_path = Path(tmp_dir)
            input_file = tmp_dir_path / "input.txt"
            with open(input_file, 'w') as f:
                for variant in hgvs_variants:
                    f.write(f"{variant}\n")

            cmd = [
                "docker", "run", "--rm",
                "-v", f"{tmp_dir}:/tmp/vep",
            ]
            
            vep_cmd = [
                "vep",
                "--input_file", "/tmp/vep/input.txt",
                "--format", "hgvs",
                "--output_file", "STDOUT",
                "--json",
                "--assembly", self.assembly,
                "--everything"
            ]

            if self.vep_data:
                cmd.extend(["-v", f"{self.vep_data}:/opt/vep/.vep"])
            else:
                vep_cmd.append("--database")

            cmd.append(self.image)
            cmd.extend(vep_cmd)

            logging.info(f"Running docker VEP: {' '.join(cmd)}")
            try:
                process = subprocess.run(cmd, capture_output=True, text=True, check=True)
                for line in process.stdout.splitlines():
                    if line.strip():
                        try:
                            results.append(json.loads(line))
                        except json.JSONDecodeError:
                            continue
            except subprocess.CalledProcessError as e:
                logging.error(f"Docker VEP failed: {e.stderr}")

        return self._structure_results(results)

from core.cache import cache

class VEPAnnotator:
    """
    Dispatcher class for VEP annotation.
    """
    def __init__(self, mode: str = "online", assembly: str = "GRCh38", vep_path: str | None = None, vep_data: str | None = None):
        self.mode = mode
        self.assembly = assembly
        self.engine: BaseVEPAnnotator
        if mode == "local":
            self.engine = LocalVEPAnnotator(assembly, vep_path or "vep", vep_data)
        elif mode == "docker":
            self.engine = DockerVEPAnnotator(assembly, vep_path or "ensemblorg/ensembl-vep", vep_data)
        else:
            self.engine = OnlineVEPAnnotator(assembly)

    def get_annotations(self, hgvs_variants: list[str]) -> dict[str, Any]:
        if not hgvs_variants:
            return {}

        results = {}
        missing = []

        # 1. Check Global Cache (using MGET for efficiency)
        cache_keys = {f"vep:{self.assembly}:{v}": v for v in hgvs_variants}
        cached_data_map = cache.get_many(list(cache_keys.keys()))

        for key, variant in cache_keys.items():
            cached_data = cached_data_map.get(key)
            if cached_data:
                results[variant] = cached_data
            else:
                missing.append(variant)

        if not missing:
            return results

        # 2. Fetch missing from the selected engine
        new_annotations = self.engine.get_annotations(missing)

        # 3. Store new results in cache (using Pipeline for efficiency)
        if new_annotations:
            new_cache_entries = {}
            for variant, data in new_annotations.items():
                cache_key = f"vep:{self.assembly}:{variant}"
                new_cache_entries[cache_key] = data
                results[variant] = data
            
            cache.set_many(new_cache_entries)

        return results
