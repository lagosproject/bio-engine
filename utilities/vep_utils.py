# SPDX-License-Identifier: MIT
import json
import logging
import os
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

from core.proxy_manager import proxy_manager
from services.reference import get_lrg_mapping


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
    def __init__(self, assembly: str = "GRCh38", timeout: float = 120.0):
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

class OpenCRAVATAnnotator(BaseVEPAnnotator):
    """
    Annotator using OpenCRAVAT CLI for local variant annotation.

    HGVS variants are resolved to genomic coordinates via Ensembl variant_recoder
    (lightweight coordinate lookup only), then annotated locally by all installed
    OC modules. Results include core schema fields plus an `oc_data` dict containing
    every field from every installed annotator — new databases are picked up
    automatically without code changes.

    Requires: pip install open-cravat && oc module install-base
    """

    _OC_IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}

    def __init__(self, assembly: str = "GRCh38", oc_path: str = "oc"):
        super().__init__(assembly)
        self.oc_path = oc_path
        self.genome = "hg38" if assembly.upper() == "GRCH38" else "hg19"
        self.ensembl_url = (
            "https://grch37.rest.ensembl.org" if assembly.upper() == "GRCH37"
            else "https://rest.ensembl.org"
        )
        self.client = proxy_manager.get_client("ensembl", timeout=60.0)
        self.headers = {"Content-Type": "application/json", "Accept": "application/json"}
        self._refseq_pattern = re.compile(r'^(NM_|NP_|NC_|NG_)', re.IGNORECASE)

        # Disk-backed NG_ -> chromosome anchor cache (learned once per reference).
        from core.config import settings as _settings
        from utilities.ng_anchor import NGAnchorResolver
        self._ng_anchor = NGAnchorResolver(
            os.path.join(_settings.cache_dir, "ng_anchors.json")
        )
        # Lazy fully-offline mappers (None = not built yet, False = unavailable).
        self._ng_mapper = None   # NG_ contiguous genomic (RefSeqGene alignments)
        self._tx_mapper = None   # NM_/NR_ spliced transcripts (UCSC ncbiRefSeq)

    def _get_ng_mapper(self):
        """Lazily build the offline NG_ mapper; disable permanently on failure."""
        if self._ng_mapper is False:
            return None
        if self._ng_mapper is None:
            try:
                from utilities.local_coordinate_mapper import LocalCoordinateMapper
                self._ng_mapper = LocalCoordinateMapper(self.assembly)
            except Exception as e:
                logging.warning(f"Offline NG mapper unavailable: {e}")
                self._ng_mapper = False
                return None
        return self._ng_mapper

    def _ng_resolve_offline(self, ng_hgvs: str) -> str | None:
        """Resolve one NG_ variant via the offline mapper; disable it on error."""
        mapper = self._get_ng_mapper()
        if not mapper:
            return None
        try:
            return mapper.resolve_ng_variant(ng_hgvs)
        except Exception as e:
            logging.warning(f"Offline NG mapper failed, disabling: {e}")
            self._ng_mapper = False
            return None

    def _ng_learn(self, seed: str) -> str | None:
        """Anchor learner: offline mapper first, Ensembl recode as last resort."""
        return self._ng_resolve_offline(seed) or self._recode_to_vcf_via_ensembl([seed]).get(seed)

    def _get_tx_mapper(self):
        """Lazily build the offline transcript mapper; disable permanently on failure."""
        if self._tx_mapper is False:
            return None
        if self._tx_mapper is None:
            try:
                from utilities.transcript_coordinate_mapper import TranscriptCoordinateMapper
                self._tx_mapper = TranscriptCoordinateMapper(self.assembly)
            except Exception as e:
                logging.warning(f"Offline transcript mapper unavailable: {e}")
                self._tx_mapper = False
                return None
        return self._tx_mapper

    def _tx_resolve_offline(self, hgvs: str) -> str | None:
        """Resolve one transcript variant via the offline mapper; disable on error."""
        mapper = self._get_tx_mapper()
        if not mapper:
            return None
        try:
            return mapper.resolve_transcript_variant(hgvs)
        except Exception as e:
            logging.warning(f"Offline transcript mapper failed, disabling: {e}")
            self._tx_mapper = False
            return None

    def _nc_to_chrom(self, nc_ac: str) -> str | None:
        base = nc_ac.split(".")[0]
        if not base.startswith("NC_"):
            return None
        try:
            num = int(base.split("_")[1])
            if 1 <= num <= 22:
                return str(num)
            elif num == 23:
                return "X"
            elif num == 24:
                return "Y"
            elif num == 12920 or num == 1807:
                return "MT"
        except (IndexError, ValueError):
            pass
        return None

    def _resolve_versioned_hgvs(self, variant: str) -> str:
        if ":" not in variant:
            return variant
        tx, rest = variant.split(":", 1)
        if "." in tx or not (tx.startswith("NM_") or tx.startswith("XM_") or tx.startswith("NR_") or tx.startswith("XR_") or tx.startswith("NG_")):
            return variant

        # Try to find version in cached gb file
        import os

        from core.config import settings
        cache_dirs = [
            settings.get_cache_dir(assembly=self.assembly),
            os.path.join(settings.cache_dir, "ncbi_cache"),
            settings.cache_dir
        ]
        for cdir in cache_dirs:
            gb_path = os.path.join(cdir, f"{tx}.gb")
            if os.path.exists(gb_path):
                try:
                    with open(gb_path) as f:
                        for line in f:
                            if line.startswith("VERSION"):
                                parts = line.split()
                                if len(parts) >= 2 and parts[1].startswith(tx):
                                    return f"{parts[1]}:{rest}"
                                break
                except Exception:
                    pass
        return variant

    def _recode_to_vcf_via_opencravat(self, hgvs_variants: list[str]) -> dict[str, str]:
        """Runs OpenCRAVAT converter locally on HGVS variants to resolve genomic VCF coordinates."""
        if not hgvs_variants:
            return {}

        result = {}
        with tempfile.TemporaryDirectory() as tmp_dir:
            # Name the file with .hgvs so that OpenCRAVAT routes it to hgvs-converter
            input_path = Path(tmp_dir) / "variants.hgvs"
            input_path.write_text("\n".join(hgvs_variants))

            cmd = [
                self.oc_path, "run", str(input_path),
                "-l", self.genome,
                "-d", tmp_dir,
                "--skip", "annotator", "postaggregator",
            ]

            logging.info(f"Running OpenCRAVAT for VCF recoding: {' '.join(cmd)}")
            try:
                # 300-second timeout is plenty for coordinate conversion on larger batches
                process = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
                if process.returncode != 0:
                    logging.warning(f"OpenCRAVAT VCF recoding failed with code {process.returncode}: {process.stderr[:200]}")
                    return {}
            except Exception as e:
                logging.error(f"OpenCRAVAT VCF recoding timed out or failed to start: {e}")
                return {}

            sqlite_files = sorted(Path(tmp_dir).glob("*.sqlite"))
            if not sqlite_files:
                return {}

            db_path = sqlite_files[0]
            import sqlite3
            conn = sqlite3.connect(db_path)
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()
            try:
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variant'")
                if not cursor.fetchone():
                    return {}

                cursor.execute("SELECT base__uid, base__chrom, base__pos, base__ref_base, base__alt_base FROM variant")
                rows = cursor.fetchall()
                for row in rows:
                    uid = row["base__uid"]
                    chrom = row["base__chrom"]
                    pos = row["base__pos"]
                    ref = row["base__ref_base"]
                    alt = row["base__alt_base"]

                    if uid is not None and chrom and pos:
                        idx = int(uid) - 1
                        if 0 <= idx < len(hgvs_variants):
                            original_hgvs = hgvs_variants[idx]
                            if chrom.startswith("chr"):
                                chrom = chrom[3:]
                            vcf_string = f"{chrom}:{pos}:{ref}:{alt}"
                            result[original_hgvs] = vcf_string
                            # Cache the result
                            cache.set(f"vcf_recode:{self.assembly}:{original_hgvs}", vcf_string)
            except Exception as e:
                logging.error(f"Failed to parse OpenCRAVAT SQLite VCF coordinates: {e}")
            finally:
                conn.close()

        return result

    def _recode_to_vcf_via_ensembl(self, hgvs_variants: list[str]) -> dict[str, str]:
        """Queries Ensembl REST API variant_recoder to resolve genomic VCF coordinates."""
        if not hgvs_variants:
            return {}

        result = {}
        variants_to_send = []
        id_map = {}

        for variant in hgvs_variants:
            ac = variant.split(':')[0] if ':' in variant else variant
            if ac.startswith("NG_"):
                lrg = get_lrg_mapping(ac)
                mapped = variant.replace(ac, lrg) if lrg else variant
                id_map[mapped] = variant
                variants_to_send.append(mapped)
            else:
                id_map[variant] = variant
                variants_to_send.append(variant)

        endpoint = f"{self.ensembl_url}/variant_recoder/human"
        chunk_size = 100
        chunks = [variants_to_send[i:i + chunk_size] for i in range(0, len(variants_to_send), chunk_size)]

        def recode_chunk(chunk):
            payload = {"ids": chunk}
            retries = 3
            chunk_results = {}
            for attempt in range(retries):
                try:
                    response = self.client.post(endpoint, json=payload, headers=self.headers)
                    if response.status_code == 429:
                        time.sleep(int(response.headers.get("Retry-After", 2)))
                        continue
                    response.raise_for_status()
                    for item in response.json():
                        for allele_info in item.values():
                            if not isinstance(allele_info, dict):
                                continue
                            sent_id = allele_info.get("input")
                            if not sent_id:
                                continue
                            original = id_map.get(sent_id, sent_id)

                            vcf_string = None
                            vcf_strings = allele_info.get("vcf_string", [])
                            if vcf_strings and isinstance(vcf_strings, list):
                                vcf_string = vcf_strings[0]

                            if not vcf_string:
                                spdi_list = allele_info.get("spdi", [])
                                if spdi_list and isinstance(spdi_list, list):
                                    for s in spdi_list:
                                        if s.startswith("NC_"):
                                            parts = s.split(":")
                                            if len(parts) == 4:
                                                nc_ac, pos0, ref, alt = parts
                                                chrom = self._nc_to_chrom(nc_ac)
                                                if chrom:
                                                    try:
                                                        pos1 = int(pos0) + 1
                                                        vcf_string = f"{chrom}:{pos1}:{ref}:{alt}"
                                                        break
                                                    except ValueError:
                                                        pass

                            if not vcf_string:
                                hgvsg_list = allele_info.get("hgvsg", [])
                                if hgvsg_list and isinstance(hgvsg_list, list):
                                    for h in hgvsg_list:
                                        if h.startswith("NC_"):
                                            match = re.match(r'^(NC_\d+)(?:\.\d+)?:g\.(\d+)([ATGC]*)(?:>([ATGC]*))?$', h, re.IGNORECASE)
                                            if match:
                                                nc_ac, pos, ref, alt = match.groups()
                                                ref = ref or ""
                                                alt = alt or ""
                                                chrom = self._nc_to_chrom(nc_ac)
                                                if chrom:
                                                    vcf_string = f"{chrom}:{pos}:{ref}:{alt}"
                                                    break

                            if vcf_string:
                                chunk_results[original] = vcf_string
                            break
                    return chunk_results
                except httpx.HTTPError as e:
                    logging.error(f"VCF recoding attempt {attempt + 1} failed for chunk {chunk}: {e}")
                    if attempt < retries - 1:
                        time.sleep(2)
            return chunk_results

        from concurrent.futures import ThreadPoolExecutor, as_completed
        new_results = {}
        with ThreadPoolExecutor(max_workers=8) as executor:
            future_to_chunk = {executor.submit(recode_chunk, chunk): chunk for chunk in chunks}
            for future in as_completed(future_to_chunk):
                try:
                    chunk_res = future.result()
                    new_results.update(chunk_res)
                except Exception as e:
                    logging.error(f"Failed to process chunk: {e}")

        if new_results:
            new_cache_entries = {f"vcf_recode:{self.assembly}:{variant}": val for variant, val in new_results.items()}
            cache.set_many(new_cache_entries)
            result.update(new_results)

        return result

    def _recode_to_vcf(self, hgvs_variants: list[str]) -> dict[str, str]:
        """
        Returns {original_hgvs: vcf_string} where vcf_string is 'chrom:pos:ref:alt'.
        Tries OpenCRAVAT converter locally first, falling back to Ensembl variant_recoder.
        """
        if not hgvs_variants:
            return {}

        # Filter out invalid variants with alternate base 'N' or 'n' to avoid API hangs/errors
        valid_variants = []
        for v in hgvs_variants:
            if v.strip().upper().endswith(">N"):
                logging.info(f"Skipping variant with unknown alternate base 'N': {v}")
                continue
            valid_variants.append(v)
        hgvs_variants = valid_variants

        if not hgvs_variants:
            return {}

        result: dict[str, str] = {}

        # 1. Check Redis cache first
        from core.cache import cache
        cache_keys = {f"vcf_recode:{self.assembly}:{v}": v for v in hgvs_variants}
        cached_data_map = cache.get_many(list(cache_keys.keys()))

        variants_to_process = []
        for key, variant in cache_keys.items():
            cached_val = cached_data_map.get(key)
            if cached_val:
                result[variant] = cached_val
            else:
                variants_to_process.append(variant)

        if not variants_to_process:
            return result

        # 1.5. Local parsing for genomic HGVS (NC_) to completely bypass Ensembl API
        for variant in list(variants_to_process):
            if variant.startswith("NC_"):
                match = re.match(r'^(NC_\d+)(?:\.\d+)?:g\.(\d+)([ATGC]*)(?:>([ATGC]*))?$', variant, re.IGNORECASE)
                if match:
                    nc_ac, pos, ref, alt = match.groups()
                    ref = ref or ""
                    alt = alt or ""
                    chrom = self._nc_to_chrom(nc_ac)
                    if chrom:
                        vcf_string = f"{chrom}:{pos}:{ref}:{alt}"
                        result[variant] = vcf_string
                        # Save in cache
                        cache.set(f"vcf_recode:{self.assembly}:{variant}", vcf_string)
                        variants_to_process.remove(variant)

        if not variants_to_process:
            return result

        # 1.6. Local NG_ -> chromosome resolution.
        # NG_ is a contiguous genomic record. We combine two local mechanisms and
        # only touch Ensembl as a last resort:
        #   - substitutions: a learned-once linear anchor (pure arithmetic, cached
        #     on disk forever) -> after the first variant of a reference, zero I/O;
        #   - the anchor is *learned* from the offline RefSeqGene mapper (no Ensembl)
        #     and falls back to a single Ensembl recode only if that mapper is down;
        #   - other variant types (indels) go straight to the offline mapper.
        # Ensembl is hit at most once per reference, and never once cached.
        ng_variants = [v for v in variants_to_process if v.startswith("NG_")]
        if ng_variants:
            try:
                ng_mapped = self._ng_anchor.resolve_batch(
                    ng_variants, self.assembly, learn_fn=self._ng_learn,
                )
                # Remaining NG_ variants (indels, etc.): offline mapper, no network.
                for v in ng_variants:
                    if v not in ng_mapped:
                        coord = self._ng_resolve_offline(v)
                        if coord:
                            ng_mapped[v] = coord
                for v, coord in ng_mapped.items():
                    result[v] = coord
                    cache.set(f"vcf_recode:{self.assembly}:{v}", coord)
                variants_to_process = [v for v in variants_to_process if v not in result]
            except Exception as e:
                logging.error(f"Local NG VCF recoding failed: {e}")

        if not variants_to_process:
            return result

        # 1.7. Local transcript (NM_/NR_/XM_/XR_) -> chromosome resolution via the
        # offline UCSC ncbiRefSeq exon table (spliced, strand- and CDS-aware).
        # Substitutions only; everything else is left to the network fallback.
        tx_variants = [v for v in variants_to_process
                       if re.match(r'^(NM_|NR_|XM_|XR_)', v, re.IGNORECASE)]
        if tx_variants:
            try:
                for v in tx_variants:
                    coord = self._tx_resolve_offline(v)
                    if coord:
                        result[v] = coord
                        cache.set(f"vcf_recode:{self.assembly}:{v}", coord)
                variants_to_process = [v for v in variants_to_process if v not in result]
            except Exception as e:
                logging.error(f"Local transcript VCF recoding failed: {e}")

        if not variants_to_process:
            return result

        # 2. Try OpenCRAVAT locally first
        try:
            versioned_to_original = {self._resolve_versioned_hgvs(v): v for v in variants_to_process}
            oc_mapped = self._recode_to_vcf_via_opencravat(list(versioned_to_original.keys()))

            # Map back to original HGVS and save in result & cache
            for v_ver, coord in oc_mapped.items():
                v_orig = versioned_to_original.get(v_ver)
                if v_orig:
                    result[v_orig] = coord
                    cache.set(f"vcf_recode:{self.assembly}:{v_orig}", coord)
                    if v_ver != v_orig:
                        cache.set(f"vcf_recode:{self.assembly}:{v_ver}", coord)

            # Remove mapped variants from process queue
            variants_to_process = [v for v in variants_to_process if v not in result]
        except Exception as e:
            logging.error(f"OpenCRAVAT local VCF recoding failed: {e}")

        # 3. Fallback to Ensembl REST API for anything still missing
        if variants_to_process:
            logging.warning(f"Falling back to Ensembl VCF recoding for {len(variants_to_process)} variants...")
            try:
                ensembl_mapped = self._recode_to_vcf_via_ensembl(variants_to_process)
                result.update(ensembl_mapped)
            except Exception as e:
                logging.error(f"Ensembl fallback VCF recoding failed: {e}")

        return result

    def _run_opencravat(self, input_lines: list[str]) -> list[dict[str, Any]]:
        """Write OC input file, run oc run, return parsed variant records from SQLite database."""
        with tempfile.TemporaryDirectory() as tmp_dir:
            input_path = Path(tmp_dir) / "variants.txt"
            input_path.write_text("\n".join(input_lines))

            cmd = [
                self.oc_path, "run", str(input_path),
                "-l", self.genome,
                "-d", tmp_dir,
            ]

            # Get installed annotators to ensure they are all run (e.g. clinvar)
            try:
                from services.opencravat import list_installed
                installed = list_installed(oc_path=self.oc_path)
                annotators = [m["name"] for m in installed if m.get("type") == "annotator"]
                if annotators:
                    cmd.extend(["-a"] + annotators)
            except Exception as e:
                logging.warning(f"Could not list installed OpenCRAVAT annotators: {e}")

            logging.info(f"Running OpenCRAVAT: {' '.join(cmd)}")

            try:
                process = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=300
                )
                if process.returncode != 0:
                    logging.error(f"OpenCRAVAT error: {process.stderr[:500]}")
                    return []
            except subprocess.TimeoutExpired:
                logging.error("OpenCRAVAT timed out after 5 minutes")
                return []
            except FileNotFoundError:
                logging.error(
                    f"OpenCRAVAT binary not found at '{self.oc_path}'. "
                    "Install with: pip install open-cravat && oc module install-base"
                )
                return []

            # Locate the generated SQLite database
            sqlite_files = sorted(Path(tmp_dir).glob("*.sqlite"))
            if not sqlite_files:
                logging.error("No SQLite database produced by OpenCRAVAT")
                return []

            db_path = sqlite_files[0]
            import sqlite3
            conn = sqlite3.connect(db_path)
            conn.row_factory = sqlite3.Row
            cursor = conn.cursor()
            try:
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variant'")
                if not cursor.fetchone():
                    logging.error("No 'variant' table found in OpenCRAVAT SQLite database")
                    return []

                cursor.execute("SELECT * FROM variant")
                rows = cursor.fetchall()
                if not rows:
                    return []

                colnames = [description[0] for description in cursor.description]
                records = []
                for row in rows:
                    rec = {}
                    for colname in colnames:
                        val = row[colname]
                        rec[colname] = val
                        if "__" in colname:
                            short = colname.split("__")[1]
                            rec[short] = val
                    records.append(rec)
                return records
            except Exception as e:
                logging.error(f"Failed to read variants from OpenCRAVAT SQLite database: {e}")
                return []
            finally:
                conn.close()

    def _map_oc_record(self, record: dict[str, Any]) -> dict[str, Any]:
        """
        Map an OC record to the standard annotation schema.
        All raw OC fields (from every installed annotator) are preserved in `oc_data`.
        """
        def first(*keys: str) -> str:
            for k in keys:
                v = record.get(k)
                if v is not None and v != "":
                    return str(v)
            return ""

        annotation: dict[str, Any] = {
            "gene_symbol": first("hugo"),
            "consequence": first("so"),
            "impact": first("impact"),
            "hgvs_c": first("hgvs_coding", "coding"),
            "hgvs_p": first("hgvs_protein", "achange"),
            "sift": "",
            "polyphen": "",
            "clin_sig": [],
            "phenotype": [],
            # All fields from all installed OC annotators — expands automatically
            # when new databases are added via `oc install <module>`
            "oc_data": {k: v for k, v in record.items() if v is not None and v != ""},
            "retrieved_at": datetime.now().isoformat(),
        }

        sift_pred = record.get("sift__prediction") or record.get("sift__pred")
        sift_score = record.get("sift__score")
        if sift_pred:
            annotation["sift"] = f"{sift_pred} ({sift_score})" if sift_score is not None else sift_pred

        pp_pred = record.get("polyphen2__hdiv_pred") or record.get("polyphen2__hvar_pred")
        pp_score = record.get("polyphen2__hdiv_score") or record.get("polyphen2__hvar_score")
        if pp_pred:
            annotation["polyphen"] = f"{pp_pred} ({pp_score})" if pp_score is not None else pp_pred

        clin_sig = record.get("clinvar__sig")
        if clin_sig:
            sigs = clin_sig if isinstance(clin_sig, list) else [clin_sig]
            annotation["clin_sig"] = sorted({str(s).replace("_", " ") for s in sigs})

        diseases = record.get("clinvar__disease_names") or record.get("omim__disease_names")
        if diseases:
            annotation["phenotype"] = (diseases if isinstance(diseases, list) else [diseases])[:10]

        return annotation

    def get_annotations(self, hgvs_variants: list[str]) -> dict[str, Any]:
        if not hgvs_variants:
            return {}

        # Step 1: resolve HGVS → VCF coordinates (Ensembl, coordinate lookup only)
        vcf_map = self._recode_to_vcf(hgvs_variants)
        if not vcf_map:
            logging.error("OpenCRAVAT: failed to resolve any HGVS to VCF coordinates")
            return {}

        # Step 2: build OC input lines; track uid (1-based line number) → original HGVS
        uid_to_hgvs: dict[int, str] = {}
        input_lines: list[str] = []

        for hgvs in hgvs_variants:
            vcf_string = vcf_map.get(hgvs)
            if not vcf_string:
                continue
            parts = vcf_string.split(":")
            if len(parts) != 4:
                continue
            chrom, pos, ref, alt = parts
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
            uid_to_hgvs[len(input_lines) + 1] = hgvs
            input_lines.append(f"{chrom}\t{pos}\t+\t{ref}\t{alt}")

        if not input_lines:
            return {}

        # Step 3: run OC locally
        oc_records = self._run_opencravat(input_lines)

        # Step 4: group by uid — OC may emit multiple transcript rows per variant
        uid_groups: dict[int, list[dict[str, Any]]] = {}
        for rec in oc_records:
            uid = rec.get("uid")
            if uid is not None:
                uid_groups.setdefault(int(uid), []).append(rec)

        # Step 5: pick most severe transcript, map to schema
        results: dict[str, Any] = {}
        for uid, records in uid_groups.items():
            hgvs = uid_to_hgvs.get(uid)
            if not hgvs:
                continue
            best = max(
                records,
                key=lambda r: self._OC_IMPACT_RANK.get(str(r.get("impact", "")).upper(), 0)
            )
            results[hgvs] = self._map_oc_record(best)

        return results

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
    def __init__(self, mode: str = "opencravat", assembly: str = "GRCh38", vep_path: str | None = None, vep_data: str | None = None):
        self.mode = mode
        self.assembly = assembly
        self.vep_path = vep_path

        # Always initialize both engines
        self.oc_engine = OpenCRAVATAnnotator(assembly, vep_path or "oc")
        self.vep_engine = OnlineVEPAnnotator(assembly)

    def get_annotations(self, hgvs_variants: list[str]) -> dict[str, Any]:
        if not hgvs_variants:
            return {}

        results = {}

        # 1. Check cache for existing annotations
        oc_cache_keys = {f"vep:opencravat:{self.assembly}:{v}": v for v in hgvs_variants}
        vep_cache_keys = {}
        if self.mode == "online":
            vep_cache_keys = {f"vep:online:{self.assembly}:{v}": v for v in hgvs_variants}

        all_keys = list(oc_cache_keys.keys()) + list(vep_cache_keys.keys())
        cached_data = cache.get_many(all_keys) if all_keys else {}

        oc_missing = []
        vep_missing = []

        # Determine what is missing from each cache
        for variant in hgvs_variants:
            oc_key = f"vep:opencravat:{self.assembly}:{variant}"
            if oc_key not in cached_data:
                oc_missing.append(variant)

            if self.mode == "online":
                vep_key = f"vep:online:{self.assembly}:{variant}"
                if vep_key not in cached_data:
                    vep_missing.append(variant)

        # 2. Fetch missing OpenCRAVAT annotations (always run)
        oc_new = {}
        if oc_missing:
            try:
                oc_new = self.oc_engine.get_annotations(oc_missing)
                # Store in cache
                if oc_new:
                    oc_cache_entries = {f"vep:opencravat:{self.assembly}:{v}": data for v, data in oc_new.items()}
                    cache.set_many(oc_cache_entries)
            except Exception as e:
                logging.error(f"OpenCRAVAT annotation failed: {e}")

        # 3. Fetch missing VEP annotations (only if self.mode == "online")
        vep_new = {}
        if self.mode == "online" and vep_missing:
            try:
                vep_new = self.vep_engine.get_annotations(vep_missing)
                # Store in cache
                if vep_new:
                    vep_cache_entries = {f"vep:online:{self.assembly}:{v}": data for v, data in vep_new.items()}
                    cache.set_many(vep_cache_entries)
            except Exception as e:
                logging.error(f"Ensembl VEP annotation failed: {e}")

        # 4. Merge results for each variant
        for variant in hgvs_variants:
            # Get OpenCRAVAT data
            oc_key = f"vep:opencravat:{self.assembly}:{variant}"
            oc_data = oc_new.get(variant) or cached_data.get(oc_key) or {}

            # Get VEP data
            vep_data = {}
            if self.mode == "online":
                vep_key = f"vep:online:{self.assembly}:{variant}"
                vep_data = vep_new.get(variant) or cached_data.get(vep_key) or {}

            if oc_data and vep_data:
                # Merge VEP and OpenCRAVAT data
                merged = {**oc_data}
                for k, v in vep_data.items():
                    if k == "vep_raw":
                        merged["vep_raw"] = v
                    elif k in ("clin_sig", "phenotype"):
                        merged[k] = sorted(list(set(merged.get(k, []) + (v or []))))
                    elif v:
                        merged[k] = v
                results[variant] = merged
            elif oc_data:
                results[variant] = oc_data
            elif vep_data:
                results[variant] = vep_data

        return results
