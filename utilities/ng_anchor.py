# SPDX-License-Identifier: MIT
"""
NG_ -> chromosome coordinate resolution via a learned-once linear anchor.
==========================================================================

A RefSeqGene (``NG_``) record is a *contiguous* genomic sequence, so the mapping
from an ``NG_``-relative position to its chromosome (``NC_``) position is a pure
linear relation::

    plus strand :  chrom_pos = ng_pos + const          (ref/alt unchanged)
    minus strand:  chrom_pos = const - ng_pos          (ref/alt complemented)

The single unknown (``chrom``, ``strand``, ``const``) is *learned once* per
``(accession, assembly)`` from a single network recode of one substitution and
then cached on disk forever. Every subsequent variant on that reference is
resolved locally with zero network calls.

This lets ``_recode_to_vcf`` answer ``NG_`` variants locally before falling back
to the OpenCRAVAT HGVS web API or Ensembl ``variant_recoder``.
"""
from __future__ import annotations

import json
import logging
import os
import re
import threading
from typing import Callable

logger = logging.getLogger(__name__)

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

# NG_ substitution: NG_008866.1:g.65631C>T
_NG_SUB_RE = re.compile(r'^(NG_\d+)(?:\.\d+)?:g\.(\d+)([ACGT]+)>([ACGT]+)$', re.IGNORECASE)
# chrom:pos:ref:alt as emitted by the existing recoders
_VCF_RE = re.compile(r'^([^:]+):(\d+):([ACGTN]*):([ACGTN]*)$', re.IGNORECASE)


def _complement(seq: str) -> str:
    return seq.translate(_COMP)


class NGAnchorResolver:
    """
    Disk-backed resolver of NG_ substitutions to ``chrom:pos:ref:alt`` strings.

    The learner callback (``learn_fn``) is injected so this class has no network
    dependency and is fully unit-testable offline. ``learn_fn(ng_hgvs)`` must
    return a ``chrom:pos:ref:alt`` string (chromosome-oriented) or ``None``.
    """

    def __init__(self, cache_path: str):
        self.cache_path = cache_path
        self._lock = threading.Lock()
        self._mem: dict[str, dict] = {}
        self._load()

    # ----- persistence -----------------------------------------------------
    def _load(self) -> None:
        try:
            if os.path.exists(self.cache_path):
                with open(self.cache_path) as f:
                    self._mem = json.load(f)
        except Exception as e:
            logger.warning(f"Could not load NG anchor cache {self.cache_path}: {e}")
            self._mem = {}

    def _save(self) -> None:
        try:
            os.makedirs(os.path.dirname(self.cache_path), exist_ok=True)
            tmp = f"{self.cache_path}.tmp"
            with open(tmp, "w") as f:
                json.dump(self._mem, f, indent=2)
            os.replace(tmp, self.cache_path)
        except Exception as e:
            logger.warning(f"Could not persist NG anchor cache {self.cache_path}: {e}")

    @staticmethod
    def _key(ng_ac: str, assembly: str) -> str:
        # Strip version so NG_008866 and NG_008866.1 share one anchor.
        base = ng_ac.split(".")[0]
        return f"{base}:{(assembly or 'GRCh38').upper()}"

    # ----- anchor learning -------------------------------------------------
    def _learn(self, ng_ac: str, assembly: str, seed_variant: str,
               learn_fn: Callable[[str], str | None]) -> dict | None:
        """Derive (chrom, strand, const) from one recoded substitution."""
        m = _NG_SUB_RE.match(seed_variant)
        if not m:
            return None
        ng_pos, ng_ref, ng_alt = int(m.group(2)), m.group(3).upper(), m.group(4).upper()
        if len(ng_ref) != 1 or len(ng_alt) != 1:
            return None  # need a clean SNV to read strand unambiguously

        vcf = learn_fn(seed_variant)
        if not vcf:
            return None
        vm = _VCF_RE.match(vcf)
        if not vm:
            return None
        chrom, pos, ref, alt = vm.group(1), int(vm.group(2)), vm.group(3).upper(), vm.group(4).upper()

        if (ref, alt) == (ng_ref, ng_alt):
            anchor = {"chrom": chrom, "strand": "+", "const": pos - ng_pos}
        elif (ref, alt) == (_complement(ng_ref), _complement(ng_alt)):
            anchor = {"chrom": chrom, "strand": "-", "const": pos + ng_pos}
        else:
            logger.warning(
                f"NG anchor learn: ref/alt mismatch for {seed_variant} -> {vcf}; not caching")
            return None

        with self._lock:
            self._mem[self._key(ng_ac, assembly)] = anchor
            self._save()
        logger.info(f"Learned NG anchor {ng_ac} ({assembly}): {anchor}")
        return anchor

    def get_anchor(self, ng_ac: str, assembly: str) -> dict | None:
        return self._mem.get(self._key(ng_ac, assembly))

    # ----- resolution ------------------------------------------------------
    def _resolve_one(self, ng_hgvs: str, anchor: dict) -> str | None:
        m = _NG_SUB_RE.match(ng_hgvs)
        if not m:
            return None  # only substitutions handled locally; rest -> network
        ng_pos, ng_ref, ng_alt = int(m.group(2)), m.group(3).upper(), m.group(4).upper()
        if len(ng_ref) != 1 or len(ng_alt) != 1:
            return None
        if anchor["strand"] == "+":
            pos = ng_pos + anchor["const"]
            ref, alt = ng_ref, ng_alt
        else:
            pos = anchor["const"] - ng_pos
            ref, alt = _complement(ng_ref), _complement(ng_alt)
        if pos <= 0:
            return None
        return f"{anchor['chrom']}:{pos}:{ref}:{alt}"

    def resolve_batch(self, ng_variants: list[str], assembly: str,
                      learn_fn: Callable[[str], str | None]) -> dict[str, str]:
        """
        Resolve as many NG_ substitutions as possible locally.

        Returns ``{ng_hgvs: 'chrom:pos:ref:alt'}`` only for the variants it could
        resolve. Variants it cannot resolve are omitted so the caller falls back
        to the existing network paths. Performs at most one network call (the
        anchor) per (accession, assembly) and zero once cached.
        """
        out: dict[str, str] = {}

        # Group variants by NG accession.
        by_ac: dict[str, list[str]] = {}
        for v in ng_variants:
            m = _NG_SUB_RE.match(v)
            if m:
                by_ac.setdefault(m.group(1), []).append(v)

        for ng_ac, variants in by_ac.items():
            anchor = self.get_anchor(ng_ac, assembly)
            if anchor is None:
                seed = variants[0]  # any clean SNV on this reference
                anchor = self._learn(ng_ac, assembly, seed, learn_fn)
                if anchor is None:
                    continue  # leave the whole group to the network fallback
            for v in variants:
                vcf = self._resolve_one(v, anchor)
                if vcf:
                    out[v] = vcf
        return out
