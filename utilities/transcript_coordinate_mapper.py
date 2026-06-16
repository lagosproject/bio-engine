# SPDX-License-Identifier: MIT
"""
Offline transcript (NM_/NR_/XM_/XR_) -> chromosome coordinate resolution.
=========================================================================

A spliced transcript maps to the genome exon-by-exon, so converting a coding
(``c.``) or non-coding (``n.``) position to a chromosome position needs the
transcript's exon genomic coordinates, its CDS boundaries and its strand.

This is the spliced counterpart of :mod:`utilities.local_coordinate_mapper`
(which handles the *contiguous* ``NG_`` case with a simple linear offset). The
exon structure is sourced once from the UCSC ``ncbiRefSeq`` table (~7 MB), built
into a small SQLite database, and then every variant resolves locally with zero
network calls.

Only clean single-nucleotide substitutions with a plain integer ``c.``/``n.``
position are resolved locally; intronic offsets (``c.123+4``), UTR (``c.-5``,
``c.*7``), ranges and indels return ``None`` so the caller falls back to the
existing network paths. This mirrors the conservative policy of the NG_ path.
"""
from __future__ import annotations

import gzip
import logging
import re
import sqlite3
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

# UCSC ncbiRefSeq table (genePredExt + bin + name2 ...): one transcript per row.
_URLS = {
    "HG38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz",
    "HG19": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz",
}

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")

# NM_000540.3:c.7300G>A  /  NR_...:n.123A>G   (clean SNV only)
_TX_SUB_RE = re.compile(
    r'^((?:NM|NR|XM|XR)_\d+)(?:\.\d+)?:([cn])\.(\d+)([ACGT])>([ACGT])$', re.IGNORECASE)

_PRIMARY_CHROM_RE = re.compile(r'^chr(\d+|X|Y|M)$')


class TranscriptCoordinateMapper:
    def __init__(self, assembly: str = "GRCh38"):
        self.assembly = assembly.upper()
        self.assembly_key = "HG19" if ("37" in self.assembly or "19" in self.assembly) else "HG38"

        from core.config import settings
        storage_dir = Path(settings.cache_dir).parent
        storage_dir.mkdir(parents=True, exist_ok=True)
        self.db_path = storage_dir / f"ncbi_refseq_tx_{self.assembly_key.lower()}.sqlite"
        self._conn = None

    # ----- db build/connect ------------------------------------------------
    def _get_connection(self):
        if self._conn is None:
            self._ensure_db()
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn

    def _ensure_db(self):
        if self.db_path.exists():
            return
        logger.info(f"TranscriptCoordinateMapper: building {self.assembly_key} DB (one-time)...")
        url = _URLS[self.assembly_key]
        req = urllib.request.Request(url, headers={"User-Agent": "curl/8"})
        raw = urllib.request.urlopen(req, timeout=120).read()
        text = gzip.decompress(raw).decode(errors="replace")

        conn = sqlite3.connect(str(self.db_path))
        try:
            conn.execute("""CREATE TABLE transcripts (
                name TEXT PRIMARY KEY, chrom TEXT, strand TEXT,
                cds_start INTEGER, cds_end INTEGER,
                exon_starts TEXT, exon_ends TEXT)""")
            seen = set()
            for line in text.splitlines():
                f = line.split("\t")
                if len(f) < 11:
                    continue
                name, chrom, strand = f[1], f[2], f[3]
                if not _PRIMARY_CHROM_RE.match(chrom):
                    continue  # skip _alt/_fix/_random/chrUn
                base = name.split(".")[0]
                if base in seen:
                    continue  # keep first primary-assembly row per accession
                seen.add(base)
                conn.execute(
                    "INSERT OR IGNORE INTO transcripts VALUES (?,?,?,?,?,?,?)",
                    (base, chrom[3:], strand, int(f[6]), int(f[7]),
                     f[9].rstrip(","), f[10].rstrip(",")))
            conn.commit()
            logger.info("TranscriptCoordinateMapper: DB built.")
        except Exception:
            conn.close()
            if self.db_path.exists():
                self.db_path.unlink()
            raise
        finally:
            conn.close()

    # ----- coordinate math -------------------------------------------------
    @staticmethod
    def _cds_start_exonic_index(exons, strand, cds_start, cds_end):
        """Count exonic bases 5' of the first coding base, in transcript order."""
        idx = 0
        if strand == "+":
            for s, e in exons:                 # ascending
                if e <= cds_start:
                    idx += e - s
                elif s <= cds_start < e:
                    return idx + (cds_start - s)
                else:
                    return idx
        else:
            for s, e in reversed(exons):       # descending = transcript 5'->3'
                if s >= cds_end:
                    idx += e - s
                elif s < cds_end <= e:
                    return idx + (e - cds_end)
                else:
                    return idx
        return idx

    @staticmethod
    def _exonic_to_genomic(exons, strand, k):
        """0-based genomic position of the k-th exonic base in transcript order."""
        if strand == "+":
            for s, e in exons:                 # ascending
                if k < e - s:
                    return s + k
                k -= e - s
        else:
            for s, e in reversed(exons):       # descending
                if k < e - s:
                    return e - 1 - k
                k -= e - s
        return None

    def resolve_transcript_variant(self, variant: str) -> str | None:
        m = _TX_SUB_RE.match(variant)
        if not m:
            return None
        name, kind, pos, ref, alt = (
            m.group(1), m.group(2).lower(), int(m.group(3)), m.group(4).upper(), m.group(5).upper())

        try:
            conn = self._get_connection()
            row = conn.execute(
                "SELECT chrom, strand, cds_start, cds_end, exon_starts, exon_ends "
                "FROM transcripts WHERE name = ?", (name.split(".")[0],)).fetchone()
            if not row:
                return None
            chrom, strand = row["chrom"], row["strand"]
            cds_start, cds_end = row["cds_start"], row["cds_end"]
            starts = [int(x) for x in row["exon_starts"].split(",") if x != ""]
            ends = [int(x) for x in row["exon_ends"].split(",") if x != ""]
            exons = sorted(zip(starts, ends))   # ascending genomic

            if kind == "c":
                base_idx = self._cds_start_exonic_index(exons, strand, cds_start, cds_end)
                exonic_index = base_idx + (pos - 1)
            else:  # n.
                exonic_index = pos - 1

            g0 = self._exonic_to_genomic(exons, strand, exonic_index)
            if g0 is None:
                return None
            if strand == "-":
                ref, alt = ref.translate(_COMP), alt.translate(_COMP)
            return f"{chrom}:{g0 + 1}:{ref}:{alt}"
        except Exception as e:
            logger.error(f"Transcript coordinate resolution failed for {variant}: {e}")
            return None

    def close(self):
        if self._conn is not None:
            self._conn.close()
            self._conn = None
