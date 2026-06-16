# SPDX-License-Identifier: MIT
import logging
import re
import sqlite3
import urllib.request
from pathlib import Path

logger = logging.getLogger(__name__)

# Alignment file URLs on NCBI FTP
URLS = {
    "GRCH38": "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/alignments/ARCHIVE/all/GCF_000001405.39_109.20211119_refseqgene_alignments.gff3",
    "GRCH37": "ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/alignments/ARCHIVE/all/GCF_000001405.25_105.20201022_refseqgene_alignments.gff3"
}

class LocalCoordinateMapper:
    def __init__(self, assembly: str = "GRCh38"):
        self.assembly = assembly.upper()
        if "37" in self.assembly or "19" in self.assembly:
            self.assembly_key = "GRCH37"
        else:
            self.assembly_key = "GRCH38"

        # Place the database file in storage (parent of cache_dir)
        from core.config import settings
        storage_dir = Path(settings.cache_dir).parent
        storage_dir.mkdir(parents=True, exist_ok=True)
        self.db_path = storage_dir / f"refseqgene_mappings_{self.assembly_key.lower()}.sqlite"
        self._conn = None

    def _get_connection(self):
        if self._conn is None:
            self._ensure_db()
            self._conn = sqlite3.connect(str(self.db_path))
            self._conn.row_factory = sqlite3.Row
        return self._conn

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
            elif num == 12920:
                return "MT"
        except (IndexError, ValueError):
            pass
        return None

    def _ensure_db(self):
        if self.db_path.exists():
            return

        logger.info(f"LocalCoordinateMapper: Database for {self.assembly_key} not found. Downloading and building...")
        url = URLS[self.assembly_key]

        # Download the file
        temp_gff = self.db_path.with_suffix(".gff3")
        try:
            urllib.request.urlretrieve(url, temp_gff)
            logger.info(f"Downloaded GFF3 alignment file to {temp_gff}")
        except Exception as e:
            logger.error(f"Failed to download RefSeqGene alignments: {e}")
            if temp_gff.exists():
                temp_gff.unlink()
            raise e

        # Parse the GFF3 and create SQLite database
        conn = sqlite3.connect(str(self.db_path))
        cursor = conn.cursor()
        cursor.execute("""
            CREATE TABLE IF NOT EXISTS mappings (
                ng_accession TEXT PRIMARY KEY,
                chromosome TEXT,
                chrom_start INTEGER,
                chrom_end INTEGER,
                strand TEXT,
                target_strand TEXT
            )
        """)

        try:
            with open(temp_gff) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) < 9:
                        continue

                    nc_ac = parts[0]
                    feature_type = parts[2]
                    if feature_type != "match":
                        continue

                    chrom_start = int(parts[3])
                    chrom_end = int(parts[4])
                    strand = parts[6]
                    attributes = parts[8]

                    # Parse Target
                    target_match = re.search(r"Target=(NG_\d+)(?:\.\d+)?\s+(\d+)\s+(\d+)\s+([\+\-])", attributes)
                    if target_match:
                        ng_ac = target_match.group(1)
                        target_strand = target_match.group(4)
                        chrom = self._nc_to_chrom(nc_ac)
                        if chrom:
                            cursor.execute("""
                                INSERT OR REPLACE INTO mappings
                                (ng_accession, chromosome, chrom_start, chrom_end, strand, target_strand)
                                VALUES (?, ?, ?, ?, ?, ?)
                            """, (ng_ac, chrom, chrom_start, chrom_end, strand, target_strand))
            conn.commit()
            logger.info("Successfully populated local RefSeqGene mappings SQLite database.")
        except Exception as e:
            logger.error(f"Failed to parse and store RefSeqGene mappings: {e}")
            conn.close()
            if self.db_path.exists():
                self.db_path.unlink()
            raise e
        finally:
            conn.close()
            if temp_gff.exists():
                temp_gff.unlink()

    def resolve_ng_variant(self, variant: str) -> str | None:
        """
        Translates an NG_ variant (e.g. NG_008866.1:g.6023G>A)
        to a VCF string (e.g. 19:38434722:G:A) locally.
        """
        match = re.match(r'^(NG_\d+)(?:\.\d+)?:g\.(\d+)([ATGC]*)(?:>([ATGC]*))?$', variant, re.IGNORECASE)
        if not match:
            return None

        ng_ac, local_pos_str, ref, alt = match.groups()
        local_pos = int(local_pos_str)
        ref = ref or ""
        alt = alt or ""

        try:
            conn = self._get_connection()
            cursor = conn.cursor()
            cursor.execute("""
                SELECT chromosome, chrom_start, chrom_end, strand, target_strand
                FROM mappings WHERE ng_accession = ?
            """, (ng_ac,))
            row = cursor.fetchone()

            if not row:
                # Try matching without the version suffix just in case
                cursor.execute("""
                    SELECT chromosome, chrom_start, chrom_end, strand, target_strand
                    FROM mappings WHERE ng_accession LIKE ?
                """, (f"{ng_ac}%",))
                row = cursor.fetchone()

            if not row:
                return None

            chrom, chrom_start, chrom_end, strand, target_strand = row

            # Decide strand mapping. If strand of target matches chromosome strand:
            # genomic_pos = chrom_start + local_pos - 1
            # If it is reverse:
            # genomic_pos = chrom_end - local_pos + 1
            is_reverse = (strand != target_strand)
            if not is_reverse:
                genomic_pos = chrom_start + local_pos - 1
            else:
                genomic_pos = chrom_end - local_pos + 1

            if strand == "-":
                comp = {"A": "T", "T": "A", "G": "C", "C": "G", "a": "t", "t": "a", "g": "c", "c": "g"}
                ref = "".join(comp.get(base, base) for base in reversed(ref))
                alt = "".join(comp.get(base, base) for base in reversed(alt))

            return f"{chrom}:{genomic_pos}:{ref}:{alt}"
        except Exception as e:
            logger.error(f"Error resolving NG_ variant locally: {e}")
            return None

    def close(self):
        if self._conn is not None:
            self._conn.close()
            self._conn = None
