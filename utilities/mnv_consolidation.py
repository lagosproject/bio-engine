# SPDX-License-Identifier: MIT
"""
Conservative MNV (multi-nucleotide variant) consolidation.
==========================================================

Adjacent single-nucleotide changes that lie on the *same allele* (cis) form a
single multi-nucleotide variant and must be annotated together: annotating them
as independent SNVs makes VEP/OpenCRAVAT evaluate each against the reference in
isolation, yielding the wrong protein consequence (e.g. two benign-looking
missense calls that are really a single nonsense/stop). HGVS requires such cis
changes to be written as one ``delins``.

Phase policy (conservative, standards-compliant):
  * only **strictly contiguous** SNVs (genomic distance 1) are merged;
  * only when every position in the run carries a single SNV call whose
    ``genotype`` is in ``mergeable_genotypes`` (default: homozygous ALT, where
    cis is trivially guaranteed);
  * heterozygous / mixed / non-contiguous runs are left untouched, because
    Sanger does not phase them — merging would fabricate a haplotype.

Non-contiguous same-codon merges (e.g. c.1 + c.3, which HGVS would still join
into a span ``delins`` including the unchanged middle base) need the reference
sequence and are intentionally out of scope here.
"""
from __future__ import annotations

import logging
from typing import Any

from utilities.ensembl_hgvs import EnsemblHGVS

logger = logging.getLogger(__name__)

DEFAULT_MERGEABLE = ("hom. ALT",)


def _is_single_base_snv(row: list[Any], idx: dict[str, int]) -> bool:
    if idx.get("type") is not None and str(row[idx["type"]]).upper() != "SNV":
        return False
    ref, alt = row[idx["ref"]], row[idx["alt"]]
    return isinstance(ref, str) and isinstance(alt, str) and len(ref) == 1 and len(alt) == 1


def consolidate_snv_runs(
    columns: list[str],
    rows: list[list[Any]],
    source_ac: str,
    h_type: str,
    cds_start: int = 0,
    cds_end: int | None = None,
    mergeable_genotypes: tuple[str, ...] = DEFAULT_MERGEABLE,
) -> tuple[list[str], list[list[Any]]]:
    """
    Merge contiguous, in-cis SNV runs into single ``delins`` rows.

    Returns a new ``(columns, rows)`` pair; the input is not mutated. Rows that
    are not part of a mergeable run (het, non-contiguous, indels, multi-allelic
    positions) are passed through unchanged and in their original order.
    """
    idx = {c: i for i, c in enumerate(columns)}
    if not {"pos", "ref", "alt"} <= idx.keys():
        return columns, rows

    out_columns = list(columns)
    if "hgvs" not in idx:
        out_columns.append("hgvs")
        idx = {c: i for i, c in enumerate(out_columns)}
    hgvs_i = idx["hgvs"]

    def _padded(row: list[Any]) -> list[Any]:
        row = list(row)
        while len(row) < len(out_columns):
            row.append("")
        return row

    geno_i = idx.get("genotype")

    # Positions with more than one SNV call are multi-allelic -> never merge.
    pos_snv_count: dict[int, int] = {}
    for row in rows:
        if _is_single_base_snv(row, idx):
            try:
                p = int(row[idx["pos"]])
            except (ValueError, TypeError):
                continue
            pos_snv_count[p] = pos_snv_count.get(p, 0) + 1

    def _eligible(row: list[Any]) -> bool:
        if not _is_single_base_snv(row, idx):
            return False
        if geno_i is None or row[geno_i] not in mergeable_genotypes:
            return False
        try:
            return pos_snv_count.get(int(row[idx["pos"]]), 0) == 1
        except (ValueError, TypeError):
            return False

    # Order eligible rows by position to detect contiguous runs.
    eligible = sorted(
        ((int(r[idx["pos"]]), ri) for ri, r in enumerate(rows) if _eligible(r)),
        key=lambda t: t[0],
    )

    # Build runs of strictly contiguous positions with identical genotype.
    runs: list[list[int]] = []
    cur: list[int] = []
    for pos, ri in eligible:
        if cur:
            prev_pos = int(rows[cur[-1]][idx["pos"]])
            same_geno = rows[ri][geno_i] == rows[cur[-1]][geno_i]
            if pos == prev_pos + 1 and same_geno:
                cur.append(ri)
                continue
            if len(cur) > 1:
                runs.append(cur)
            cur = []
        cur = cur or []
        cur.append(ri)
    if len(cur) > 1:
        runs.append(cur)

    # Map: first-row-index -> merged row ; member-row-index -> drop
    merged_at: dict[int, list[Any]] = {}
    drop: set[int] = set()
    for run in runs:
        run_sorted = sorted(run, key=lambda ri: int(rows[ri][idx["pos"]]))
        first = run_sorted[0]
        start = int(rows[first][idx["pos"]])
        ref = "".join(str(rows[ri][idx["ref"]]) for ri in run_sorted)
        alt = "".join(str(rows[ri][idx["alt"]]) for ri in run_sorted)

        merged = _padded(rows[first])
        merged[idx["ref"]] = ref
        merged[idx["alt"]] = alt
        if idx.get("type") is not None:
            merged[idx["type"]] = "MNV"
        try:
            merged[hgvs_i] = EnsemblHGVS.format_hgvs(
                source_ac, h_type, start, ref, alt, cds_start=cds_start, cds_end=cds_end)
        except Exception as e:
            logger.error(f"Failed to format merged MNV HGVS at {source_ac}:{start}: {e}")
            merged[hgvs_i] = ""
        merged_at[first] = merged
        drop.update(run_sorted[1:])

    # Rebuild rows preserving original order.
    new_rows: list[list[Any]] = []
    for ri, row in enumerate(rows):
        if ri in drop:
            continue
        new_rows.append(merged_at[ri] if ri in merged_at else _padded(row))
    return out_columns, new_rows
