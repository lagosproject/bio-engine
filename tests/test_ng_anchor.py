# SPDX-License-Identifier: MIT
"""
Battery for NGAnchorResolver: local NG_ -> chromosome resolution.

Ground truth (NG_ -> NC_) was produced by the live pipeline (Ensembl
variant_recoder) and is cached in real RYR1 jobs:

    NG_008866.1:g.65631C>T -> NC_000019.10:g.38494330C>T  (19:38494330:C:T)
    NG_008866.1:g.65880G>A -> NC_000019.10:g.38494579G>A  (19:38494579:G:A)

The decisive assertions are (1) local output == network ground truth, and
(2) network calls drop from one-per-job to one-ever-per-reference.
"""
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from utilities.ng_anchor import NGAnchorResolver

# What the CURRENT network path (Ensembl) returns, as chrom:pos:ref:alt.
GROUND_TRUTH = {
    "NG_008866.1:g.65631C>T": "19:38494330:C:T",
    "NG_008866.1:g.65880G>A": "19:38494579:G:A",
}


class CountingLearner:
    """Stand-in for the one-shot Ensembl recode; counts how often it's hit."""
    def __init__(self, table):
        self.table = table
        self.calls = 0

    def __call__(self, ng_hgvs):
        self.calls += 1
        return self.table.get(ng_hgvs)


class TestNGAnchorResolver(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.cache_path = os.path.join(self.tmp.name, "ng_anchors.json")

    def tearDown(self):
        self.tmp.cleanup()

    def _resolver(self):
        return NGAnchorResolver(self.cache_path)

    # 1) Correctness: local result matches the network ground truth exactly.
    def test_local_matches_network_ground_truth(self):
        learner = CountingLearner(GROUND_TRUTH)
        res = self._resolver().resolve_batch(
            list(GROUND_TRUTH.keys()), "GRCh38", learner)
        self.assertEqual(res, GROUND_TRUTH)

    # 2) Only ONE network call is needed to resolve the whole batch (the anchor).
    def test_single_anchor_call_for_whole_batch(self):
        learner = CountingLearner(GROUND_TRUTH)
        batch = list(GROUND_TRUTH.keys()) + [
            "NG_008866.1:g.5131A>G", "NG_008866.1:g.158000T>A",
        ]
        res = self._resolver().resolve_batch(batch, "GRCh38", learner)
        self.assertEqual(learner.calls, 1)
        self.assertEqual(len(res), len(batch))
        # Linear law from the proven anchor (offset 38428699).
        self.assertEqual(res["NG_008866.1:g.5131A>G"], "19:38433830:A:G")

    # 3) Persistence: a second resolver reuses the disk cache, ZERO calls.
    def test_anchor_persists_to_disk(self):
        learner1 = CountingLearner(GROUND_TRUTH)
        self._resolver().resolve_batch(list(GROUND_TRUTH.keys()), "GRCh38", learner1)
        self.assertEqual(learner1.calls, 1)

        learner2 = CountingLearner(GROUND_TRUTH)
        res = self._resolver().resolve_batch(list(GROUND_TRUTH.keys()), "GRCh38", learner2)
        self.assertEqual(learner2.calls, 0)   # served from disk
        self.assertEqual(res, GROUND_TRUTH)

    # 4) Network-call reduction across many jobs on the same reference.
    def test_calls_drop_across_jobs(self):
        learner = CountingLearner(GROUND_TRUTH)
        for _ in range(5):  # five jobs, shared disk cache
            self._resolver().resolve_batch(
                list(GROUND_TRUTH.keys()), "GRCh38", learner)
        self.assertEqual(learner.calls, 1)   # vs 5 for the current pipeline

    # 5) Minus-strand reference: positions mirror and bases complement.
    def test_minus_strand(self):
        # Synthetic minus-strand gene: NG pos 100 (C>T) maps to chrom 5000 (G>A).
        # const = pos + ng_pos = 5100 ; another pos 110 -> 4990, complemented.
        gt = {"NG_999999.1:g.100C>T": "7:5000:G:A"}
        learner = CountingLearner(gt)
        r = self._resolver()
        out = r.resolve_batch(
            ["NG_999999.1:g.100C>T", "NG_999999.1:g.110A>G"], "GRCh38", learner)
        self.assertEqual(out["NG_999999.1:g.100C>T"], "7:5000:G:A")
        self.assertEqual(out["NG_999999.1:g.110A>G"], "7:4990:T:C")

    # 6) Conservative: non-substitutions are left for the network (omitted).
    def test_indels_left_to_network(self):
        learner = CountingLearner(GROUND_TRUTH)
        batch = ["NG_008866.1:g.65631C>T", "NG_008866.1:g.100_102del"]
        out = self._resolver().resolve_batch(batch, "GRCh38", learner)
        self.assertIn("NG_008866.1:g.65631C>T", out)
        self.assertNotIn("NG_008866.1:g.100_102del", out)

    # 7) Anchor learn aborts (no cache) if the recode is inconsistent.
    def test_inconsistent_recode_not_cached(self):
        bad = {"NG_008866.1:g.65631C>T": "19:38494330:A:G"}  # ref/alt don't match
        learner = CountingLearner(bad)
        r = self._resolver()
        out = r.resolve_batch(["NG_008866.1:g.65631C>T"], "GRCh38", learner)
        self.assertEqual(out, {})
        self.assertIsNone(r.get_anchor("NG_008866.1", "GRCh38"))

    # 8) Assembly isolation: GRCh37 and GRCh38 keep separate anchors.
    def test_assembly_isolation(self):
        gt37 = {"NG_008866.1:g.65631C>T": "19:38985232:C:T"}
        r = self._resolver()
        r.resolve_batch(list(GROUND_TRUTH.keys()), "GRCh38", CountingLearner(GROUND_TRUTH))
        r.resolve_batch(list(gt37.keys()), "GRCh37", CountingLearner(gt37))
        self.assertEqual(r.get_anchor("NG_008866", "GRCh38")["const"], 38428699)
        self.assertEqual(r.get_anchor("NG_008866", "GRCh37")["const"], 38919601)


if __name__ == "__main__":
    unittest.main(verbosity=2)
