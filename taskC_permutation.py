# -*- coding: utf-8 -*-
"""TASK C — permutation tests for the chromosome-specificity contrast.

T_obs = mean(matched delta) - mean(mismatched delta) = +0.479
Simple permutation:    randomly shuffle the MATCH/MISMATCH label across all 25 rows
Stratified permutation: within each gene_set (C), shuffle which trisomy group is
                        designated MATCH (5 groups, each with 1 MATCH + 4 MISMATCH)
10,000 draws each.
"""
import io
import sys
from collections import defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd

WS = Path(r"F:\xinlunwen")
df = pd.read_csv(WS / "results" / "dosage_specificity.tsv", sep="\t")

# parse
gsets = sorted(df.gene_set.unique(), key=lambda s: int(s))
match_delta = {}
mism_delta = defaultdict(dict)   # gene_set -> {trisomy: mean_delta}
for r in df.itertuples():
    if r.test_type == "MATCH":
        match_delta[r.gene_set] = r.mean_delta
    else:
        mism_delta[r.gene_set][r.trisomy] = r.mean_delta

# T_obs
matched_vals = list(match_delta.values())
mismatched_vals = [v for gs in mism_delta.values() for v in gs.values()]
T_obs = np.mean(matched_vals) - np.mean(mismatched_vals)
print("=" * 80)
print("TASK C  permutation tests")
print("=" * 80)
print("T_obs = %.4f" % T_obs)
print("  matched    n=%d  mean=%+.4f" % (len(matched_vals), np.mean(matched_vals)))
print("  mismatched n=%d  mean=%+.4f" % (len(mismatched_vals), np.mean(mismatched_vals)))

# ---------------------------------------------------------------- simple permutation
# shuffle the MATCH/MISMATCH label across all 25 rows
all_vals = matched_vals + mismatched_vals  # 5 + 20 = 25 values
rng = np.random.default_rng(42)
N = 10000
null_simple = np.zeros(N)
for i in range(N):
    perm = rng.permutation(all_vals)
    null_simple[i] = np.mean(perm[:5]) - np.mean(perm[5:])
p_simple = (np.sum(np.abs(null_simple) >= abs(T_obs)) + 1) / (N + 1)
print("\nSimple permutation (%d draws):" % N)
print("  null mean %+.4f  SD %.4f  95%% range [%.4f, %.4f]"
      % (np.mean(null_simple), np.std(null_simple),
         np.percentile(null_simple, 2.5), np.percentile(null_simple, 97.5)))
print("  p (two-sided) = %.6f" % p_simple)

# ---------------------------------------------------------------- stratified permutation
# within each gene_set C, shuffle which trisomy group is designated MATCH
# (the 4 MISMATCH groups follow from the remaining 4)
# this preserves the (embryo x lineage) pairing and the group sizes
null_strat = np.zeros(N)
for i in range(N):
    dms, dmm = [], []
    for C in gsets:
        # available trisomy groups for this gene set = the 5 test chromosomes
        # (1 is MATCH, 4 are MISMATCH; we shuffle which one is designated MATCH)
        available = list(mism_delta[C].keys()) + [C]  # all 5 trisomy values for this C
        all_d = {**mism_delta[C], C: match_delta[C]}
        chosen_match = rng.choice(available)
        dms.append(all_d[chosen_match])
        for D in available:
            if D != chosen_match:
                dmm.append(all_d[D])
    null_strat[i] = np.mean(dms) - np.mean(dmm)
p_strat = (np.sum(np.abs(null_strat) >= abs(T_obs)) + 1) / (N + 1)
print("\nStratified permutation (%d draws):" % N)
print("  null mean %+.4f  SD %.4f  95%% range [%.4f, %.4f]"
      % (np.mean(null_strat), np.std(null_strat),
         np.percentile(null_strat, 2.5), np.percentile(null_strat, 97.5)))
print("  p (two-sided) = %.6f" % p_strat)

# ---------------------------------------------------------------- summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("  T_obs                    = %+.4f" % T_obs)
print("  simple permutation    p  = %.6f" % p_simple)
print("  stratified permutation p = %.6f" % p_strat)

with open(WS / "results" / "permutation_results.tsv", "w", encoding="utf-8") as fh:
    fh.write("metric\tvalue\n")
    fh.write("T_obs\t%.6f\n" % T_obs)
    fh.write("simple_p\t%.6f\n" % p_simple)
    fh.write("stratified_p\t%.6f\n" % p_strat)
    fh.write("n_draws\t%d\n" % N)
    fh.write("null_simple_mean\t%.6f\n" % np.mean(null_simple))
    fh.write("null_simple_sd\t%.6f\n" % np.std(null_simple))
    fh.write("null_strat_mean\t%.6f\n" % np.mean(null_strat))
    fh.write("null_strat_sd\t%.6f\n" % np.std(null_strat))

np.save(WS / "results" / "permutation_null_simple.npy", null_simple)
np.save(WS / "results" / "permutation_null_strat.npy", null_strat)
print("\nwrote results/permutation_results.tsv and .npy files")
