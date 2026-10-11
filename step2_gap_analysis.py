# -*- coding: utf-8 -*-
"""STEP 2 — dosage compensation gap analysis.

gap = log2(1.5) - delta = 0.585 - delta
gap > 0 means the trisomic chromosome is under-expressed relative to the
uncompensated full-trisomy expectation (i.e. partial dosage compensation).
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
import pandas as pd
from scipy import stats

WS = Path(r"F:\xinlunwen")
LOG15 = 0.5849625
dr = pd.read_csv(WS / "results" / "dosage_results.tsv", sep="\t")
dr["gap"] = LOG15 - dr["delta"]

# chromosome lengths (hg19, Mb)
CHR_LEN = {1: 249.25, 2: 243.19, 3: 198.02, 4: 190.21, 5: 181.53, 6: 170.80,
           7: 159.34, 8: 145.14, 9: 138.39, 10: 133.80, 11: 135.09, 12: 133.27,
           13: 114.36, 14: 107.35, 15: 102.53, 16: 90.34, 17: 81.20, 18: 78.07,
           20: 64.44, 21: 46.94, 22: 50.82}

# gene counts per chromosome from the cleaned set
dg = pd.read_csv(WS / "results" / "dosage_genes_clean.tsv", sep="\t")
genes_per_chrom = dg.groupby("chrom").size().to_dict()
dr["n_genes_on_chrom"] = dr["chrom"].map(genes_per_chrom)
dr["chrom_length_Mb"] = dr["chrom"].map(CHR_LEN)

# ---------------------------------------------------------------- per-strata output
out = dr[["chrom", "lineage", "n_trisomy_emb", "n_disomy_emb", "delta", "gap",
          "t_stat", "p_value", "n_genes_on_chrom", "chrom_length_Mb"]].copy()
out.to_csv(WS / "results" / "dosage_compensation_gap.tsv", sep="\t", index=False)
print("wrote dosage_compensation_gap.tsv (%d rows)" % len(out))

# ---------------------------------------------------------------- by lineage
print("\n" + "=" * 80)
print("GAP BY LINEAGE")
print("=" * 80)
lineages = ["Mural.TE", "Polar.TE", "EPI", "PE"]
groups = {}
for L in lineages:
    g = dr[dr.lineage == L]["gap"].values
    groups[L] = g
    n = len(g)
    mean = np.mean(g)
    sd = np.std(g, ddof=1) if n > 1 else 0
    se = sd / np.sqrt(n) if n > 1 else 0
    lo = mean - 1.96 * se
    hi = mean + 1.96 * se
    print("  %-10s  n=%2d  mean=%+.4f  SD=%.4f  95%% CI [%+.4f, %+.4f]"
          % (L, n, mean, sd, lo, hi))

# ANOVA
f_stat, p_anova = stats.f_oneway(*[groups[L] for L in lineages])
# Kruskal-Wallis
h_stat, p_kw = stats.kruskal(*[groups[L] for L in lineages])
print("\n  One-way ANOVA: F = %.3f, p = %.4g" % (f_stat, p_anova))
print("  Kruskal-Wallis: H = %.3f, p = %.4g" % (h_stat, p_kw))

# ---------------------------------------------------------------- correlations
print("\n" + "=" * 80)
print("CORRELATIONS (across all 49 strata)")
print("=" * 80)
for label, x in [("n_trisomic", dr.n_trisomy_emb),
                 ("chromosome length (Mb)", dr.chrom_length_Mb),
                 ("n_genes on chromosome", dr.n_genes_on_chrom)]:
    rho, p = stats.spearmanr(x, dr.gap)
    print("  gap vs %-28s  rho = %+.4f  p = %.4g" % (label, rho, p))

# ---------------------------------------------------------------- save JSON
summary = {
    "log2_1_5": LOG15,
    "gap_range": [float(dr.gap.min()), float(dr.gap.max())],
    "gap_mean": float(dr.gap.mean()),
    "lineages": {L: {"n": int(len(groups[L])),
                     "mean": float(np.mean(groups[L])),
                     "sd": float(np.std(groups[L], ddof=1)),
                     "ci_lo": float(np.mean(groups[L]) - 1.96 * np.std(groups[L], ddof=1) / np.sqrt(len(groups[L]))),
                     "ci_hi": float(np.mean(groups[L]) + 1.96 * np.std(groups[L], ddof=1) / np.sqrt(len(groups[L])))}
                 for L in lineages},
    "anova_F": float(f_stat), "anova_p": float(p_anova),
    "kruskal_H": float(h_stat), "kruskal_p": float(p_kw),
    "spearman_n_trisomic": {"rho": float(stats.spearmanr(dr.n_trisomy_emb, dr.gap)[0]),
                            "p": float(stats.spearmanr(dr.n_trisomy_emb, dr.gap)[1])},
    "spearman_chrom_length": {"rho": float(stats.spearmanr(dr.chrom_length_Mb, dr.gap)[0]),
                              "p": float(stats.spearmanr(dr.chrom_length_Mb, dr.gap)[1])},
    "spearman_n_genes": {"rho": float(stats.spearmanr(dr.n_genes_on_chrom, dr.gap)[0]),
                         "p": float(stats.spearmanr(dr.n_genes_on_chrom, dr.gap)[1])},
}
with open(WS / "results" / "compensation_stats.json", "w") as fh:
    json.dump(summary, fh, indent=2)
print("\nwrote compensation_stats.json")
