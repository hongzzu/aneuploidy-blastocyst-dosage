# -*- coding: utf-8 -*-
"""TASK B alternative — rank-based GSEA (no background needed).

Runs GSEA-prerank on a per-gene dosage ranking built from the gene x stratum
Welch t statistics in gene_level_dosage_t.tsv.  Because ranking uses ALL genes
with a computed t (no expression-floor filter, no background definition), this
is immune to the background-matching artefact that collapsed the ORA result.
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd
import gseapy as gp

WS = Path(r"F:\xinlunwen")
RES = WS / "results"

# ---------------------------------------------------------------- ranking
t = pd.read_csv(RES / "gene_level_dosage_t.tsv", sep="\t")
print("gene x stratum tests: %d" % len(t))

# rank on the mural.TE stratum (the dominant layer; matches Fig 3a-c)
mte = t[t.lineage == "Mural.TE"].copy()
mte = mte.sort_values("t", ascending=False).drop_duplicates("gene", keep="first")
print("unique genes with a mural.TE t statistic: %d" % len(mte))

rnk = mte[["gene", "t"]].rename(columns={"gene": "gene", "t": "t_muralTE"})
rnk = rnk.sort_values("t_muralTE", ascending=False)
print("ranking range: t = %.2f .. %.2f" % (rnk.t_muralTE.min(), rnk.t_muralTE.max()))

# ---------------------------------------------------------------- GSEA-prerank
LIB = ["MSigDB_Hallmark_2020", "KEGG_2021_Human", "Reactome_2022"]
print("\nrunning GSEA-prerank ...", flush=True)
results = {}
for lib in LIB:
    print("  %s ..." % lib, flush=True)
    try:
        r = gp.prerank(rnk=rnk, gene_sets=lib, outdir=None,
                       min_size=15, max_size=500, permutation_num=1000,
                       seed=42, no_plot=True, threads=4)
        df = r.res2d
        df["NES"] = df["NES"].astype(float)
        df["FDR q-val"] = df["FDR q-val"].astype(float)
        df = df.sort_values("NES", ascending=False)
        results[lib] = df
        sig = df[df["FDR q-val"] < 0.05]
        print("     sets tested: %d   NES>0 & FDR<0.05: %d   NES<0 & FDR<0.05: %d"
              % (len(df),
                 len(df[(df.NES > 0) & (df["FDR q-val"] < 0.05)]),
                 len(df[(df.NES < 0) & (df["FDR q-val"] < 0.05)])))
        for _, r_ in df.head(10).iterrows():
            print("       %-50s NES=%+6.2f  FDR=%.4g"
                  % (str(r_["Term"])[:50], r_["NES"], r_["FDR q-val"]))
    except Exception as e:
        print("     FAILED: %s: %s" % (type(e).__name__, e))

# ---------------------------------------------------------------- key pathways
print("\n" + "=" * 100)
print("KEY PATHWAYS  (rank-based GSEA)")
print("=" * 100)
KEY = ["MYC Targets V1", "MYC Targets V2", "E2F Targets",
       "Oxidative Phosphorylation", "G2-M Checkpoint", "DNA Repair",
       "Ribosome", "Cell Cycle", "Unfolded Protein Response"]
for k in KEY:
    best = None
    for lib, df in results.items():
        m = df[df["Term"].str.contains(k, case=False, na=False)]
        if len(m):
            r_ = m.iloc[0]
            if best is None or abs(r_["NES"]) > abs(best[1]):
                best = (lib, r_["NES"], r_["FDR q-val"])
    if best:
        print("  %-28s %-24s NES=%+6.2f  FDR=%.4g  %s"
              % (k, best[0][:24], best[1], best[2],
                 "SIG" if best[2] < 0.05 else "not sig"))
    else:
        print("  %-28s not tested (set size out of range)" % k)

# ---------------------------------------------------------------- write
out = RES / "enrichment_rank_gsea.tsv"
with open(out, "w", encoding="utf-8") as fh:
    fh.write("library\tterm\tNES\tFDR_q_val\tFDR_p_val\n")
    for lib, df in results.items():
        for _, r in df.iterrows():
            fh.write("%s\t%s\t%.4f\t%.4g\t%.4g\n"
                     % (lib, r["Term"], r["NES"], r["FDR q-val"], r["NOM p-val"]))
print("\nwrote %s" % out)
