"""TASK 2a: functional enrichment of the cleaned dosage-sensitive gene set.

Uses gseapy with a local background so the test is not biased by the fact that
only genes present on the expression matrix could ever be detected.

Two enrichment runs:
  (1) Enrichr online libraries (GO BP/CC/MF, KEGG, Reactome) with the cleaned
      gene set as foreground and the tested-gene universe as background.
  (2) A focused check of the specific themes claimed in the Discussion:
      DNA repair, chromatin, cell cycle, mitochondrial/ribosomal.
"""
import json
import sys
from collections import Counter
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import gseapy as gp

WS = Path(r"F:\xinlunwen")
RES = WS / "results"

# ---------------- foreground / background ----------------
fg = []
with open(RES / "dosage_genes_final.tsv", "r", encoding="utf-8") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    c = {x: i for i, x in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= len(hdr):
            fg.append(f[c["gene"]])
fg = sorted(set(fg))
print(f"foreground (cleaned dosage genes): {len(fg):,}")

bg = set()
with open(RES / "gene_level_dosage_t.tsv", "r", encoding="utf-8") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    c = {x: i for i, x in enumerate(hdr)}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= len(hdr):
            bg.add(f[c["gene"]])
bg = sorted(bg)
print(f"background (all genes testable): {len(bg):,}")

# ---------------- enrichr ----------------
LIB = ["GO_Biological_Process_2023", "GO_Molecular_Function_2023",
       "GO_Cellular_Component_2023", "KEGG_2021_Human",
       "Reactome_2022", "MSigDB_Hallmark_2020"]
results = {}
for lib in LIB:
    print(f"\nquerying {lib} ...", flush=True)
    try:
        enr = gp.enrichr(gene_list=fg, gene_sets=lib, background=bg,
                         organism="human", outdir=None, cutoff=0.5)
        df = enr.results
        print(f"  columns: {list(df.columns)}")
        df = df.sort_values("Adjusted P-value")
        results[lib] = df
        sig = df[df["Adjusted P-value"] < 0.05]
        print(f"  terms tested: {len(df):,}   FDR<0.05: {len(sig):,}")
        ovcol = "Overlap" if "Overlap" in df.columns else None
        for _, r in sig.head(8).iterrows():
            ov = r[ovcol] if ovcol else ""
            print(f"    {str(r['Term'])[:62]:<64} FDR={r['Adjusted P-value']:.2e}  {ov}")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")

# ---------------- save ----------------
out = RES / "enrichment_results.tsv"
with open(out, "w", encoding="utf-8") as fh:
    fh.write("library\tterm\toverlap\tp_value\tadjusted_p_value\tgenes\n")
    for lib, df in results.items():
        has_ov = "Overlap" in df.columns
        for _, r in df.iterrows():
            ov = r["Overlap"] if has_ov else ""
            fh.write(f"{lib}\t{r['Term']}\t{ov}\t{r['P-value']:.3e}\t"
                     f"{r['Adjusted P-value']:.3e}\t{r.get('Genes', '')}\n")
print(f"\nwrote {out}")

# ---------------- focused theme check ----------------
print()
print("=" * 100)
print("FOCUSED THEME CHECK (Discussion claims)")
print("=" * 100)
THEMES = {
    "DNA repair": ["XRCC4", "RAD51", "MSH4", "MSH5", "MLH3", "BAP1", "MNAT1",
                   "RAD50", "MRE11", "NBN", "BRCA1", "BRCA2", "PARP1", "FANCM"],
    "Chromatin / histone": ["SETD3", "CHD8", "EZH2", "KDM5C", "H2AFZ", "HIST1H1C",
                            "KAT6A", "KMT2A", "SMARCA4", "ARID1A"],
    "Cell cycle": ["ANAPC13", "PLK1", "AURKB", "BUB1", "BUB1B", "MAD2L1", "TTK",
                   "NDC80", "MIS12", "CDC20", "CDK1", "CCNB1"],
    "Mitochondrial / ribosomal": ["NDUFB6", "NDUFB4", "COX17", "MRPL22", "MRPL52",
                                  "MRPL50", "MRPS6", "MRPL51", "SDHB", "ATP5F1A"],
    "RNA processing": ["DICER1", "GEMIN2", "DHX29", "DBR1", "TSEN54", "SNRNP200"],
}
fgs = set(fg)
print(f"{'theme':<28} {'in foreground':>14} {'theme size':>11} {'%':>7}")
for name, genes in THEMES.items():
    hit = [g for g in genes if g in fgs]
    tot = len(genes)
    print(f"{name:<28} {len(hit):>6}/{tot:<7} {'':>11} {100*len(hit)/tot:>6.1f}%")
    if hit:
        print(f"    hits: {', '.join(hit)}")

# overall summary
summary = {}
for lib, df in results.items():
    sig = df[df["Adjusted P-value"] < 0.05]
    summary[lib] = dict(tested=int(len(df)), significant=int(len(sig)),
                        top=sig.head(10)["Term"].tolist() if len(sig) else [])
json.dump(summary, open(RES / "enrichment_summary.json", "w"), indent=2)
print(f"\nwrote {RES / 'enrichment_summary.json'}")
