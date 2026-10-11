"""SPECIFICITY CONTROL: cross-chromosome test.

The strongest possible validation. If the dosage effect is real and specific:
    delta(genes of C) in embryos trisomic for C   >>   delta(genes of C) in embryos trisomic for D
i.e. the elevation must follow the chromosome that is actually trisomic.

A non-specific artefact (e.g. aneuploid embryos generally over-expressing
whatever set we pick) would show elevation regardless of the match.
"""
import csv
import gzip
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
WS = Path(r"F:\xinlunwen")
META = WS / "02_SingleCell_Data" / "metadata_for_R.csv.gz"
MAT = WS / "embryo.14082cells.raw.counts.normalized.edgeR.ERCC.pseudo.counts.txt.gz"
OUT = WS / "results" / "dosage_specificity.tsv"
GENE_RE = re.compile(r"^(.*?)-(ENSG\d+)(?:\.\d+)?$")

cell_meta = {}
with gzip.open(META, "rt", encoding="utf-8", errors="replace") as fh:
    rd = csv.reader(fh)
    header = next(rd)
    ci = {c: i for i, c in enumerate(header)}
    for f in rd:
        if len(f) == len(header):
            cell_meta[f[ci["Cell.Id"]]] = (f[ci["Embryo.Id"]], f[ci["Lineage"]],
                                           f[ci["Karyotype.By.PGT.A"]])
with gzip.open(MAT, "rt", encoding="utf-8", errors="replace") as fh:
    rd = csv.reader(fh)
    head = next(rd)
    cid2col = {h.strip().strip('"'): i for i, h in enumerate(head[1:])}

grp = defaultdict(list)
for cid, (emb, lin, kar) in cell_meta.items():
    j = cid2col.get(cid)
    if j is not None:
        grp[(emb, lin)].append(j)
grp = {k: v for k, v in grp.items() if len(v) >= 3}
gkeys = list(grp.keys())

ov_epi = WS / "data" / "Ovary" / "Ovary.epi"
ensg2chrom = {}
with open(ov_epi, "r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 5:
            ensg2chrom[f[1].strip().split(".")[0]] = f[0].strip().replace("chr", "")

print("building matrix ...")
gene_means, gene_chroms = [], []
with gzip.open(MAT, "rt", encoding="utf-8", errors="replace") as fh:
    rd = csv.reader(fh)
    next(rd)
    for row in rd:
        if len(row) < 2:
            continue
        m = GENE_RE.match(row[0].strip().strip('"'))
        if not m:
            continue
        ch = ensg2chrom.get(m.group(2))
        if not ch:
            continue
        vals = row[1:]
        means = [0.0] * len(grp)
        for gi, cols in enumerate(grp.values()):
            s = 0.0; k = 0
            for j in cols:
                if j < len(vals):
                    try:
                        s += float(vals[j]); k += 1
                    except ValueError:
                        pass
            means[gi] = s / k if k else 0.0
        gene_means.append(means)
        gene_chroms.append(ch)
G = len(gene_means)
print(f"  genes: {G:,}")

pos_sum = [0.0] * len(grp)
pos_cnt = [0] * len(grp)
for gi in range(len(grp)):
    s = 0.0; k = 0
    for gene in gene_means:
        v = gene[gi]
        if v > 0:
            s += v; k += 1
    pos_sum[gi] = s
    pos_cnt[gi] = k

by_chrom = defaultdict(list)
for i, ch in enumerate(gene_chroms):
    by_chrom[ch].append(i)

emb_karyo = {}
for cid, (emb, lin, kar) in cell_meta.items():
    emb_karyo.setdefault(emb, kar)


def tri_chrom(k):
    if not k.endswith("T"):
        return None
    c = k[:-1]
    if c == "X":
        return "X"
    if c.isdigit() and 1 <= int(c) <= 22:
        return str(int(c))
    return None


def delta_for(idx_on, gi):
    so = 0.0; ko = 0
    idx_on = list(idx_on)
    for i in idx_on:
        v = gene_means[i][gi]
        if v > 0:
            so += v; ko += 1
    if ko == 0:
        return None
    off_s = pos_sum[gi] - so
    off_c = pos_cnt[gi] - ko
    if off_c <= 0 or off_s <= 0:
        return None
    return math.log2(so / ko) - math.log2(off_s / off_c)


LINEAGE = "Mural.TE"
lg = [gi for gi, g in enumerate(gkeys) if g[1] == LINEAGE]
print(f"lineage {LINEAGE}: {len(lg)} groups")

# precompute per-group trisomy chromosome
g_tri = {gi: tri_chrom(emb_karyo.get(gkeys[gi][0], "")) for gi in lg}

TESTS = ["1", "7", "13", "16", "22"]
print()
print("=" * 100)
print("SPECIFICITY: delta of chromosome-C genes in embryos trisomic for C vs for D")
print("=" * 100)
print(f"{'gene set':<10} {'trisomic grp':<14} {'nEmb':>5} {'mean delta':>12} {'vs MATCH':>10}")
print("-" * 100)

rows = []
for C in TESTS:
    if C not in by_chrom:
        continue
    idx = by_chrom[C]
    # delta per group
    d_by_g = {}
    for gi in lg:
        d = delta_for(idx, gi)
        if d is not None:
            d_by_g[gi] = d
    # matched group
    match = [d for gi, d in d_by_g.items() if g_tri.get(gi) == C]
    if len(match) < 3:
        continue
    mm = statistics.mean(match)
    print(f"{'chr' + C:<10} {'MATCH (' + C + 'T)':<14} {len(match):>5} {mm:>12.4f} {'—':>10}")
    rows.append((C, "MATCH", C, len(match), mm))
    # mismatched groups
    for D in TESTS:
        if D == C:
            continue
        mm2 = [d for gi, d in d_by_g.items() if g_tri.get(gi) == D]
        if len(mm2) < 3:
            continue
        m2 = statistics.mean(mm2)
        diff = mm - m2
        print(f"{'chr' + C:<10} {'MISMATCH (' + D + 'T)':<14} {len(mm2):>5} {m2:>12.4f} {diff:>+10.4f}")
        rows.append((C, "MISMATCH", D, len(mm2), m2, diff))

with open(OUT, "w", encoding="utf-8") as fh:
    fh.write("gene_set\ttest_type\ttrisomy\tn_embryos\tmean_delta\tdiff_vs_match\n")
    for r in rows:
        fh.write("\t".join(str(x) for x in r) + "\n")

print("-" * 100)
m = [r[4] for r in rows if r[1] == "MATCH"]
x = [r[4] for r in rows if r[1] == "MISMATCH"]
print()
print("=" * 100)
print("VERDICT")
print("=" * 100)
if m and x:
    print(f"  MATCH    (genes of C in C-trisomic embryos): mean delta = {statistics.mean(m):+.4f}  (n={len(m)})")
    print(f"  MISMATCH (genes of C in D-trisomic embryos): mean delta = {statistics.mean(x):+.4f}  (n={len(x)})")
    print(f"  specificity gap = {statistics.mean(m) - statistics.mean(x):+.4f}")
    print()
    if statistics.mean(m) > statistics.mean(x) + 0.2:
        print("  -> STRONG SPECIFICITY: chromosome-gene elevation tracks the trisomic")
        print("     chromosome. The dosage effect is chromosome-specific, not a global")
        print("     artefact of aneuploid embryos.")
    else:
        print("  -> WEAK/NO SPECIFICITY: elevation appears regardless of which chromosome")
        print("     is trisomic; the effect is not chromosome-specific.")
print(f"\nwrote {OUT}")
