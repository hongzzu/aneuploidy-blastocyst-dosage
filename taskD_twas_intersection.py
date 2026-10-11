# -*- coding: utf-8 -*-
"""Write the TWAS intersection and the full core-gene TWAS table."""
import gzip, csv, json
from pathlib import Path

WS = Path(r"F:\xinlunwen")
RES = WS / "results"
p = WS / "work" / "agwas" / "twas_sumstats.total.csv.gz"

rows_by_gene = {}
with gzip.open(p, 'rt', encoding='utf-8', errors='replace') as fh:
    rd = csv.reader(fh)
    head = next(rd)
    pi, gi, ci, si = head.index('PHENO'), head.index('GENE'), head.index('COMBINED_PVALUE'), head.index('is_sig')
    for f in rd:
        if len(f) < 9:
            continue
        gene = f[gi].strip()
        pheno = f[pi].strip()
        try:
            pval = float(f[ci])
        except ValueError:
            continue
        is_sig = f[si].strip().upper() == 'TRUE'
        rows_by_gene.setdefault(gene, {})[pheno] = (pval, is_sig)

# define significance
THRESH = 3e-6
sig_by_pheno = {}
all_phenos = set()
for gene, d in rows_by_gene.items():
    for pheno, (pval, is_sig) in d.items():
        all_phenos.add(pheno)
        if is_sig or pval < THRESH:
            sig_by_pheno.setdefault(pheno, set()).add(gene)

# intersection aneuploidy x female_MeanCO
inter = sig_by_pheno.get('aneuploidy', set()) & sig_by_pheno.get('female_MeanCO', set())
inter = {g.strip() for g in inter}

# write intersection
with open(RES / "twas_intersection.tsv", "w", encoding="utf-8") as fh:
    fh.write("gene\taneuploidy_p\taneuploidy_sig\tfemale_MeanCO_p\tfemale_MeanCO_sig\n")
    for g in sorted(inter):
        a = rows_by_gene.get(g, {}).get('aneuploidy', (None, False))
        m = rows_by_gene.get(g, {}).get('female_MeanCO', (None, False))
        fh.write("%s\t%.3g\t%s\t%.3g\t%s\n"
                 % (g.strip(), a[0] if a[0] else 'NA', a[1], m[0] if m[0] else 'NA', m[1]))

# write full core-gene table
core = ['C14orf39', 'RIBC2', 'SMC1B', 'RNF212', 'CCNB1IP1', 'SIX6', 'MRPL51']
with open(RES / "twas_core_genes_all_phenotypes.tsv", "w", encoding="utf-8") as fh:
    pheno_order = sorted(all_phenos)
    fh.write("gene\t" + "\t".join(pheno_order) + "\n")
    for g in core:
        d = rows_by_gene.get(g, {})
        row = []
        for ph in pheno_order:
            v = d.get(ph, (None, False))
            row.append("%.3g%s" % (v[0], "*" if v[1] else ""))
        fh.write(g + "\t" + "\t".join(row) + "\n")

# check core candidates
core_in_inter = [g for g in core if g in inter]
print("aneuploidy significant: %d" % len(sig_by_pheno.get('aneuploidy', set())))
print("female_MeanCO significant: %d" % len(sig_by_pheno.get('female_MeanCO', set())))
print("intersection: %d" % len(inter))
print("core candidates in intersection: %d / %d" % (len(core_in_inter), len(core)))
print("  %s" % core_in_inter)

# also check all significant phenotypes
print("\nAll phenotypes with significant genes:")
for ph in sorted(sig_by_pheno):
    print("  %-28s %d genes" % (ph, len(sig_by_pheno[ph])))
    if len(sig_by_pheno[ph]) <= 10:
        print("      %s" % sorted(g.strip() for g in sig_by_pheno[ph]))
