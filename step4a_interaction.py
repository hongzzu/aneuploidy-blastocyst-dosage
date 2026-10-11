# -*- coding: utf-8 -*-
"""Task 1: interaction term on the 5 test chromosomes (complete 4-lineage data)."""
import sys, pandas as pd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import statsmodels.formula.api as smf
from statsmodels.stats.anova import anova_lm

df = pd.read_csv(r"F:\xinlunwen\results\dosage_results.tsv", sep="\t")
sub = df[df.chrom.isin([1, 7, 13, 16, 22])].copy()
print("subset: %d rows, chromosomes %s, lineages %s"
      % (len(sub), sorted(sub.chrom.unique()), sorted(sub.lineage.unique())))
print("complete design:", len(sub) == 20)

try:
    model = smf.ols("delta ~ C(chrom) * C(lineage)", data=sub)
    res = model.fit()
    table = anova_lm(res, typ=2)
    print(table)
    print("\nR-squared: %.4f" % res.rsquared)
except Exception as e:
    print("FAILED: %s: %s" % (type(e).__name__, e))
