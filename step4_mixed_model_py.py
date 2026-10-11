# -*- coding: utf-8 -*-
"""STEP 4 — mixed-effects model using statsmodels MixedLM.

delta ~ C(chromosome) + C(lineage) + C(chromosome):C(lineage)
groups = chromosome (random intercept per chromosome)
"""
import sys
import warnings
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd
import statsmodels.formula.api as smf

df = pd.read_csv(r"F:\xinlunwen\results\dosage_results.tsv", sep="\t")
print("rows: %d  chromosomes: %d  lineages: %d"
      % (len(df), df.chrom.nunique(), df.lineage.nunique()))

# Model: fixed effects = chromosome + lineage + interaction
# Random effect = chromosome (random intercept per chromosome)
# Note: chromosome is both a fixed effect and the grouping variable is not
# ideal, so we use chromosome as the group and test lineage + interaction as fixed.

# Model without interaction first (simpler, avoids singular fit)
print("\n=== Model A: delta ~ C(lineage), groups = chromosome ===")
try:
    md = smf.mixedlm("delta ~ C(lineage)", df, groups=df["chrom"])
    fit = md.fit(reml=True)
    print(fit.summary())
except Exception as e:
    print("Model A FAILED: %s: %s" % (type(e).__name__, e))

# Model with interaction using lineage as group (chromosome as fixed effect)
print("\n=== Model B: delta ~ C(chrom), groups = lineage ===")
try:
    md2 = smf.mixedlm("delta ~ C(chrom)", df, groups=df["lineage"])
    fit2 = md2.fit(reml=True)
    print(fit2.summary())
except Exception as e:
    print("Model B FAILED: %s: %s" % (type(e).__name__, e))

# One-way ANOVA equivalent: delta ~ C(lineage) + C(chrom) + C(lineage):C(chrom)
print("\n=== Type II ANOVA: delta ~ C(lineage) + C(chrom) + C(lineage):C(chrom) ===")
import statsmodels.api as sm
model = smf.ols("delta ~ C(lineage) * C(chrom)", data=df)
res = model.fit()
from statsmodels.stats.anova import anova_lm
table = anova_lm(res, typ=2)
print(table)
print("\nR-squared: %.4f" % res.rsquared)
print("Adj R-squared: %.4f" % res.rsquared_adj)

# Extract key p-values
print("\n=== KEY P-VALUES ===")
for idx in table.index:
    print("  %-30s F = %8.3f  p = %.4g" % (idx, table.loc[idx, "F"], table.loc[idx, "PR(>F)"]))
