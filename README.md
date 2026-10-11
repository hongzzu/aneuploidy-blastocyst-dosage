# Chromosome-specific dosage effects in human aneuploid blastocysts

Analysis code for:

> **Within-embryo pseudobulk reveals chromosome-specific dosage effects in human aneuploid blastocysts**

## Contents

| Script | Purpose |
|---|---|
| `dosage_specificity.py` | Layer 3 core: chromosome-mismatch specificity contrast (matched vs mismatched delta) |
| `taskC_permutation.py` | Permutation confirmation of the specificity contrast (10,000 draws; simple + stratified) |
| `taskB_rank_gsea.py` | Rank-based GSEA of the per-gene dosage ranking (MSigDB Hallmark / KEGG / Reactome) |
| `task2_enrichment.py` | Over-representation analysis (superseded by rank-GSEA; kept for provenance) |
| `step2_gap_analysis.py` | Per-stratum dosage compensation gap and lineage comparisons |
| `step4_mixed_model_py.py` | Mixed-effects model (delta ~ chromosome + lineage) |
| `step4a_interaction.py` | Interaction-term subset analysis |
| `mr_power.py` | Post-hoc MR power calculation (mRnd model, Brion et al. 2013) |
| `taskD_twas_intersection.py` | TWAS intersection between aneuploidy and crossover-count phenotypes |
| `audit_all.py` | Whole-manuscript consistency audit (word counts, panel citations, key numbers) |
| `check_markers.py` | Markdown-artifact detector for rendered DOCX files |
| `finalise_authoryear.py` | Numeric-to-author-year citation conversion and alphabetical reference sort |
| `merge_tables.py` | Assembly of Tables 1 and 2 into a single Word document |

## Data

All input data are public:

- Maternal meiotic aneuploidy GWAS and TWAS summary statistics — Carioscia et al., *Nature* (2026)
- GTEx v8 ovary cis-eQTL — GTEx Consortium (2020)
- eQTLGen whole-blood cis-eQTL — Vosa et al. (2021)
- Blastocyst single-cell RNA-seq count matrix and metadata — Wang et al., *Nat Genet* (2024)

Accession numbers are listed in Table 1 of the manuscript.

## Environment

Python 3.12.14 with NumPy, pandas, SciPy, statsmodels, GSEApy 1.3.1.
Full software versions are in Supplementary Methods S6 of the manuscript.

## Note

`dosage_specificity.py` and `taskB_rank_gsea.py` read the blastocyst count matrix
(`embryo.14082cells.raw.counts.normalized.edgeR.ERCC.pseudo.counts.txt.gz`) from the
path encoded in each script; edit the `WS` / `MAT` constants to point at your local copy.
