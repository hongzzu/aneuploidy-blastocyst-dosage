# How to publish (I cannot run network commands in this session)

## GitHub
1. Create a repository named `aneuploidy-blastocyst-dosage` under your account.
2. From the `zenodo_stage` parent folder:
       cd aneuploidy-blastocyst-dosage
       git init
       git add .
       git commit -m "Analysis code for within-embryo dosage study"
       git branch -M main
       git remote add origin https://github.com/<username>/aneuploidy-blastocyst-dosage.git
       git push -u origin main

## Zenodo
3. Sign in at zenodo.org with your GitHub account and enable the
   `aneuploidy-blastocyst-dosage` repository in Account > GitHub.
4. Cut a GitHub release (e.g. v1.0.0). Zenodo archives it automatically and mints a DOI.
5. Copy the DOI (10.5281/zenodo.XXXXXX) into the Code availability statement below.

## Then update the manuscript
Replace the Code availability statement in `论文2_终稿_v3.md` (3 places are consistent)
with:

    Analysis scripts are available at https://github.com/<username>/aneuploidy-blastocyst-dosage
    and archived at https://doi.org/10.5281/zenodo.<DOI>.

`<username>` and `<DOI>` are the only two placeholders.
