# P20 TNBC Analysis

This repository contains historical RNA-seq differential-expression analyses
for P20 triple-negative breast cancer (TNBC) cohorts, including race, obesity,
stage, age, and survival-related comparisons.

It combines three closely related directories extracted from the former
`workcodes` monorepo. Private reports and rendered notebook output are
intentionally excluded from this public repository.

## Repository layout

```text
.
├── analysis/
│   ├── p20/             # P20 cohort analyses and the source R Markdown report
│   ├── blia_vs_blis/    # BLIA/BLIS and adipose-related comparisons
│   └── obesity/         # obesity-related differential-expression analyses
├── environment.yml
└── .gitignore
```

## Environment setup

Create the Conda environment:

```bash
conda env create -f environment.yml
conda activate p20-tnbc-analysis
```

The environment provides R, Pandoc, and the CRAN/Bioconductor packages imported
by the current scripts.

## Analysis organization

The `analysis/p20/DEedgeRLimma*.R` files represent related analyses performed
for different cohort sizes or clinical comparisons. Their filenames preserve
the historical sample counts and comparison labels.

The primary methods are:

- count filtering and normalization with `edgeR`;
- linear modeling and contrasts with `limma`;
- exploratory and result visualization with `ggplot2`, `gplots`, and
  `pheatmap`;
- report source in `analysis/p20/v4.Rmd`.

These scripts are retained as separate historical analyses. They have not yet
been consolidated into a single parameterized workflow.

## Data and report policy

Raw sequencing data, count matrices containing study-level observations,
controlled clinical data, generated results, rendered HTML, and private DOCX
reports are not stored in Git. Keep controlled data in institutionally approved
storage.

Public manifests, checksums, source code, and small non-sensitive aggregate
results may be committed when appropriate.

## Portability status

This is an initial history-preserving extraction. Most scripts still contain
absolute paths under `/private/path/` and expect external count matrices,
clinical annotations, reference genes, or helper AWK scripts. Those inputs must
be represented by project-relative configuration before the analyses are
considered portable or reproducible.

Some `.shed` files are historical command records rather than directly
executable shell scripts.
