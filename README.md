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

The public history has been anonymized. Personal paths, host names, commit
emails, internal sample identifiers, individual-level examples, and historical
command-record files have been removed.

The analysis scripts obtain private input locations from environment variables:

- `P20_DATA_DIR`: directory containing controlled P20 input tables;
- `BLIA_BLIS_DATA_DIR`: directory containing BLIA/BLIS inputs;
- `OBESITY_DATA_DIR`: directory containing obesity-analysis inputs;
- `CYTOBAND_FILE`: cytoband annotation file;
- `P20_BILATERAL_EXCLUSIONS_FILE`: one controlled sample ID per line;
- `P20_DISCORDANT_SAMPLES_FILE`: one controlled sample ID per line.

Unset optional exclusion-list variables are treated as empty lists. These
private files must remain outside the repository and in institutionally
approved storage.
