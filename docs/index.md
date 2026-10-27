---
layout: full
homepage: true
disable_anchors: true
description: Fast genomic mixed-model analysis for GWAS and GP of complex traits
---

## fastGMA

`fastGMA` is a C++ toolkit for genomic mixed-model analysis of complex traits. It provides efficient methods for:

- genome-wide association studies (GWAS),
- variance component and heritability estimation, and
- genomic prediction.

## Modules

| Mode | Purpose | Guide |
| --- | --- | --- |
| `--make-grm`, `--process-grm` | Build additive, dominance and epistatic GRMs; merge, convert, invert and decompose GRMs | [GRM](documentation/03_GRM/) |
| `--uvlmm` | Univariate linear mixed model with one or more GRMs | [uvLMM](documentation/04_UVLMM/) |
| `--remma` | Rapid linear mixed model for non-additive analysis | [REMMA](documentation/05_REMMA/) |
| `--glmm` | Generalized linear mixed model for binary traits | [GLMM](documentation/06_GLMM/) |
| `--mvlmm` | Multivariate linear mixed model for several correlated traits | [mvLMM](documentation/07_MVLMM/) |
| `--longwas` | Longitudinal data analysis with random regression models | [lonGWAS](documentation/08_LONGWAS/) |
| `--fast-reml` | REML estimation of genetic, GxE and heteroscedastic noise variances | [fastREML](documentation/09_FASTREML/) |

`fastgma -h` lists all modes, and `fastgma -h <n>` prints the full option list of mode `n`.

## Inputs and outputs

- **Genotypes:** PLINK binary files (`--bfile <prefix>` for `.bed`, `.bim`, `.fam`) or a feature matrix (`--ffile <prefix>` for `.fmat`, `.bim`, `.fam`), for every module except fastREML, which reads only GRMs.
- **Phenotypes and covariates:** a whitespace-separated table with a header row and the individual ID in the first column. Traits, covariates and class variables are selected by column name.
- **GRMs:** three files that share one prefix, given as `--grm <prefix>`:
  - `<prefix>.grm.bin`: the lower triangle of the GRM including the diagonal, stored row by row ($$A_{11}, A_{21}, A_{22}, A_{31}, \ldots$$) as little-endian float32.
  - `<prefix>.grm.id`: one line per individual (`FID IID`), in the row order of the matrix.
  - `<prefix>.grm.N.bin`: the number of SNPs used for each pair of individuals, in the same layout as `.grm.bin`.

  GRMs are built with `fastgma --make-grm` (see [GRM](documentation/03_GRM/)). This is the GCTA binary GRM format, so a GRM built by `fastgma` can be used directly in GCTA, and a GRM built by GCTA can be used directly in every `fastgma` module.

## Quick start

```bash
# 1. Additive GRM
fastgma --make-grm --bfile test --out test

# 2. Mixed-model GWAS with heritability estimation
fastgma --uvlmm --grm test --bfile test \
  --data pheno.txt --trait height --covar age --class sex \
  --out height
```

See [Installation](documentation/02_installation/) to download `fastgma`.

## Contact

If you have questions, please open an issue on [GitHub](https://github.com/chaoning/fastGMA/issues) or email ningchao91@gmail.com.
