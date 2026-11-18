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

## Citation

If you use fastGMA, please cite:

- Dan Wang, Hui Tang, Jian-Feng Liu, Shizhong Xu, Qin Zhang, Chao Ning. Rapid epistatic mixed-model association studies by controlling multiple polygenic effects. *Bioinformatics* 36(19):4833–4837, 2020. [doi:10.1093/bioinformatics/btaa610](https://doi.org/10.1093/bioinformatics/btaa610)

If you use the modules below, please also cite the papers on their methods:

- `REMMA`:
  - Dan Wang, Hui Tang, Jian-Feng Liu, Shizhong Xu, Qin Zhang, Chao Ning. Rapid epistatic mixed-model association studies by controlling multiple polygenic effects. *Bioinformatics* 36(19):4833–4837, 2020. [doi:10.1093/bioinformatics/btaa610](https://doi.org/10.1093/bioinformatics/btaa610)
  - Chao Ning, Dan Wang, Huimin Kang, Raphael Mrode, Lei Zhou, Shizhong Xu, Jian-Feng Liu. A rapid epistatic mixed-model association analysis by linear retransformations of genomic estimated values. *Bioinformatics* 34(11):1817–1825, 2018. [doi:10.1093/bioinformatics/bty017](https://doi.org/10.1093/bioinformatics/bty017)
- `lonGWAS`:
  - Chao Ning, Dan Wang, Lei Zhou, Julong Wei, Yuanxin Liu, Huimin Kang, Shengli Zhang, Xiang Zhou, Shizhong Xu, Jian-Feng Liu. Efficient multivariate analysis algorithms for longitudinal genome-wide association studies. *Bioinformatics* 35(23):4879–4885, 2019. [doi:10.1093/bioinformatics/btz304](https://doi.org/10.1093/bioinformatics/btz304)
  - Chao Ning, Dan Wang, Xianrui Zheng, Qin Zhang, Shengli Zhang, Raphael Mrode, Jian-Feng Liu. Eigen decomposition expedites longitudinal genome-wide association studies for milk production traits in Chinese Holstein. *Genetics Selection Evolution* 50:12, 2018. [doi:10.1186/s12711-018-0383-0](https://doi.org/10.1186/s12711-018-0383-0)
  - Chao Ning, Huimin Kang, Lei Zhou, Dan Wang, Haifei Wang, Aiguo Wang, Jinluan Fu, Shengli Zhang, Jianfeng Liu. Performance gains in genome-wide association studies for longitudinal traits via modeling time-varied effects. *Scientific Reports* 7:590, 2017. [doi:10.1038/s41598-017-00638-2](https://doi.org/10.1038/s41598-017-00638-2)

## Contact

If you have questions, please open an issue on [GitHub](https://github.com/chaoning/fastGMA/issues) or email ningchao91@gmail.com.
