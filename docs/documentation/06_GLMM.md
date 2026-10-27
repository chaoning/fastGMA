---
layout: page
permalink: /documentation/06_GLMM/
title: GLMM
description: Generalized linear mixed model for binary traits
order: 6
---

`GLMM` (generalized linear mixed model) is the GWAS engine of `fastGMA` for **binary (case/control) traits**. It first fits a logistic **null model** with one GRM as the polygenic random effect, then tests every SNP for an **additive and a dominance effect** jointly.

It is invoked through the top-level dispatcher with `--glmm`:

```bash
fastgma --glmm [options]
```

Run `fastgma -h 5` to print the full option list.

## 1. What it does

The trait $$y$$ must be coded `0`/`1`, and both cases and controls must be present. The null model is a logistic mixed model

$$
\operatorname{logit}(\mu) = \eta = Xb + Zg, \qquad g \sim N(0, \sigma^2_g G)
$$

where $$\mu = P(y = 1)$$, $$Xb$$ are the fixed effects (intercept, covariates, class variables), $$G$$ is the GRM and $$Z$$ maps records to GRM individuals. An individual may have several records.

> **Repeated records.** The model has no permanent environmental effect: the records of one individual are treated as conditionally independent given its genetic effect $$g$$. Non-genetic causes shared by the records of an individual are therefore absorbed into $$g$$ or ignored, which can overstate the evidence. The log warns (Data section) when the data contain repeated records.

The null model is fitted by penalized quasi-likelihood (PQL). Each outer iteration forms the working response and weights

$$
\begin{aligned}
y^* &= \eta + \frac{y - \mu}{\mu(1 - \mu)} \\
W &= \operatorname{diag}\{\mu(1 - \mu)\} \\
V &= W^{-1} + \sigma^2_g Z G Z'
\end{aligned}
$$

and estimates $$\sigma^2_g$$ from the linear mixed model on $$y^*$$ with AI-REML (WEMAI: the AI matrix is blended with the EM matrix when needed to keep the variance positive). The fixed effects $$b$$ and the random effects $$g$$ are then updated and the next iteration starts. The outer loop stops when both the change in $$\sigma^2_g$$ and the change in $$b$$ are below `--cc-par0`. There is no residual variance parameter: the residual covariance is $$W^{-1}$$.

The SNP test reuses $$y^*$$ and $$V$$ from the final linearization of the null model. For each SNP, the model

$$
y^* = Xb + Z x_{\text{add}} a + Z x_{\text{dom}} d + Zg + e
$$

is fitted by generalized least squares. Both the additive effect $$a$$ and the dominance effect $$d$$ get a Wald test (chi-square, 1 df). Each effect is adjusted for the other, so `p_add` is not the p-value of an additive-only scan. Effects are on the logit (log-odds) scale. `GLMM` has no separate additive-only test.

Let $$x$$ be the count of the A1 allele, $$p$$ the A1 frequency and $$q = 1 - p$$. The additive code is $$x - 2p$$, with missing genotypes set to $$2p$$. The dominance code is the statistical coding $$0 \to -2p^2$$, $$1 \to 2pq$$, $$2 \to -2q^2$$, the same as in the dominance GRM (`--grm-type dom`). A missing genotype gets the mean dominance code of the observed genotypes. With this coding, $$a$$ is the allele substitution effect and $$d$$ the dominance deviation, as in [REMMA](../05_REMMA/). The flag `--dom` names this test; it is optional, because `GLMM` has no other test.

**Non-identifiable SNPs.** If one of the three genotype classes is absent, or the two codes are collinear after adjusting for the fixed effects, the dominance effect cannot be estimated. `beta_dom`, `se_dom` and `p_dom` are then `nan`, and the `_add` columns hold the additive-only test. The log counts these SNPs. These rules are the same as in [REMMA](../05_REMMA/).

**Small genotype classes (`--min-geno-count N`).** A dominance effect estimated from a very small genotype class can give an extreme `p_dom`. With `--min-geno-count N` (`N > 0`), the joint fit is only done when each genotype class `0`/`1`/`2` has at least `N` analysed individuals with a non-missing genotype; otherwise the SNP is treated like a non-identifiable SNP (`_dom` columns `nan`, additive-only `_add` columns). The default `0` turns the rule off. The log reports the setting and the number of SNPs affected.

The null model stores and inverts dense $$n \times n$$ matrices ($$n$$ = number of records), so memory grows with $$n^2$$ and the null-model time with $$n^3$$.

## 2. Inputs

### Phenotype / covariate table (`--data`)

A whitespace- or tab-separated table with a header row. The first column is the individual ID. It is matched against the IID column of the `.grm.id` file and of the `.fam` file. Columns are selected by name.

| Option | Meaning |
| --- | --- |
| `--data` | Phenotype / covariate table (**required**) |
| `--trait` | Column name of the binary trait (**required**). Give exactly one column. Values must be `0` or `1`. |
| `--covar` | Column name(s) of continuous covariates |
| `--class` | Column name(s) of categorical covariates |
| `--missing-data` | Tokens treated as missing in the table (default `NA NaN nan N/A ...`) |

Records are kept when the ID is in the GRM and none of the used columns is missing. Every kept ID must also be in the `.fam` file. Linearly dependent fixed-effect columns are removed automatically, and the log lists them. A trait with only cases or only controls is rejected.

### GRM input (`--grm`)

The prefix of one GRM (required), for example `--grm output/test` for `output/test.grm.bin` and `output/test.grm.id`. Build it with `fastgma --make-grm`; GCTA GRMs are also accepted. See [GRM](../03_GRM/).

### Genotype input (`--bfile`)

The PLINK binary prefix (`.bed/.bim/.fam`), required.

## 3. Outputs

All files share the `--out` prefix and are tab-separated with a header row. With `--snp-set` the prefix becomes `<out>.<start>_<end>` (0-based start, exclusive end); with `--parallel m i` it becomes `<out>.<m>_<i>`.

| File | Contents |
| --- | --- |
| `<out>.var` | Polygenic variance $$\sigma^2_g$$ with its standard error (logit scale) |
| `<out>.sln` | Fixed-effect estimates of the null model (logit scale) |
| `<out>.assoc` | One row per SNP with the association results |
| `<out>.log` | Log of the run (appended; `<out>.<m>_<i>.log` with `--parallel`, `<out>.<start>_<end>.log` with `--snp-set`) |

### `<out>.var`

The header is `component variance se`, followed by one row for the GRM component, named `G`, with the variance and its SE on the logit scale.

The SE comes from the inverse AI matrix of the working model at the final linearization. There is no `residual` row, because the working residual covariance is $$W^{-1}$$ (dispersion fixed at 1). No `.h2` file is written.

### `<out>.sln`

The columns are `type term levels effect_logit se`, as in [REMMA](../05_REMMA/) but with the estimate column named for its scale: `type` is `0` for the intercept (`mu`), `1` for a covariate and `2` for a class level; a dependent column is listed with effect `0` and se `NA`. The estimates and their SEs both come from the final linearization (the $$y^*$$, $$W$$ and $$V$$ used by the SNP test): $$\hat b = (X' V^{-1} X)^{-1} X' V^{-1} y^*$$ with SEs from $$(X' V^{-1} X)^{-1}$$.

### `<out>.assoc`

```
chr	snp	pos	a1	a2	af	N	beta_add	se_add	p_add	beta_dom	se_dom	p_dom
```

| Column | Meaning |
| --- | --- |
| `chr`, `snp`, `pos` | Chromosome, SNP ID and position from the `.bim` file |
| `a1`, `a2` | Counted (effect) allele and the other allele |
| `af` | Frequency of `a1` in the analyzed individuals |
| `N` | Number of analyzed individuals with a non-missing genotype |
| `beta_add`, `se_add`, `p_add` | Additive effect, standard error and Wald p-value |
| `beta_dom`, `se_dom`, `p_dom` | Dominance effect, standard error and Wald p-value |

The effect, SE and p-value columns are `nan` when the SNP fails the `--maf` or `--missing-rate` filter. For non-identifiable SNPs only the `_dom` columns are `nan` (see Section 1). With `--p-cut` below `1`, a row is written when `p_add` or `p_dom` is at or below the cut.

### Convergence

The log reports each PQL iteration and whether the null model converged. Two warnings can appear:

- **Outer loop not converged.** The PQL loop reached `--maxiter0` iterations. The results use the last estimates and may be unreliable; increase `--maxiter0`.
- **Inner loop not converged.** The inner AI-REML loop also stops after `--maxiter0` iterations. The log reports in how many outer iterations this happened.

## 4. Options

### Main options

| Option | Default | Meaning |
| --- | --- | --- |
| `--data` | — | Phenotype / covariate table (**required**) |
| `--trait` | — | Binary trait column, exactly one (**required**) |
| `--grm` | — | GRM prefix (**required**) |
| `--bfile` | — | PLINK binary prefix (**required**) |
| `--out` | — | Output prefix (**required**) |
| `--dom` | — | The additive + dominance test (optional; it is the only test) |
| `--covar`, `--class` | none | Covariate and class columns |
| `--missing-data` | `NA NaN nan N/A ...` | Missing-value tokens in the data file |
| `--maf` | `0.01` | Keep SNPs with $$\min(p, 1 - p)$$ at or above this, in `[0, 0.5]` |
| `--missing-rate` | `0.05` | Keep SNPs with missing rate at or below this |
| `--min-geno-count` | `0` | Minimum number of individuals in each genotype class for the dominance test; otherwise `_dom` is `nan` and the additive-only test is reported (`0` = off) |
| `--snp-set <first> <last>` | — | Test SNPs from the first to the last SNP ID (inclusive) |
| `--parallel <m> <i>` | — | Split the SNPs into `m` parts and test part `i` (1-based) |
| `--p-cut` | `1` | Write SNPs whose `p_add` or `p_dom` is at or below this; `1` writes all |
| `--threads` | `10` | Number of MKL/OpenMP threads |

`--snp-set`, `--parallel` and `--start-pos`/`--end-pos` are mutually exclusive.

### Advanced

| Option | Default | Meaning |
| --- | --- | --- |
| `--npart-snp` | `20` | Number of chunks the tested SNPs are read in (memory only; does not change results) |
| `--start-pos` | `0` | First SNP to test (0-based `.bim` index) |
| `--end-pos` | `-1` | Stop before this SNP index; `-1` (or a value past the last SNP) means the end of the file |
| `--maxiter0` | `100` | Max iterations, used for both the outer PQL loop and the inner AI-REML loop |
| `--cc-par0` | `1e-7` | Convergence criterion for the change in variance and in fixed effects |
| `--cc-gra0` | `1e-6` | Convergence criterion for the REML score in the inner loop |

Set the environment variable `FASTGMA_LOG_LEVEL=debug` to log the inner AI-REML iterations.

## 5. Examples

### Build the GRM and run the scan

```bash
fastgma --make-grm --bfile test --out output/test

fastgma --glmm --dom \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait disease \
  --covar age \
  --class sex \
  --threads 4 \
  --out output/test_glmm
```

### Split a scan into parallel jobs

Every job refits the same null model, so the results match a single full run.

```bash
for i in 1 2 3 4; do
  fastgma --glmm --dom --grm output/test --bfile test \
    --data test_simu_pheno.txt --trait disease \
    --parallel 4 $i --out output/test_glmm
done

(head -1 output/test_glmm.4_1.assoc
 tail -n +2 -q output/test_glmm.4_{1,2,3,4}.assoc) > output/test_glmm.assoc
```
