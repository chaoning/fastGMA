---
layout: page
permalink: /documentation/06_GLMM/
title: GLMM
description: Generalized linear mixed model for binary traits
order: 6
---

`GLMM` (generalized linear mixed model) is the GWAS engine of `fastGMA` for **binary (case/control) traits**. It fits a logistic **null model** with one or more GRMs as polygenic random effects, then tests every SNP for an **additive effect** with a score test. SNPs with a small p-value get a **saddlepoint approximation (SPA)**, which keeps the test calibrated when cases and controls are unbalanced.

It is invoked through the top-level dispatcher with `--glmm`:

```bash
fastgma --glmm [options]
```

Run `fastgma -h 5` to print the full option list.

## 1. What it does

### Null model

The trait $$y$$ must be coded `0`/`1`, and both cases and controls must be present. The null model is the logistic mixed model

$$
\operatorname{logit}(\mu) = \eta = Xb + Z g_1 + \cdots + Z g_K, \qquad g_k \sim N(0, \sigma^2_k G_k)
$$

where $$\mu = P(y = 1)$$, $$Xb$$ are the fixed effects (intercept, covariates, class variables), $$G_k$$ are the GRMs and $$Z$$ maps records to individuals. An individual may have several records.

> **Repeated records.** The model has no permanent environmental effect: the records of one individual are treated as conditionally independent given its genetic effects. Non-genetic causes shared by the records of an individual are absorbed into $$g$$ or ignored, which can overstate the evidence. The log warns when the data contain repeated records.

The fit starts, as in SAIGE and GMMAT, from the logistic regression without random effects (IRLS), with every $$\sigma^2_k = 1$$. It then uses penalized quasi-likelihood (PQL). Each outer iteration forms the working response and weights

$$
\begin{aligned}
y^* &= \eta + \frac{y - \mu}{\mu(1 - \mu)} \\
W &= \operatorname{diag}\{\mu(1 - \mu)\} \\
V &= W^{-1} + \textstyle\sum_k \sigma^2_k Z G_k Z'
\end{aligned}
$$

and estimates the $$\sigma^2_k$$ from the linear mixed model on $$y^*$$ with AI-REML. As in SAIGE and GMMAT, a variance that an AI step takes below `1e-5` is set to `0` (a boundary estimate, logged); a later step can take it back above `0`, which is also logged. Only if the AI step is not finite is the AI matrix blended with the EM matrix (WEMAI) to keep the variances positive. The fixed effects $$b$$ and the random effects $$g_k$$ are then updated for the next iteration. The outer loop stops when the changes in the variances and in $$b$$ are both below `--cc-par0`. There is no residual variance: the working residual covariance is $$W^{-1}$$.

The null model stores and inverts dense $$n \times n$$ matrices ($$n$$ = number of records), so memory grows with $$n^2$$ and time with $$n^3$$.

### SNP test

The test uses $$y^*$$, $$\mu$$ and $$V$$ from the final linearization of the null model. Let $$P = V^{-1} - V^{-1}X(X'V^{-1}X)^{-1}X'V^{-1}$$, $$S = Z'PZ$$ and $$s = Z'Py^*$$, and let $$x$$ be the centred SNP code of the individuals (the A1 count minus $$2p$$, with missing genotypes set to $$2p$$). The score of the SNP is $$T = x's$$ with variance $$x'Sx$$; at convergence $$Py^* = y - \mu$$, so $$T = x'Z'(y - \mu)$$ is the score test of GMMAT (`glmm.score`). `GLMM` reports

$$
\hat\beta = \frac{x's}{x'Sx}, \qquad \operatorname{se}(\hat\beta) = \frac{1}{\sqrt{x'Sx}}, \qquad p_\text{normal} = P\left(\chi^2_1 > \frac{T^2}{x'Sx}\right)
$$

$$\hat\beta$$ is the one-step (score) estimate of the effect on the logit (log-odds) scale. It is close to the maximum-likelihood estimate for small effects.

**Saddlepoint approximation.** With few cases (or few controls), the score of a rare variant is far from normal, and $$p_\text{normal}$$ can be much too small. For SNPs with $$p_\text{normal}$$ below `--spa-cut` (default `0.05`), `GLMM` computes $$p$$ by the saddlepoint approximation of SAIGE (Zhou et al. 2018). The score is written as $$T = \sum_r \tilde g_r (y_r - \mu_r)$$ over the records, with the record codes adjusted for the covariates, $$\tilde g = g - X(X'WX)^{-1}X'Wg$$, and the $$y_r$$ treated as independent Bernoulli($$\mu_r$$). $$T$$ is rescaled from its variance $$x'Sx$$ (which includes the GRMs) to the Bernoulli variance $$\tilde g' W \tilde g$$; SAIGE uses one variance ratio estimated from random markers, while `GLMM` uses the exact ratio of each SNP. For these SNPs, `se` is recomputed from the SPA p-value, $$\operatorname{se} = |\hat\beta| / \sqrt{\chi^2_{1, p}}$$, so that `beta`, `se` and `p` agree. For all other SNPs, $$p = p_\text{normal}$$. `--spa-cut 0` turns the approximation off and `--spa-cut 1` applies it to every SNP.

In a simulation with 3,000 individuals, 2% cases and 4,000 null rare variants (MAF 0.002–0.02), the normal approximation gave 21 SNPs with p < 1e-3 and 9 with p < 1e-4 (expected 4 and 0.4); SPA gave 5 and 0. With 50% cases the two agreed.

## 2. Inputs

### Phenotype / covariate table (`--data`)

A whitespace- or tab-separated table with a header row. The first column is the individual ID. It is matched against the IID column of the `.grm.id` files and of the `.fam` file. Columns are selected by name.

| Option | Meaning |
| --- | --- |
| `--data` | Phenotype / covariate table (**required**) |
| `--trait` | Column name of the binary trait (**required**). Give exactly one column. Values must be `0` or `1`. |
| `--covar` | Column name(s) of continuous covariates |
| `--class` | Column name(s) of categorical covariates |
| `--missing-data` | Tokens treated as missing in the table (default `NA NaN nan N/A ...`) |

Records are kept when the ID is in the GRMs and none of the used columns is missing. Every kept ID must also be in the `.fam` file. Linearly dependent fixed-effect columns are removed automatically, and the log lists them. A trait with only cases or only controls is rejected.

### GRM input (`--grm`, `--mgrm`)

Give one of:

| Option | Meaning |
| --- | --- |
| `--grm` | One GRM prefix, for example `--grm output/test` for `output/test.grm.bin` and `output/test.grm.id`. The component is named `G` |
| `--mgrm` | One or more GRM prefixes, one random effect each (for example additive and dominance GRMs). The components are named after the GRMs, as in [REMMA](../05_REMMA/) |

All GRMs must list the same individuals in the same order. Build them with `fastgma --make-grm`; GCTA GRMs are also accepted. See [GRM](../03_GRM/).

### Genotype input (`--bfile`, `--ffile`)

| Option | Meaning |
| --- | --- |
| `--bfile` | PLINK binary genotype prefix (`.bed/.bim/.fam`) |
| `--ffile` | Feature-matrix prefix (`.fam/.bim/.fmat`). Each `.fmat` row is one feature and each column one sample, for imputed dosages or any other continuous per-sample matrix |
| `--missing-geno` | Tokens treated as missing in the `.fmat` file (default `NA NaN nan N/A ...`) |

Give exactly one of `--bfile` and `--ffile`. A `--ffile` feature is centred by its mean, missing values get the mean, `--maf` is not applied, and `af`, `af_case` and `af_ctrl` are feature means.

## 3. Outputs

All files share the `--out` prefix and are tab-separated with a header row. With `--snp-set` the prefix becomes `<out>.<start>_<end>` (0-based start, exclusive end); with `--parallel m i` it becomes `<out>.<m>_<i>`.

| File | Contents |
| --- | --- |
| `<out>.var` | Polygenic variance(s) with standard errors (logit scale) |
| `<out>.sln` | Fixed-effect estimates of the null model (logit scale); with `--predict`, also the predicted genetic values (BLUP) and their reliabilities |
| `<out>.res` | With `--predict`: the linear predictor of each record split into its fixed and random parts, and the fitted probability |
| `<out>.assoc` | One row per SNP with the association results |
| `<out>.log` | Log of the run (appended; `<out>.<m>_<i>.log` with `--parallel`, `<out>.<start>_<end>.log` with `--snp-set`) |

### `<out>.var`

The header is `component variance se`, followed by one row per GRM component (`G` with `--grm`) with the variance and its SE on the logit scale. The SE comes from the inverse AI matrix of the working model at the final linearization; a variance on the boundary (`0`) has the SE `NA`. There is no `residual` row, because the working residual covariance is $$W^{-1}$$ (dispersion fixed at 1). No `.h2` file is written.

### `<out>.sln`

The columns are `type term levels effect_logit se`, and with `--predict` also `rel`, as in [uvLMM](../04_UVLMM/) but with the estimate column named for its scale: `type` is `0` for the intercept (`mu`), `1` for a covariate and `2` for a class level; a dependent column is listed with effect `0` and se `NA`. The estimates and their SEs come from the final linearization (the $$y^*$$, $$W$$ and $$V$$ of the SNP test): $$\hat b = (X' V^{-1} X)^{-1} X' V^{-1} y^*$$ with SEs from $$(X' V^{-1} X)^{-1}$$.

With `--predict`, type `3` rows hold the predicted genetic values of every GRM individual (including those without records, predicted through their relationships), one block per GRM. `term` is the component as named in `<out>.var` and `levels` the individual ID. As in uvLMM, they come from the working model of the final linearization: $$\hat g_k = \sigma^2_k G_k Z' P y^*$$ on the logit scale, `se` is the square root of the prediction error variance $$\mathrm{PEV}_k = \sigma^2_k \operatorname{diag}(G_k) - \sigma^4_k \operatorname{diag}(G_k Z'PZ G_k)$$, and `rel` is the reliability $$1 - \mathrm{PEV} / (\sigma^2_k G_{jj})$$ (`NA` on the fixed-effect rows and for a variance on the boundary). Numbers keep 10 significant digits.

### `<out>.res`

With `--predict`, one row per record: `id observed fixed random mu`. `observed` is the `0`/`1` trait, `fixed` is $$x_i'\hat b$$ and `random` the sum of the predicted genetic values of the record, both on the logit scale, and `mu` $$= 1 / (1 + e^{-(\text{fixed} + \text{random})})$$ is the fitted probability. At convergence the mean of `mu` equals the fraction of cases (the score equation of the intercept).

### `<out>.assoc`

```
chr	snp	pos	a1	a2	af	N	N_case	N_ctrl	af_case	af_ctrl	beta	se	p_normal	p
```

| Column | Meaning |
| --- | --- |
| `chr`, `snp`, `pos` | Chromosome, SNP ID and position from the `.bim` file |
| `a1`, `a2` | Counted (effect) allele and the other allele |
| `af` | Frequency of `a1` in the analysed individuals |
| `N` | Number of analysed individuals with a non-missing genotype |
| `N_case`, `N_ctrl` | Number of case and control records with a non-missing genotype |
| `af_case`, `af_ctrl` | Frequency of `a1` over those case and control records |
| `beta`, `se` | Additive effect (logit scale) and its standard error (from the SPA p-value when SPA was applied) |
| `p_normal` | p-value of the score test by the normal (chi-square) approximation |
| `p` | Final p-value: SPA when `p_normal` is below `--spa-cut`, otherwise `p_normal` |

The effect, SE and p-value columns are `nan` when the SNP fails the `--maf` or `--missing-rate` filter or gives a non-finite statistic. If the saddlepoint equation cannot be solved for a SNP, `p` keeps `p_normal`; the log counts these SNPs. With `--p-cut` below `1`, a row is written when `p` is at or below the cut.

### Convergence

The log reports each PQL iteration and whether the null model converged. Two warnings can appear:

- **Outer loop not converged.** The PQL loop reached `--maxiter0` iterations. The results use the last estimates and may be unreliable; increase `--maxiter0`.
- **Inner loop not converged.** The inner AI-REML loop also stops after `--maxiter0` iterations. The log reports in how many outer iterations this happened.

With few cases the polygenic variance is often estimated at the boundary: in a simulation with 3,000 individuals and 2% cases, `G` was set to `0` in the first iteration and the null model converged after 2 iterations. With `G = 0` the null model is the logistic regression, and the SNP test is its score test (with SPA).

## 4. Options

### Main options

| Option | Default | Meaning |
| --- | --- | --- |
| `--data` | — | Phenotype / covariate table (**required**) |
| `--trait` | — | Binary trait column, exactly one (**required**) |
| `--grm` / `--mgrm` | — | One GRM prefix, or one or more GRM prefixes; exactly one is **required** |
| `--bfile` / `--ffile` | — | PLINK binary prefix, or feature-matrix prefix; exactly one is **required** |
| `--out` | — | Output prefix (**required**) |
| `--covar`, `--class` | none | Covariate and class columns |
| `--missing-data` | `NA NaN nan N/A ...` | Missing-value tokens in the data file |
| `--missing-geno` | `NA NaN nan N/A ...` | Missing-value tokens in the `.fmat` file |
| `--maf` | `0.01` | Keep SNPs with $$\min(p, 1 - p)$$ at or above this, in `[0, 0.5]` (`--bfile` only) |
| `--missing-rate` | `0.05` | Keep SNPs with missing rate at or below this |
| `--spa-cut` | `0.05` | Saddlepoint approximation for SNPs with `p_normal` below this; `0` = off, `1` = every SNP |
| `--predict` | off | Also write the predicted genetic values and their reliabilities to `<out>.sln`, and `<out>.res` |
| `--snp-set <first> <last>` | — | Test SNPs from the first to the last SNP ID (inclusive) |
| `--parallel <m> <i>` | — | Split the SNPs into `m` parts and test part `i` (1-based) |
| `--p-cut` | `1` | Write SNPs whose `p` is at or below this; `1` writes all |
| `--threads` | `10` | Number of MKL/OpenMP threads |

`--snp-set`, `--parallel` and `--start-pos`/`--end-pos` are mutually exclusive.

### Advanced

| Option | Default | Meaning |
| --- | --- | --- |
| `--npart-snp` | `20` | Minimum number of chunks the tested SNPs are read in; chunks are kept to about 128 MB of genotype codes (memory only; does not change results) |
| `--start-pos` | `0` | First SNP to test (0-based `.bim` index) |
| `--end-pos` | `-1` | Stop before this SNP index; `-1` (or a value past the last SNP) means the end of the file |
| `--maxiter0` | `100` | Max iterations, used for both the outer PQL loop and the inner AI-REML loop |
| `--cc-par0` | `1e-7` | Convergence criterion for the change in variances and in fixed effects |
| `--cc-gra0` | `1e-6` | Convergence criterion for the REML score in the inner loop |

Set the environment variable `FASTGMA_LOG_LEVEL=debug` to log the inner AI-REML iterations.

## 5. Examples

### Build the GRM and run the scan

```bash
fastgma --make-grm --bfile test --out output/test

fastgma --glmm \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait disease \
  --covar age \
  --class sex \
  --threads 4 \
  --out output/test_glmm
```

### Additive and dominance GRMs

```bash
fastgma --make-grm --grm-type dom --bfile test --out output/test_dom

fastgma --glmm \
  --mgrm output/test output/test_dom \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait disease \
  --out output/test_glmm_mgrm
```

### Split a scan into parallel jobs

Every job refits the same null model, so the results match a single full run.

```bash
for i in 1 2 3 4; do
  fastgma --glmm --grm output/test --bfile test \
    --data test_simu_pheno.txt --trait disease \
    --parallel 4 $i --out output/test_glmm
done

(head -1 output/test_glmm.4_1.assoc
 tail -n +2 -q output/test_glmm.4_{1,2,3,4}.assoc) > output/test_glmm.assoc
```
