---
layout: page
permalink: /documentation/05_REMMA/
title: REMMA
description: Additive, dominance and epistasis GWAS with a multi-GRM random effect model
order: 5
---

`REMMA` (Random Effect Model for Multiple Associations) tests SNPs for additive, dominance and additive-by-additive epistatic effects. It first fits a **null model** with one or more GRMs as random effects (for example additive, dominance and epistatic GRMs). It then runs one per-SNP or per-SNP-pair test using the variance components from that null model.

It is invoked through the top-level dispatcher with `--remma` (there is no short form):

```bash
fastgma --remma <test> [options]
```

For the full option list, run:

```bash
fastgma -h 4
```

## 1. What it does

The null model is

$$
\begin{aligned}
y &= Xb + Z u_1 + \cdots + Z u_K + e \\
u_k &\sim N(0, \sigma^2_k G_k), \qquad e \sim N(0, \sigma^2_e I) \\
V &= \textstyle\sum_k \sigma^2_k Z G_k Z' + \sigma^2_e I
\end{aligned}
$$

$$Xb$$ holds the fixed effects (intercept, covariates, class variables). Each $$u_k$$ is a random effect with covariance matrix $$G_k$$, one per GRM given with `--mgrm`. $$Z$$ maps records to individuals, so repeated records per individual are allowed and are used by every test. The variance components are estimated once by AI-REML (WEMAI: the AI matrix is blended with the EM matrix when needed to keep the variances positive). The start values are scale-aware: the phenotypic variance is split evenly over the $$K$$ GRMs and the residual. `--solver general` (default) fits on the dense $$n \times n$$ $$V$$ of the $$n$$ records; `--solver mme` fits the same model on the mixed model equations, whose size does not grow with the number of records, and is much faster with many repeated records (see [uvLMM](../04_UVLMM/), Section 2). Both give the same estimates and tests. `mme` adds `--val` (default `0.001`) to each GRM diagonal so that the GRM can be inverted; give `general` the same `--val` to fit the identical model. All tests below reuse $$V$$ from this fit, and the variances are not re-estimated for each SNP.

Let $$P = V^{-1} - V^{-1} X (X' V^{-1} X)^{-1} X' V^{-1}$$, $$S = Z'PZ$$ and $$s = Z'Py$$. Genotypes are coded per individual as the number of copies of the A1 allele (`0/1/2`). Missing genotypes are replaced by the SNP mean $$2p$$ for the additive code.

## 2. Tests

Choose exactly one test. The run stops with an error if none or more than one is given.

| Flag | Test | Model per SNP / pair |
| --- | --- | --- |
| `--add` | Additive | Exact GLS statistic with the null-model `V` |
| `--dom` | Dominance | Joint GLS fit: fixed effects + additive + dominance term |
| `--epis` | Additive-by-additive epistasis | Approximate test for every SNP pair, exact re-test below `--p-approx-cut` |
| `--epis --marginal` | Additive-by-additive epistasis | Full GLS model per pair with both main effects (slower) |

`--marginal` is only valid together with `--epis`.

### Additive test (`--add`)

Let $$x$$ be the centered SNP code (the dosage minus $$2p$$). With $$V$$ fixed, the GLS estimate of the SNP effect and its standard error are

$$
\hat\beta = \frac{x's}{x'Sx}, \qquad \operatorname{se}(\hat\beta) = \sqrt{\frac{1}{x'Sx}}, \qquad \frac{\hat\beta^2}{\operatorname{se}(\hat\beta)^2} \sim \chi^2_1
$$

This is the exact GLS statistic for $$y = Xb + Zxa + \cdots + e$$ at the null $$V$$. It needs only $$s$$ and $$S$$, so it is very fast.

### Dominance test (`--dom`)

For each SNP, `REMMA` fits

$$
y = Xb + Z x_a a + Z x_d d + Z u_1 + \cdots + Z u_K + e
$$

by generalized least squares with the fixed $$V$$ from the null model. $$x_a$$ is the centered additive code $$x - 2p$$. $$x_d$$ is the statistical dominance code $$0 \to -2p^2$$, $$1 \to 2pq$$, $$2 \to -2q^2$$, the same as in the dominance GRM (`--grm-type dom`) (see [GRM](../03_GRM/)). Here $$p$$ is the A1 frequency and $$q = 1 - p$$. A missing genotype gets the mean dominance code of the observed genotypes, so it carries no information (like the mean-imputed additive code). Wald tests (chi-square, 1 df) are reported for both $$a$$ and $$d$$ from the same fit.

With this coding, $$a$$ is the allele substitution effect, the quantity estimated by an ordinary additive GWAS, and $$d$$ is the dominance deviation. The genotypic additive effect (half the difference between the two homozygotes) is $$a - d(q - p)$$.

**Non-identifiable SNPs.** If one of the three genotype classes is absent among the analysed individuals, or the two codes are collinear after adjusting for the fixed effects, the dominance effect cannot be estimated. `REMMA` then reports `beta_dom`, `se_dom` and `p_dom` as `nan`, and the `_add` columns hold the additive-only test (the same statistic as `--add`). The log counts these SNPs.

**Small genotype classes (`--min-geno-count N`).** A dominance effect estimated from a very small genotype class (for example a few minor-allele homozygotes) can give an extreme `p_dom`. With `--min-geno-count N` (`N > 0`), the joint fit is only done when each genotype class `0`/`1`/`2` has at least `N` analysed individuals with a non-missing genotype; otherwise the SNP is treated exactly like a non-identifiable SNP (`_dom` columns `nan`, `_add` columns from the additive-only test). The default `0` turns the rule off (only an absent class triggers the fallback). The log reports the setting and how many tested SNPs it affected. The rule applies to `--dom` only: `--epis` uses additive codes and has no dominance term, so it ignores the option (with a warning).

### Epistasis test (`--epis`)

For a SNP pair $$i, j$$ the interaction code is $$z = x_i \circ x_j$$, the elementwise product of the two centered additive codes, centered again by its mean. The test uses

$$
\hat\beta = \frac{z's}{z'Sz}, \qquad \operatorname{se}(\hat\beta) = \sqrt{\frac{1}{z'Sz}}
$$

Computing $$z'Sz$$ for every pair is slow, so `REMMA` approximates it by $$z'Sz \approx \gamma\, z'z$$:

1. It draws `--num-random-snp` random SNP pairs from all SNPs in `--bfile`, with `--seed` (default `42`), and computes $$z'Sz / z'z$$ for each.
2. $$\gamma$$ is the mean of this ratio over the pairs that pass the MAF and missing-rate filters and give a finite, positive ratio. The ratio does not depend on the phenotype, so no p-value filter is applied (as in `UVLMM` and `MVLMM`). At least 10% of the random pairs must be usable, or the run stops.
3. Every pair is tested with the approximate statistic. Pairs with approximate p below `--p-approx-cut` are recomputed with the exact $$z'Sz$$.

Because the random pairs come from the whole genotype file and the seed is fixed, every job of a split scan uses the same $$\gamma$$, and reruns give identical results.

With `--marginal`, each pair is instead fitted by GLS as

$$
y = Xb + Z x_i a_i + Z x_j a_j + Z z c + Z u_1 + \cdots + Z u_K + e
$$

and the Wald test of $$c$$ is reported. If the two SNPs are in complete LD, so their main effects are collinear, the redundant main effect $$x_j$$ is dropped and $$[X, x_i, z]$$ is fitted. No $$\gamma$$ is estimated, and `--p-approx-cut`, `--num-random-snp` and `--seed` are not used. This is much slower.

**Which pairs are tested.** A run (or a `--snp-set` / `--parallel` / `--start-pos`/`--end-pos` job) tests every pair $$(i, j)$$ with $$j > i$$ and SNP $$i$$ in the tested range; $$j$$ runs up to the last SNP in `--bfile`. Jobs that cover the SNPs once therefore test every pair exactly once. `--parallel m i` splits the SNPs so that each part holds about the same number of pairs (early parts hold fewer SNPs). The tested range is held in memory.

## 3. Inputs

### Phenotype / covariate table (`--data`)

A whitespace- or tab-separated table with a header row. The first column must be the individual ID; it is matched to the IID of the `.grm.id` files and of the `.fam` file. Columns are selected by name. Only individuals present in both the table and the GRMs are used, and they must all be in the `.fam` file. An individual may have several records.

| Option | Meaning |
| --- | --- |
| `--data` | Phenotype/data file (**required**) |
| `--trait` | Column name of the trait — exactly one (**required**) |
| `--covar` | Column name(s) of continuous covariates |
| `--class` | Column name(s) of categorical covariates |
| `--missing-data` | Tokens treated as missing in the table (default `NA NaN nan N/A ...`) |

Records with a missing value in any used column are removed.

### GRM input

| Option | Meaning |
| --- | --- |
| `--mgrm` | One or more GRM prefixes (**required**). Each GRM is one random effect in the null model. |

All GRMs must list the same individuals in the same order. GRMs are referenced by prefix (`--mgrm output/add` reads `output/add.grm.bin` and `output/add.grm.id`). GRMs built by `fastgma --make-grm` and by GCTA are both accepted. A typical set is:

| Random effect | How to build it |
| --- | --- |
| Additive | `fastgma --make-grm --grm-type add ` |
| Dominance | `fastgma --make-grm --grm-type dom` |
| Additive × additive | `fastgma --process-grm --make-epis --prefix1 <add> --prefix2 <add>` |

See [GRM](../03_GRM/) for details. A single additive GRM is also allowed.

### Genotype input

| Option | Meaning |
| --- | --- |
| `--bfile` | PLINK binary genotype prefix (`.bed/.bim/.fam`) (**required**) |

Only PLINK binary input is supported.

## 4. Outputs

All files share the `--out` prefix and are tab-separated with a header row. With `--snp-set` the prefix becomes `<out>.<start>_<end>` (0-based start, exclusive end); with `--parallel m i` it becomes `<out>.<m>_<i>`. The suffix applies to every file of that job.

| File | Written by | Contents |
| --- | --- | --- |
| `<out>.var` | all tests | Null-model variance components with standard errors (inverse AI matrix at the estimates) |
| `<out>.h2` | all tests | Per-component proportion of variance $$\sigma^2_k d_k / V_p$$ ($$d_k$$ = mean GRM diagonal), SE by the delta method |
| `<out>.sln` | all tests | Fixed-effect estimates of the null model |
| `<out>.assoc` | all tests | Association results |
| `<out>.random.assoc` | `--epis` without `--marginal` | Results for the random pairs used to estimate $$\gamma$$ |
| `<out>.log` | all tests | Log of the run (appended; `<out>.<m>_<i>.log` with `--parallel`, `<out>.<start>_<end>.log` with `--snp-set`) |

### Null-model files

`<out>.var` and `<out>.h2` have one row per GRM, named by the basename of its `--mgrm` prefix, and `residual` last. `<out>.h2` adds a `total` row: the summed proportion of all GRM components.

The proportions put each GRM variance on the phenotypic scale by the mean diagonal $$d_k$$ of its kernel $$Z G_k Z'$$ over the analysed records (the mean of the diagonal elements $$G_k(j,j)$$ over the individuals with records, each counted once per record):

$$
h^2_k = \frac{\sigma^2_k d_k}{V_p}, \qquad V_p = \sum_k \sigma^2_k d_k + \sigma^2_e
$$

so GRMs of different scale (for example an $$A \circ A$$ epistatic GRM, whose mean diagonal can be far from 1) are comparable. This is the definition used by [fastREML](../09_FASTREML/) and [UVLMM](../04_UVLMM/). The $$d_k$$ are written to the log. SEs are from the delta method on the inverse AI matrix. `<out>.var` holds the unscaled variances, so $$\sigma^2_k / \sum \sigma^2$$ (the convention of GCTA `--reml`) can be computed from it.

The example output below is from the mouse data (`example/mouse/geno/plink`, `example/mouse/pheno/pheno_simu.txt`, `--trait trait --covar e1 --class c1`) with the GRMs `add` (`--grm-type add`), `dom` (`--grm-type dom`) and `axa` ($$\text{add} \circ \text{add}$$):

```
component	variance	se
add	1.15076	0.214824
dom	0.00432849	0.0642098
axa	0.0225108	0.186052
residual	3.20265	0.197495
```

```
component	h2	se
add	0.265346	0.0396146
dom	0.000985636	0.014621
axa	0.00531368	0.0439219
residual	0.728355	0.0492958
total	0.271645	0.0492958
```

(Mean diagonals in this example: `add` 1.0139, `dom` 1.00126, `axa` 1.03794.)

`<out>.sln` has the columns `type term levels effect se`. `type` is `0` for the intercept (`mu`), `1` for a covariate and `2` for a class level. A class level removed because it is linearly dependent (the reference level) is shown with effect `0` and se `NA`. Numbers keep 10 significant digits.

```
type	term	levels	effect	se
0	mu	1	-0.353782	0.0912941
1	e1	1	-0.677201	0.0512313
2	c1	0	0.643054	0.128113
2	c1	1	1.34628	0.128405
2	c1	2	0	0
```

### `--add`

```
chr snp pos a1 a2 af N beta se p
```

### `--dom`

```
chr snp pos a1 a2 af N beta_add se_add p_add beta_dom se_dom p_dom
```

Columns: `chr`, `snp`, `pos`, `a1` and `a2` come from the `.bim` file; `af` is the A1 allele frequency in the analysed individuals and `N` the number of analysed individuals with a non-missing genotype. `beta`/`se`/`p` are the additive results. `beta_add`/`se_add`/`p_add` and `beta_dom`/`se_dom`/`p_dom` are the additive and dominance terms from the joint fit (see Section 2 for non-identifiable SNPs).

With the default `--p-cut 1`, every SNP in the tested range gets one row. SNPs that fail the MAF or missing-rate filter, or give a non-finite statistic, report `nan` for the effect, SE and p-value.

### `--epis`

```
chr_1 snp_1 pos_1 a1_1 a2_1 af_1 N_1 chr_2 snp_2 pos_2 a1_2 a2_2 af_2 N_2 beta se p
```

The first seven columns describe the first SNP of the pair, the next seven the second SNP, and `beta`/`se`/`p` the interaction effect. `<out>.random.assoc` has the same columns and lists every random pair, with `nan` for pairs that fail the filters.

Only pairs in which both SNPs pass the filters are written. Without `--marginal`, pairs with approximate p below `--p-approx-cut` carry the exact `beta`, `se` and `p`; the others carry the approximate values. With the default `--p-cut 1`, every tested pair is written, which is about $$M^2/2$$ rows for $$M$$ SNPs. Use `--p-cut` to keep the file small.

## 5. Options

### Main options

| Option | Default | Meaning |
| --- | --- | --- |
| `--data` | — | Phenotype/covariate table (**required**) |
| `--mgrm` | — | GRM prefix(es) of the null model (**required**) |
| `--bfile` | — | PLINK binary prefix (**required**) |
| `--out` | — | Output prefix (**required**) |
| `--trait` | — | Trait column, exactly one (**required**) |
| `--covar`, `--class` | none | Covariate and class columns |
| `--missing-data` | `NA NaN nan N/A ...` | Missing-value tokens in the data file |
| `--add` / `--dom` / `--epis` | — | The test (exactly one) |
| `--marginal` | off | With `--epis`: full GLS model per pair |
| `--maf` | `0.01` | Keep SNPs with $$\min(p, 1 - p)$$ at or above this, in `[0, 0.5]` |
| `--missing-rate` | `0.05` | Keep SNPs with missing rate at or below this, in `[0, 1]` |
| `--min-geno-count` | `0` | `--dom`: minimum number of individuals in each genotype class for the dominance test; otherwise `_dom` is `nan` and the additive-only test is reported (`0` = off) |
| `--snp-set <first> <last>` | — | Test SNPs from the first to the last SNP ID (inclusive); with `--epis`, pairs whose first SNP is in this set |
| `--parallel <m> <i>` | — | Split the scan into `m` parts and run part `i` (1-based); with `--epis`, parts hold about equal numbers of pairs |
| `--p-cut` | `1` | Write SNPs/pairs whose final p is at or below this. For dominance a row is kept if `p_add` or `p_dom` passes. `1` writes all |
| `--p-approx-cut` | `1e-3` | `--epis`: re-test pairs exactly when the approximate p is below this; `1` = exact test for every pair |
| `--solver` | `general` | Null-model algorithm: `general` (dense $$V$$) or `mme` (mixed model equations; for many repeated records) |
| `--val` | `0.001` for `mme`, `0` for `general` | Value added to the GRM diagonal when it is read |
| `--threads` | `10` | Number of threads |

`--snp-set`, `--parallel` and `--start-pos`/`--end-pos` are mutually exclusive.

### Advanced

| Option | Default | Meaning |
| --- | --- | --- |
| `--npart-snp` | `20` | Number of chunks the tested SNPs are read in by `--add` and `--dom` (memory only; not used by `--epis`) |
| `--start-pos` | `0` | First tested SNP (0-based index in the `.bim` file) |
| `--end-pos` | `-1` | End of the tested SNPs (0-based, exclusive). `-1`, or a value above the number of SNPs, means the last SNP |
| `--num-random-snp` | `10000` | `--epis`: random SNP pairs for the $$\gamma$$ calibration; capped at the number of SNPs |
| `--seed` | `42` | `--epis`: random seed for the calibration pairs |
| `--maxiter0` | `100` | Max iterations for the null-model fit |
| `--cc-par0` | `1e-7` | Null-model convergence: norm of the parameter change |
| `--cc-gra0` | `1e-6` | Null-model convergence: norm of the gradient |
| `--cc-logl0` | `5e-5` | Null-model convergence: change in REML log-likelihood |

Set the environment variable `FASTGMA_LOG_LEVEL=debug` to log per-chunk progress.

## 6. Examples

### Build the GRMs

```bash
fastgma --make-grm --grm-type add --bfile test --out output/add
fastgma --make-grm --grm-type dom --bfile test --out output/dom
fastgma --process-grm --make-epis --prefix1 output/add --prefix2 output/add --out output/axa
```

### Additive test with three GRMs

```bash
fastgma --remma --add \
  --bfile test \
  --mgrm output/add output/dom output/axa \
  --data test_simu_pheno.txt \
  --trait pheno \
  --covar age \
  --class sex \
  --threads 4 \
  --out output/test_remma_add
```

### Dominance test

```bash
fastgma --remma --dom \
  --bfile test \
  --mgrm output/add output/dom output/axa \
  --data test_simu_pheno.txt \
  --trait pheno \
  --covar age \
  --class sex \
  --out output/test_remma_dom
```

Add `--start-pos 0 --end-pos 200` to test only the first 200 SNPs.

### Epistasis scan, split into parallel jobs

```bash
for i in 1 2 3 4; do
  fastgma --remma --epis \
    --bfile test \
    --mgrm output/add output/dom output/axa \
    --data test_simu_pheno.txt \
    --trait pheno \
    --num-random-snp 1000 \
    --p-cut 1e-3 \
    --parallel 4 $i \
    --out output/test_remma_epis
done
```

Each job writes `output/test_remma_epis.4_<i>.*`. Each job tests the pairs whose first SNP is in its part, and the parts hold about equal numbers of pairs. Together they cover every pair once. Concatenate the `.assoc` files (keeping one header) to get the genome-wide result.

### Epistasis with the full model per pair

```bash
fastgma --remma --epis --marginal \
  --bfile test \
  --mgrm output/add output/dom output/axa \
  --data test_simu_pheno.txt \
  --trait pheno \
  --start-pos 0 --end-pos 10 \
  --out output/test_remma_epis_marginal
```

This fits one GLS model per pair, so keep the range small.
