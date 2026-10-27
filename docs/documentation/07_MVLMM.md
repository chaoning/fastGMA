---
layout: page
permalink: /documentation/07_MVLMM/
title: mvLMM
description: Multivariate linear mixed model for several correlated traits
order: 7
---

`MVLMM` (multivariate linear mixed model) is the multi-trait GWAS engine of `fastGMA`. For two or more quantitative traits it first fits a **multivariate null model** to estimate the genetic and residual covariance matrices between the traits (and from them the heritabilities and genetic correlations), then runs a **genome-wide scan** that tests each SNP jointly on all traits.

It is invoked through the top-level dispatcher with `--mvlmm`:

```bash
fastgma --mvlmm [options]
```

Run `fastgma -h 6` for the full option list.

## 1. What it does

With $$t$$ traits, the records of all traits are modelled jointly:

$$
\begin{aligned}
y &= Xb + Zu + e \\
\operatorname{Var}(u) &= \Sigma_g \otimes G \\
\operatorname{Var}(e) &= \Sigma_e \otimes I
\end{aligned}
$$

$$y$$ stacks the records of all traits, trait by trait. $$Xb$$ holds the fixed effects of each trait (an intercept plus that trait's covariates and class variables). $$u$$ is the polygenic effect with GRM $$G$$, and $$e$$ is the residual. $$\Sigma_g$$ and $$\Sigma_e$$ are $$t \times t$$ covariance matrices: their diagonals are the genetic and residual variances of the traits, and their off-diagonals the covariances between traits. Each extra GRM (`--mgrm`) or non-genetic random effect (`--random`) adds one more $$t \times t$$ covariance matrix.

Each trait uses its own non-missing records, so the traits do not need to be measured on the same individuals (with `--solver general`). Records of different traits from the same individual share the residual covariance.

**Null model.** All covariance matrices are estimated by REML with the average-information (AI) algorithm. When the AI step would leave the valid parameter space (a covariance matrix that is not positive definite), the update is chosen, by the REML log-likelihood, among an AI step blended with the EM step, the AI step shortened to 90% of the distance to the boundary, and an AI step in the Cholesky parameterisation. This lets the fit reach estimates on the boundary, such as a genetic correlation near ±1 or a variance near 0. If a covariance matrix ends close to singular, the log says so. Iteration stops when the parameter change is below `--cc-par0`, the gradient norm is below `--cc-gra0`, or the change in log-likelihood is below `--cc-logl0`, or after `--maxiter0` iterations.

**SNP test.** For each SNP, the centred count of `a1` is added as a fixed effect with one coefficient per trait, and the covariance matrices are kept at their null-model values. The $$t$$ effects are tested jointly ("no effect on any trait") with a score/Wald statistic that follows a chi-square distribution with $$t$$ degrees of freedom.

**Approximate test.** The exact test needs $$F = S'PS$$ (the $$t \times t$$ information matrix of the SNP effects) for every SNP. `MVLMM` approximates it by $$F \approx D\, \Gamma\, D$$, where $$D = \operatorname{diag}(\sqrt{x_a'x_a})$$, $$x_a$$ is the centred SNP over the records of trait $$a$$ (the same $$x$$ for every trait with complete data, so that $$F \approx (x'x)\, \Gamma$$), and $$\Gamma$$ is the mean of $$D^{-1} F D^{-1}$$ over `--num-random-snp` random SNPs that pass the MAF and missing-rate filters. $$F$$ does not depend on the phenotypes, so no SNP is excluded from the calibration for being associated. The random SNPs are drawn from the whole genotype file with `--seed` (default `42`), so $$\Gamma$$ is the same in reruns and in every `--parallel` / `--snp-set` job. SNPs with approximate p below `--p-approx-cut` are re-tested exactly. `--p-approx-cut 1` gives the exact test for every SNP; `--p-approx-cut 0` keeps the approximate test only.

> **Choosing `--p-approx-cut`.** The exact test is computed in batches and is cheap for samples below about 10,000 individuals: use `--p-approx-cut 1` there. The approximation works well in large samples of distantly related individuals and less well in small or structured samples. In the 1,062-individual mouse example, the approximate and exact `-log10 p` correlated at 0.94 and the only SNP with exact p < 1e-3 was not re-tested at the default cut (it was at `--p-approx-cut 0.01`). In a human sample of 427 individuals the correlation was 0.99996. If you use the approximate test on structured data, choose a larger cut such as `0.01`.

## 2. Solvers

The fitting algorithm is chosen with `--solver` (default `general`):

| `--solver` | Method | Notes |
| --- | --- | --- |
| `general` (default) | Dense AI-REML on the covariance matrix of all records (one row per record and trait) | Allows missing trait values, several GRMs (`--mgrm`), non-genetic random effects (`--random`), per-trait covariates and `--predict` |
| `eigen` | The GRM is eigendecomposed once ($$G = U \Lambda U'$$) and $$y$$ and $$X$$ are rotated by $$U'$$. The covariance matrix becomes block diagonal with one $$t \times t$$ block $$\lambda_i \Sigma_g + \Sigma_e$$ per individual, so REML costs $$O(n)$$ per evaluation | Needs complete (balanced) data and exactly one `--grm`. Supports `--refit` |

The two solvers maximise the same likelihood. On data that both can fit they give the same `.var`, `.sln`, `.h2`, `.rg` and `.assoc` files.

Rules checked at start-up. Each violation stops the run with an error; the eigen-solver errors suggest `--solver general`:

- `--solver eigen` needs every individual with records to be observed for every trait (with its covariates). Remove the incomplete individuals from the data file, or use `--solver general`.
- `--solver eigen` needs exactly one GRM given with `--grm` (no `--mgrm`).
- `--solver eigen` does not allow non-genetic random effects. `--random` must be omitted or name only the ID column.
- `--solver eigen` needs the same `--covar` and `--class` columns for all traits.
- `--predict` is not available with `--solver eigen`.
- `--refit` is accepted only with `--solver eigen`.

**Per-SNP refit (`--solver eigen --refit`).** Each SNP that is tested exactly (every SNP with `--p-approx-cut 1`) gets its own REML fit of $$\Sigma_g$$ and $$\Sigma_e$$ with the SNP in the model, starting from the null estimates. The `p` column is then the joint Wald test with the re-estimated covariance matrices. An approximate likelihood-ratio p-value `p_lrt` (ML log-likelihoods evaluated at the REML estimates of the two models) and the convergence diagnostics `cc_par` / `cc_gra` are added. Without `--refit` the covariance matrices stay at the null estimates for every SNP.

## 3. Inputs

### Data table (`--data`)

A whitespace- or tab-separated table with a header row. The first column must be the individual ID. It must match the IID (second column) of the `.grm.id` file and of the PLINK `.fam` file. Rows whose ID is not in the GRM(s) are dropped.

Columns are selected by **header name**. A token that is not a header name but a non-negative integer is read as a 0-based column index (column `0` is the ID).

| Option | Meaning |
| --- | --- |
| `--data` | Data file (**required**) |
| `--trait` | Column names of the traits (**required**; two or more, one is allowed) |
| `--covar` | Continuous covariates. Give once for all traits, or once per trait in `--trait` order |
| `--class` | Categorical covariates. Give once for all traits, or once per trait in `--trait` order |
| `--random` | ID column of each random effect (see below). Default: the ID column for every GRM |
| `--missing-data` | Tokens treated as missing in the data file (default `NA NaN nan N/A ...`; giving the option replaces the list) |

An intercept is always fitted for each trait. Per-trait covariates: `--trait y1 y2 --covar E1 E2 --covar E1` fits `E1` and `E2` for `y1` and `E1` for `y2`. `--covar E1` alone fits `E1` for both traits. `--covar` or `--class` given a number of times other than 1 or the number of traits is an error. Fixed-effect columns that are linearly dependent on earlier ones (for example one level of each class variable) are dropped and reported in the log.

A record of a trait is dropped when the trait, one of that trait's covariates or class variables, or one of the `--random` columns is missing. The other traits of the same individual are kept.

Each individual may have **at most one record per trait**. A repeated ID within a trait stops the run with an error; average or remove repeated records first.

### GRMs and random effects (`--grm`, `--mgrm`, `--random`)

| Option | Meaning |
| --- | --- |
| `--grm` | Prefix of a GRM. Its covariance block is named `G` |
| `--mgrm` | Prefixes of further GRMs (for example a dominance GRM). Each block is named after the prefix basename. `--solver general` only |
| `--random` | One ID column per GRM (in the order: `--mgrm` prefixes first, then `--grm`), then the columns of non-genetic random effects |

At least one of `--grm` / `--mgrm` is required. The first entries of `--random` link the GRMs to a data column (normally the ID column). Further entries are non-genetic random effects, such as a cage or batch column, with independent levels. A non-genetic block is named after its column.

| Model | Options |
| --- | --- |
| One GRM | `--grm output/test` |
| Dominance + additive GRM | `--mgrm output/test_dom --grm output/test` |
| GRM + non-genetic effect in column `cage` | `--grm output/test --random 0 cage` |

In the last row, `0` is the ID column (its header name works too). When `--random` is given it must list at least one column per GRM.

GRMs are referenced by prefix and must be in GCTA format (`.grm.bin`, `.grm.id`). GRMs built by `fastgma --make-grm` and by GCTA are both accepted; see [GRM](../03_GRM/).

### Genotype input

| Option | Meaning |
| --- | --- |
| `--bfile` | PLINK binary genotype prefix (`.bed/.bim/.fam`) |
| `--ffile` | Feature-matrix prefix (`.fam/.bim/.fmat`). Each `.fmat` row is one feature and each column one sample, for imputed dosages, expression or any other continuous per-sample matrix |
| `--missing-geno` | Tokens treated as missing in the `.fmat` file (default `NA NaN nan N/A ...`) |

Give exactly one of `--bfile` and `--ffile`. Every individual with records must be in the `.fam` file. A `--ffile` feature is centred by its mean, missing values get the mean, `--maf` is not applied and `af` is the feature mean. The SNP scan always runs; there is no variance-components-only mode.

## 4. Outputs

All files share the `--out` prefix and are tab-separated with a header row. Values that are not finite print as `nan` (`NA` in `.sln`).

| File | Contents |
| --- | --- |
| `<out>.var` | Estimated (co)variance components with standard errors |
| `<out>.var_cov` | Sampling covariance matrix of the components (inverse AI matrix) |
| `<out>.h2` | Heritability of each trait, SE by the delta method |
| `<out>.rg` | Genetic correlation of each pair of traits, SE by the delta method |
| `<out>.sln` | Fixed-effect estimates; with `--predict`, also the predicted random effects |
| `<out>.assoc` | One row per SNP with the joint association test |
| `<out>.log` | Log of the run (appended) |

### `<out>.var`

One row per element of the lower triangle of each covariance matrix:

```
component	row	col	variance	se
G	1	1	1.14819	0.16287
G	2	1	0.674167	0.113107
G	2	2	0.687829	0.11156
residual	1	1	0.961924	0.0498522
residual	2	1	0.212369	0.0385441
residual	2	2	1.15062	0.0563373
```

| Column | Meaning |
| --- | --- |
| `component` | `G` for the `--grm` effect, the prefix basename for each `--mgrm` GRM, the column name for each non-genetic `--random` effect, and `residual` last |
| `row`, `col` | Trait numbers (1-based, in `--trait` order). `row` = `col` = $$i$$ is the variance of trait $$i$$; `row` = $$i$$, `col` = $$j$$ is the covariance between traits $$i$$ and $$j$$ |
| `variance`, `se` | Estimate and standard error (square root of the diagonal of the inverse AI matrix at the estimates) |

The blocks appear in the order `--mgrm` GRMs, `--grm`, non-genetic effects, `residual`.

### `<out>.var_cov`

The full sampling covariance matrix of the `.var` estimates, in the same order. The header row is `parameter` followed by one label per parameter, and each row starts with its label. Labels are `<component>_<row>_<col>`, for example `G_1_1`, `G_2_1`, `residual_2_2`. Values are written with 10 significant digits, as in `.sln` (`.var`, `.h2`, `.rg` and `.assoc` use 6).

### `<out>.h2` and `<out>.rg`

```
trait	h2	se
y1	0.548342	0.0404541
y2	0.378117	0.0420733
```

```
trait1	trait2	rg	se
y2	y1	0.758612	0.0555244
```

Each GRM block $$m$$ is weighted by $$d_m$$, the mean diagonal element of its GRM over the analysed records (pooled over all traits; with complete data, the mean diagonal of the GRM over the individuals), as in `UVLMM`, `REMMA` and `--fast-reml`; non-genetic effects and the residual have $$d = 1$$. The `h2` of trait $$i$$ is $$\sum_{m \in \text{GRM}} d_m S_{m,ii} \big/ \sum_{m} d_m S_{m,ii}$$, where $$S_{m,ii}$$ is its variance in block $$m$$ and the denominator runs over all blocks. A non-genetic random effect such as `cage` therefore counts in the denominator only. The `rg` between traits $$i$$ and $$j$$ is computed from the weighted sum of the GRM blocks, $$S_g = \sum_{m \in \text{GRM}} d_m S_m$$: $$r_g = S_{g,ij} / \sqrt{S_{g,ii}\, S_{g,jj}}$$; with a single GRM, $$d$$ cancels and `rg` is unchanged. The $$d_m$$ are written to the log. Both SEs use the delta method with `<out>.var_cov`. `.rg` has one row per pair (the later trait in `--trait` order first) and only a header for a single trait.

### `<out>.sln`

```
trait	type	term	levels	effect	se
y1	0	mu	1	1.44219	0.093429
y1	1	cov1	1	0.414169	0.0498185
y1	2	sex	1	-0.347646	0.110603
y1	2	sex	2	0	NA
```

| Column | Meaning |
| --- | --- |
| `trait` | Trait name |
| `type` | `0` intercept, `1` covariate, `2` class level, `3` GRM random effect, `4` non-genetic random effect |
| `term` | Fixed effects: column name (or SNP ID). Random effects: the variance component as named in `<out>.var` |
| `levels` | Class level, or individual ID / level for random effects (`1` for intercept and covariates) |
| `effect`, `se` | Estimate and standard error (for random effects, the square root of the prediction error variance, PEV) |
| `rel` | With `--predict`: reliability $$1 - \mathrm{PEV} / \mathrm{prior\ variance}$$ of a random-effect level; `NA` on fixed-effect rows |

A fixed-effect column dropped as linearly dependent has `effect` `0` and `se` `NA` (here `sex` level `2`). Numbers keep 10 significant digits (the rows above are from an earlier version that printed 6). With `--predict`, one row per trait and level is added for each random effect (`type` 3 or 4). For GRM random effects, `--predict` also keeps the GRM individuals without records, so their genetic values are predicted too.

### `<out>.assoc`

One row per SNP. For the traits `y1` and `y2`:

```
chr	snp	pos	a1	a2	af	N	beta_y1	se_y1	beta_y2	se_y2	cov_y2_y1	p
```

With `--solver eigen --refit`:

```
chr	snp	pos	a1	a2	af	N	cc_par	cc_gra	beta_y1	se_y1	beta_y2	se_y2	cov_y2_y1	p_lrt	p
```

| Column | Meaning |
| --- | --- |
| `chr`, `snp`, `pos`, `a1`, `a2` | Chromosome, SNP ID, position, effect allele and other allele from the `.bim` file |
| `af` | Frequency of `a1` among the individuals with records (the feature mean with `--ffile`) |
| `N` | Number of individuals with records and a non-missing genotype |
| `cc_par`, `cc_gra` | `--refit` only: norm of the last parameter change and of the gradient of the per-SNP REML fit (`nan` for SNPs that were not refitted) |
| `beta_<trait>`, `se_<trait>` | Effect of one copy of `a1` on the trait, and its standard error |
| `cov_<ti>_<tj>` | Sampling covariance between the effect estimates of traits `ti` and `tj`, for every pair with `ti` after `tj` in `--trait` order (for three traits: `cov_y2_y1 cov_y3_y1 cov_y3_y2`) |
| `p_lrt` | `--refit` only: approximate likelihood-ratio p-value (chi-square, `t` df) |
| `p` | Joint p-value (chi-square, `t` df) |

For a SNP that was not re-tested exactly, `beta`, `se`, `cov` and `p` come from the approximate test. SNPs that fail `--maf` or `--missing-rate` and SNPs whose test is singular (for example monomorphic among the records of one trait) have `nan` in the effect and p-value columns.

With `--p-cut` below `1`, only SNPs with a finite `p` at or below the cut are written. The default `1` writes every SNP, including the `nan` rows.

## 5. Key options

### Null model

| Option | Default | Meaning |
| --- | --- | --- |
| `--solver` | `general` | `general` or `eigen` (see Section 2) |
| `--maxiter0` | `100` | Max iterations for the null-model fit |
| `--continue` | off | Start from the estimates in an existing `<out>.var` (ignored with a warning if they do not fit the model or are not positive definite) |
| `--predict` | off | Also write the predicted random effects to `<out>.sln` (`--solver general`) |
| `--threads` | `10` | Number of MKL/OpenMP threads |

### SNP scan

| Option | Default | Meaning |
| --- | --- | --- |
| `--maf` | `0.01` | Minimum minor allele frequency `min(af, 1 − af)` (`--bfile` only) |
| `--missing-rate` | `0.05` | Keep a SNP if its missing rate is at or below this |
| `--p-approx-cut` | `1e-3` | SNPs with approximate p below this are re-tested exactly. `0` = approximate test only; `1` = exact test for every SNP |
| `--refit` | off | `--solver eigen`: re-estimate the covariance matrices for every exactly tested SNP |
| `--p-cut` | `1` | In `[0, 1]`. `1` writes every SNP (including `nan` rows); below `1`, only SNPs with a finite `p` at or below the cut |

### Restricting / parallelizing the scan

| Option | Meaning |
| --- | --- |
| `--snp-set <first> <last>` | Test only the SNPs from the first to the last SNP ID (inclusive). The output prefix becomes `<out>.<start>_<end>` (0-based start, exclusive end) |
| `--parallel <m> <i>` | Split the SNPs into `m` equal parts and test part `i` (1-based). The output prefix becomes `<out>.<m>_<i>` |
| `--start-pos` / `--end-pos` | 0-based index of the first SNP and one past the last SNP (`-1` = up to the last SNP). Advanced; output prefix unchanged |

`--snp-set` and `--parallel` cannot be used together or with `--start-pos` / `--end-pos`. Every job fits the same null model and writes all output files (`.var`, `.var_cov`, `.h2`, `.rg`, `.sln`, `.assoc`) to its own prefix. The log is `<out>.<m>_<i>.log` for `--parallel`, `<out>.<start>_<end>.log` for `--snp-set` and `<out>.log` otherwise. Because `--seed` is fixed, all jobs use the same $$\Gamma$$, and the merged `.assoc` files match a single full run.

### Advanced (rarely needed)

Grouped under `Advanced` in `--help`.

| Option | Default | Meaning |
| --- | --- | --- |
| `--num-random-snp` | `2000` | Random SNPs for the approximate-test calibration (capped at the number of SNPs) |
| `--seed` | `42` | Random seed for drawing the calibration SNPs |
| `--npart-snp` | `20` | Minimum number of chunks the tested SNPs are read in; chunks are kept to about 128 MB of genotype codes (memory only; does not change results) |
| `--cc-par0` / `--cc-gra0` | `1e-7` / `1e-6` | Null model: parameter-change / gradient-norm convergence criteria |
| `--cc-logl0` | `5e-5` | Null model: log-likelihood-change convergence criterion |
| `--maxiter` | `30` | `--refit` only: max iterations per SNP |
| `--cc-par` / `--cc-gra` / `--cc-logl` | `1e-7` / `1e-6` / `5e-5` | `--refit` only: per-SNP convergence criteria |

Set the environment variable `FASTGMA_LOG_LEVEL=debug` to log per-block progress.

## 6. Examples

Build the GRM first (see [GRM](../03_GRM/)):

```bash
fastgma --make-grm --bfile test --out output/test
```

In the examples, `test_simu_pheno.txt` has the ID in its first column, the traits `y1`, `y2` and `y3`, and the covariates `E1` and `E2`.

### Two-trait GWAS with one GRM

```bash
fastgma --mvlmm \
  --data test_simu_pheno.txt \
  --grm output/test \
  --bfile test \
  --trait y1 y2 \
  --covar E1 \
  --out output/test_mvlmm
```

### Exact test for every SNP, write only SNPs with p ≤ 1e-3

```bash
fastgma --mvlmm \
  --data test_simu_pheno.txt \
  --grm output/test \
  --bfile test \
  --trait y1 y2 \
  --covar E1 \
  --p-approx-cut 1 \
  --p-cut 1e-3 \
  --out output/test_mvlmm_exact
```

### Eigen solver with a per-SNP refit (complete data)

```bash
fastgma --mvlmm --solver eigen --refit \
  --data test_simu_pheno.txt \
  --grm output/test \
  --bfile test \
  --trait y1 y2 \
  --covar E1 \
  --p-approx-cut 1 \
  --out output/test_mvlmm_eigen
```

Every individual with a record must have both traits and `E1`.

### Per-trait covariates, three traits

`E1` and `E2` for `y1`, `E1` for `y2` and `E2` for `y3`:

```bash
fastgma --mvlmm \
  --data test_simu_pheno.txt \
  --grm output/test \
  --bfile test \
  --trait y1 y2 y3 \
  --covar E1 E2 --covar E1 --covar E2 \
  --out output/test_mvlmm_cov
```

`--class` is given the same way.

### GRM plus a non-genetic random effect

```bash
fastgma --mvlmm \
  --data test_simu_pheno.txt \
  --grm output/test \
  --bfile test \
  --trait y1 y2 \
  --covar E1 \
  --random 0 cage \
  --out output/test_mvlmm_cage
```

`<out>.var` then has the blocks `G`, `cage` and `residual`.

### Dominance and additive GRMs, with predicted random effects

```bash
fastgma --make-grm --bfile test --grm-type dom --out output/test_dom

fastgma --mvlmm \
  --data test_simu_pheno.txt \
  --mgrm output/test_dom \
  --grm output/test \
  --bfile test \
  --trait y1 y2 \
  --covar E1 \
  --predict \
  --out output/test_mvlmm_mgrm
```

`<out>.var` then has the blocks `test_dom`, `G` and `residual`, and `.h2` / `.rg` combine `test_dom` and `G`.

### Split a scan into 10 parallel jobs

```bash
fastgma --mvlmm \
  --data test_simu_pheno.txt \
  --grm output/test \
  --bfile test \
  --trait y1 y2 \
  --covar E1 \
  --parallel 10 1 \
  --out output/test_mvlmm
```

Repeat with part indices `2` through `10`, then concatenate `output/test_mvlmm.10_*.assoc` (keeping one header).

## 7. Run time and memory

`--solver general` builds dense matrices with one row per record, $$N$$ rows in total (up to $$tn$$ for $$t$$ traits and $$n$$ individuals). It holds about two $$N \times N$$ matrices plus the GRMs themselves (8 bytes per element), and each REML iteration costs $$O(N^3)$$. It is suited to data sets with up to a few tens of thousands of records. `--solver eigen` costs one $$O(n^3)$$ eigendecomposition and then $$O(n)$$ per REML iteration; the exact SNP test costs $$O(n^2)$$ per SNP for the rotation.

Wall time with 4 threads on a shared node (exact test for every SNP unless noted):

| Data | `--solver general` | `--solver eigen` |
| --- | --- | --- |
| Mouse, 2 traits with missing values (2,358 records, 1,295 individuals, 1,407 SNPs) | 5.7 s | — |
| Mouse, 2 traits, complete cases (1,062 individuals, 1,407 SNPs) | 6.4 s | 0.97 s |
| Human, 2 traits (427 individuals, 256,868 SNPs) | 22.0 s | 12.3 s |
| Human, as above, default `--p-approx-cut 1e-3` | 6.9 s | 5.4 s |

For the univariate counterpart, see [uvLMM](../04_UVLMM/).
