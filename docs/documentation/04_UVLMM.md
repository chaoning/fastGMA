---
layout: page
permalink: /documentation/04_UVLMM/
title: uvLMM
description: Univariate linear mixed model for single trait
order: 4
---

`UVLMM` (univariate linear mixed model) is the single-trait GWAS engine of `fastGMA`. For one quantitative phenotype it first fits a **null model** to estimate the variance components (and heritability), then performs a **genome-wide single-SNP association scan** while controlling for the polygenic background and any additional random effects.

It is invoked through the top-level dispatcher with `--uvlmm` (short form `-u`):

```bash
fastgma --uvlmm [options]
```

Run `fastgma -h 3` for the full option list.

## 1. What it does

For each SNP, `UVLMM` tests the additive effect under a mixed model

$$
y = Xb + s\,g + u + e
$$

where $$Xb$$ are the fixed effects (intercept, covariates, class variables), $$s\,g$$ is the SNP being tested ($$g$$ its genotypes, $$s$$ its effect), $$u$$ is the polygenic random effect captured by one or more GRMs, and $$e$$ is the residual. The null model (no SNP) is fitted once to estimate the variance components. The genome-wide scan then uses a fast approximate test for every SNP and an exact test only for promising SNPs.

**Approximate test.** The exact test needs $$g'Pg$$ for every SNP. `UVLMM` approximates it by $$g'Pg \approx \gamma\, g'g$$. $$\gamma$$ is the mean of $$g'Pg / g'g$$ over `--num-random-snp` random SNPs that pass the MAF and missing-rate filters. The random SNPs are drawn from the whole genotype file with `--seed` (default `42`), so $$\gamma$$ is the same in reruns and in every `--parallel` / `--snp-set` job. SNPs with approximate p below `--p-approx-cut` are re-tested exactly. `--p-approx-cut 1` gives the exact test for every SNP; `--p-approx-cut 0` keeps the approximate test only.

## 2. Model and solver

The model is chosen with `--model` (default `lmm`). The fitting algorithm is chosen with `--solver`:

| `--model` | Model | `--solver` | Requires |
| --- | --- | --- | --- |
| `lmm` (default) | One GRM as random effect | `eigen` (default), `mme` or `general` | `--grm` |
| `lm` | Linear model, no random effect | not used | — |
| `multi-grm` | Several GRMs as random effects | `general` (default) or `mme` | `--mgrm` |
| `repeat` | GRMs plus a permanent (individual) effect for repeated records | `mme` (default) or `general` | `--mgrm` |

The three solvers:

| `--solver` | Method | Notes |
| --- | --- | --- |
| `eigen` | The GRM is eigendecomposed once and the data are rotated. REML and the per-SNP tests then cost $$O(n)$$ per evaluation. | `--model lmm` only. Needs one record per individual. `--random` is not allowed. Supports `--refit` / `--no-refit`. |
| `mme` | AI-REML on Henderson's mixed model equations (WEMAI as for `general`). The equations have one row per fixed effect and per random-effect level, so their size does not grow with the number of records; each GRM is inverted once. | Allows repeated records and `--random`. Needs positive-definite GRMs (see `--val` below). Same estimates and tests as `general`. Exact re-tests use the null variances. |
| `general` | Dense AI-REML on the $$n \times n$$ covariance matrix (WEMAI: the AI matrix is blended with the EM matrix when needed to keep variances positive). | Allows repeated records and `--random`. Exact re-tests use the null variances. |

**`mme` or `general`.** Both fit the same model; they differ in cost. Each `general` iteration factorizes the $$n \times n$$ covariance matrix of the $$n$$ records, while `mme` works with a matrix of size $$p + \sum_k q_k$$ ($$p$$ fixed effects, $$q_k$$ levels of random effect $$k$$). With many records per individual, as for repeated measures, `mme` is much faster and needs far less memory: on `example/mouse_long` (19,392 records of 1,212 mice, `--model repeat`, 8 threads) the null model and the scan of 11,833 SNPs take about 4 s and 0.2 GB in total. Each iteration costs about the cube of the matrix size in both, so `mme` is faster when $$p + \sum_k q_k$$ is smaller than the number of records. Each GRM adds a block for all individuals, so with three GRMs and two records per individual `mme` is slower than `general`. With one record per individual and one GRM, `eigen` is the fastest. When `--solver general` is used and the mixed model equations are smaller than the number of records, the log suggests `--solver mme`.

**GRM diagonal (`--val`).** `mme` needs the inverse of each GRM, and a centred GRM is singular (its rows sum to zero). `mme` therefore adds `--val` (default `0.001`) to the GRM diagonal when it reads the GRM, as `--longwas` does. The other solvers add nothing unless `--val` is given. To compare `mme` with `general` on the same model, give both the same `--val`. If a GRM is still not positive definite (for example, duplicated individuals), the run stops and suggests a larger `--val`.

Rules checked at start-up:

- `--solver eigen` or `--solver mme` with `--model lm` is an error. `--solver general` with `--model lm` is ignored with a warning.
- `--solver eigen` with `--model multi-grm` or `--model repeat` is an error.
- `--random` with `--solver eigen` is an error. Use `--model lmm --solver mme` (or `general`) to add non-genetic random effects.
- An individual with more than one record is an error with `--solver eigen`. Use `--solver mme` (or `general`), or `--model repeat` for a permanent effect.
- `--model repeat` cannot be combined with `--random`.
- `--refit` is accepted only with `--solver eigen`.

**Per-SNP refit (`--solver eigen`).** By default (`--refit`), each SNP that is re-tested exactly gets its own REML fit of the variance components with the SNP in the model. Score, likelihood-ratio and Wald p-values are then reported. With `--no-refit` the exact test uses the null variances, so the Wald p equals the score p. `--refit --p-approx-cut 1` refits every SNP.

## 3. Inputs

### Phenotype / covariate table (`--data`)

A whitespace- or tab-separated table with a header row. The first column must be the individual ID. Columns are selected by name. Only individuals present in both the data table and the GRM are retained.

| Option | Meaning |
| --- | --- |
| `--data` | Data file (**required**) |
| `--trait` | Column name of the trait — exactly one (**required**; UVLMM is univariate) |
| `--covar` | Column name(s) of continuous covariates |
| `--class` | Column name(s) of categorical covariates |
| `--random` | Column name(s) of additional non-genetic random effects (`--solver mme` or `general`) |
| `--missing-data` | Tokens treated as missing in the data table (default `NA NaN nan N/A ...`) |

### GRM input

| Option | Meaning |
| --- | --- |
| `--grm` | Prefix of a single GRM (`--model lmm`) |
| `--mgrm` | Prefixes of one or more GRMs (`--model multi-grm` / `--model repeat`) |
| `--val` | Value added to the GRM diagonal when it is read (default `0.001` for `--solver mme`, `0` otherwise; see Section 2) |

All `--mgrm` GRMs must list the same individuals in the same order. As elsewhere in `fastGMA`, GRMs are referenced by prefix (e.g. `--grm output/test`), not by a single file name. GRMs built by `fastgma --make-grm` and by GCTA are both accepted; see [GRM](../03_GRM/) for how to build them.

### Genotype input

| Option | Meaning |
| --- | --- |
| `--bfile` | PLINK binary genotype prefix (`.bed/.bim/.fam`) |
| `--ffile` | Feature-matrix prefix (`.fam/.bim/.fmat`). Each `.fmat` row is one feature and each column one sample — suitable for imputed dosages, expression, or any continuous per-sample matrix. |
| `--missing-geno` | Tokens treated as missing in the `.fmat` file (default `NA NaN nan N/A ...`) |

Give only one of `--bfile` / `--ffile`. A `--ffile` feature is centred by its mean, missing values get the mean, `--maf` is not applied and `af` is the feature mean. If no genotype file is supplied, the command fits the null model, writes `.var`, `.h2` and `.sln`, and stops before the SNP scan.

## 4. Outputs

All files share the `--out` prefix and are tab-separated with a header row:

| File | Contents |
| --- | --- |
| `<out>.var` | Estimated variance components with standard errors (inverse AI matrix at the estimates) |
| `<out>.h2` | Per-component proportion $$h^2_i = \sigma^2_i d_i / \sum_k \sigma^2_k d_k$$ ($$d_i$$ = mean GRM diagonal, 1 for non-GRM effects), SE by the delta method (not written for `--model lm`) |
| `<out>.sln` | Fixed-effect estimates; with `--predict`, also the predicted random effects (BLUP) and their reliabilities |
| `<out>.res` | With `--predict`: each record split into its fixed, random and residual parts |
| `<out>.assoc` | One row per SNP with the association results |
| `<out>.log` | Log of the run (appended; `<out>.<m>_<i>.log` with `--parallel`, `<out>.<start>_<end>.log` with `--snp-set`) |

`<out>.var` and `<out>.h2` are keyed by component name rather than a positional
index: `G` for the `--grm` effect, the GRM prefix basename for each
`--mgrm` entry, `permanent` for the `--model repeat` permanent-environment
effect, the column name for each `--random` effect, and `residual` last.
`--model lm` writes only the `residual` row.

```
component	variance	se
G	1.26776	0.206283
residual	2.39709	0.105419
```

`<out>.h2` adds a trailing `total` row holding the combined proportion of the
GRM components — the number usually quoted as the heritability. A
`permanent` or `--random` effect contributes variance but is not counted here:

```
component	h2	se
G	0.349279	0.0403661
residual	0.650721	0.0403661
total	0.349279	0.0403661
```

Each GRM variance is put on the phenotypic scale by the mean diagonal $$d_k$$ of its
kernel $$Z G_k Z'$$ over the analysed records (for one record per individual, the mean
diagonal of the GRM over the analysed individuals); a `permanent` / `--random` effect
and the residual have $$d = 1$$:

$$
h^2_k = \frac{\sigma^2_k d_k}{V_p}, \qquad V_p = \sum_k \sigma^2_k d_k \quad \bigl(\text{the mean of } \operatorname{diag}(V)\bigr)
$$

This is the definition used by `--fast-reml` (`variance_prop`) and the PVE of GEMMA,
and the same for every model and solver, so `--grm A` and `--mgrm A` give the same
`.h2`. It matters for GRMs whose mean diagonal is far from 1 (for example an $$A \circ A$$
epistatic GRM). For a standard additive GRM $$d$$ is $$1 + F$$ averaged over the sample
(close to 1); GCTA `--reml` reports the unscaled `V(G)/Vp`, which can be computed from
`<out>.var` as $$\sigma^2_G / (\sigma^2_G + \sigma^2_e)$$. The $$d_k$$ are written to the log
("mean diagonal over the records"). SEs use the delta method on the inverse AI matrix.

Values that are not finite (for example an SE that could not be derived from the
AI matrix) print as `NA`.

`<out>.sln` has the columns `type term levels effect se`, and with `--predict` also `rel`. `type` is `0` for the intercept (`mu`), `1` for a covariate and `2` for a class level. A fixed-effect column removed as linearly dependent (for example the reference level of a class variable) is listed with effect `0` and se `NA`. With `--predict`, type `3` rows hold the predicted genetic values (BLUP) and type `4` rows the predicted `--random` / `permanent` effects. For these rows `term` is the variance component as named in `<out>.var` (`G`, the GRM prefix basename, `permanent` or the `--random` column), `levels` is the individual ID or level, `se` is the square root of the prediction error variance (PEV), and `rel` is the reliability $$1 - \mathrm{PEV} / (\sigma^2_k K_{jj})$$ ($$K_{jj}$$ = GRM diagonal, `1` for a non-genetic effect); `rel` is `NA` on the fixed-effect rows. GRM rows cover every GRM individual, including those without records, whose values are predicted through their relationships. Numbers keep 10 significant digits. From the [repeated-records](#repeated-records) example with `--predict` (first rows of each kind):

```
type	term	levels	effect	se	rel
0	mu	1	19.7737604	0.07620444974	NA
2	Sex	1	-2.65796099	0.09581998468	NA
2	Sex	2	0	NA	NA
2	week	1	-13.23082676	0.05535989006	NA
...
3	mouse_long	1419	0.4080401921	0.3347456694	0.71612076
...
4	permanent	1000	-0.7212183419	0.4809213132	0.9077624809
...
```

`<out>.res` (with `--predict`) has one row per record: `id observed fixed random residual`, where `fixed` is $$x_i'\hat b$$, `random` is the sum of the predicted random effects of the record and `residual` is $$\hat e_i$$, so that `observed` = `fixed` + `random` + `residual`. `observed` - `fixed` is the phenotype corrected for the fixed effects. From the same run:

```
id	observed	fixed	random	residual
1419	4.9	3.884972647	0.4058137606	0.6092135924
1419	8.4	7.348924696	0.4058137606	0.6452615437
...
```

The `<out>.assoc` columns depend on the solver.

`--solver eigen` (`--model lmm` default):

```
chr snp pos a1 a2 af N cc_par cc_gra beta se p_score p_lrt p_wald
```

`cc_par` and `cc_gra` are the convergence diagnostics of the per-SNP refit (norm of the last parameter change and of the gradient). They are written only with `--refit` (the default) and are `nan` for SNPs that were not re-tested exactly. With `--no-refit` the header is `chr snp pos a1 a2 af N beta se p_score p_lrt p_wald`. For a SNP that was not re-tested exactly, `beta` and `se` are the approximate values and the three p-value columns all hold the approximate p.

`--model lm` (exact least-squares test for every SNP):

```
chr snp pos a1 a2 af N beta se p_score p_lrt p_wald
```

`--solver general` or `mme` (`multi-grm`, `repeat`, or `--model lmm` with either solver): only the Wald test is computed, so a single `p` column is written:

```
chr snp pos a1 a2 af N beta se p
```

Columns: chromosome `chr`, SNP ID `snp`, position `pos`, effect allele `a1` and other allele `a2`, allele frequency of `a1` (or feature mean) `af`, `N` (number of individuals with a non-missing genotype), the effect estimate `beta` with its standard error `se`, and the p-values. SNPs that fail the MAF/missing-rate filters report `nan`. With `--p-cut` below `1`, only SNPs whose final p (`p_wald`, or `p` for the general and mme solvers) is at or below the cut are written.

## 5. Key options

### Null model

| Option | Default | Meaning |
| --- | --- | --- |
| `--model` | `lmm` | `lmm`, `lm`, `multi-grm` or `repeat` |
| `--solver` | `eigen` for `lmm`; `general` for `multi-grm`; `mme` for `repeat` | `eigen`, `mme` or `general` (see Section 2) |
| `--val` | `0.001` for `mme`, `0` otherwise | Value added to the GRM diagonal (see Section 2) |
| `--maxiter0` | `100` | Max iterations for the null-model fit |
| `--continue` | off | Start from the variances in an existing `<out>.var` |
| `--predict` | off | Also write the predicted random effects and their reliabilities to `<out>.sln`, and `<out>.res` |
| `--threads` | `10` | Number of MKL/OpenMP threads |

### SNP scan

| Option | Default | Meaning |
| --- | --- | --- |
| `--maf` | `0.01` | Minimum minor allele frequency $$\min(p, 1 - p)$$ (PLINK input only) |
| `--missing-rate` | `0.05` | Keep a SNP if its missing rate is at or below this |
| `--p-approx-cut` | `1e-3` | SNPs with approximate p below this are re-tested exactly. `0` = approximate test only; `1` = exact test for every SNP |
| `--refit` / `--no-refit` | `--refit` | `--solver eigen`: re-estimate the variance components for every exactly re-tested SNP |
| `--maxiter` | `30` | Max iterations for the per-SNP refit |
| `--p-cut` | `1` | Write only SNPs whose final p-value is at or below this (`1` = all, including `nan` rows) |

### Restricting / parallelizing the scan

| Option | Meaning |
| --- | --- |
| `--snp-set <first> <last>` | Test only the SNPs from the first to the last SNP ID (inclusive). The output prefix becomes `<out>.<start>_<end>` (0-based start, exclusive end), for example `<out>.<start>_<end>.assoc` and the log `<out>.<start>_<end>.log`. |
| `--parallel <m> <i>` | Split the SNPs into `m` equal parts and test part `i` (1-based). The output prefix becomes `<out>.<m>_<i>`, for example `<out>.<m>_<i>.assoc` and the log `<out>.<m>_<i>.log`. |

`--snp-set` and `--parallel` cannot be used together; `--start-pos` / `--end-pos` apply only without them. Each job fits the same null model and writes all its files (`.var`, `.h2`, `.sln`, `.assoc`, `.log`) to its own prefix. Because `--seed` is fixed, all jobs use the same $$\gamma$$, and the merged `.assoc` files match a single full run.

### Advanced (rarely needed)

Grouped under `Advanced` in `--help`; defaults are almost always fine.

| Option | Default | Meaning |
| --- | --- | --- |
| `--cc-par0` / `--cc-gra0` | `1e-7` / `1e-6` | Null-model parameter / gradient convergence criteria |
| `--cc-logl0` | `5e-5` | Null-model convergence criterion for the change in log-likelihood |
| `--cc-par` / `--cc-gra` | `1e-7` / `1e-6` | Per-SNP refit parameter / gradient convergence criteria |
| `--start-pos` / `--end-pos` | `0` / `-1` | 0-based index of the first SNP and one past the last SNP to test (`-1` = up to the last SNP); output prefix unchanged |
| `--npart-snp` | `20` | Minimum number of chunks the tested SNPs are read in; chunks are kept to about 128 MB of genotype codes (memory only; does not change results) |
| `--num-random-snp` | `2000` | Random SNPs for the approximate-test calibration (capped at the number of SNPs) |
| `--seed` | `42` | Random seed for drawing the calibration SNPs |

Set the environment variable `FASTGMA_LOG_LEVEL=debug` to log per-block progress.

## 6. Examples

### Standard single-GRM LMM scan

```bash
fastgma --uvlmm --model lmm \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --covar age \
  --class sex \
  --out output/test_uvlmm
```

`--model lmm` and `--solver eigen` are the defaults, so they can be omitted.

### Exact test with a per-SNP refit for every SNP

```bash
fastgma --uvlmm \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --p-approx-cut 1 \
  --out output/test_uvlmm_exact
```

### Variance components only (no SNP scan)

```bash
fastgma --uvlmm --model lmm \
  --grm output/test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --out output/test_uvlmm
```

### Single GRM with a non-genetic random effect

```bash
fastgma --uvlmm --model lmm --solver general \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --random litter \
  --out output/test_uvlmm_litter
```

### Repeated records

```bash
fastgma --make-grm --bfile example/mouse_long/plink --out output/mouse_long
fastgma --uvlmm --model repeat \
  --mgrm output/mouse_long \
  --bfile example/mouse_long/plink \
  --data example/mouse_long/phe.balance.txt \
  --trait trait \
  --class Sex week \
  --out output/mouse_long_repeat
```

16 weekly body weights of each of 1,212 mice: the GRM, a permanent effect per mouse and the residual. `--solver mme` is the default for `--model repeat`.

### Multiple GRMs

```bash
fastgma --uvlmm --model multi-grm \
  --mgrm output/grm_add output/grm_dom \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --out output/test_uvlmm_mgrm
```

### Split a scan into 10 parallel jobs

```bash
fastgma --uvlmm --model lmm \
  --grm output/test \
  --bfile test \
  --data test_simu_pheno.txt \
  --trait pheno \
  --parallel 10 1 \
  --out output/test_uvlmm
```

Repeat with part indices `2` through `10`, then concatenate `output/test_uvlmm.10_*.assoc` (keeping one header).
