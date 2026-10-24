---
layout: page
permalink: /documentation/08_LONGWAS/
title: lonGWAS
description: Longitudinal GWAS with random regression models
order: 8
---

`LonGWAS` is the longitudinal GWAS module of `fastGMA`. It analyses a trait measured repeatedly on the same individuals over time, such as body weight recorded every week. It first fits a **random regression null model** to estimate the genetic and permanent environmental covariance functions. It then tests every SNP for an effect on the **trajectory** of the trait, not only on its mean.

It is invoked through the top-level dispatcher with `--longwas`:

```bash
fastgma --longwas [options]
```

Run `fastgma -h 7` (or `fastgma --longwas --help`) for the full option list.

## 1. What it does

Let `y_ij` be record `j` of individual `i`, taken at time `t_ij`. Time is mapped to `[-1, 1]` and expanded in normalized Legendre polynomials `φ_k`. `LonGWAS` fits

```
y_ij = Σ_{k=0..f} b_k φ_k(t_ij)          fixed mean curve (plus --time-covar / --time-class curves)
     + c                                 time-independent fixed effects (--covar / --class)
     + Σ_{k=0..f} s_k x_i φ_k(t_ij)      SNP effect on the curve (tested)
     + Σ_{k=0..a} u_ik φ_k(t_ij)         additive genetic random regression
     + Σ_{k=0..p} pe_ik φ_k(t_ij)        permanent environmental random regression
     + e_ij

u  ~ N(0, K_a ⊗ G)      K_a: (a+1) × (a+1) genetic covariance matrix of the coefficients
pe ~ N(0, K_p ⊗ I)      K_p: (p+1) × (p+1) permanent environmental covariance matrix
e  ~ N(0, σ²_e I)       one residual variance for all records
```

where `x_i` is the centered genotype of individual `i`, `G` is the additive GRM, and `f`, `a` and `p` are the orders set by `--forder`, `--aorder` and `--porder`.

The Legendre basis is

```
φ_k(x) = sqrt((2k + 1) / 2) · P_k(x),      x = 2 (t − t_min) / (t_max − t_min) − 1
```

with `P_k` the Legendre polynomial of degree `k`, and `t_min`, `t_max` the smallest and largest time point among the analysed records.

The null model (no SNP) is fitted once by REML with a weighted AI/EM algorithm (WEMAI). At each iteration it uses the smallest EM weight that keeps `K_a` and `K_p` positive definite and `σ²_e` positive; if no blend does, the EM step is halved until the update is valid. The SNP scan reuses these variance components (unless `--refit` is given). For each SNP the `f + 1` regression coefficients `s_0 ... s_f` are estimated jointly, and `H0: s_0 = ... = s_f = 0` is tested with a Wald chi-square test with `f + 1` degrees of freedom. A SNP that changes the level, the slope or the shape of the curve can therefore be detected.

The covariance between two time points follows from the coefficient matrices. For example, the additive genetic variance at time `t` is `φ(t)ᵀ K_a φ(t)`, with `φ(t) = (φ_0(t), ..., φ_a(t))ᵀ`. `LonGWAS` writes `K_a`, `K_p` and `σ²_e` to `<out>.var`, and the heritability at each observed time point,

```
h2(t) = d φ_a(t)ᵀ K_a φ_a(t) / (d φ_a(t)ᵀ K_a φ_a(t) + φ_p(t)ᵀ K_p φ_p(t) + σ²_e),
```

to `<out>.h2`. Here `d` is the mean diagonal element of the GRM as fitted (including `--val`) over the records, as in the other modules; it puts the additive genetic variance on the scale of the phenotypic variance of a record. It is written to the log.

### Solvers

| `--solver` | Method | Data |
| --- | --- | --- |
| `general` (default) | REML on the mixed model equations, with the permanent-environment equations absorbed | Any data: individuals may have different numbers of records and different time points |
| `eigen` | The GRM is eigendecomposed once and the data rotated, so `V` becomes block diagonal with `T × T` blocks | Balanced data only: every individual has exactly one record at each of the same `T` time points (after records with missing values are removed). Otherwise the run stops and suggests `--solver general` |

Both solvers read the same long-format input and give the same estimates on balanced data. `--solver eigen` does not allow `--predict` or `--use-vmat`. It is the only solver that supports `--refit`.

### SNP test

For each SNP the test uses the null-model variance components. The exact test needs `F' P F` for the SNP design `F`; `LonGWAS` approximates it by `Γ · s' s`, where `s` is the centered genotype and the `(f+1) × (f+1)` matrix `Γ` is calibrated on `--num-random-snp` random SNPs (drawn with `--seed`, default `42`) that pass the filters and give a well-defined test. `F' P F / s' s` does not depend on the phenotype, so no p-value filter is applied (as in `UVLMM` and `MVLMM`). At least 10% of the random SNPs must qualify, or the run stops. SNPs with approximate p below `--p-approx-cut` are re-tested exactly. `--p-approx-cut 1` gives the exact test for every SNP. With `--solver eigen --refit`, each exactly tested SNP also gets its own REML fit of the variance components with the SNP in the model.


## 2. Model terms

| Term | Set by | Default | Description |
| --- | --- | --- | --- |
| Mean curve | `--forder` | `2` | Fixed regression of order `f` on time (always included) |
| Time-varying fixed effects | `--time-covar`, `--time-class` | none | Each covariate, or each level of a class variable, gets its own fixed curve of order `f` |
| Time-independent fixed effects | `--covar`, `--class` | none | One effect per covariate or class level, constant over time |
| SNP effect | `--forder` | `2` | Fixed regression of order `f`, tested with `f + 1` degrees of freedom |
| Additive genetic effect | `--aorder` | `2` | Random regression of order `a`, covariance `K_a ⊗ G` |
| Permanent environmental effect | `--porder` | `2` | Random regression of order `p`, covariance `K_p ⊗ I` |
| Residual | — | — | One homogeneous residual variance |

Order `0` means a constant (intercept-only) term; order `1` adds a linear term, order `2` a quadratic term, and so on. The same `--forder` is used for the fixed curves and for the SNP effect. The highest order must be smaller than the number of distinct time points.

## 3. Inputs

### Long-format data file (`--data`)

One row per record, so an individual measured at 16 time points has 16 rows. The file is whitespace- or tab-separated and its first line is a header. **Columns are selected by name.** A 0-based column number is still accepted when no header column has that name, with a deprecation warning.

```
ID      Sex  week  trait
1419    1    1     4.9
1419    1    2     8.4
1419    1    3     10.5
...
1422    2    1     4.6
```

| Option | Meaning |
| --- | --- |
| `--data` | Long-format data file (**required**) |
| `--id-col` | Column of the individual ID (default: the first column) |
| `--trait` | Trait column, exactly one (**required**) |
| `--tpoint` | Time point column (**required**); must be numeric |
| `--covar` | Time-independent continuous covariates |
| `--class` | Time-independent class variables |
| `--time-covar` | Continuous covariates with a time-varying effect |
| `--time-class` | Class variables with a time-varying effect |
| `--missing-data` | Tokens treated as missing in the data file (default `NA NaN nan N/A ...`) |

Notes:

- Records whose ID is not in the GRM are dropped. IDs are matched to the IID of `.grm.id`.
- With `--solver general` the records do not need to be sorted, and individuals may have different numbers of records and different time points. `--solver eigen` needs balanced data (Section 1).
- A record with a missing value in any used column is dropped. A non-numeric time point that is not a missing token stops the run.
- Dependent fixed-effect columns (for example the last level of a class variable, which is confounded with the mean curve) are removed automatically and reported in the log.

### GRM input (`--grm`)

| Option | Meaning |
| --- | --- |
| `--grm` | Prefix of the additive GRM (**required**) |
| `--val` | Value added to the GRM diagonal to keep it positive definite (default `0.001`) |

As elsewhere in `fastGMA`, the GRM is referenced by prefix (`--grm output/test` reads `output/test.grm.bin` and `output/test.grm.id`). Build it with `fastgma --make-grm`; GCTA GRMs are also accepted. See [GRM](../03_GRM/). `LonGWAS` adds `--val` (default `0.001`, the same option name as in `--make-grm`/`--process-grm`) to the diagonal of the GRM when it reads it, for both solvers. If the GRM is still not positive definite (for example duplicated or near-identical individuals), the run stops and suggests a larger `--val`.

### Genotype input (`--bfile`)

| Option | Meaning |
| --- | --- |
| `--bfile` | PLINK binary genotype prefix (`.bed/.bim/.fam`). Optional: without it only the null model is fitted and no `.assoc` is written. |

The SNP effect is estimated per copy of the `a1` allele (column 5 of the `.bim` file).

## 4. Outputs

All files share the `--out` prefix and are tab-separated with a header row.

| File | Contents |
| --- | --- |
| `<out>.var` | Estimated (co)variance components of the null model with standard errors |
| `<out>.h2` | Heritability at each observed time point, with standard errors |
| `<out>.sln` | Fixed-effect estimates; with `--predict`, also the predicted random regression coefficients |
| `<out>.assoc` | One row per SNP with the association results (only with `--bfile`) |
| `<out>.log` | Log of the run (appended; see below for `--snp-set` / `--parallel`) |

With `--snp-set` the output prefix becomes `<out>.<start>_<end>` (0-based start, exclusive end); with `--parallel m i` it becomes `<out>.<m>_<i>`. The suffix applies to every file of that job (`.var`, `.h2`, `.sln`, `.assoc` and `.log`).

### `<out>.var`

The lower triangle of `K_a` (`G`), then the lower triangle of `K_p` (`permanent`), then the residual variance. `row` and `col` are 1-based indices of the Legendre coefficients (`1` = `φ_0`, `2` = `φ_1`, ...). Rows with `row = col` are variances; the others are covariances. With the default orders (`--aorder 2 --porder 2`), from the [standard scan](#standard-scan) example:

```
component	row	col	variance	se
G	1	1	1.50282	0.378203
G	2	1	0.519918	0.148037
G	2	2	0.241275	0.0679534
G	3	1	-0.202555	0.0586241
G	3	2	-0.060573	0.0226511
G	3	3	0.039713	0.0117012
permanent	1	1	7.89254	0.36572
permanent	2	1	3.07149	0.150174
permanent	2	2	1.56857	0.0756402
permanent	3	1	-0.923918	0.0588253
permanent	3	2	-0.315397	0.0253611
permanent	3	3	0.155919	0.0147514
residual	1	1	1.46686	0.0190305
```

Standard errors come from the AI matrix. Values that are not finite print as `NA`.

**Starting values.** By default, the diagonals of `K_a` and `K_p` and the residual variance start at the phenotypic variance divided by `aorder + porder + 3`, and the covariances at `0`. With `--continue`, the `variance` column of an existing `<out>.var` is used instead, provided the file has the expected rows for the chosen `--aorder`/`--porder` and gives positive definite matrices; otherwise the defaults are used with a warning. Without `--continue`, an existing `<out>.var` is not read and is overwritten.

### `<out>.h2`

```
time	h2	se
```

One row per distinct observed time point, in increasing order. `h2` is defined in Section 1; the SE is from the delta method with the inverse AI matrix.

### `<out>.sln`

```
type	term	levels	effect	se
```

| `type` | `term` | `levels` |
| --- | --- | --- |
| `0` | Mean curve: `mu_t0`, `mu_t1`, ... (one per Legendre coefficient) | `1` |
| `1` | Covariate name; for `--time-covar`, `<name>_t0`, `<name>_t1`, ... | `1` |
| `2` | Class variable name; for `--time-class`, `<name>_t0`, `<name>_t1`, ... | Class level |
| `3` | Additive genetic coefficient `G_t0`, `G_t1`, ... (with `--predict`) | Individual ID |
| `4` | Permanent environmental coefficient `permanent_t0`, `permanent_t1`, ... (with `--predict`) | Individual ID |

A fixed-effect column removed as dependent is listed with `0 0` in the `effect` and `se` columns. With `--predict`, the `G_tk` rows cover every individual in the GRM, including those without records. Example (default orders, `--class Sex`, from the [prediction](#predict-the-random-regression-coefficients) example):

```
type	term	levels	effect	se
0	mu_t0	1	24.221	0.0925801
0	mu_t1	1	5.59215	0.0383282
0	mu_t2	1	-2.07213	0.0168577
2	Sex	1	-0.374745	0.0623099
2	Sex	2	0	0
3	G_t0	1419	0.625807	0.552273
...
```

### `<out>.assoc`

For the default `--forder 2`:

```
chr	snp	pos	a1	a2	af	N	beta_t0	se_t0	beta_t1	se_t1	beta_t2	se_t2	cov_t1_t0	cov_t2_t0	cov_t2_t1	p
```

With `--solver eigen --refit`, the columns `cc_par cc_gra` follow `N`. The number of `beta_t*`, `se_t*` and `cov_t*` columns depends on `--forder`.

| Column | Meaning |
| --- | --- |
| `chr`, `snp`, `pos`, `a1`, `a2` | Chromosome, SNP ID, position and alleles from the `.bim` file; effects are per copy of `a1` |
| `af` | Frequency of `a1` in the analysed individuals |
| `N` | Number of analysed individuals with a non-missing genotype |
| `cc_par`, `cc_gra` | Convergence diagnostics of the per-SNP refit (norm of the last parameter change and of the gradient); `nan` for SNPs not tested exactly |
| `beta_tk`, `se_tk` | Estimated SNP effect on Legendre coefficient `k` and its standard error |
| `cov_tk_tl` | Sampling covariance of `beta_tk` and `beta_tl` |
| `p` | Wald chi-square p-value with `forder + 1` degrees of freedom |

SNPs that fail the `--maf` or `--missing-rate` filter are written with `nan` results. For SNPs that were not re-tested exactly, the effects, SEs, covariances and `p` are the approximate values. With the default `--p-cut 1` every SNP is written, including those with `nan` p; with `--p-cut` below `1`, only SNPs with a finite `p` at or below the cut are written.

## 5. Key options

### Model

| Option | Default | Meaning |
| --- | --- | --- |
| `--solver` | `general` | `general` or `eigen` (Section 1) |
| `--forder` | `2` | Order of the fixed mean curve, the time-varying fixed effects and the SNP effect |
| `--aorder` | `2` | Order of the additive genetic random regression |
| `--porder` | `2` | Order of the permanent environmental random regression |

The orders must be `0` or larger. The number of (co)variance parameters grows quadratically with `--aorder` and `--porder`.

### Null model

| Option | Default | Meaning |
| --- | --- | --- |
| `--maxiter0` | `100` | Maximum number of iterations |
| `--continue` | off | Start from the estimates in an existing `<out>.var` |
| `--predict` | off | Also write the predicted random regression coefficients to `<out>.sln` (`--solver general` only) |
| `--threads` | `10` | Number of MKL/OpenMP threads |

The fit stops when any one of the three convergence criteria (Advanced options) is met. If it does not converge within `--maxiter0` iterations, the log warns and the last estimates are still written and used.

### SNP scan

| Option | Default | Meaning |
| --- | --- | --- |
| `--maf` | `0.01` | Minimum minor allele frequency |
| `--missing-rate` | `0.05` | Maximum SNP missing rate |
| `--p-approx-cut` | `1e-3` | SNPs with approximate p below this are re-tested exactly; `1` = exact test for every SNP |
| `--refit` | off | `--solver eigen` only: re-estimate the variance components for every exactly tested SNP |
| `--p-cut` | `1` | In `[0, 1]`. Write only SNPs with a finite `p` at or below this; `1` writes every SNP |
| `--seed` | `42` | Random seed for sampling the calibration SNPs |
| `--snp-set <first> <last>` | — | Test SNPs from the first to the last SNP ID (inclusive) |
| `--parallel <m> <i>` | — | Split the SNPs into `m` parts and test part `i` (1-based) |

`--snp-set`, `--parallel` and `--start-pos`/`--end-pos` are mutually exclusive.

### Advanced

| Option | Default | Meaning |
| --- | --- | --- |
| `--start-pos` | `0` | 0-based index of the first SNP to test |
| `--end-pos` | `-1` | 0-based index one past the last SNP to test; `-1` (or a value beyond the last SNP) tests to the end |
| `--npart-snp` | `20` | Number of chunks the SNP range is read in (memory/IO; does not change results) |
| `--num-random-snp` | `2000` | Random SNPs for the approximate-test calibration; capped at the number of SNPs |
| `--use-vmat` | off | `--solver general`: build the full phenotypic covariance matrix `V` (records × records) for the test instead of using the mixed-model equations. Needs much more memory |
| `--cc-par0` | `1e-7` | Null model: convergence when the norm of the parameter change is below this |
| `--cc-gra0` | `1e-6` | Null model: convergence when the norm of the gradient is below this |
| `--cc-logl0` | `5e-5` | Null model: convergence when the change in log-likelihood is below this |
| `--maxiter` | `30` | Per-SNP refit (`--solver eigen --refit`): max iterations |
| `--cc-par` / `--cc-gra` / `--cc-logl` | `1e-7` / `1e-6` / `5e-5` | Per-SNP refit convergence criteria |

Set the environment variable `FASTGMA_LOG_LEVEL=debug` to log per-chunk progress.

## 6. Examples

The examples use the mouse body-weight data in `example/mouse_long/` (1,212 mice weighed in weeks 1–16; Gray et al. 2015, *Genetics* 201: 213–228) and are run from the repository root. Both phenotype files are in long format with the columns `ID Sex week trait`:

| File | Records | Design |
| --- | --- | --- |
| `phe.balance.txt` | 19,392 | Balanced: every mouse has one record in each of the 16 weeks |
| `phe.unbalance.txt` | 15,499 | Unbalanced: about 20% of the records of `phe.balance.txt` removed at random, leaving 6–16 records per mouse |

First build the GRM:

```bash
fastgma --make-grm --bfile example/mouse_long/plink --out output/mouse_long
```

### Standard scan

```bash
fastgma --longwas \
  --bfile example/mouse_long/plink \
  --grm output/mouse_long \
  --data example/mouse_long/phe.unbalance.txt \
  --trait trait \
  --tpoint week \
  --class Sex \
  --threads 4 \
  --out output/mouse_longwas
```

This fits quadratic curves (`--forder 2 --aorder 2 --porder 2`, the defaults) with a sex effect that is constant over time.

### Variance components only

```bash
fastgma --longwas \
  --grm output/mouse_long \
  --data example/mouse_long/phe.unbalance.txt \
  --trait trait \
  --tpoint week \
  --class Sex \
  --out output/mouse_longwas_vc
```

Without `--bfile`, `LonGWAS` writes `.var`, `.h2` and `.sln` and skips the SNP scan.

### Balanced data with the eigen solver

```bash
fastgma --longwas --solver eigen \
  --bfile example/mouse_long/plink \
  --grm output/mouse_long \
  --data example/mouse_long/phe.balance.txt \
  --trait trait \
  --tpoint week \
  --class Sex \
  --out output/mouse_longwas_eigen
```

Every individual must have one record in each week. On these data the `general` solver (the [standard scan](#standard-scan) with `--data example/mouse_long/phe.balance.txt`) gives the same variance components and p-values, but the eigen solver is faster and needs less memory (here 7.6 s and 123 MB against 16 s and 573 MB). Add `--refit` to re-estimate the variance components for the exactly tested SNPs.

The eigen solver stops on unbalanced data and names a missing record. With `--data example/mouse_long/phe.unbalance.txt`:

```
[error] The data are not balanced: 3893 of the 19392 individual x time-point cells (1212 individuals x 16 time points) have no record; e.g. individual 1419 has no record at time 5. --solver eigen needs every individual to have exactly one record at each of the same time points (after removing records with missing values). Use --solver general for these data.
```

### Linear curves with a sex-specific mean curve

```bash
fastgma --longwas \
  --bfile example/mouse_long/plink \
  --grm output/mouse_long \
  --data example/mouse_long/phe.unbalance.txt \
  --trait trait \
  --tpoint week \
  --time-class Sex \
  --forder 1 --aorder 1 --porder 1 \
  --maxiter0 400 \
  --threads 4 \
  --out output/mouse_longwas_linear
```

`--time-class Sex` gives each sex its own fixed linear curve. The SNP test then has 2 degrees of freedom. Linear curves fit these growth data much worse than quadratic curves (null-model logL −21,072 against −13,769), and the fit converges slowly: it needs 314 iterations here, so `--maxiter0` is raised from the default `100`. Stopping at 100 iterations gives a warning, but the variance components already agree with the converged values to three significant digits.

### Predict the random regression coefficients

```bash
fastgma --longwas \
  --bfile example/mouse_long/plink \
  --grm output/mouse_long \
  --data example/mouse_long/phe.unbalance.txt \
  --trait trait \
  --tpoint week \
  --class Sex \
  --predict \
  --continue \
  --threads 4 \
  --out output/mouse_longwas
```

`<out>.sln` then also holds `G_tk` and `permanent_tk` for every individual. `--continue` starts the fit from `output/mouse_longwas.var` of the first example.

### Split a scan into parallel jobs

```bash
for i in 1 2 3 4; do
  fastgma --longwas \
    --bfile example/mouse_long/plink \
    --grm output/mouse_long \
    --data example/mouse_long/phe.unbalance.txt \
    --trait trait \
    --tpoint week \
    --class Sex \
    --parallel 4 $i \
    --p-cut 0.05 \
    --out output/mouse_longwas
done
```

Each job writes `output/mouse_longwas.4_<i>.assoc` with only the SNPs with `p <= 0.05`. Because `--seed` is fixed, all jobs use the same calibration.
