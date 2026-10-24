---
layout: page
permalink: /documentation/09_FASTREML/
title: fastREML
description: REML variance components and heritability for G, GxE and noise-by-environment effects
order: 9
---

`fastREML` estimates genetic (G), gene-by-environment (GxE), noise-by-environment (NxE) and residual variance components, and the corresponding heritabilities, directly from a dense GRM. Small samples are fitted exactly; large samples (up to UK Biobank scale) are fitted without ever forming the GRM in memory, by streaming it from disk.

It is invoked through the top-level dispatcher with `--fast-reml`:

```bash
fastgma --fast-reml [options]
```

## 1. What it does

For one quantitative trait and `q` interacting environments `W` (an `n x q` matrix), `fastREML` fits

```
y = X b + g + u_gxe + u_nxe + e
var(y) = V = sg2 K + sgxe2 K_gxe + snxe2 N_E + se2 I
```

where

| Term | Covariance | Meaning |
| --- | --- | --- |
| `g` | `sg2 K` | additive genetic effect; `K` is the GRM |
| `u_gxe` | `sgxe2 K_gxe`, `K_gxe = K o (W W' / q)` | gene-by-environment interaction (`o` = element-wise product) |
| `u_nxe` | `snxe2 N_E`, `N_E = diag(W_i W_i' / q)` | noise-by-environment: residual variance that changes with the environment |
| `e` | `se2 I` | residual |

`X` contains the intercept, the `--covar` and `--class` columns and the environment main effects (`--env-int` columns). The environment columns are standardized (mean 0, variance 1) unless `--no-standardize-env` is given. `--no-noisebye` drops the NxE term, giving a three-component model.

The variance components are estimated by **restricted maximum likelihood (REML)** with the average-information (AI) algorithm. By default, the fixed effects stay in the model:

```
-2 logL = log|V| + log|X' V^-1 X| + y' P y,   P = V^-1 - V^-1 X (X' V^-1 X)^-1 X' V^-1
```

With `--residualize`, the trait is first regressed on `X`, and the residuals are then fitted with `P = V^-1`. This is the fastGxE-compatible variant. Because it ignores the degrees of freedom used by `X`, it gives somewhat smaller variance components. On the mouse example below, with 41 fixed-effect columns, G is 0.873 instead of 0.913.

Heritabilities are the proportions of the phenotypic variance, `Vp = sg2 d_g + sgxe2 d_gxe + snxe2 d_nxe + se2`, where `d_*` is the mean diagonal element of each kernel. The SEs use the delta method with the inverse AI matrix.

## 2. Exact and implicit computation

| Path | Used when | How |
| --- | --- | --- |
| exact | `n <= --exact-threshold` (default 8,000) and no `--force-implicit` | `K` is read into memory; `V` is factorized (Cholesky) and inverted every iteration; traces are exact. |
| implicit | otherwise | `K` is never formed. Products `K B` and `K_gxe B` stream the lower triangle of the GRM from disk (or from RAM with `--in-memory`). `V^-1 B` is solved by block preconditioned conjugate gradients (PCG), traces are estimated with `--num-rand-vec` random probe vectors (Hutchinson estimator), and `log|V|` is estimated by stochastic Lanczos quadrature. |

The implicit path is stochastic: its estimates vary with `--seed` and `--num-rand-vec` within the probe noise. With a fixed seed and thread count it is exactly reproducible. With the default strict REML it also solves for `V^-1 X`, which adds `p` columns (the number of fixed-effect columns) to each PCG solve.

Starting values come from a randomized method-of-moments (Haseman-Elston) estimate using the same number of probes (`<out>.mom`).

**Boundary estimates.** When the REML optimum of a component is zero (for example NxE when there is no noise heterogeneity), the Newton step would make it negative. That component is then held at a tiny positive value (`1e-6 x` the sum of the components), and the other components take the Newton step. The log shows `boundary: NxE = ...` on those iterations.

## 3. Inputs

### Data table (`--data`)

A whitespace- or tab-separated table with a header row. The first column must be the individual ID. Only individuals in both the data table and the GRM, with no missing value in the used columns, are analysed, and each individual must appear once.

| Option | Meaning |
| --- | --- |
| `--trait` | Column name of the trait (exactly one) |
| `--covar` | Continuous covariates (column names; ranges such as `pc1:pc10` are allowed) |
| `--class` | Categorical covariates (one-hot coded) |
| `--env-int` | Interacting environment columns (**required**; ranges allowed, e.g. `e1:e40`). These are also fitted as fixed main effects. They must not also be given in `--covar` or `--class`. |
| `--missing-data` | Tokens treated as missing (default: `NA`, `NaN`, `nan`, `N/A`, …) |

### GRM (`--grm`)

A dense GRM prefix in GCTA format: `<grm>.grm.bin` (float32, or legacy float64; the file size tells which) and `<grm>.grm.id`. GRMs built by `fastgma --make-grm` and by GCTA are both accepted; see [GRM](../03_GRM/).

## 4. Outputs

All files use the `--out` prefix and are tab-separated with a header row. Values that are not defined are written as `nan`.

| File | Columns |
| --- | --- |
| `<out>.var` | `component variance se diag_mean variance_prop h2 h2_se`: one row each for `G`, `GxE`, `NxE` (absent with `--no-noisebye`) and `residual`, then a `total` row whose `h2`/`h2_se` hold the total genetic heritability (G + GxE). `diag_mean` is the mean kernel diagonal and `variance_prop` the proportion of `Vp`. `h2` is defined for the genetic components only. |
| `<out>.mom` | `component mom_variance mom_used diag_mean variance_prop`: the method-of-moments estimate and the positive starting values actually used, plus an `h2_total` row |
| `<out>.iter` | `iter sigma_g2 sigma_gxe2 sigma_nxe2 sigma_e2 neg2logL score_norm delta_norm ai_condition pcg_iters pcg_max_rel_resid`: one row per REML iteration (values at the start of the iteration; `pcg_iters` is 0 and `pcg_max_rel_resid` is `nan` on the exact path) |
| `<out>.diag` | `iid K_diag Kgxe_diag Nxe_diag V_diag`: per-individual kernel diagonals and the diagonal of `V` at the estimates |
| `<out>.samples` | `iid grm_index analysis_index`: the analysed individuals with their GRM row |
| `<out>.log` | the run log |

Example `<out>.var` (mouse example, default strict REML):

```
component	variance	se	diag_mean	variance_prop	h2	h2_se
G	0.912913	0.169888	1.0139	0.241628	0.241628	0.0365331
GxE	0.697459	0.20997	1.01353	0.184534	0.184534	0.0546111
NxE	3.80856e-06	0.619244	1	9.94221e-07	nan	nan
residual	2.1982	0.590571	1	0.573837	nan	nan
total	nan	nan	nan	nan	0.426162	0.0607458
```

## 5. Key options

### Model and likelihood

| Option | Default | Meaning |
| --- | --- | --- |
| `--residualize` | off | fastGxE-compatible: residualize the trait on `X`, then use `P = V^-1` |
| `--no-noisebye` | off | Drop the NxE component |
| `--no-standardize-env` | off | Use the environment columns as given |
| `--log-trait` | off | Log-transform the trait (must be positive) |
| `--int-trait` | off | Rank-based inverse-normal transform of the residualized trait |

### Convergence

| Option | Default | Meaning |
| --- | --- | --- |
| `--maxiter0` | 20 | Maximum REML iterations |
| `--cc-par0` | 1e-4 | Stop when the relative parameter change is below this |
| `--cc-gra0` | 1e-2 | … or the score norm is below this |
| `--cc-logl0` | 1e-2 | … or the change in `-2logL` is below this |

### Computation

| Option | Default | Meaning |
| --- | --- | --- |
| `--threads` | 10 | Number of threads |
| `--exact-threshold` | 8000 | Use the exact path when `n` is at most this |
| `--force-implicit` | off | Always use the implicit path |
| `--num-rand-vec` | 32 | Random probe vectors for the trace and method-of-moments estimates |
| `--seed` | 42 | Seed of the probe vectors |
| `--in-memory` | off | Implicit path: load the kept GRM triangle into RAM once |
| `--grm-single` | off | With `--in-memory`: store the GRM in single precision (half the RAM) |
| `--pcg-tol` / `--pcg-maxiter` | 1e-5 / 500 | PCG relative residual tolerance / iteration limit. An iteration whose PCG solves stop above `--pcg-tol` logs a warning and cannot end REML as converged. |
| `--precond-rank` / `--precond-oversample` | 0 / 10 | Rank (0 = diagonal only) and oversampling of a randomized low-rank preconditioner for `K` |
| `--no-pcg-warmstart` | off | Do not start each iteration's PCG from the previous solution |
| `--logl-lanczos` | 15 | Lanczos steps for the `log|V|` estimate (0 = do not compute `-2logL` on the implicit path) |
| `--matvec-block` | 0 (auto) | Right-hand-side columns per GRM pass; caps implicit-path memory |
| `--dry-run` | off | Read the data, cache the kernel diagonals, time one GRM product, then exit |

## 6. Examples

All examples use the simulated mouse data (`example/mouse/pheno/pheno_gxe_simu.txt`: 1,304 individuals, 40 environments `e1`–`e40`, trait `trait`) and a GRM made with `fastgma --make-grm --bfile example/mouse/geno/plink --out mouse`.

### Default: strict REML, exact path

```bash
fastgma --fast-reml --grm mouse \
  --data example/mouse/pheno/pheno_gxe_simu.txt \
  --trait trait --env-int e1:e40 \
  --maxiter0 50 --threads 2 \
  --out mouse_gxe
```

It converges in 4 iterations (less than 1 s). NxE sits on the boundary, h2_G = 0.242 and h2_total = 0.426.

### fastGxE-compatible residualized fit

```bash
fastgma --fast-reml --grm mouse \
  --data example/mouse/pheno/pheno_gxe_simu.txt \
  --trait trait --env-int e1:e40 \
  --maxiter0 50 --threads 2 --residualize \
  --out mouse_gxe_res
```

### Implicit (streaming) path, as used for large samples

```bash
fastgma --fast-reml --grm mouse \
  --data example/mouse/pheno/pheno_gxe_simu.txt \
  --trait trait --env-int e1:e40 \
  --maxiter0 50 --threads 2 \
  --force-implicit --num-rand-vec 100 \
  --out mouse_gxe_implicit
```

Add `--in-memory` (and `--grm-single`) when the kept GRM triangle fits in RAM. Otherwise every product re-streams the GRM file.
