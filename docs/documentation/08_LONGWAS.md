---
layout: page
permalink: /documentation/08_LONGWAS/
title: lonGWAS
description: Longitudinal data analysis with random regression models
order: 8
---

`LonGWAS` is the longitudinal GWAS module of `fastGMA`. It analyses a trait measured repeatedly on the same individuals over time, such as body weight recorded every week. It first fits a **random regression null model** to estimate the genetic and permanent environmental covariance functions. It then tests every SNP for an effect on the **trajectory** of the trait, not only on its mean.

It is invoked through the top-level dispatcher with `--longwas`:

```bash
fastgma --longwas [options]
```

Run `fastgma -h 7` (or `fastgma --longwas --help`) for the full option list.

**Quick start.** With a GRM and a long-format data file (one row per record):

```bash
fastgma --make-grm --bfile example/mouse_long/plink --out output/mouse_long
fastgma --longwas --bfile example/mouse_long/plink --grm output/mouse_long \
  --data example/mouse_long/phe.unbalance.txt --trait trait --tpoint week --class Sex \
  --out output/mouse_longwas
```

This writes the variance components (`.var`), the heritability at each time point (`.h2`), the fixed effects (`.sln`) and one association test per SNP (`.assoc`); `--h2-cum` adds the heritability of the trait summed over a period (`.h2cum`). Section 7 has more examples.

## 1. Model

Let $$y_{ij}$$ be record $$j$$ of individual $$i$$, taken at time $$t_{ij}$$. Time is mapped to $$[-1, 1]$$ and expanded in normalized Legendre polynomials $$\phi_k$$. `LonGWAS` fits

$$
\begin{aligned}
y_{ij} = {} & \textstyle\sum_{k=0}^{f} b_k \phi_k(t_{ij}) && \text{fixed mean curve} \\
& + c && \text{time-independent fixed effects} \\
& + \textstyle\sum_{k=0}^{f} s_k x_i \phi_k(t_{ij}) && \text{SNP effect on the curve (tested)} \\
& + \textstyle\sum_{k=0}^{a} u_{ik} \phi_k(t_{ij}) && \text{additive genetic random regression} \\
& + \textstyle\sum_{k=0}^{p} \mathit{pe}_{ik} \phi_k(t_{ij}) && \text{permanent environmental random regression} \\
& + e_{ij}
\end{aligned}
$$

$$
u \sim N(0, K_a \otimes G), \qquad \mathit{pe} \sim N(0, K_p \otimes I), \qquad e \sim N(0, \sigma^2_e I)
$$

where $$x_i$$ is the centred genotype of individual $$i$$, $$G$$ is the additive GRM, and $$f$$, $$a$$ and $$p$$ are the orders set by `--forder`, `--aorder` and `--porder`. $$K_a$$ and $$K_p$$ are the $$(a+1) \times (a+1)$$ and $$(p+1) \times (p+1)$$ covariance matrices of the genetic and permanent environmental coefficients, and $$\sigma^2_e$$ is one residual variance for all records. `--time-covar` and `--time-class` add further fixed curves to the mean curve, and $$c$$ holds the time-independent fixed effects of `--covar` and `--class`.

The Legendre basis is

$$
\phi_k(x) = \sqrt{\frac{2k + 1}{2}}\, P_k(x), \qquad x = \frac{2(t - t_{\min})}{t_{\max} - t_{\min}} - 1
$$

with $$P_k$$ the Legendre polynomial of degree $$k$$, and $$t_{\min}$$, $$t_{\max}$$ the smallest and largest time point among the analysed records. Because time is rescaled to $$[-1, 1]$$, only the relative spacing of the time points matters: weeks 1–16 and days 7–112 give the same fit.

### Covariance functions and heritability

With $$\phi_a(t) = (\phi_0(t), \ldots, \phi_a(t))^\top$$, the additive genetic variance at time $$t$$ is $$\phi_a(t)^\top K_a \phi_a(t)$$, and the genetic covariance between two time points is $$\phi_a(t_1)^\top K_a \phi_a(t_2)$$ (likewise for the permanent environment with $$K_p$$). `LonGWAS` writes $$K_a$$, $$K_p$$ and $$\sigma^2_e$$ to `<out>.var`, and the heritability at each observed time point,

$$
h^2(t) = \frac{d\, \phi_a(t)^\top K_a \phi_a(t)}{d\, \phi_a(t)^\top K_a \phi_a(t) + \phi_p(t)^\top K_p \phi_p(t) + \sigma^2_e},
$$

to `<out>.h2`. Here $$d$$ is the mean diagonal element of the GRM as fitted (including `--val`) over the records, as in the other modules; it puts the additive genetic variance on the scale of the phenotypic variance of a record. It is written to the log. The genetic correlation between two time points follows from `<out>.var` as

$$
r_g(t_1, t_2) = \frac{\phi_a(t_1)^\top K_a \phi_a(t_2)}{\sqrt{\phi_a(t_1)^\top K_a \phi_a(t_1)\; \phi_a(t_2)^\top K_a \phi_a(t_2)}}.
$$

**Cumulative heritability (`--h2-cum`).** For a trait summed over a period, as 305-day milk yield is the sum of the daily yields of a test-day model, let $$S$$ be the $$m$$ observed time points in $$[t_1, t_2]$$, $$w_a = \sum_{t \in S} \phi_a(t)$$ and $$w_p = \sum_{t \in S} \phi_p(t)$$. The cumulative breeding value of an individual is $$w_a^\top u_i$$ (the sum of its genetic values at the time points of $$S$$), and

$$
h^2_{\text{cum}} = \frac{d\, w_a^\top K_a w_a}{d\, w_a^\top K_a w_a + w_p^\top K_p w_p + m\, \sigma^2_e},
$$

with independent residuals at the time points. The mean over the period gives the same value, and $$m = 1$$ gives $$h^2(t)$$. The standard errors of $$h^2(t)$$ and $$h^2_{\text{cum}}$$ are both from the delta method with the inverse AI matrix.

## 2. Fitting and testing

### Null model

The null model (no SNP) is fitted once by REML with a weighted AI/EM algorithm (WEMAI). At each iteration it uses the smallest EM weight that keeps $$K_a$$ and $$K_p$$ positive definite and $$\sigma^2_e$$ positive; if no blend does, the EM step is halved until the update is valid. The fit stops when one of the three convergence criteria (`--cc-par0`, `--cc-gra0`, `--cc-logl0`) is met; if it does not converge within `--maxiter0` iterations, the log warns and the last estimates are still written and used.

**Starting values.** By default, the diagonals of $$K_a$$ and $$K_p$$ and the residual variance start at the phenotypic variance divided by `aorder + porder + 3`, and the covariances at `0`. With `--continue`, the `variance` column of an existing `<out>.var` is used instead, provided the file has the expected rows for the chosen `--aorder`/`--porder` and gives positive definite matrices; otherwise the defaults are used with a warning. Without `--continue`, an existing `<out>.var` is not read and is overwritten.

### Solvers

The two solvers read different data layouts (Section 4):

| `--solver` | Data file | Method |
| --- | --- | --- |
| `general` (default) | **Long format**: one row per record (`--trait`, `--tpoint`). Individuals may have different numbers of records and different time points | REML on the mixed model equations, with the permanent-environment equations absorbed |
| `eigen` | **Wide format**: one row per individual and one column per time point (`--trait w1:w16`, `--tpoint-values`). Every individual needs a value in every time column | The GRM is eigendecomposed once and the data rotated, so the covariance matrix becomes block diagonal with $$T \times T$$ blocks ($$T$$ time points) |

On the same balanced data both solvers give the same estimates and p-values. `eigen` is faster and needs less memory: on the example data (1,212 mice, 16 weeks, 4 threads) it takes 2.8 s and 124 MB against 15.9 s and 585 MB for `general` on the long version of the same data. `--solver eigen` does not allow `--predict` or `--use-vmat`; it is the only solver that supports `--refit`.

**Data size.** `general` keeps two dense matrices whose size grows with the square of the number of individuals: the absorbed mixed-model equations, of dimension $$p_X + (a+1)N$$ ($$p_X$$ fixed-effect columns, $$N$$ individuals in the GRM), and $$F'PF$$ for the SNP test, of dimension $$(f+1)n$$ ($$n$$ individuals with records). Each takes $$8 d^2$$ bytes for dimension $$d$$; with the default orders this is about 7 GB each for 10,000 individuals. The number of records hardly matters. `eigen` needs $$O(n^2)$$ memory for the GRM and its eigenvectors only.

### SNP test

The SNP scan uses the null-model variance components (unless `--refit` is given). For each SNP the $$f + 1$$ regression coefficients $$s_0, \ldots, s_f$$ are estimated jointly, and $$H_0\colon s_0 = \cdots = s_f = 0$$ is tested with a Wald chi-square test with $$f + 1$$ degrees of freedom, so a SNP that changes the level, the slope or the shape of the curve can be detected.

The exact test needs $$F'PF$$ for the SNP design $$F$$; `LonGWAS` approximates it by $$\Gamma\, x'x$$, where $$x$$ is the centred genotype and the $$(f+1) \times (f+1)$$ matrix $$\Gamma$$ is calibrated on `--num-random-snp` random SNPs (drawn with `--seed`, default `42`) that pass the filters and give a well-defined test. $$F'PF / x'x$$ does not depend on the phenotype, so no p-value filter is applied (as in `UVLMM` and `MVLMM`). At least 10% of the random SNPs must qualify, or the run stops. SNPs with approximate p below `--p-approx-cut` are re-tested exactly; `--p-approx-cut 1` gives the exact test for every SNP. With `--solver eigen --refit`, each exactly tested SNP also gets its own REML fit of the variance components with the SNP in the model.

## 3. Model terms

| Term | Set by | Default | Description |
| --- | --- | --- | --- |
| Mean curve | `--forder` | `2` | Fixed regression of order $$f$$ on time (always included) |
| Time-varying fixed effects | `--time-covar`, `--time-class` | none | Each covariate, or each level of a class variable, gets its own fixed curve of order $$f$$ |
| Time-independent fixed effects | `--covar`, `--class` | none | One effect per covariate or class level, constant over time |
| SNP effect | `--forder` | `2` | Fixed regression of order $$f$$, tested with $$f + 1$$ degrees of freedom |
| Additive genetic effect | `--aorder` | `2` | Random regression of order $$a$$, covariance $$K_a \otimes G$$ |
| Permanent environmental effect | `--porder` | `2` | Random regression of order $$p$$, covariance $$K_p \otimes I$$ |
| Residual | — | — | One homogeneous residual variance |

Order `0` means a constant (intercept-only) term; order `1` adds a linear term, order `2` a quadratic term, and so on. The same `--forder` is used for the fixed curves and for the SNP effect. The highest order must be smaller than the number of distinct time points. The number of (co)variance parameters, $$(a+1)(a+2)/2 + (p+1)(p+2)/2 + 1$$, grows quadratically with `--aorder` and `--porder`.

### Choosing the orders

- **`--aorder`, `--porder`.** Fits with the same `--forder` (and the same fixed effects) can be compared by their REML log-likelihood, printed in the log as `logL = ... at the estimates`, for example with a likelihood-ratio test or AIC. On the example data, `--aorder 1 --porder 1` (7 parameters) gives logL −14,192.5 and the default `--aorder 2 --porder 2` (13 parameters) −13,768.6, so the quadratic random regressions fit far better. REML log-likelihoods with different `--forder` are not comparable, because the fixed effects differ.
- **High orders need much information.** With `--forder 4 --aorder 4 --porder 4` (31 parameters), neither solver converges in 100 iterations on the example data: the estimate of $$K_a$$ has eigenvalues 1.82, 0.078, 0.0099, 0.0019 and 0.0006, so the highest-order coefficients are hardly identified and every iteration needs the EM blend to stay positive definite. Raise the orders step by step and keep those the data support.
- **`--forder`** sets both the shape of the mean curve and the degrees of freedom of the SNP test; more coefficients detect more complex SNP effects but spread the test over more degrees of freedom.

## 4. Inputs

### Data file (`--data`)

The file is whitespace- or tab-separated and its first line is a header. **Columns are selected by name.** A 0-based column number is still accepted when no header column has that name, with a deprecation warning.

**Long format (`--solver general`).** One row per record, so an individual measured at 16 time points has 16 rows. `--trait` names the trait column and `--tpoint` the time column.

```
ID      Sex  week  trait
1419    1    1     4.9
1419    1    2     8.4
1419    1    3     10.5
...
1422    2    1     4.6
```

**Wide format (`--solver eigen`).** One row per individual and one column per time point. `--trait` lists the time columns in time order, either by name or as a range of adjacent columns (`w1:w16`, as in `fastREML`), and `--tpoint-values` gives a file with the time of each of these columns: one number per line, no header, in the same order.

```
ID    Sex  w1   w2   w3    ...  w16          weeks.txt:  1
1419  1    4.9  8.4  10.5  ...  ...                      2
1422  2    4.6  8.3  10.9  ...  ...                      ...
...                                                      16
```

The number of time values must equal the number of `--trait` columns, and the values must be distinct numbers. Every individual needs a value in every time column: a missing value stops the run, naming the individual and the column. Individuals whose time-independent covariates are missing are dropped before this check.

| Option | Meaning |
| --- | --- |
| `--data` | Data file (**required**) |
| `--id-col` | Column of the individual ID (default: the first column) |
| `--trait` | Long format: the trait column, exactly one (**required**). Wide format: the time columns in time order, e.g. `w1:w16` |
| `--tpoint` | Long format: the time point column (**required**); must be numeric |
| `--tpoint-values` | Wide format: file with the time of each `--trait` column, one value per line, no header (**required**) |
| `--covar` | Time-independent continuous covariates |
| `--class` | Time-independent class variables |
| `--time-covar` | Continuous covariates with a time-varying effect |
| `--time-class` | Class variables with a time-varying effect |
| `--missing-data` | Tokens treated as missing in the data file (default `NA NaN nan N/A ...`) |

Notes:

- Records whose ID is not in the GRM are dropped. IDs are matched to the IID of `.grm.id`.
- Long format: the records do not need to be sorted. A record with a missing value in any used column is dropped. A non-numeric time point that is not a missing token stops the run.
- Long format: time points are compared exactly, so values that differ only by rounding (for example `0.5` and `0.5000000001` from an export) are different time points. `LonGWAS` warns when two time points differ by less than $$10^{-6}$$ of the time range; round the time column in that case.
- Dependent fixed-effect columns (for example the last level of a class variable, which is confounded with the mean curve) are removed automatically and reported in the log.

### GRM input (`--grm`)

| Option | Meaning |
| --- | --- |
| `--grm` | Prefix of the additive GRM (**required**) |
| `--val` | Value added to the GRM diagonal to keep it positive definite (default `0.001`) |

As elsewhere in `fastGMA`, the GRM is referenced by prefix (`--grm output/test` reads `output/test.grm.bin` and `output/test.grm.id`). Build it with `fastgma --make-grm`; GCTA GRMs are also accepted. See [GRM](../03_GRM/). `--val` is added to the GRM diagonal when the GRM is read, for both solvers. If the GRM is still not positive definite (for example duplicated or near-identical individuals), the run stops and suggests a larger `--val`.

### Genotype input (`--bfile` / `--ffile`)

| Option | Meaning |
| --- | --- |
| `--bfile` | PLINK binary genotype prefix (`.bed/.bim/.fam`). The SNP effect is estimated per copy of the `a1` allele (column 5 of the `.bim` file) |
| `--ffile` | Feature-matrix prefix (`.fam/.bim/.fmat`). Each `.fmat` row is one feature and each column one sample, whitespace-separated, for example imputed dosages. Features are centred by their mean and tested like SNPs; the `--maf` filter does not apply, and `af` in `<out>.assoc` is the feature mean |
| `--missing-geno` | Tokens treated as missing in the `.fmat` file (default `NA NaN nan N/A ...`; PLINK BED encodes missing calls itself) |

Give at most one of `--bfile` and `--ffile`. Without either, only the null model is fitted and no `.assoc` is written.

## 5. Outputs

All files share the `--out` prefix and are tab-separated with a header row.

| File | Contents |
| --- | --- |
| `<out>.var` | Estimated (co)variance components of the null model with standard errors |
| `<out>.h2` | Heritability at each observed time point, with standard errors |
| `<out>.h2cum` | With `--h2-cum`: cumulative heritability over each requested period, with standard errors |
| `<out>.sln` | Fixed-effect estimates; with `--predict`, also the predicted random regression coefficients |
| `<out>.assoc` | One row per SNP with the association results (only with `--bfile` or `--ffile`) |
| `<out>.log` | Log of the run (appended; see below for `--snp-set` / `--parallel`) |

With `--snp-set` the output prefix becomes `<out>.<start>_<end>` (0-based start, exclusive end); with `--parallel m i` it becomes `<out>.<m>_<i>`. The suffix applies to every file of that job (`.var`, `.h2`, `.h2cum`, `.sln`, `.assoc` and `.log`).

### `<out>.var`

The lower triangle of $$K_a$$ (`G`), then the lower triangle of $$K_p$$ (`permanent`), then the residual variance. `row` and `col` are 1-based indices of the Legendre coefficients (`1` = $$\phi_0$$, `2` = $$\phi_1$$, ...). Rows with `row = col` are variances; the others are covariances. Standard errors come from the inverse AI matrix; values that are not finite print as `NA`. With the default orders, from the [standard scan](#standard-scan) example:

```
component	row	col	variance	se
G	1	1	1.50281	0.378199
G	2	1	0.519915	0.148035
G	2	2	0.241273	0.0679531
G	3	1	-0.202554	0.0586236
G	3	2	-0.0605725	0.022651
G	3	3	0.0397128	0.0117012
permanent	1	1	7.89243	0.365713
permanent	2	1	3.07146	0.150173
permanent	2	2	1.56857	0.0756402
permanent	3	1	-0.923908	0.0588248
permanent	3	2	-0.315397	0.0253611
permanent	3	3	0.155919	0.0147515
residual	1	1	1.46686	0.0190305
```

### `<out>.h2`

One row per distinct observed time point, in increasing order; `h2` is defined in Section 1, and the SE is from the delta method with the inverse AI matrix. From the same run:

```
time	h2	se
1	0.044509	0.0156255
2	0.0559842	0.0156134
3	0.0776433	0.0188936
...
16	0.121616	0.0295892
```

### `<out>.h2cum`

One row per `--h2-cum` interval: `from to n_times h2 se var_g var_pe var_e`, where `n_times` is the number of observed time points in the interval and `var_g` ($$d\, w_a^\top K_a w_a$$), `var_pe` and `var_e` ($$m \sigma^2_e$$) are the variance components of the summed trait. From the [cumulative heritability](#cumulative-heritability) example:

```
from	to	n_times	h2	se	var_g	var_pe	var_e
1	16	16	0.154694	0.0342144	182.8	975.416	23.4697
9	16	8	0.146791	0.033394	82.1053	465.494	11.7349
```

### `<out>.sln`

```
type	term	levels	effect	se	rel
```

`rel` is written with `--predict` only.

| `type` | `term` | `levels` |
| --- | --- | --- |
| `0` | Mean curve: `mu_t0`, `mu_t1`, ... (one per Legendre coefficient) | `1` |
| `1` | Covariate name; for `--time-covar`, `<name>_t0`, `<name>_t1`, ... | `1` |
| `2` | Class variable name; for `--time-class`, `<name>_t0`, `<name>_t1`, ... | Class level |
| `3` | Additive genetic coefficient `G_t0`, `G_t1`, ... (with `--predict`) | Individual ID |
| `4` | Permanent environmental coefficient `permanent_t0`, `permanent_t1`, ... (with `--predict`) | Individual ID |

A fixed-effect column removed as dependent is listed with effect `0` and se `NA`. With `--predict`, the `G_tk` rows cover every individual in the GRM, including those without records, and `rel` is the reliability of each coefficient, $$1 - \mathrm{PEV} / (K_a[k,k]\, G_{jj})$$ for `G_tk` and $$1 - \mathrm{PEV} / K_p[k,k]$$ for `permanent_tk`; it is `NA` on the fixed-effect rows. Numbers keep 10 significant digits. The genetic value of an individual at time $$t$$ is $$\sum_k \hat u_{ik} \phi_k(t)$$. Example (default orders, `--class Sex`, from the [prediction](#predict-the-random-regression-coefficients) example):

```
type	term	levels	effect	se	rel
0	mu_t0	1	24.22101506	0.09258007857	NA
0	mu_t1	1	5.592146947	0.03832815573	NA
0	mu_t2	1	-2.072129875	0.01685774373	NA
2	Sex	1	-0.3747445202	0.06230991615	NA
2	Sex	2	0	NA	NA
3	G_t0	1419	0.6258066228	0.5522730352	0.6848795603
...
```

### `<out>.assoc`

For the default `--forder 2`, from the [standard scan](#standard-scan) example:

```
chr	snp	pos	a1	a2	af	N	beta_t0	se_t0	beta_t1	se_t1	beta_t2	se_t2	cov_t1_t0	cov_t2_t0	cov_t2_t1	p
1	UNC010515443	3658709	G	W	0.495462	1212	0.378564	0.228033	0.224081	0.098128	-0.0403298	0.0408635	0.0187006	-0.00664458	-0.00205035	0.143166
...
```

With `--solver eigen --refit`, the columns `cc_par cc_gra` follow `N`. The number of `beta_t*`, `se_t*` and `cov_t*` columns depends on `--forder`.

| Column | Meaning |
| --- | --- |
| `chr`, `snp`, `pos`, `a1`, `a2` | Chromosome, SNP ID, position and alleles from the `.bim` file; effects are per copy of `a1` |
| `af` | Frequency of `a1` in the analysed individuals (with `--ffile`: the feature mean) |
| `N` | Number of analysed individuals with a non-missing genotype |
| `cc_par`, `cc_gra` | Convergence diagnostics of the per-SNP refit (norm of the last parameter change and of the gradient); `nan` for SNPs not tested exactly |
| `beta_tk`, `se_tk` | Estimated SNP effect on Legendre coefficient $$k$$ and its standard error |
| `cov_tk_tl` | Sampling covariance of `beta_tk` and `beta_tl` |
| `p` | Wald chi-square p-value with `forder + 1` degrees of freedom |

**The SNP effect at a time point.** The effect of one copy of `a1` at time $$t$$ is $$\beta(t) = \sum_k \beta_k \phi_k(t)$$, with $$\phi_k$$ the Legendre basis of Section 1 (using $$t_{\min}$$, $$t_{\max}$$ of the analysed records, which are written to the log). Its variance is $$\phi(t)^\top \Sigma_\beta\, \phi(t)$$, where $$\Sigma_\beta$$ has `se_tk`$$^2$$ on the diagonal and `cov_tk_tl` off the diagonal, so the SNP effect curve and a pointwise confidence band follow from one `.assoc` row.

SNPs that fail the `--maf` or `--missing-rate` filter are written with `nan` results. For SNPs that were not re-tested exactly, the effects, SEs, covariances and `p` are the approximate values. With the default `--p-cut 1` every SNP is written, including those with `nan` p; with `--p-cut` below `1`, only SNPs with a finite `p` at or below the cut are written.

## 6. Key options

### Model

| Option | Default | Meaning |
| --- | --- | --- |
| `--solver` | `general` | `general` (long format) or `eigen` (wide format) (Section 2) |
| `--forder` | `2` | Order of the fixed mean curve, the time-varying fixed effects and the SNP effect |
| `--aorder` | `2` | Order of the additive genetic random regression |
| `--porder` | `2` | Order of the permanent environmental random regression |

The orders must be `0` or larger (see [Choosing the orders](#choosing-the-orders)).

### Null model

| Option | Default | Meaning |
| --- | --- | --- |
| `--maxiter0` | `100` | Maximum number of iterations |
| `--continue` | off | Start from the estimates in an existing `<out>.var` |
| `--predict` | off | Also write the predicted random regression coefficients to `<out>.sln` (`--solver general` only) |
| `--h2-cum <from> <to>` | — | Heritability of the trait summed over the observed time points in `[from, to]`, with SE, to `<out>.h2cum`; may be given several times |
| `--threads` | `10` | Number of MKL/OpenMP threads |

### SNP scan

| Option | Default | Meaning |
| --- | --- | --- |
| `--maf` | `0.01` | Minimum minor allele frequency (`--bfile` only) |
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
| `--npart-snp` | `20` | Minimum number of chunks the tested SNPs are read in; chunks are kept to about 128 MB of genotype codes (memory only; does not change results) |
| `--num-random-snp` | `2000` | Random SNPs for the approximate-test calibration; capped at the number of SNPs |
| `--use-vmat` | off | `--solver general`: build the full phenotypic covariance matrix `V` (records × records) for the test instead of using the mixed-model equations. Needs much more memory |
| `--cc-par0` | `1e-7` | Null model: convergence when the norm of the parameter change is below this |
| `--cc-gra0` | `1e-6` | Null model: convergence when the norm of the gradient is below this |
| `--cc-logl0` | `5e-5` | Null model: convergence when the change in log-likelihood is below this |
| `--maxiter` | `30` | Per-SNP refit (`--solver eigen --refit`): max iterations |
| `--cc-par` / `--cc-gra` / `--cc-logl` | `1e-7` / `1e-6` / `5e-5` | Per-SNP refit convergence criteria |

Set the environment variable `FASTGMA_LOG_LEVEL=debug` to log per-chunk progress.

## 7. Examples

The examples use the mouse body-weight data in `example/mouse_long/` (1,212 mice weighed in weeks 1–16; Gray et al. 2015, *Genetics* 201: 213–228) and are run from the repository root:

| File | Format | Contents |
| --- | --- | --- |
| `phe.balance.txt` | Long, `ID Sex week trait` | 19,392 records: every mouse has one record in each of the 16 weeks |
| `phe.unbalance.txt` | Long, `ID Sex week trait` | 15,499 records: a subset of `phe.balance.txt` with about 20% of the records removed, leaving 6–16 records per mouse |
| `phe.wide.txt` | Wide, `ID Sex w1 ... w16` | The data of `phe.balance.txt`, one row per mouse |
| `weeks.txt` | One value per line | The weeks of the columns `w1` ... `w16` (1 to 16), for `--tpoint-values` |

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

This fits quadratic curves (`--forder 2 --aorder 2 --porder 2`, the defaults) with a sex effect that is constant over time. The null model converges in 11 iterations (logL −13,768.6).

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

Without `--bfile` or `--ffile`, `LonGWAS` writes `.var`, `.h2` and `.sln` and skips the SNP scan.

### Balanced data with the eigen solver

```bash
fastgma --longwas --solver eigen \
  --bfile example/mouse_long/plink \
  --grm output/mouse_long \
  --data example/mouse_long/phe.wide.txt \
  --trait w1:w16 \
  --tpoint-values example/mouse_long/weeks.txt \
  --class Sex \
  --threads 4 \
  --out output/mouse_longwas_eigen
```

The wide file holds the same data as `phe.balance.txt`, and the [standard scan](#standard-scan) with `--data example/mouse_long/phe.balance.txt` gives the same variance components and p-values with `--solver general` (Section 2 compares time and memory). Add `--refit` to re-estimate the variance components for the exactly tested SNPs.

A missing value in a time column stops the run and names the first one, for example:

```
[error] The data are not balanced: 1 of the 19392 individual x time-point cells (1212 individuals x 16 time points) have no value; e.g. individual 1422 has none in column w5 (time 5). --solver eigen needs a value in every time column for every individual (after removing individuals with a missing covariate). Fill in or remove the missing values, or use the long format with --solver general.
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

`--time-class Sex` gives each sex its own fixed linear curve, and the SNP test has 2 degrees of freedom. The fit converges slowly: it needs 313 iterations here, so `--maxiter0` is raised from the default `100`. Stopping at 100 iterations gives a warning, but the variance components are then already within 0.3% of the converged values.

### Cumulative heritability

```bash
fastgma --longwas \
  --grm output/mouse_long \
  --data example/mouse_long/phe.unbalance.txt \
  --trait trait \
  --tpoint week \
  --class Sex \
  --h2-cum 1 16 --h2-cum 9 16 \
  --out output/mouse_longwas_cum
```

This adds `<out>.h2cum` with the heritability of the body weight summed over weeks 1–16 and over weeks 9–16 (0.155 and 0.147, against 0.04–0.13 for the single weeks in `<out>.h2`): summing averages out the residual of the single records.

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
