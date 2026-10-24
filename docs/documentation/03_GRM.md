---
layout: page
permalink: /documentation/03_GRM/
title: GRM
description: Build and process genomic relationship matrices
order: 3
---

The GRM module builds the genomic relationship matrices (GRMs) used as random-effect covariance matrices by `UVLMM`, `REMMA`, `fastREML` and the other mixed-model modules. It also post-processes existing GRMs. The files it writes use the GCTA format, so a GRM built by `fastGMA` can be used in GCTA, and a GRM built by GCTA can be used in `fastGMA`.

There are two entry points:

| Command | Short form | Purpose |
| --- | --- | --- |
| `fastgma --make-grm [options]` | `-m` | Compute a GRM from genotypes or a feature matrix |
| `fastgma --process-grm [options]` | `-d` | Merge, convert, invert, decompose or filter an existing GRM |

As elsewhere in `fastGMA`, a GRM is referred to by its prefix: `--grm output/test` means the files `output/test.grm.bin`, `output/test.grm.id` and so on.

## 1. Building a GRM (`--make-grm`)

### GRM types

Choose the matrix with `--grm-type`. Let `x` be the genotype dosage (0/1/2 copies of the A1 allele), `p` the A1 frequency and `q = 1 − p`.

| `--grm-type` | Matrix | Genotype coding for each SNP | GCTA equivalent |
| --- | --- | --- | --- |
| `agrm` (default) | Additive | `x − 2p` | `--make-grm` |
| `dgrm_gs` | Dominance, breeding-value ("statistical") coding | `0 → −2p²`, `1 → 2pq`, `2 → −2q²` | `--make-grm-d` |

Under Hardy–Weinberg equilibrium the dominance code is uncorrelated with the additive code, so fitting `agrm` and `dgrm_gs` together splits the genetic variance into the classical additive and dominance variances.

The scaling is set with `--code-type`:

| `--code-type` | Additive GRM | GCTA equivalent |
| --- | --- | --- |
| `2` (default) | each SNP standardized by its own SD, then averaged: `A_jk = (1/N_jk) Σ_i z_ij z_ik`, with `z = (x − 2p)/√(2pq)` | `--make-grm` |
| `1` | VanRaden-type: `A_jk = Σ_i (x_ij − 2p_i)(x_ik − 2p_i) / Σ_i 2p_i q_i` | `--make-grm-alg 1` |

Both additive scalings are special cases of `--alpha`, which assumes that the variance of a SNP effect is proportional to `(2pq)^α`. With `v_i = 2p_i q_i`, the additive GRM is

`A_jk = Σ_i v_i^α (x_ij − 2p_i)(x_ik − 2p_i) / Σ_i v_i^(1+α)`,

so that the expected diagonal is 1 (with `--ffile`, `v_i` is the observed variance of feature `i`). `--alpha -1` gives `--code-type 2` and `--alpha 0` gives `--code-type 1`, with identical output files. Under `α = −1` every SNP explains the same share of the genetic variance; values above −1 shift the variance toward common SNPs and values below −1 toward rare SNPs. Speed et al. (2017, *Nature Genetics* 49:986–992) recommended `α = −0.25` for human complex traits. To choose α for your data, build GRMs over a grid of values (for example −1, −0.75, −0.5, −0.25, 0) and compare the REML log-likelihoods. `--alpha` cannot be combined with `--code-type` and applies only to `agrm`. With missing genotypes, each element is rescaled by `M / N_jk` as for the other scalings (see "Quality control and missing genotypes" below).

For the dominance GRM the two `--code-type` choices apply, using the variance of the dominance coding instead of `2pq`. An epistatic GRM (for example additive-by-additive) is built from existing GRMs with `--process-grm --make-epis` (see Section 3).

### Inputs

Provide exactly one of:

| Option | Input | Files read |
| --- | --- | --- |
| `--bfile <prefix>` | PLINK binary genotypes | `<prefix>.bed`, `<prefix>.bim`, `<prefix>.fam` |
| `--ffile <prefix>` | Any numeric feature-by-sample matrix (gene expression, methylation, dosages, ...) | `<prefix>.fmat`, `<prefix>.bim`, `<prefix>.fam` |

With `--ffile`:

- `<prefix>.fmat` has one row per feature and one whitespace-separated column per sample, in `.fam` order. It has no header and no row names.
- `<prefix>.bim` needs one line per feature; only the line count is used.
- Each feature is centered by its mean and scaled by its standard deviation, both computed over the observed values. Negative values are fine.
- Missing values are written as any of the `--missing-geno` tokens (default `NA`, `NaN`, `N/A` and similar).
- `--maf` does not apply, and only `--grm-type agrm` is supported.

### Selecting individuals and SNPs

| Option | File or value | Effect |
| --- | --- | --- |
| `--keep <file>` | one individual per line: `IID`, or `FID IID` | use only these individuals |
| `--remove <file>` | same format | leave these individuals out |
| `--extract <file>` | SNP IDs in the first column | use only these SNPs |
| `--exclude <file>` | SNP IDs in the first column | leave these SNPs out |
| `--chr <list>` | chromosomes, e.g. `1-22`, `1,3,5` or `2 7 19` | use only these chromosomes |
| `--not-chr <list>` | same format | leave these chromosomes out |

A line with one column is matched to the IID; a line with two or more columns is matched to FID and IID, as in PLINK. A leading `chr` in chromosome names is ignored. The options can be combined: individuals are kept and then removed, and SNPs are extracted, excluded and then restricted by chromosome. The output keeps the `.fam` order of the individuals. Allele frequencies and missing rates are computed in the selected individuals, so the result is identical to running on a file subset with `plink --keep-allele-order --make-bed`. The SNP options require `--bfile`.

A leave-one-chromosome-out (LOCO) GRM for chromosome `c` is built with `--not-chr c`, and the GRM of one chromosome with `--chr c`.

### Allele frequencies

By default `p` is estimated from the selected individuals. `--read-freq <file>` uses given A1 frequencies instead, for example from a reference population, so that GRMs built on different subsets of individuals are on the same scale. The file needs a header line with the columns `SNP`, `A1` and `FREQ` (other columns are ignored), which is the format written by `--write-freq`:

```text
CHR	SNP	A1	A2	FREQ	NOBS
0	rs1155740	T	C	0.27383863	409
```

`FREQ` is the frequency of the allele in the `A1` column. If that allele is the `.bim` A2, the frequency is converted to `1 − FREQ`. SNPs that are absent from the file, or whose `A1` matches neither `.bim` allele, are excluded and counted in the log. The frequencies are used for centering, scaling and the `--maf` filter; the missing-rate filter still uses the analysed individuals. With the same frequencies, the GRM is identical to GCTA `--update-freq`.

`--write-freq` writes `<out>.freq` with one line per selected SNP: the A1 frequency used (estimated or read) and `NOBS`, the number of individuals with an observed genotype. With `--npart`, only part 1 writes the file. Both options require `--bfile`.

### Quality control and missing genotypes

A SNP (or feature) is excluded when any of the following holds:

- its minor allele frequency `min(p, 1 − p)` is below `--maf` (BED input only)
- its missing rate is above `--missing-rate`
- it has zero variance

Missing genotypes are handled as in GCTA. A missing call contributes nothing to the cross-product, and each pair of samples is scaled by `N_jk`, the number of SNPs observed in both samples. These counts are written to `.grm.N.bin`.

> **Note.** Because each element has its own denominator, a GRM built from data with missing genotypes is not guaranteed to be positive semi-definite (GCTA has the same property). With complete data it always is. If you need to invert such a GRM, use a larger `--val` with `--process-grm --make-inv`.

### Options

| Option | Default | Meaning |
| --- | --- | --- |
| `--bfile` / `--ffile` | — | Input prefix (exactly one is required) |
| `--out` | — | Output prefix (required) |
| `--grm-type` | `agrm` | `agrm` or `dgrm_gs` |
| `--code-type` | `2` | `2`: per-SNP standardization; `1`: VanRaden-type scaling |
| `--alpha` | `-1` | Additive GRM with SNP-effect variance proportional to `(2pq)^α`; excludes `--code-type` |
| `--maf` | `0.01` | Minimum minor allele frequency (BED input) |
| `--missing-rate` | `0.05` | Maximum missing rate per SNP/feature |
| `--keep` / `--remove` | — | Individuals to use / leave out (see above) |
| `--extract` / `--exclude` | — | SNPs to use / leave out (BED input) |
| `--chr` / `--not-chr` | — | Chromosomes to use / leave out (BED input) |
| `--read-freq` | — | A1 frequencies to use instead of estimating them (BED input) |
| `--write-freq` | off | Write `<out>.freq` (BED input) |
| `--val` | `0` | Value added to the diagonal of the output |
| `--npart <n> <i>` | `1 1` | Compute only part `i` of `n` of the lower triangle (see below) |
| `--npart-snp` | `20` | Number of SNP chunks read from disk in turn; raise it to lower memory use |
| `--missing-geno` | `NA NaN N/A ...` | Tokens treated as missing in `--ffile` input |
| `--threads` | `10` | Number of threads |

The full genotype matrix is never held in memory: the SNPs are read in `--npart-snp` chunks. Memory is dominated by the GRM block (about `8 x n x n` bytes for all rows) plus one chunk (about `8 x n x m / npart-snp` bytes for `m` SNPs), so raise `--npart-snp` when `m` is large. For very large samples, split the lower triangle into row blocks with `--npart` and merge them afterwards:

```bash
for i in 1 2 3 4; do
  fastgma --make-grm --bfile test --npart 4 $i --out output/test
done
fastgma --process-grm --merge --grm output/test --npart 4
```

Each job writes `output/test.4_<i>.grm.*` and its own log `output/test.4_<i>.log`. The merge step concatenates the `.grm.id`, `.grm.bin` and `.grm.N.bin` blocks into `output/test.grm.*` and copies `.grm.N` from the first block. The result is identical to a single full run.

## 2. Output files

| File | Content |
| --- | --- |
| `<out>.grm.bin` | Lower triangle including the diagonal, row by row (`A_11, A_21, A_22, A_31, ...`), little-endian float32 |
| `<out>.grm.id` | One line per sample: `FID<TAB>IID`, in matrix row order |
| `<out>.grm.N.bin` | Per-pair SNP counts `N_jk` in the same layout as `.grm.bin` (float32). GCTA needs this file for `--mgrm` merging, `--grm-cutoff`, `--make-bK-sparse` and `--adj-grm` |
| `<out>.grm.N` | Text file: number of SNPs used and the scale factor |
| `<out>.freq` | With `--write-freq`: `CHR SNP A1 A2 FREQ NOBS` for every selected SNP |

These are exactly the files GCTA writes (plus the text `.grm.N`).

Every module that reads a GRM accepts both value widths of `.grm.bin` (float32 as written now, or float64 from older fastGMA versions; the width is detected from the file size) and a `.grm.id` with one column (IID) or two columns (FID, IID). Samples are matched by IID (the second column, or the only column).

Each run also writes a log. It is appended to `<out>.log`, so several steps that share one `--out` prefix collect in one file. A partitioned run (`--npart n i`) logs to `<out>.<n>_<i>.log`. Set the environment variable `FASTGMA_LOG_LEVEL=debug` for more detail.

Reading a GRM in Python:

```python
import numpy as np

ids = [line.split()[-1] for line in open("output/test.grm.id") if line.strip()]
n = len(ids)
tri = np.fromfile("output/test.grm.bin", dtype=np.float32)
assert tri.size == n * (n + 1) // 2

A = np.zeros((n, n))
A[np.tril_indices(n)] = tri
A = A + A.T - np.diag(np.diag(A))
```

## 3. Processing an existing GRM (`--process-grm`)

Pass `--process-grm` together with exactly one operation:

| Operation | What it does | Main output |
| --- | --- | --- |
| `--merge` | Concatenate the `--npart` blocks of a partitioned run | `<grm>.grm.bin/.id/.N.bin/.N` |
| `--pca` | Leading eigenvectors of the GRM | `<out>.eigenval`, `<out>.eigenvec` |
| `--make-epis` | Elementwise product of two GRMs (e.g. additive × additive) | `<out>.grm.bin/.id` |
| `--reformat` | Convert to another storage format (`--out-fmt`) | see the format table below |
| `--make-inv` | Inverse of the GRM after adding `--val` to the diagonal (LAPACK Cholesky) | inverse in the `--out-fmt` format, plus `<out>.grm.id` |
| `--out-close` | Greedily select individuals to drop until no pair exceeds `--cut-value` | `<out>.grm.pruned_iid` |
| `--group` | Split samples into connected groups of relatives (pairs above `--cut-value`) | `<grm>.grm.group`, `<grm>.grm.group.size`, `<grm>.grm.sp.bin` |

| Option | Default | Meaning |
| --- | --- | --- |
| `--grm` | — | Input GRM prefix |
| `--out` | — | Output prefix |
| `--npart` | — | Number of partitions (`--merge`); must be greater than 1 |
| `--prefix1`, `--prefix2` | — | The two input GRMs (`--make-epis`); their `.grm.id` must have the same order |
| `--out-fmt` | `2` | Output format for `--reformat` / `--make-inv` (table below) |
| `--val` | `0.001` | Diagonal value added before `--make-inv` (this default is separate from the `--make-grm` default of `0`) |
| `--num-pc` | `20` | Number of leading principal components written by `--pca`; `0` writes all |
| `--cut-value` | `0.05` | Relatedness threshold for `--out-close` / `--group` |
| `--sparse` | off | The input is a sparse `.grm.sp.bin` (`--reformat`, `--pca`, `--out-close`; not `--make-epis`, `--make-inv` or `--group`) |
| `--threads` | `10` | Number of threads |

Output formats for `--out-fmt`:

| `--out-fmt` | File | Layout |
| --- | --- | --- |
| `0` | `<out>.grm.matrix` | Full `n × n` text matrix |
| `1` | `<out>.grm.index_triplet` | `row col value` per line, lower triangle, 1-based indices |
| `2` | `<out>.grm.iid_triplet` | `IID_row IID_col value` per line, lower triangle |
| `3` | `<out>.grm.bin` (dense) or `<out>.grm.sp.bin` (sparse) | Dense: GCTA float32 lower triangle. Sparse: records of `int64 row, int64 col, float64 value`, 0-based |

Notes:

- **`--merge`.** The blocks are read from `<grm>.<n>_<i>.grm.*` for `i = 1..n` and the merged GRM is written to `<grm>.grm.*`. `.grm.N.bin` is merged when the blocks have it. `--out` is not used, so the log goes to the console only unless `--out` is given.
- **`--pca`.** The files follow the GCTA layout. `<out>.eigenval` has one eigenvalue per line, largest first. `<out>.eigenvec` has one line per sample in `.grm.id` order with the columns `FID IID PC1 PC2 ...` (tab-separated, no header line). With a one-column `.grm.id`, the IID is written in both ID columns. `--num-pc 20` (default) computes only the 20 leading eigenpairs, which is faster than a full decomposition; `--num-pc 0` writes all of them.
- **`--make-epis`.** The two GRMs must have the same `.grm.id` order. The output is float32 `<out>.grm.bin` plus `<out>.grm.id`; `.grm.N` is the concatenation of the two input `.grm.N` files when both exist. No `.grm.N.bin` is written.
- **`--group`.** The output goes next to the input GRM, not to `--out`. Besides the group tables, it writes a sparse GRM (`.grm.sp.bin`) that keeps only within-group pairs.
- **`--make-inv`.** The inverse is computed with a LAPACK Cholesky factorization (`dpotrf`, `dpotri`). A GRM centered by sample allele frequencies has rank at most `n − 1`, so it cannot be inverted without a diagonal value. That is why `--val` defaults to `0.001` here, while `--make-grm` itself adds nothing. With missing genotypes the GRM can have slightly negative eigenvalues (Section 1). If the matrix is not positive definite, `--make-inv` stops, reports the smallest eigenvalue and asks to rerun with `--val` above a given value. With `--out-fmt 3` the inverse is stored as float32, like every `.grm.bin`.

## 4. Compatibility with GCTA

| fastGMA | GCTA |
| --- | --- |
| `--make-grm` | `--make-grm` |
| `--make-grm --code-type 1` | `--make-grm-alg 1` |
| `--make-grm --grm-type dgrm_gs` | `--make-grm-d` |
| `--keep`, `--remove`, `--extract`, `--exclude`, `--chr` | same options |
| `--read-freq` | `--update-freq` |
| `--process-grm --merge` | concatenating `--make-grm-part` outputs |
| `--process-grm --pca` | `--pca` (same `.eigenval` / `.eigenvec` layout) |

With the same SNPs, the values agree with GCTA to float32 precision (largest difference about 1e-7), with or without missing genotypes.

- **GCTA GRMs in fastGMA.** Pass the GCTA prefix directly, for example `--grm gcta_out`.
- **fastGMA GRMs in GCTA.** `--reml`, `--mlma`, `--pca`, `--grm-cutoff`, `--make-bK-sparse` and `--mgrm` all read the files directly.
- **Files from older fastGMA versions** (float64 `.grm.bin`, one-column `.grm.id`) are still read by every module (see Section 2).

Differences to keep in mind when comparing the two tools:

| Item | fastGMA | GCTA |
| --- | --- | --- |
| Default `--maf` | `0.01` | `0` (only monomorphic SNPs removed) |
| Default missing-rate filter | `0.05` | none |
| Chromosomes used | all SNPs in the input | autosomes only by default |
| Partition file names | `<out>.<n>_<i>.grm.*` | `<out>.part_<n>_<i>.grm.*` |
| X chromosome GRM | not supported | `--make-grm-xchr` |
| Feature-matrix input (`--ffile`) | supported | not available |

To reproduce a GCTA GRM exactly, use `--maf 0` (or the same `--maf` in both tools), `--missing-rate 1`, and only autosomal SNPs (`--chr 1-22` for human data).

**LDAK.** With the same SNPs and complete genotypes, `--alpha α` gives the kinship matrix of LDAK `--calc-kins-direct --ignore-weights YES --power α` up to a constant factor: LDAK divides by the observed mean diagonal, so its matrix has mean diagonal exactly 1, whereas fastGMA (like GCTA) divides by the expected value `Σ_i v_i^(1+α)`. After dividing the fastGMA matrix by its mean diagonal, the two agree to float32 precision. With missing genotypes they also differ in the treatment of missing calls (LDAK uses one denominator for all pairs instead of `N_jk`).

## 5. Examples

### Additive GRM for UVLMM

```bash
fastgma --make-grm --bfile test --out output/test
fastgma --uvlmm --grm output/test --bfile test \
  --data test_simu_pheno.txt --trait pheno --out output/test_uvlmm
```

### Leave-one-chromosome-out GRMs

```bash
for c in $(seq 1 22); do
  fastgma --make-grm --bfile test --not-chr $c --out output/test_loco$c
done
```

`output/test_loco<c>` is used when testing the SNPs of chromosome `c`.

### Choosing α

Build a GRM for each value on a grid and compare the REML log-likelihoods of the null model. Without `--bfile`, `--uvlmm` fits only the null model:

```bash
for a in -1 -0.75 -0.5 -0.25 0; do
  fastgma --make-grm --bfile test --alpha $a --out output/alpha_$a
  fastgma --uvlmm --grm output/alpha_$a --data test_simu_pheno.txt --trait pheno --out output/alpha_$a
done
grep "logL (REML" output/alpha_*.log
```

All models have the same fixed effects and the same number of variance components, so the log-likelihoods can be compared directly; the highest one gives the best-supported α. For the example data:

| α | `G` | `residual` | Heritability | REML logL |
| --- | --- | --- | --- | --- |
| −1 | 0.2204 | 0.6551 | 0.2533 | −187.9204 |
| −0.75 | 0.2221 | 0.6534 | 0.2553 | −187.9097 |
| −0.5 | 0.2237 | 0.6518 | 0.2571 | −187.8994 |
| −0.25 | 0.2251 | 0.6504 | 0.2587 | −187.8895 |
| 0 | 0.2263 | 0.6491 | 0.2601 | −187.8802 |

With 427 individuals the log-likelihood changes by only 0.04 over the grid, so these data cannot distinguish the values of α; estimating α needs large samples.

### Additive, dominance and epistatic GRMs for a multi-GRM model

```bash
fastgma --make-grm --grm-type agrm    --bfile test --out output/add
fastgma --make-grm --grm-type dgrm_gs --bfile test --out output/dom
fastgma --process-grm --make-epis --prefix1 output/add --prefix2 output/add --out output/axa

fastgma --uvlmm --model multi-grm --mgrm output/add output/dom output/axa \
  --bfile test --data test_simu_pheno.txt --trait pheno --out output/test_mgrm
```

### GRM from gene expression

```bash
# expr.fmat: genes x samples; expr.fam: samples; expr.bim: one line per gene
fastgma --make-grm --ffile expr --missing-rate 0.2 --out output/expr
```

### PCA, related-sample pruning and inverse

```bash
fastgma --process-grm --pca --grm output/test --out output/test             # top 20 PCs
fastgma --process-grm --pca --num-pc 0 --grm output/test --out output/test  # all PCs
fastgma --process-grm --out-close --grm output/test --cut-value 0.125 --out output/test
fastgma --process-grm --make-inv --grm output/test --out-fmt 3 --out output/test.inv
```

If the genotypes contain missing calls, `--make-inv` may stop with "not positive definite". The message gives the smallest eigenvalue and the `--val` needed; rerun with a larger value, for example `--val 0.02`.

### Use the GRM in GCTA

```bash
gcta64 --grm output/test --pheno test.phen --reml --out output/test_gcta
gcta64 --grm output/test --make-bK-sparse 0.05 --out output/test_sp
```
