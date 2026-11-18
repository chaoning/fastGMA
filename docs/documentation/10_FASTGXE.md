---
layout: page
permalink: /documentation/10_FASTGXE/
title: fastGxE
description: Genome-wide genotype-by-environment (GxE) association analysis at biobank scale, a separate program
order: 10
---

[fastGxE](https://chaoning.github.io/fastGxE/) is a separate program for genome-wide genotype-by-environment (GxE) association analysis at biobank scale. It tests one SNP at a time for interactions with many environmental factors together. It controls for polygenic effects, polygenic-by-environment effects and noise that differs between environments, so its *p*-values stay calibrated. For each detected locus, it uses mmSuSiE to identify the environmental factors that drive the interaction.

fastGxE is not part of `fastgma`; it has its own executable and documentation at <https://chaoning.github.io/fastGxE/>. Within fastGMA, [fastREML](../09_FASTREML/) estimates the variance components of the same model (G, GxE, noise-by-environment and residual).

Reference: Chao Ning, Xiang Zhou. Genome-wide detection of genotype–environment interactions at biobank scale with fastGxE. *Nature Genetics*, 2026. [doi:10.1038/s41588-026-02809-2](https://doi.org/10.1038/s41588-026-02809-2)
