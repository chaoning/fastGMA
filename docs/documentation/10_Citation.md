---
layout: page
permalink: /documentation/10_Citation/
title: Citation
description: How to cite fastGMA and the method papers behind its modules
order: 10
---

## 1. Citing fastGMA

There is no paper on fastGMA as a whole yet. If you use fastGMA in published work, please cite:

1. The software, with the version you ran:

   > Chao Ning. fastGMA: Fast genomic mixed-model analysis for genome-wide association studies and genomic prediction of complex traits. Version `<version>`. <https://github.com/chaoning/fastGMA>

   `<version>` is the build hash printed by `fastgma --version` (for example, `fastGMA 27c5b10` gives `27c5b10`).

2. As the paper for fastGMA, Wang et al. (2020), which describes the mixed model with multiple polygenic effects (several GRMs) that the REMMA module is based on:

   > Dan Wang, Hui Tang, Jian-Feng Liu, Shizhong Xu, Qin Zhang, Chao Ning. Rapid epistatic mixed-model association studies by controlling multiple polygenic effects. *Bioinformatics* 36(19):4833–4837, 2020. [doi:10.1093/bioinformatics/btaa610](https://doi.org/10.1093/bioinformatics/btaa610)

When a paper on fastGMA is published, it will replace the second item. **Cite this repository** on the [GitHub repository page](https://github.com/chaoning/fastGMA) gives the second item in APA and BibTeX. BibTeX for the software:

```bibtex
@misc{fastGMA,
  author       = {Ning, Chao},
  title        = "{fastGMA}: Fast genomic mixed-model analysis for genome-wide association studies and genomic prediction of complex traits",
  year         = {2026},
  note         = {Version <version>},
  howpublished = {\url{https://github.com/chaoning/fastGMA}}
}
```

BibTeX for the papers is at the end of this page.

## 2. Method papers

If you use `REMMA` or `lonGWAS`, please also cite the papers that describe their methods. The other modules (GRM, uvLMM, GLMM, mvLMM and fastREML) have no method paper of their own, so cite fastGMA as in section 1.

### REMMA

- Dan Wang, Hui Tang, Jian-Feng Liu, Shizhong Xu, Qin Zhang, Chao Ning. Rapid epistatic mixed-model association studies by controlling multiple polygenic effects. *Bioinformatics* 36(19):4833–4837, 2020. [doi:10.1093/bioinformatics/btaa610](https://doi.org/10.1093/bioinformatics/btaa610)
- Chao Ning, Dan Wang, Huimin Kang, Raphael Mrode, Lei Zhou, Shizhong Xu, Jian-Feng Liu. A rapid epistatic mixed-model association analysis by linear retransformations of genomic estimated values. *Bioinformatics* 34(11):1817–1825, 2018. [doi:10.1093/bioinformatics/bty017](https://doi.org/10.1093/bioinformatics/bty017)

### lonGWAS

- Chao Ning, Dan Wang, Lei Zhou, Julong Wei, Yuanxin Liu, Huimin Kang, Shengli Zhang, Xiang Zhou, Shizhong Xu, Jian-Feng Liu. Efficient multivariate analysis algorithms for longitudinal genome-wide association studies. *Bioinformatics* 35(23):4879–4885, 2019. [doi:10.1093/bioinformatics/btz304](https://doi.org/10.1093/bioinformatics/btz304)
- Chao Ning, Dan Wang, Xianrui Zheng, Qin Zhang, Shengli Zhang, Raphael Mrode, Jian-Feng Liu. Eigen decomposition expedites longitudinal genome-wide association studies for milk production traits in Chinese Holstein. *Genetics Selection Evolution* 50:12, 2018. [doi:10.1186/s12711-018-0383-0](https://doi.org/10.1186/s12711-018-0383-0)
- Chao Ning, Huimin Kang, Lei Zhou, Dan Wang, Haifei Wang, Aiguo Wang, Jinluan Fu, Shengli Zhang, Jianfeng Liu. Performance gains in genome-wide association studies for longitudinal traits via modeling time-varied effects. *Scientific Reports* 7:590, 2017. [doi:10.1038/s41598-017-00638-2](https://doi.org/10.1038/s41598-017-00638-2)

### BibTeX

```bibtex
@article{Wang2020REMMA,
  author  = {Wang, Dan and Tang, Hui and Liu, Jian-Feng and Xu, Shizhong and Zhang, Qin and Ning, Chao},
  title   = {Rapid epistatic mixed-model association studies by controlling multiple polygenic effects},
  journal = {Bioinformatics},
  volume  = {36},
  number  = {19},
  pages   = {4833--4837},
  year    = {2020},
  doi     = {10.1093/bioinformatics/btaa610}
}

@article{Ning2018REMMA,
  author  = {Ning, Chao and Wang, Dan and Kang, Huimin and Mrode, Raphael and Zhou, Lei and Xu, Shizhong and Liu, Jian-Feng},
  title   = {A rapid epistatic mixed-model association analysis by linear retransformations of genomic estimated values},
  journal = {Bioinformatics},
  volume  = {34},
  number  = {11},
  pages   = {1817--1825},
  year    = {2018},
  doi     = {10.1093/bioinformatics/bty017}
}

@article{Ning2019lonGWAS,
  author  = {Ning, Chao and Wang, Dan and Zhou, Lei and Wei, Julong and Liu, Yuanxin and Kang, Huimin and Zhang, Shengli and Zhou, Xiang and Xu, Shizhong and Liu, Jian-Feng},
  title   = {Efficient multivariate analysis algorithms for longitudinal genome-wide association studies},
  journal = {Bioinformatics},
  volume  = {35},
  number  = {23},
  pages   = {4879--4885},
  year    = {2019},
  doi     = {10.1093/bioinformatics/btz304}
}

@article{Ning2018lonGWAS,
  author  = {Ning, Chao and Wang, Dan and Zheng, Xianrui and Zhang, Qin and Zhang, Shengli and Mrode, Raphael and Liu, Jian-Feng},
  title   = {Eigen decomposition expedites longitudinal genome-wide association studies for milk production traits in {Chinese Holstein}},
  journal = {Genetics Selection Evolution},
  volume  = {50},
  pages   = {12},
  year    = {2018},
  doi     = {10.1186/s12711-018-0383-0}
}

@article{Ning2017lonGWAS,
  author  = {Ning, Chao and Kang, Huimin and Zhou, Lei and Wang, Dan and Wang, Haifei and Wang, Aiguo and Fu, Jinluan and Zhang, Shengli and Liu, Jianfeng},
  title   = {Performance gains in genome-wide association studies for longitudinal traits via modeling time-varied effects},
  journal = {Scientific Reports},
  volume  = {7},
  pages   = {590},
  year    = {2017},
  doi     = {10.1038/s41598-017-00638-2}
}
```
