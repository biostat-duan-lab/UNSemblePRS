This repository provides R code for UNSemblePRS, an unsupervised ensemble learning framework for efficiently integrating pre-trained polygenic risk scores (PRS). For methodological details, see our medRxiv preprint: https://www.medrxiv.org/content/10.1101/2025.01.06.25320058v2

[utils_UNSemblePRS.R] contains the core implementation of UNSemblePRS and includes functions to compute partial R2 for PRS after adjusting for covariates (e.g., sex and genetic principal components) in both continuous and binary trait models.

[example_code.R] provides code for generating example pre-trained models and demonstrates how to use UNSemblePRS to compute the final aggregated PRS.

[eval_UNSemblePRS.R] contains code for the All of Us analyses, including comparisons with competing methods as well as sex- and ancestry-stratified analyses. To protect data privacy, this script includes analysis code only and does not contain any individual-level data.

[AoU_compute_PRS_PGScatalog.ipynb] provides R code to compute polygenic risk scores (PRSs) within the All of Us (AoU) Research Program using pre-trained models obtained from the PGS Catalog.

For more details on computing PRSs using the PGS Catalog, we recommend that users explore "pgsc_calc" (https://github.com/PGScatalog/pgsc_calc), a pipeline for calculating PRSs or PGss using scoring files published in the PGS Catalog and/or custom scoring files.

Reference
Lambert, Wingfield, et al. (2024). Enhancing the Polygenic Score Catalog with tools for score calculation and ancestry normalization. Nature Genetics. https://doi.org/10.1038/s41588-024-01937-x
