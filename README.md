This repository provides R code for UNSemblePRS, an unsupervised ensemble learning framework for efficiently integrating pre-trained polygenic risk scores (PRS). For methodological details, see our medRxiv preprint: https://www.medrxiv.org/content/10.1101/2025.01.06.25320058v2

The core implementation is available in "utils_UNSemblePRS.R", which also includes functions to compute partial R2 for PRS after adjusting for covariates (e.g., sex and genetic principal components) in both continuous and binary trait models.

"example_code.R" provides code for generating example pre-trained models and demonstrates how to use UNSemblePRS to compute the final aggregated PRS.

"eval_UNSemblePRS.R" contains code for the All of Us analyses, including comparisons with competing methods as well as sex- and ancestry-stratified analyses. To protect data privacy, this script includes analysis code only and does not contain any individual-level data.
