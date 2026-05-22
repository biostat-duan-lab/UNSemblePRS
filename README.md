# UNSemblePRS

UNSemblePRS is an R-based framework for unsupervised ensemble learning of polygenic risk scores (PRS). This repository provides utility functions for integrating and evaluating multiple pre-trained PRS models across diverse prediction settings and populations.

For methodological details, please see our medRxiv preprint:

> Chenyin Gao, Yu-Jyun Huang, et al.  
> *Unsupervised Ensemble Learning for Efficient Integration of Pre-trained Polygenic Risk Scores.*  
> https://www.medrxiv.org/content/10.1101/2025.01.06.25320058v2

---

## Requirements

UNSemblePRS is implemented in R and requires the following packages.

### Core dependencies

```r
install.packages(c("kernlab", "sparsepca"))
```

### Additional packages for evaluation utilities

```r
install.packages(c("dplyr", "RISCA"))
```

---

## Installation

UNSemblePRS is currently not available on CRAN. Users should clone the GitHub repository and manually source the main utility file. The primary function for implementing UNSemblePRS is `UNSemblePRS()`.

### Clone the repository

```bash
git clone https://github.com/biostat-duan-lab/UNSemblePRS.git
cd UNSemblePRS
```

### Load the main functions in R

```r
source("utils_UNSemblePRS.R")
```

### Run the example workflow

```r
source("example_code.R")
```

---

## Repository Structure

```text
UNSemblePRS/
├── utils_UNSemblePRS.R
├── example_code.R
├── eval_UNSemblePRS.R
├── AoU_compute_PRS_PGScatalog.ipynb
├── PGS_catalog_info/
├── README.md
└── data/
```

---

## File Descriptions

### `example_code.R`

Provides a complete example workflow for generating example pre-trained PRS models and demonstrates how to use UNSemblePRS to compute the final aggregated PRS.

### `utils_UNSemblePRS.R`

Contains the core implementation of UNSemblePRS, including functions for:

- unsupervised ensemble learning of PRS models,
- PRS aggregation,
- partial R² computation after covariate adjustment,
- evaluation for both continuous and binary traits.

Covariates may include sex, age, and genetic principal components.

### `eval_UNSemblePRS.R`

Contains the analysis code used for the All of Us (AoU) evaluation experiments, including:

- comparisons with competing PRS integration methods,
- sex-stratified analyses,
- ancestry-stratified analyses.

To protect participant privacy, this script contains analysis workflows only and does not include any individual-level data.

### `AoU_compute_PRS_PGScatalog.ipynb`

Provides code for computing polygenic risk scores (PRSs) within the All of Us Research Program using pre-trained scoring files obtained from the PGS Catalog.

### `PGS_catalog_info/`

Contains supplementary metadata files related to the PGS Catalog resources used in this study:

- `pgs_all_metadata_v2024.xlsx`  
  PGS Catalog metadata corresponding to the version used in this study.

- `PGSID_by_phenotype.xlsx`  
  Complete list of PGS IDs evaluated for each phenotype.

- `AoU_PGSID_to_remove.csv`  
  List of PGS IDs excluded from analysis because the corresponding models were trained using All of Us genetic and phenotypic data, which could introduce potential data leakage.

---

## Recommended External Tool

We highly recommend that users explore `pgsc_calc`:

https://github.com/PGScatalog/pgsc_calc

`pgsc_calc` is a pipeline for calculating PRSs/PGSs using scoring files published in the PGS Catalog and/or custom scoring files.

Reference:

> Lambert SA, Wingfield B, et al. (2024).  
> *Enhancing the Polygenic Score Catalog with tools for score calculation and ancestry normalization.*  
> Nature Genetics.  
> https://doi.org/10.1038/s41588-024-01937-x

---


## Citation

If you use UNSemblePRS in your research, please cite the corresponding manuscript when available.

---

## Contact

For questions, bug reports, or feature requests, please open an issue on the GitHub repository:

https://github.com/biostat-duan-lab/UNSemblePRS
