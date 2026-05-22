# UNSemblePRS

UNSemblePRS is an R-based framework for ensemble polygenic risk score (PRS) modeling and evaluation. The repository provides utility functions for integrating and evaluating PRS methods across multiple prediction settings.

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

UNSemblePRS is currently not available on CRAN. Users should clone the GitHub repository and source the main utility file manually.

### Clone the repository

```bash
git clone https://github.com/biostat-duan-lab/UNSemblePRS.git
cd UNSemblePRS
```

### Load the main functions in R

```r
source("utils_UNSemblePRS.R")
```

### Run the example script

```r
source("example_code.R")
```

---

## Repository Structure

```text
UNSemblePRS/
├── utils_UNSemblePRS.R   # Main utility functions
├── example_code.R        # Example workflow
├── README.md             # Documentation
└── data/                 # Example or supporting datasets
```

---

## Notes for Version Control

System-generated files such as `.DS_Store` and `.Rhistory` are not necessary for reproducibility or version control. These files have been removed from the repository and added to `.gitignore`.

### Recommended `.gitignore`

```gitignore
.DS_Store
.Rhistory
.RData
.Rproj.user/
```

---

## Citation

If you use UNSemblePRS in your research, please cite the corresponding manuscript when available.

---

## Contact

For questions or bug reports, please open an issue on the GitHub repository:

https://github.com/biostat-duan-lab/UNSemblePRS





This repository provides R code for UNSemblePRS, an unsupervised ensemble learning framework for efficiently integrating pre-trained polygenic risk scores (PRS). For methodological details, see our medRxiv preprint: https://www.medrxiv.org/content/10.1101/2025.01.06.25320058v2

[example_code.R] provides code for generating example pre-trained models and demonstrates how to use UNSemblePRS to compute the final aggregated PRS.

[utils_UNSemblePRS.R] contains the core implementation of UNSemblePRS and includes functions to compute partial R2 for PRS after adjusting for covariates (e.g., sex and genetic principal components) in both continuous and binary trait models.

[eval_UNSemblePRS.R] contains code for the All of Us analyses, including comparisons with competing methods as well as sex- and ancestry-stratified analyses. To protect data privacy, this script includes analysis code only and does not contain any individual-level data.

[AoU_compute_PRS_PGScatalog.ipynb] provides R code to compute polygenic risk scores (PRSs) within the All of Us (AoU) Research Program using pre-trained models obtained from the PGS Catalog.

[PGS_catalog_info/ folder] contains: (i) PGS Catalog metadata corresponding to the version used in this study (pgs_all_metadata_v2024.xlsx); (ii) the complete list of PGS IDs for all pre-trained models evaluated for each trait (PGSID_by_phenotype.xlsx); and (iii) the PGS IDs excluded from analysis due to potential data leakage, as the corresponding models were trained using All of Us (AoU) genetic and phenotypic data (AoU_PGSID_to_remove.csv)."

We highly recommend that users explore "pgsc_calc" (https://github.com/PGScatalog/pgsc_calc), a pipeline for calculating PRSs or PGSs using scoring files published in the PGS Catalog and/or custom scoring files.

Reference: Lambert, Wingfield, et al. (2024). Enhancing the Polygenic Score Catalog with tools for score calculation and ancestry normalization. Nature Genetics. https://doi.org/10.1038/s41588-024-01937-x
