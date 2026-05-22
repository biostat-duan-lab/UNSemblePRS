# UNSemblePRS

**UNSemblePRS** is an R-based framework for **UN**supervised en**Semble** learning of **P**olygenic **R**isk **S**cores (PRS). It provides utility functions for integrating and evaluating multiple pre-trained PRS models across diverse prediction settings and populations.

For methodological details, please see our medRxiv preprint:

> Chenyin Gao, Yu-Jyun Huang, et al.
> *Unsupervised Ensemble Learning for Efficient Integration of Pre-trained Polygenic Risk Scores.*
> medRxiv (2025). https://www.medrxiv.org/content/10.1101/2025.01.06.25320058v2

---

## Requirements

UNSemblePRS is implemented in **R** (version ≥ 4.0.0 recommended) and depends on the following CRAN packages.

| Purpose | Packages |
| --- | --- |
| Core functionality | `kernlab`, `sparsepca` |
| Evaluation utilities | `dplyr`, `RISCA` |

You can install all required packages in a single step:

```r
install.packages(c("kernlab", "sparsepca", "dplyr", "RISCA"))
```

No additional system dependencies are required beyond a standard R installation.

---

## Installation

UNSemblePRS is **not** currently available on CRAN. To use it, clone this repository and source the main utility file directly.

**1. Clone the repository**

```bash
git clone https://github.com/biostat-duan-lab/UNSemblePRS.git
cd UNSemblePRS
```

**2. Load the core functions in R**

```r
source("utils_UNSemblePRS.R")
```

The primary function is `UNSemblePRS()`.

**3. (Optional) Run the example workflow**

```r
source("example_code.R")
```

---

## Quick Start

```r
# Load the core functions
source("utils_UNSemblePRS.R")

# Run UNSemblePRS by using simulated pre-trained PRS scores (prepared in matrix format)
# (see example_code.R for a fully reproducible example)
ensemble_prs <- UNSemblePRS(prs_matrix)
```

See [`example_code.R`](example_code.R) for a complete, runnable demonstration that generates example pre-trained PRS models and computes the final aggregated PRS.

---

## Repository Structure

```text
UNSemblePRS/
├── README.md
├── utils_UNSemblePRS.R              # Core implementation
├── example_code.R                   # Reproducible example workflow
├── eval_UNSemblePRS.R               # All of Us evaluation analyses
├── AoU_compute_PRS_PGScatalog.ipynb # PRS computation in All of Us
├── PGS_catalog_info/                # PGS Catalog metadata
```

---

## File Descriptions

### `utils_UNSemblePRS.R`
Contains the core implementation of UNSemblePRS, including functions for:
- unsupervised ensemble learning of PRS models,
- PRS aggregation,
- partial R² computation after covariate adjustment,
- evaluation for both continuous and binary traits.

Covariates may include sex, age, and genetic principal components.

### `example_code.R`
Provides a complete example workflow that generates example pre-trained PRS models and demonstrates how to use UNSemblePRS to compute the final aggregated PRS.

### `eval_UNSemblePRS.R`
Contains the analysis code used for the All of Us (AoU) evaluation experiments, including:
- comparisons with competing PRS integration methods,
- sex-stratified analyses,
- ancestry-stratified analyses.

To protect participant privacy, this script contains analysis workflows only and does not include any individual-level data.

### `AoU_compute_PRS_PGScatalog.ipynb`
Provides code for computing polygenic risk scores within the All of Us Research Program using pre-trained scoring files obtained from the PGS Catalog. 

**Note:** This notebook cannot be run directly, as it requires access to individual-level genetic data that is only available within the secure AoU Researcher Workbench. To reproduce this workflow, users must register for an All of Us Researcher Workbench account, set up a workspace, and preprocess the individual-level genetic data within that controlled tier environment. The code is provided as a reference/template; individual-level data are not included in order to protect participant privacy.


### `PGS_catalog_info/`
Contains supplementary metadata files related to the PGS Catalog resources used in this study:
- `pgs_all_metadata_v2024.xlsx` — PGS Catalog metadata corresponding to the version used in this study.
- `PGSID_by_phenotype.xlsx` — Complete list of PGS IDs evaluated for each phenotype.
- `AoU_PGSID_to_remove.csv` — List of PGS IDs excluded from analysis because the corresponding models were trained using All of Us genetic and phenotypic data, which could introduce potential data leakage.

---

## Recommended External Tool

We recommend that users explore [`pgsc_calc`](https://github.com/PGScatalog/pgsc_calc), a pipeline for calculating PRSs/PGSs using scoring files published in the PGS Catalog and/or custom scoring files.

> Lambert, S. A., et al. (2024). *Enhancing the Polygenic Score Catalog with tools for score calculation and ancestry normalization.* **Nature Genetics.** https://doi.org/10.1038/s41588-024-01937-x

---

## Citation

If you use UNSemblePRS in your research, please cite:

> Chenyin Gao, Yu-Jyun Huang, et al. *Unsupervised Ensemble Learning for Efficient Integration of Pre-trained Polygenic Risk Scores.* medRxiv (2025). https://www.medrxiv.org/content/10.1101/2025.01.06.25320058v2

---

## Contact

For questions, bug reports, or feature requests, please [open an issue](https://github.com/biostat-duan-lab/UNSemblePRS/issues) on the GitHub repository.
