## Repository for: "Effect of Statewide Prenatal Substance Exposure Policy on Infant Maltreatment Report Rates" - Health Affairs Ms. No. 2024-01160

This repository contains the code and data necessary to reproduce the results
presented in the manuscript "Effect of Statewide Prenatal Substance Exposure
Policy on Infant Maltreatment Report Rates" (Health Affairs Ms. No. 2024-01160).
It provides all materials to allow for full computational reproducibility of the
study. Note, however, that the dataset used for the analysis is not public and,
therefore, is not included in the repository.

**Repository Structure:**

*   **`reproducible-code.R`:** The primary R script containing the code for data
    processing, model fitting, and generation of the tables and figures
    presented in the paper.
*   **`stan/`:** This directory contains the Stan model (`.stan` file) used for
    the analysis, along with any necessary helper files.
*   **`data/`:** This directory contains the datasets used in the analysis,
    stored as comma-separated value (`.csv`) files.  These datasets are
    sufficient to reproduce all figures and tables.
*   **`figs/`:** This directory contains high-resolution versions of the figures
    included in the manuscript, stored as Portable Document Format (`.pdf`)
    files.
*   **`renv.lock`:** This file captures the R environment used for the analysis,
    including the R version and the specific versions of all R packages. This
    ensures that the code can be run with the same dependencies used in the
    original analysis. Use the `renv` package to restore the environment
    (`renv::restore()`).

**Software and Package Versions:**

*   **R:** The specific R version is recorded in the `renv.lock` file.
*   **R Packages:** All R package dependencies and their versions are listed in
    the `renv.lock` file.
*   **Stan:** Version 2.32.2 was used.  Instructions for installing Stan can be
    found at
    [https://mc-stan.org/users/interfaces/](https://mc-stan.org/users/interfaces/).
    It is recommended to use CmdStanR
    ([https://mc-stan.org/cmdstanr/](https://mc-stan.org/cmdstanr/)) to
    interface with Stan from R.

**Instructions for Reproduction:**

1.  **Clone this repository:**
    ```bash
    git clone <repository_url>
    ```
    (Replace `<repository_url>` with the actual URL of the repository.)

2.  **Install R and Stan:** Ensure you have `R` and `Stan` (version 2.32.2 or
    compatible) installed.

3.  **Restore the R environment:**  Open the project in `R`. Install the `renv` package if you don't have it (`install.packages("renv")`). Then, restore the project's R environment using:
    ```R
    renv::restore()
    ```

4.  **Run the analysis:** Execute the `reproducible-code.R` script. This will
    fit the models, and generate the figures and tables.
