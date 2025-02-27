# Version 1.4-0 (CRAN)

#### ***Expected Release Date:*** 2025-03-15

### New Features

-   Implemented a validation check to identify dates with completely missing data in BD_insitu. This feature allows users to detect and visualize dates where all recorded values are NA, preventing the model from processing them. If such dates are found, the system will trigger a warning, ensuring data completeness before running the Random Forest predictions.

<!-- -->

-   Two additional categorical metrics have been added when “training” has a value other than 1. The added metrics are: success ratio (SR), bias score (bscore).

<!-- -->

-   An update of the vignettes was made to address the improvements introduced in the previous versions.
-   Added a workflow within GitHub Actions to calculate code coverage on Windows, ubuntu and macOS.

# Version 1.3-0 (CRAN)

### New Features

-   Removed dependency on 'dplyr' and migrated all code to 'data.table' to ensure efficiency and speed for large data sets.

-   Added functionality to apply point-to-pixel validation. The metrics analyzed are: Pearson correlation coefficient (CC), root mean square error (RMSE), modified Kling-Gupta efficiency (KGE), relative bias (PBIAS), probability of detection (POD), false alarm rate (FAR), critical success index (CSI).

-   Removed dependencies on external libraries for FAR, POD, CSI calculations. Calculations are now performed using R base functions.

-   A complete refactoring of the code has been carried out to improve its efficiency and ease of maintenance.

### Bug Fixed

-   Fixed a bug in 'wet.day' when set to False, rounding was still performed.
-   Fixed 'pboptions' slash bug that caused a new slash to be created at each iteration when setting char = “=”.

# Version 1.2-2 (release-CRAN)

### Bug Fixed

-   We modified the description of the package to meet the corrections suggested by CRAN.

-   Replaced by \dontrun with \donttest. due to the time of execution of the example. (\> 5 seconds)

# Version 1.2-1 (release-CRAN)

### Bug Fixed

-   The word quantile mapping was changed to Quantile Mapping due to CRAN's comment of “ Words possibly misspelled in DESCRIPTION”.

-   A test optimization was performed to address the problem of “the error Executing ‘testthat.R’ [421s/114s] the execution of the R code in ‘testthat.R’ had a CPU time 3.7 times higher than the elapsed time” reported by CRAN.

# Version 1.2-0 (release)

### New Features

-   The versioning used has been modified. The semantic versioning has been migrated to the versioning used by CRAN.
-   Changed the format for compatibility with the S3 method.
-   The input data verification methods were refactored.
-   Updated the data for the examples and internal tests

# Version 1.1.3

### New Features

-   Added support for adjusting the simulated distribution by quantile mapping and nonparametric quantile mapping.
-   We have changed the ranger library to Random Forest to make it compatible with terra predict.
-   Improvements have been made to the code to improve interpretation.

### Bug Fixed

-   Fixed an ID mess error when transforming the .csv with coordinates to a vector.
