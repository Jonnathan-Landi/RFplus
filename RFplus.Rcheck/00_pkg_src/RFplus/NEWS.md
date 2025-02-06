# Version 1.2-1 (release)

### New Features

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
