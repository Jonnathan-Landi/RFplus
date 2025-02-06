# RFplus: Multi-Ensemble Algorithm for Reducing Satellite Precipitation Bias <img src="./inst/logos/logo_RFplus.png" align="right" width="150" />

# Overview
RFplus is an advanced multi-ensemble algorithm developed by Landi et al. (2025) to reduce bias in satellite precipitation products. The primary goal of RFplus is to enhance the accuracy of satellite-based precipitation estimates by leveraging in situ station data and the Random Forest algorithm. Although RFplus is currently designed for precipitation bias correction, ongoing research is investigating its potential applications for other meteorological variables, such as temperature and wind speed.
How RFplus Works

RFplus operates through a structured multi-step process to train and apply models for satellite data bias adjustment:

1. Double Ensemble Approach

First Model (Model 1):

Trains a Random Forest model using covariates (e.g., satellite products and Digital Elevation Model [DEM]) and in situ station data.

Incorporates additional features such as Euclidean distance and altitude differences, derived from the DEM, as covariates.

Second Model (Model 2):

Trains on the residuals of Model 1, predicting the differences between observed and predicted values.

Aims to correct residual errors from Model 1 and minimize remaining biases.

2. Quantile Mapping Correction (QDM)

Combines predictions from both models.

Applies quantile mapping correction (QDM) to adjust the final distribution of predicted values using the nearest in situ stations.

3. Output Generation

Produces bias-corrected satellite precipitation maps.

Optionally saves the final corrected maps as NetCDF files.

Key Requirements for RFplus

To ensure the optimal performance of RFplus, the following requirements must be met:

Data Quality:

In situ station data must undergo rigorous quality control and homogenization before use.

Stations should have less than 10% missing data to avoid skewing the quantile correction.

Input Data Consistency:

Covariates must share the same spatial extent and coordinate reference system (CRS).

All covariates should have the same number of layers.

Mandatory Use of DEM:

A DEM is required as an input parameter.

The DEM should have a single layer, which will be replicated to match the number of layers in other covariates.

Covariate Structure:

Covariates should be provided as a list of raster layers.

The classes of all covariates must be consistent.

## [Installation](https://github.com/Jonnathan-Landi/RFplus)

If you want to test the development version, you can do it from [GitHub](https://github.com/Jonnathan-Landi/RFplus):

```R
if (!require(devtools)) install.packages("devtools") 
library(devtools) 
install_github("Jonnathan-Landi/RFplus")
```

## Installation and Usage
RFplus can be implemented in R with the following function:
```R
RFplus = function(Covariates, BD_Insitu, Cords_Insitu, ntree = 2000, threshold = NULL,
                  n_round = NULL, save_model = FALSE, name_save = NULL, seed = 123, ratio = 15) {
  # Implementation details (refer to full documentation for function specifics).
}
```
Example Workflow:

Preprocess your in situ data to ensure quality control.

Prepare covariates, ensuring they share the same extent, CRS, and layer count.

Provide a DEM as part of the covariates list.

Call the RFplus function with your data.

Outputs

The RFplus function produces:

Bias-corrected raster maps of precipitation.

Optional NetCDF files if save_model = TRUE.

## Limitations and Future Research

While RFplus currently focuses on precipitation bias correction, research is underway to evaluate its capability to reduce bias in other meteorological variables, such as temperature and wind speed. Users should ensure proper preprocessing and adhere to input requirements for optimal results.

Citation

Landi et al. (2025). RFplus: A Multi-Ensemble Model for Satellite Precipitation Bias Correction
