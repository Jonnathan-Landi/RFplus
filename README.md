# **RFplus:** *Progressive Bias Correction of Satellite Environmental Data* <img src="man/figures/logo_RFplus.png" align="right" width="250"/>

<!-- CRAN:Check -->

[![R-CMD-check](https://github.com/Jonnathan-Landi/RFplus/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Jonnathan-Landi/RFplus/actions/workflows/R-CMD-check.yaml) [![CRAN status](https://www.r-pkg.org/badges/version/RFplus)](https://cran.r-project.org/package=RFplus) [![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/RFplus)](https://www.r-pkg.org/pkg/RFplus)

# Overview

RFplus is an advanced multi-ensemble algorithm developed by Landi et al. (2025) to reduce bias in satellite precipitation products. The primary goal of RFplus is to enhance the accuracy of satellite-based precipitation estimates by leveraging in situ station data and the Random Forest algorithm. Although RFplus is currently designed for precipitation bias correction, ongoing research is investigating its potential applications for other meteorological variables, such as temperature and wind speed. How RFplus Works

RFplus operates through a structured multi-step process to train and apply models for satellite data bias adjustment:

## **Double Ensemble Approach**

### First Model (Model 1):

Trains a Random Forest model using covariates (e.g., satellite products and Digital Elevation Model [DEM]) and in situ station data.

Incorporates additional features such as Euclidean distance and altitude differences, derived from the DEM, as covariates.

### Second Model (Model 2):

Trains on the residuals of Model 1, predicting the differences between observed and predicted values.

Aims to correct residual errors from Model 1 and minimize remaining biases.

### Quantile Mapping Correction (QDM)

Combines predictions from both models.

Applies quantile mapping correction (QDM) to adjust the final distribution of predicted values using the nearest in situ stations.

3.  Output Generation

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

## Installation

`RFplus` is available on CRAN, so, to get the current version you can use:

``` r
install.packages("RFplus")
```

If you want to test the development version, you can do it from [GitHub](https://github.com/Jonnathan-Landi/RFplus):

``` r
if (!require(devtools)) install.packages("devtools") 
library(devtools) 
install_github("Jonnathan-Landi/RFplus")
```

## Installation and Usage

RFplus can be implemented in R an example of RFplus as follows:

``` r
library(RFplus)
library(terra)
library(data.table)

# Load the data
data("BD_Insitu", package = "RFplus")
data("Cords_Insitu", package = "RFplus")

# Load the covariates
Covariates <- list(
 MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "RFplus")),
 CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "RFplus")),
 DEM = terra::rast(system.file("extdata/DEM.nc", package = "RFplus"))
 )

# Apply the RFplus bias correction model
model = RFplus(BD_Insitu, Cords_Insitu, Covariates, n_round = 1, wet.day = 0.1, ntree = 2000, seed = 123, training = 0.8, Rain_threshold = 0.1, method = "RQUANT", ratio = 5, save_model = FALSE, name_save = NULL
)

# Visualize the results
# Precipitation results within the study area
modelo_rainfall = model$Ensamble

# Validation statistic results 
metrics = model$Validation
# Note: In the above example we used 80% of the data for training and 20% for # model validation.  
```

## Limitations and Future Research

While RFplus currently focuses on precipitation bias correction, research is underway to evaluate its capability to reduce bias in other meteorological variables, such as temperature and wind speed. Users should ensure proper preprocessing and adhere to input requirements for optimal results.
