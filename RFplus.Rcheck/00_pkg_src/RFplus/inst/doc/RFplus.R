## -----------------------------------------------------------------------------
# Load necessary libraries
library(RFplus)
library(terra)
library(data.table)

## -----------------------------------------------------------------------------
# Load the in-situ data and the coordinates
data("BD_Insitu", package = "RFplus")
data("Cords_Insitu", package = "RFplus")

# Load the covariates
MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "RFplus"))
CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "RFplus"))
DEM = terra::rast(system.file("extdata/DEM.nc", package = "RFplus"))


## -----------------------------------------------------------------------------
# Adjust individual covariates to the covariate format required by RFplus
Covariates = list(MSWEP = MSWEP, CHIRPS = CHIRPS, DEM = DEM)

# Apply the RFplus -----------------------------------------------------------
# 1. Define categories to categorize rainfall intensity
Rain_threshold = list(
  no_rain = c(0, 1),
  light_rain = c(1, 5),
  moderate_rain = c(5, 20),
  heavy_rain = c(20, 40),
  violent_rain = c(40, 100)
)
# 2. Apply de the model
model = RFplus(BD_Insitu, Cords_Insitu, Covariates, n_round = 1, wet.day = 0.1 , ntree = 2000, seed = 123, training = 0.8, Rain_threshold = Rain_threshold, method = "RQUANT", ratio = 5, save_model = FALSE, name_save = NULL
)

## -----------------------------------------------------------------------------
# Precipitation results within the study area
modelo_rainfall = model$Ensamble

# Validation statistic results 
# goodness-of-fit metrics
metrics_gof = model$Validation$gof

# Categorical Metrics
metrics_categoricxal = model$Validation$categorical_metrics
# Note: In the above example we used 80% of the data for training and 20% for model validation. 

## -----------------------------------------------------------------------------
# First layer of the QUANT method
plot(modelo_rainfall[[1]])


## -----------------------------------------------------------------------------
Rain_threshold = list(
  no_rain = c(0, 1), # No precipitation
  light_rain = c(1, 5), # Light rainfall
  moderate_rain = c(5, 20), # Moderate rainfall
  heavy_rain = c(20, 40), # Heavy rainfall
  violent_rain = c(40, 100) # violent rain
)

