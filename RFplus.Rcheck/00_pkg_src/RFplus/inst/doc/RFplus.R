## -----------------------------------------------------------------------------
# Load necessary libraries
library(RFplus)
library(terra)
library(data.table)

## -----------------------------------------------------------------------------
# Load the in-situ data and the coordinates
data("BD_Insitu")
data("Cords_Insitu")

# Load the covariates
MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "RFplus"))
CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "RFplus"))
DEM = terra::rast(system.file("extdata/DEM.nc", package = "RFplus"))


## -----------------------------------------------------------------------------
# Adjust individual covariates to the covariate format required by RFplus
Covariates_list = list(MSWEP = MSWEP, CHIRPS = CHIRPS, DEM = DEM)

# Verify the extension and the reference coordinate system of the covariates.
# extension
geom_check = sapply(Covariates_list[-1], function(r) terra::compareGeom(Covariates_list[[1]], r))

if (!all(geom_check)) {
  stop("Error: The spatial geometry of the covariates does not match.")
}

# reference coordinate system
crs_list = list(terra::crs(MSWEP), terra::crs(CHIRPS), terra::crs(DEM))

if (!all(sapply(crs_list, function(x) identical(x, crs_list[[1]])))) {
  stop("Error: The coordinate reference system (CRS) of the covariates does not match.")
}



## -----------------------------------------------------------------------------
# Apply the RFplus bias correction model (example using "QUANT" method)
model_example = RFplus(
  BD_Insitu = BD_Insitu, Cords_Insitu = Cords_Insitu, Covariates = Covariates_list, 
  n_round = 1, wet.day = 0.1, ntree = 2000, seed = 123, method = "QUANT", 
  ratio = 15, save_model = FALSE, name_save = NULL
)

# Other methods such as "RQUANT" and "none" can be used by changing the 'method' argument.


## -----------------------------------------------------------------------------
# First layer of the QUANT method
plot(model_example[[1]])




