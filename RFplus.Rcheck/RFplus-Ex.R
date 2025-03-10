pkgname <- "RFplus"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "RFplus-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('RFplus')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BD_Insitu")
### * BD_Insitu

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BD_Insitu
### Title: Precipitation Station Measurement Dataset
### Aliases: BD_Insitu
### Keywords: datasets

### ** Examples

data(BD_Insitu)
## You can use str(BD_Insitu) to get a description of the structure
## or view some of the first rows using head(BD_Insitu)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BD_Insitu", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Cords_Insitu")
### * Cords_Insitu

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Cords_Insitu
### Title: Precipitation Station Coordinates Dataset
### Aliases: Cords_Insitu
### Keywords: datasets

### ** Examples

data(Cords_Insitu)
## You can use str(Cords_Insitu) to get a description of the structure
## or view some of the first rows using head(Cords_Insitu)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Cords_Insitu", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("RFplus")
### * RFplus

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: RFplus
### Title: Machine learning algorithm for fusing ground and satellite
###   precipitation data.
### Aliases: RFplus RFplus.default RFplus.data.table

### ** Examples

## No test: 
# Load the libraries
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
model = RFplus(BD_Insitu, Cords_Insitu, Covariates, n_round = 1, wet.day = 0.1,
        ntree = 2000, seed = 123, training = 1, Rain_threshold = 0.1, method = "RQUANT",
        ratio = 5, save_model = FALSE, name_save = NULL)
# Visualize the results
# Precipitation results within the study area
modelo_rainfall = model$Ensamble
# Validation statistic results
# goodness-of-fit metrics
metrics_gof = model$Validation$gof

# categorical metrics
metrics_cat = model$Validation$categorical_metrics
# Note: In the above example we used 80% of the data for training and 20% for # model validation.
## End(No test)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("RFplus", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
