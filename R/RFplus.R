#' RFplus: Bias correction model using Random Forest for satellite-derived environmental data.
#'
#' This function applies a Random Forest-based bias correction method to satellite-derived environmental data, such as precipitation, temperature, and other relevant meteorological or environmental variables.
#' The model corrects biases in these satellite products by incorporating in-situ measurements and environmental covariates, including the Digital Elevation Model (DEM).
#'
#' @param Covariates A list of covariates used as independent variables.
#'                   The selected covariates should be provided as a list of \code{SpatRaster} objects.
#'                   Each covariate must be stored in a separate object within the list.
#'                   Additionally, make sure to include the DEM (Digital Elevation Model).
#' @param BD.Insitu A \code{data.table} file containing in-situ measurements.
#'                  It must include columns such as \code{Date}, and the observed variable(s).
#'                  The \code{Date} column should contain the measurement dates,
#'                  and the column(s) corresponding to the variable(s) should have the observed values.
#' @param BD.Insitu A \code{data.table} file containing the in-situ station coordinates.
#'                  This file must have 3 columns: the first column is \code{Cod}, which indicates the unique code
#'                  for identifying each station. This code should match the one used in the \code{data.table} of the 'BD.Insitu'.
#'                  The second and third columns represent the \code{X} and \code{Y} coordinates of each in-situ station, respectively.
#' @param ntree An integer representing the number of trees to use in the random forest model. Default is 2000.
#' @param threshold A numeric value for setting a threshold for the final predictions. Values below this threshold will be set to 0. Default is NULL.
#' @param save_model A logical value indicating whether to save the final corrected model in NetCDF format.
#'                  Default is \code{FALSE}, which does not save the model.
#'                  If \code{TRUE}, the corrected model will be saved to a NetCDF file."
#'                  It is important that, before executing the function, the directory where the file will be saved is defined.
#'                  using \code{setwd(“path/directory”)}.
#' @param name_save The name under which the model will be saved if \code{save_model = TRUE}.
#'                  By default, the model is saved with the name \code{"Model_RFplus"}.
#'                  The name should be provided without an extension, as the function will append the appropriate file extension (.nc).
#' @param seed An integer value to set the reproducibility seed in the random forest model. Default is "123".
#'
#' @return A \code{terra} raster object containing the corrected environmental data.
#'         The output will have the same spatial resolution and extent as the input covariates.
#' @author Jonnathan Augusto Landi Bermeo
#'
#' @import terra
#' @import data.table
#' @import dplyr
#' @import randomForest
#' @import pbapply
#' @import ncdf4
#'
#' @export

RFplus = function(Covariates, BD_Insitu, Cords_Insitu, ntree = 2000,
                  threshold = NULL, save_model = F, name_save = NULL, seed = 123) {
  ##############################################################################
  #                      Check the input data of the covariates                #
  ##############################################################################
  # Verify the class of covariates
  class_check = sapply(Covariates, class)
  if (length(unique(class_check)) > 1) stop("The class of the covariates are different (all classes should be similar).")

  # Verify the extent of covariates
  ext_check = sapply(Covariates, function(x) terra::ext(x))
  if (!all(sapply(ext_check, function(x) x == ext_check[[1]]))) stop("The extension of the covariates are different (all extensions should be similar).")

  # Verify the crc of covariates
  crs_check = sapply(Covariates, function(x) terra::crs(x))
  if (length(unique(crs_check)) > 1) stop("The crs of the covariates are different (all crs should be similar).")

  ##############################################################################
  #                    Check input data from on-site stations                  #
  ##############################################################################
  # Verify the BD.Insitu is data.table
  if (!is.data.table(BD_Insitu)) stop("The data of the on-site stations should be a data.table.")

  # Verify the columns of Cords_Insitu
  if (!is.data.table(Cords_Insitu)) stop("The coordinate data of the on-site stations must be a data.table.")

  # Check that the coordinate names appear in the observed data
  if (!all(Cords_Insitu$Cod %in% setdiff(colnames(BD_Insitu), "Date"))) stop("The names of the coordinates do not appear in the observed data.")

  ##############################################################################
  #                             Algorithm start                                #
  ##############################################################################
  # DEM identification
  nlyr_covs = sapply(Covariates, nlyr)
  indx_dem = which(nlyr_covs == 1)
  if (length(indx_dem) == 0) stop("The DEM should have a single layer. Check the input DEM.")
  nlyrs_tots = which(nlyr_covs != 1)

  # Equalize DEM layers
  nlyr_rep = nlyr_covs[nlyrs_tots[1]]
  DEM = Covariates[[indx_dem]]
  DEM = rast(replicate(nlyr_rep, DEM))
  Covariates[[indx_dem]] = DEM

  ##############################################################################
  #                           Check the input data on site                     #
  ##############################################################################
  # Training data for randomForest.
  data_train = data.table::melt(
    BD_Insitu,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  ) %>%
    dplyr::mutate(ID = as.numeric(factor(Cod)))

  # Extract Dates
  Dates_extracted = base::unique(data_train$Date)
  Points_Train = data.table::merge.data.table(data_train, Cords_Insitu, by = "Cod")
  Points_Train = base::unique(Points_Train[,.(ID, Cod, X, Y)])

  # Extraction of a sample to save the correction results
  Sample_lyrs = Covariates[[2]][[1]]

  # RandonForest function for bias correction
  RF_train = function(day_COV, fecha) {
    data_obs = data_train[Date == as.Date(fecha), ] # Filtering the observation 'i'
    data_obs = data_obs[!is.na(var)] # Method for training models independently of NAs
    data_obs$ID = 1:nrow(data_obs) # ID for the data.table

    points_EstTrain = data.table::merge.data.table(data_obs, Points_Train, by = "Cod")
    points_EstTrain = terra::vect(points_EstTrain, geom = c("X", "Y"), crs = terra::crs(Sample_lyrs))

    # Calculate the Euclidean distance of my training set
    distance_ED = terra::distance(Sample_lyrs, points_EstTrain, rasterize = FALSE)
    names(distance_ED) = "Distance_ED"

    # Adequacy of covariates for training.
    day_COV$distance_ED = distance_ED
    names_covs = names(day_COV)

    # Merge the covariates with the training data
    data_cov = lapply(day_COV, function(x) terra::extract(x, points_EstTrain))
    data_cov = Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_cov)
    names(data_cov)[2:length(data_cov)] = names_covs

    # Data train RF
    dt.train = merge(data_obs[,.(ID, Date, var)], data_cov, by = "ID")
    dt.train = dt.train[, !(c("ID", "Date")), with = FALSE]

    # Generate RF model to predict values
    set.seed(seed)

    RF_model = suppressWarnings(
      randomForest::randomForest(var ~ ., data = dt.train, na.action = na.omit, ntree = ntree)
    )

    # Generate RF model to predict residuals
    pred_m1_RF = predict(RF_model, dt.train[,-("var")])
    val_m1_RF = data.table(Obs = dt.train$var, sim = pred_m1_RF)
    val_m1_RF$residuals = val_m1_RF$Obs - val_m1_RF$sim

    dt.train_resi = cbind(residuals = val_m1_RF$residuals, dt.train[, setdiff(names(dt.train), "var"), with = FALSE])

    RF_model_resdls = suppressWarnings(
      randomForest::randomForest(residuals ~ ., data = dt.train_resi, na.action = na.omit, ntree = ntree)
    )

    # Make predictions over the entire satellite dataset.
    cov_Sat = terra::rast(day_COV)

    # 1. value prediction
    sat_1 = terra::predict(cov_Sat, RF_model)

    # 2. Prediction of residuals
    sat_2 = terra::predict(cov_Sat, RF_model_resdls)

    # Generate the final model
    final_model = sat_1 + sat_2
    final_model = terra::app(final_model, fun = function(x) round(x, 1))

    if (!is.null(threshold))  final_model = terra::app(final_model, fun = function(x) ifelse(x < threshold, 0, x))

    return(final_model)
  }

  # Execute analysis all layers
  pbapply::pboptions(type = "timer", use_lb = T, style = 1, char = " ")
  message("Creating corrected models. Please wait...")
  raster_prediction = pbapply::pbsapply(Dates_extracted, function(fecha) {
    day_COV = lapply(Covariates, function(x) x[[which(Dates_extracted == fecha)]])
    prediction_lyr = RF_train(day_COV, fecha)
    names(prediction_lyr) = fecha
    return(prediction_lyr)
  }, simplify = FALSE)

  raster_final = terra::rast(raster_prediction)
  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Saving model. Please wait.")
    if (is.null(name_save)) name_save = "Model_RFplus"
    name_saving = paste0(name_save, ".nc")
    terra::writeCDF(raster_final, filename= name_saving, overwrite=TRUE)
  }

  return(raster_final)
}
