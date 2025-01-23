#' RFplus: Bias Correction of Satellite Products Using Random Forest and Quantile Regression Forests
#'
#' RFplus applies a comprehensive bias correction approach that combines three modeling steps to enhance satellite-derived environmental data. This ensemble method includes:
#' 1) Initial Random Forest predictions,
#' 2) Residual corrections using a secondary Random Forest, and
#' 3) Outlier preservation through Quantile Regression Forests (QRF).
#' The final predictions are calculated as a weighted average of all steps.
#'
#' @details The RFplus methodology operates in three stages:
#' - **Base Prediction**: A Random Forest model generates predictions based on satellite data and environmental covariates.
#' - **Residual Correction**: A second Random Forest model corrects residual errors from the base predictions.
#' - **Outlier Adjustment**: A Quantile Regression Forest computes prediction quantiles (Q5, Q8, Q9) to preserve distribution tails.
#'
#' Final predictions are calculated as:
#' \deqn{(Q5 + Q8 + Q9 + (adjusted prediction - residuals prediction) / 4}
#' This formulation balances overall accuracy while maintaining outliers.
#'
#' @param Covariates A list of `terra::SpatRaster` objects representing predictive features, including:
#' - Satellite-derived data layers
#' - Digital Elevation Model (DEM)
#' - Additional environmental covariates as needed
#' @param BD_Insitu A `data.table` containing in-situ measurements with the following columns:
#' - `Date`: Measurement dates (Date format)
#' - Observed variable columns with numerical values
#' @param Cords_Insitu A `data.table` with station coordinates, including:
#' - `Cod`: Unique station ID (matching IDs in `BD_Insitu`)
#' - `X`: Longitude values
#' - `Y`: Latitude values
#' @param ntree An integer specifying the number of trees for all forest models (RF, residual RF, and QRF). Default: 2000.
#' @param threshold A numeric value defining the lower limit for predictions. Predictions below this threshold are set to 0. Default: `NULL` (retains all values, useful for temperature predictions).
#' @param n_round An integer indicating the number of decimal places for rounding final predictions. Default: `NULL` (no rounding applied).
#' @param seed An integer for setting the random seed to ensure reproducibility. Default: 123.
#' @param save_model A logical value indicating whether to save the final model output as a NetCDF file. Default: `FALSE`.
#' @param name_save A character string specifying the base file name for saved outputs (without the `.nc` extension). Default: "Model_RFplus".
#'
#' @return A `terra::SpatRaster` containing the bias-corrected predictions. Spatial properties match those of the input covariates.
#'
#' @note This approach ensures outliers are preserved while reducing overall bias. The QRF component specifically enhances the preservation of distribution tails through its quantile-based estimates.
#'
#' @author Jonnathan Augusto Landi Bermeo
#'
#' @importFrom terra ext crs nlyr rast vect extract distance writeCDF
#' @import data.table
#' @importFrom dplyr mutate
#' @importFrom randomForest tuneRF
#' @importFrom ranger ranger
#' @importFrom pbapply pbsapply pboptions
#' @importFrom ncdf4 nc_open ncvar_get nc_close
#' @importFrom methods as
#' @importFrom stats na.omit predict
#' @importFrom utils globalVariables
#' @export

RFplus = function(Covariates, BD_Insitu, Cords_Insitu, ntree = 2000, threshold = NULL,
                  n_round = NULL, save_model = F, name_save = NULL, seed = 123) {

  ##############################################################################
  #                      Check the input data of the covariates                #
  ##############################################################################
  # Verify the class of covariates
  if (length(unique(sapply(Covariates, class))) > 1) stop("The class of the covariates are different (all classes should be similar).")

  # Verify the extent of covariates
  ext_check = sapply(Covariates, function(x) terra::ext(x))
  if (!all(sapply(ext_check, function(x) x == ext_check[[1]]))) stop("The extension of the covariates are different (all extensions should be similar).")

  # Verify the crc of covariates
  if (length(unique(sapply(Covariates, function(x) terra::crs(x)))) > 1) stop("The crs of the covariates are different (all crs should be similar).")

  ##############################################################################
  #                    Check input data from on-site stations                  #
  ##############################################################################
  # Verify the BD_Insitu is data.table
  if (!is.data.table(BD_Insitu)) stop("The data of the on-site stations should be a data.table.")

  # Verify the columns of Cords_Insitu
  if (!is.data.table(Cords_Insitu)) stop("The coordinate data of the on-site stations must be a data.table.")

  # Check that the coordinate names appear in the observed data
  if (!all(Cords_Insitu$Cod %in% setdiff(colnames(BD_Insitu), "Date"))) stop("The names of the coordinates do not appear in the observed data.")

  ##############################################################################
  #                           Operations with covariables                      #
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
  #                         Prepare data for training                          #
  ##############################################################################
  # Layer to sample
  Sample_lyrs = Covariates[[2]][[1]]

  # Data for training
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
  Points_VectTrain = terra::vect(Points_Train, geom = c("X", "Y"), crs = terra::crs(Sample_lyrs))

  # Calculate the Distance Euclidean
  distance_ED = setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    terra::distance(Sample_lyrs, Points_VectTrain[i, ], rasterize = FALSE)
  }), Points_VectTrain$Cod)

  # Calculate altitude difference
  difference_altitude = setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    Covariates$DEM[[1]] - terra::extract(Covariates$DEM[[1]], Points_VectTrain[i, ])[, 2]
  }), Points_VectTrain$Cod)

  ##############################################################################
  #                              Function of RF models                         #
  ##############################################################################
  RF_Modelplus = function(day_COV, fecha) {

    for (i in seq_along(day_COV)) {
      names(day_COV)[i] = names(day_COV[[i]])
    }

    data_obs = data_train[Date == as.Date(fecha), ] # Filtering the observation 'i'
    data_obs = data_obs[!is.na(var)] # Method for training models independently of NAs
    data_obs$ID = 1:nrow(data_obs)

    points_EstTrain = data.table::merge.data.table(data_obs[,.(ID, Cod)], Points_Train[,.(Cod, X, Y)], by = "Cod")
    points_EstTrain = terra::vect(points_EstTrain, geom = c("X", "Y"), crs = terra::crs(Sample_lyrs))

    # Covariates extras
    day_COV$dist_ED = terra::rast(distance_ED[Points_Train$Cod]) %>%
      setNames(paste("dist_ED_", seq_along(Points_Train$Cod), sep = ""))

    day_COV$diff_alt = terra::rast(difference_altitude[points_EstTrain$Cod]) %>%
      setNames(paste("diff_alt_", seq_along(points_EstTrain$Cod), sep = ""))

    # Training the model
    data_cov = lapply(day_COV, function(x) terra::extract(x, points_EstTrain))
    data_cov = Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_cov)

    dt.train = merge(data_obs[,.(ID, Date, var)], data_cov, by = "ID")
    dt.train = dt.train[, !(c("ID", "Date")), with = FALSE]

    # Search for the best mtry
    set.seed(seed)
    Hyper_P1 = suppressWarnings(
      randomForest::tuneRF(
        x = dt.train[, setdiff(names(dt.train), "var"), with = FALSE],
        y = dt.train$var,
        ntreeTry = 100,
        stepFactor = 1,
        doBest = T,
        plot = F,
        seed = seed,
        trace = 0
      )
    )

    # Model 1
    Model_P1 = ranger::ranger(
      formula = var ~ .,
      data = dt.train,
      num.trees = ntree,
      seed = seed,
      mtry = Hyper_P1$mtry
    )

    # Prediction of the Model 1
    pred_P1 = predict(Model_P1, dt.train[,-("var")])
    val_P1 = data.table(Obs = dt.train$var, sim = pred_P1$predictions)
    val_P1$residuals = val_P1$Obs - val_P1$sim
    dt.train_resi = cbind(residuals = val_P1$residuals, dt.train[, setdiff(names(dt.train), "var"), with = FALSE])

    # Search for the best mtry for the residuals
    set.seed(seed)
    Hyper_P2 = suppressWarnings(
      randomForest::tuneRF(
        x = dt.train_resi[, setdiff(names(dt.train_resi), "residuals"), with = FALSE],
        y = dt.train_resi$residuals,
        ntreeTry = 2000,
        stepFactor = 1,
        doBest = T,
        plot = F,
        seed = seed,
        trace = 0
      )
    )

    # Model 2 (Residuals)
    Model_P2 = ranger::ranger(
      formula = residuals ~ .,
      data = dt.train_resi,
      num.trees = ntree,
      seed = seed,
      mtry = Hyper_P2$mtry
    )

    # Predictions using the two models
    cov_Sat = terra::rast(day_COV)
    data_SatTrain = terra::as.data.frame(cov_Sat, xy = TRUE, cells = FALSE, na.rm = FALSE)
    data_SatTrain = as.data.table(data_SatTrain)
    cords = data_SatTrain[, .(x, y)]
    data_SatTrain = data_SatTrain[, setdiff(names(data_SatTrain), c("x", "y", "ID")), with = FALSE]

    pred_1 = predict(Model_P1, data_SatTrain)
    pred_2 = predict(Model_P2, data_SatTrain)
    Phase_1_2 = pred_1$predictions + pred_2$predictions

    if (!is.null(n_round)) Phase_1_2 = round(Phase_1_2, n_round)
    if (!is.null(threshold))  Phase_1_2 = ifelse(Phase_1_2 < threshold, 0, Phase_1_2)

    Phase_1_2 = cbind(cords[, c("x", "y")], Sim = Phase_1_2)
    ############################################################################
    #                         Phase 3 (Model of QRF)                           #
    ############################################################################
    # Model QRF
    set.seed(seed)
    model_QRF = ranger::ranger(
      formula = var ~ .,
      data = dt.train,
      num.trees = ntree,
      seed = seed,
      quantreg = TRUE
    )

    pred_completas = predict(model_QRF, data_SatTrain, type = "quantiles", quantiles = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9), na.action = na.omit, seed = seed)
    pred_completas = cbind(cords[, c("x", "y")], pred_completas$predictions)
    names(pred_completas) = c("x", "y","Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9")

    # Merge Phase 1 and 2 with QRF
    final_model = merge(Phase_1_2, pred_completas, by = c("x", "y"))
    final_model$mean = rowMeans(final_model[, c("Q5", "Q8","Q9", "Sim")], na.rm = TRUE)
    final_model = final_model[,.(x, y, mean)]

    # Ensamble of the models
    Ensamble = terra::rast(data.frame(
      x = final_model$x,
      y = final_model$y,
      value = final_model$mean
    ), crs = terra::crs(Sample_lyrs))

    return(Ensamble)
  }

  # Execute analysis
  pbapply::pboptions(type = "timer", use_lb = T, style = 1, char = " ")
  message("Analysis in progress. Please wait....")

  raster_Model = pbapply::pbsapply(Dates_extracted, function(fecha) {
    day_COV = lapply(Covariates, function(x) x[[which(Dates_extracted == fecha)]])
    prediction_lyr = RF_Modelplus(day_COV, fecha)
    return(prediction_lyr)
  }, simplify = FALSE)

  Ensamble = terra::rast(raster_Model)

  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Saving model. Please wait.")
    if (is.null(name_save)) name_save = "Model_RFplus"
    name_saving = paste0(name_save, ".nc")
    terra::writeCDF(Ensamble, filename = name_saving, overwrite=TRUE)
  }
  return(Ensamble)
}
