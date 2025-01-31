# Declare global variables to prevent R CMD check warnings
utils::globalVariables(c("Cod", "ID", "X", "Y", "Date", "Obs", "sim", "residuals", "var"))
#' Bias Correction of Satellite Products Using Hybrid Random Forest and Quantile Mapping
#'
#' @description
#' Applies a three-step bias correction combining Random Forest predictions, residual correction,
#' and distribution adjustment using quantile mapping methods.
#'
#' @details
#' `RFplus` implements a hybrid approach:
#' 1. **Base Prediction**: Initial Random Forest model using satellite data and covariates.
#' 2. **Residual Correction**: Secondary Random Forest to correct residuals from step 1.
#' 3. **Distribution Adjustment**: Quantile mapping (QUANT/RQUANT) to match observed distribution.
#'
#' The final prediction combines all three steps, preserving outliers while reducing bias.
#'
#' @param BD_Insitu A `data.table` with in-situ measurements containing:
#'  - `Date`: Date of measurements (Date format)
#'  - Columns with observed values (numeric)
#' @param Cords_Insitu A `data.table` with station coordinates:
#'  - `Cod`: Station ID (matching `BD_Insitu`)
#'  - `X`: Longitude
#'  - `Y`: Latitude
#' @param Covariates List of `terra::SpatRaster` objects containing:
#'  - Satellite data layers
#'  - Digital Elevation Model (DEM)
#'  - Other environmental covariates
#' @param ntree Number of trees for Random Forest (default: 2000)
#' @param method Distribution adjustment method: "RQUANT", "QUANT", or "none" (default: "none")
#' @param ratio Maximum search radius (km) for nearest station adjustment (default: 15)
#' @param wet.day Logical threshold for wet-day correction or numeric value (see Note)
#' @param seed Random seed for reproducibility (default: 123)
#' @param n_round Number of decimal places for predictions (default: NULL)
#' @param save_model Save output as NetCDF file? (default: FALSE)
#' @param name_save Base name for output file (default: "Model_RFplus")
#' @param ... Additional arguments passed to methods
#'
#' @note The `wet.day` parameter can be:
#' - `FALSE`: No wet-day correction
#' - `TRUE`: Apply correction with default threshold (0.1 mm)
#' - Numeric value: Custom threshold
#'
#' @rdname RFplus
#' @export
RFplus <- function(BD_Insitu, ...) {
  UseMethod("RFplus")
}

#' @rdname RFplus
#' @export
RFplus.default <- function(BD_Insitu, Cords_Insitu, Covariates, n_round = NULL, wet.day = FALSE,
                           ntree = 2000, seed = 123, method=c("RQUANT","QUANT","none"), ratio = 15,
                           save_model = FALSE, name_save = NULL, ...) {

  if (!inherits(BD_Insitu, "data.table")) stop("BD_Insitu must be a data.table.")

  RFplus.data.table(BD_Insitu = BD_Insitu, Cords_Insitu = Cords_Insitu, Covariates = Covariates,
                    n_round = n_round, wet.day = wet.day, ntree = ntree, seed = seed, method = method,
                    ratio = ratio, save_model = save_model, name_save = name_save)
}

#' @rdname RFplus
#' @export
RFplus.data.table <- function(BD_Insitu, Cords_Insitu, Covariates, n_round = NULL, wet.day = FALSE,
                              ntree = 2000, seed = 123, method=c("RQUANT","QUANT","none"), ratio = 15,
                              save_model = FALSE, name_save = NULL) {

  ##############################################################################
  #                      Check the input data of the covariates                #
  ##############################################################################
  # Verify that covariates is a list
  if (!inherits(Covariates, "list")) stop("Covariates must be a list.")

  # Verify that the covariates are type SpatRaster
  if (!all(sapply(Covariates, function(x) inherits(x, "SpatRaster")))) stop("The covariates must be of type SpatRaster.")

  # Verify the extent of covariates
  ext_list = lapply(Covariates, terra::ext)
  if (!all(sapply(ext_list, function(x) x == ext_list[[1]]))) stop("The extension of the covariates are different (all extensions should be similar).")

  # Verify the crc of covariates
  if (length(unique(vapply(Covariates, terra::crs, character(1)))) > 1) stop("The crs of the covariates are different (all crs should be similar).")

  ##############################################################################
  #                    Check input data from on-site stations                  #
  ##############################################################################
  # Verify the BD_Insitu is data.table
  if (!inherits(BD_Insitu, "data.table")) stop("The data of the on-site stations should be a data.table.")

  # Verify the columns of Cords_Insitu
  if (!inherits(Cords_Insitu, "data.table")) stop("The coordinate data of the on-site stations must be a data.table.")

  # Check that the coordinate names appear in the observed data
  if (!all(Cords_Insitu$Cod %chin% dplyr::setdiff(names(BD_Insitu), "Date"))) stop("The names of the coordinates do not appear in the observed data.")

  ##############################################################################
  #                Checking the input parameters for quantile mapping          #
  ##############################################################################
  # Verify the method
  method = match.arg(method, choices = c("RQUANT", "QUANT", "none"))

  ##############################################################################
  #               Verify that there is a DEM and manage DEM layers.            #
  ##############################################################################
  # Check if there is a DEM layer
  nlyr_covs = sapply(Covariates, function(x) terra::nlyr(x))
  index_dem = which(nlyr_covs == 1)
  if (length(index_dem) == 0) stop("A single layer covariate was not found. Possibly the DEM was not entered.")

  # Replicating the DEM at covariate scale
  nlyrs_tots = which(nlyr_covs != 1)
  nlyr_rep = nlyr_covs[nlyrs_tots[1]]
  DEM = Covariates[[index_dem]]
  DEM = terra::rast(replicate(nlyr_rep, DEM))
  Covariates[[index_dem]] = DEM

  # Verify the layers of the covariates.
  if (length(unique(sapply(Covariates, function(x) terra::nlyr(x)))) > 1) stop("The number of covariate layers does not match. Check the input data.")

  ##############################################################################
  #                         Prepare data for training                          #
  ##############################################################################
  # Layer to sample
  Sample_lyrs = DEM[[1]]

  # Data for training
  data_train = data.table::melt(
    BD_Insitu,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  ) %>%
    dplyr::mutate(ID = as.numeric(factor(Cod)))

  # Date of the data
  Dates_extracted = base::unique(data_train$Date)
  Points_Train = merge(data_train, Cords_Insitu, by = "Cod")
  setDT(Points_Train)

  Points_Train = base::unique(Points_Train[,.(ID, Cod, X, Y)])
  setorder(Points_Train, ID)
  Points_VectTrain = terra::vect(Points_Train, geom = c("X", "Y"), crs = crs(Sample_lyrs))

  # Calculate the Distance Euclidean
  distance_ED = setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    terra::distance(Sample_lyrs, Points_VectTrain[i, ], rasterize = FALSE)
  }), Points_VectTrain$Cod)

  # Calculate altitude difference
  difference_altitude = setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    Covariates$DEM[[1]] - terra::extract(Covariates$DEM[[1]], Points_VectTrain[i, ])[, 2]
  }), Points_VectTrain$Cod)

  ##############################################################################
  ##############################################################################
  #                    Progressive correction methodology                      #
  ##############################################################################
  ##############################################################################

  # Model of the Random Forest for the progressive correction 1 y 2 ------------
  RF_Modelplus = function(day_COV, fecha) {
    for (i in seq_along(day_COV)) {
      names(day_COV)[i] = names(day_COV[[i]])
    }

    data_obs = data_train[Date == as.Date(fecha), ]
    points_EstTrain = merge(data_obs[,.(ID, Cod)], Points_Train[,.(Cod, X, Y)], by = "Cod")
    setorder(points_EstTrain, ID)
    points_EstTrain = terra::vect(points_EstTrain, geom = c("X", "Y"), crs = terra::crs(Sample_lyrs))

    # Covariates extras
    day_COV$dist_ED = terra::rast(distance_ED[Points_Train$Cod]) %>%
      setNames(paste("dist_ED_", seq_along(Points_Train$Cod), sep = ""))

    day_COV$diff_alt = terra::rast(difference_altitude[points_EstTrain$Cod]) %>%
      setNames(paste("diff_alt_", seq_along(points_EstTrain$Cod), sep = ""))

    #                   Training the model Random Forest (1)                   #
    data_cov = lapply(day_COV, function(x) terra::extract(x, points_EstTrain))
    data_cov = Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_cov)

    dt.train = merge(data_obs[,.(ID, Date, var)], data_cov, by = "ID")
    dt.train = dt.train[, setdiff(names(dt.train), "Date"), with = FALSE]

    # Train Model 1
    set.seed(seed)
    Model_P1 = suppressWarnings(
      randomForest::randomForest(
        formula = var ~ .,
        data = dt.train[, setdiff(names(dt.train), "ID"), with = FALSE],
        ntree = ntree,
        na.action = stats::na.omit
      )
    )

    #                   Training the model Random Forest (2)                   #
    # Prediction of the Model 1
    pred_P1 = predict(Model_P1, dt.train[, setdiff(names(dt.train), c("ID", "var")), with = FALSE], na.rm = TRUE)
    # val_P1 = cbind(ID = dt.train$ID, Obs = dt.train$var, sim = pred_P1)
    val_P1 = data.table(ID = dt.train$ID,
                        Obs = dt.train$var, sim = pred_P1)

    val_P1[, residuals := Obs - sim]
    dt.train_resi = cbind(residuals = val_P1$residuals, dt.train[, setdiff(names(dt.train), c("ID", "var")), with = FALSE])

    # Train Model 2
    Model_P2 = suppressWarnings(randomForest::randomForest(
      formula = residuals ~ .,
      data = dt.train_resi,
      ntree = ntree,
      na.action = stats::na.omit
    )
    )

    # Create the corrected model 1
    cov_Sat = terra::rast(day_COV)
    pred_1 = predict(object = cov_Sat, model = Model_P1, na.rm = TRUE, fun = predict)
    pred_2 = predict(object = cov_Sat, model = Model_P2, na.rm = TRUE, fun = predict)
    Ensamble = pred_1 + pred_2

    # Extra operations in case they have been established
    if (!is.null(n_round)) Ensamble = terra::app(Ensamble, function(x) round(x, n_round))
    if (wet.day) Ensamble = terra::app(Ensamble, function(x) ifelse(x < wet.day, 0, x))
    return(Ensamble)
  }

  pbapply::pboptions(type = "timer", use_lb = T, style = 1, char = " ")
  message("Analysis in progress: Stage 1 of 2. Please wait...")

  raster_Model = pbapply::pbsapply(Dates_extracted, function(fecha) {
    day_COV = lapply(Covariates, function(x) x[[which(Dates_extracted == fecha)]])
    prediction_lyr = RF_Modelplus(day_COV, fecha)
    return(prediction_lyr)
  }, simplify = FALSE)

  Ensamble = terra::rast(raster_Model)

  # Model of the QM or QDM correction ------------------------------------------
  if (method == "none") {
    message("Analysis completed, QUANT or RQUANT correction phase not applied.")
  } else if (method %in% c("RQUANT","QUANT")) {
    message(paste0("Analysis in progress: Stage 2 of 2. Correction by: ", method, ". Please wait..."))

    data_CM = terra::extract(Ensamble, Points_VectTrain)
    names_train = data_train %>%
      dplyr::select(ID, Cod) %>%
      dplyr::distinct() %>%
      data.table::data.table()

    data_CM$ID = names_train$Cod[match(data_CM$ID, names_train$ID)]
    data_CM = na.omit(data_CM)
    names = as.character(data_CM$ID)

    data_CM = data.table(t(data_CM[, -1]))
    colnames(data_CM) = names
    data_CM = cbind(Date = BD_Insitu[, Date], data_CM)

    common_columns = dplyr::intersect(colnames(data_CM), colnames(BD_Insitu))
    common_columns = dplyr::setdiff(common_columns, "Date")

    res_interpolation = lapply(common_columns, function(col) {
      data_obs = BD_Insitu[, .(Date, Obs = get(col))]
      data_sim = data_CM[, .(Date, Sim = get(col))]
      merge(data_obs, data_sim, by = "Date")
    })

    names(res_interpolation) = common_columns
    res_interpolation = lapply(res_interpolation, function(x) na.omit(x))

    data_complete = data.table(terra::as.data.frame(Ensamble, xy = TRUE))
    colnames(data_complete) = c("x", "y", as.character(Dates_extracted))

    points =Cords_Insitu[Cod %chin% names, ]
    points = terra::vect(points, geom = c("X", "Y"), crs = crs(Sample_lyrs))
    dat_final = data.table()

    if (method == "QUANT") {
      cuantiles = lapply(res_interpolation, function(x) fitQmapQUANT(x$Obs, x$Sim, method = method, wet.day = wet.day))
      for (i in 1:nrow(data_complete)) {
        x = data_complete[i, x]
        y = data_complete[i, y]

        distances = terra::distance(terra::vect(data.table(cbind(x = x,y = y)), geom = c("x", "y"), crs = crs(Sample_lyrs)),
                                    points,
                                    unit = "km")

        distances = data.table(t(distances))
        distances = cbind(distances, Cod = points$Cod)
        names(distances) = c("dist", "Cod")

        if (any(distances$dist <= ratio)){
          inds = base::which.min(distances$dist)
          name = distances[inds, .(Cod)]

          data = data.table(Sim = t(data_complete[i, -c("x", "y")]))
          data_corregido = doQmapQUANT(data$Sim.V1, cuantiles[[name$Cod]])
          data_sat = cbind(data_complete[i, c("x", "y")], t(data_corregido))
          colnames(data_sat) = c("x", "y", as.character(Dates_extracted))
          dat_final = rbind(dat_final, data_sat)
        } else {
          dat_final = rbind(dat_final, data_complete[i, ])
        }
      }
    } else {
      cuantiles = lapply(res_interpolation, function(x) fitQmapRQUANT(x$Obs, x$Sim, method = method))
      for (i in 1:nrow(data_complete)) {
        x = data_complete[i, x]
        y = data_complete[i, y]

        distances = terra::distance(terra::vect(data.table(cbind(x = x,y = y)), geom = c("x", "y"), crs = crs(Sample_lyrs)),
                                    points,
                                    unit = "km")

        distances = data.table(t(distances))
        distances = cbind(distances, Cod = points$Cod)
        names(distances) = c("dist", "Cod")

        if (any(distances$dist <= ratio)){
          inds = base::which.min(distances$dist)
          name = distances[inds, .(Cod)]

          data = data.table(Sim = t(data_complete[i, -c("x", "y")]))
          data_corregido = doQmapRQUANT(data$Sim.V1, cuantiles[[name$Cod]])
          data_sat = cbind(data_complete[i, c("x", "y")], t(data_corregido))
          colnames(data_sat) = c("x", "y", as.character(Dates_extracted))
          dat_final = rbind(dat_final, data_sat)
        } else {
          dat_final = rbind(dat_final, data_complete[i, ])
        }
      }
    }
    Ensamble = terra::rast(dat_final, crs = crs(Sample_lyrs))

  }

  message("Analysis completed.")
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
