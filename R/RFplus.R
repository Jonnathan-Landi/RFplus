# Declare global variables to prevent R CMD check warnings
utils::globalVariables(c("Cod", "ID", "X", "Y", "Date", "Obs", "sim", "residuals", "var", ".", ".SD"))
#' Machine learning algorithm for fusing ground and satellite precipitation data.
#'
#' @description
#' MS-GOP is a machine learning algorithm for merging satellite-based and ground precipitation data.
#' It combines Random Forest for spatial prediction, residual modeling for bias correction, and quantile mapping for final adjustment, ensuring accurate precipitation estimates across different temporal scales
#'
#' @details
#' The `RFplus` method implements a three-step approach:
#' 1. **Base Prediction**: Random Forest model is trained using satellite data and covariates.
#' 2. **Residual Correction**: A second Random Forest model is used to correct the residuals from the base prediction.
#' 3. **Distribution Adjustment**: Quantile mapping (QUANT or RQUANT) is applied to adjust the distribution of satellite data to match the observed data distribution.
#'
#' The final result combines all three steps, correcting the biases while preserving the outliers, and improving the accuracy of satellite-derived data such as precipitation and temperature.
#'
#' @param BD_Insitu `data.table` containing the ground truth measurements (dependent variable) used to train the RFplus model.
#'   Each column represents a ground station, and station identifiers are stored as column names. The class of `BD_Insitu`
#'   must be `data.table`. Each row represents a time step with measurements of the corresponding station.
#' @param Cords_Insitu `data.table` containing metadata for the ground stations. Must include the following columns:
#'   - `Cod`: Unique identifier for each ground station.
#'   - `X`: Latitude of the station in UTM format.
#'   - `Y`: Longitude of the station in UTM format.
#' @param Covariates A list of covariates used as independent variables in the RFplus model. Each covariate should be a
#'   `SpatRaster` object (from the `terra` package) and can represent satellite-derived weather variables or a Digital
#'    Elevation Model (DEM). All covariates should have the same number of layers (bands), except for the DEM, which must have only one layer.
#' @param n_round Numeric indicating the number of decimal places to round the corrected values. If `n_round` is set to `NULL`, no rounding is applied.
#' @param wet.day Numeric value indicating the threshold for wet day correction. Values below this threshold will be set to zero.
#'   - `wet.day = FALSE`: No correction is applied (default).
#'   - For wet day correction, provide a numeric threshold (e.g., `wet.day = 0.1`).
#' @param ntree Numeric indicating the maximum number trees to grow in the Random Forest algorithm. The default value is set to 2000. This should not be set to too small a number, to ensure that every input row gets predicted at least a few times. If this value is too low, the prediction may be biased.
#' @param seed Integer for setting the random seed to ensure reproducibility of results (default: 123).
#' @param training Numerical value between 0 and 1 indicating the proportion of data used for model training. The remaining data are used for validation. Note that if you enter, for example, 0.8 it means that 80 % of the data will be used for training and 20 % for validation.
#' If you do not want to perform validation, set training = 1. (Default training = 1).
#' @param Rain_threshold Numeric value that defines the precipitation threshold for classifying rainfall events.
#' Precipitation values above this threshold will be considered as rainfall events, while values below it will be considered as no-rain events.
#' This parameter is used to calculate key performance metrics such as the Probability of Detection (POD), False Alarm Rate (FAR), and Critical Success Index (CSI),
#' which help assess the accuracy of rainfall event predictions.
#' Note: This parameter should only be set when `training` is not equal to 1, as it is needed to calculate the POD, FAR, and CSI metrics.
#' The default value for this parameter is 0.1.
#' @param method A character string specifying the quantile mapping method used for distribution adjustment. Options are:
#'   - `"RQUANT"`: Robust quantile mapping to adjust satellite data distribution to observed data.
#'   - `"QUANT"`: Standard quantile mapping.
#'   - `"none"`: No distribution adjustment is applied. Only Random Forest-based bias correction and residual correction are performed.
#' @param ratio integer Maximum search radius (in kilometers) for the quantile mapping setting using the nearest station. (default = 15 km)
#' @param save_model Logical value indicating whether the corrected raster layers should be saved to disk. The default is `FALSE`.
#'    If set to `TRUE`, make sure to set the working directory beforehand using `setwd(path)` to specify where the files should be saved.
#' @param name_save Character string. Base name for output file (default: NULL). The output file will be saved as "Model_RFplus.nc".
#' If you set a different name, make sure you do not set the ".nc" format,
#' as the code will internally assign it.
#' @param ... Additional arguments to pass to the underlying methods (e.g., for model tuning or future extensions).
#' @examples
#' \donttest{
#' # Load the libraries
#' library(terra)
#' library(data.table)
#'
#' # Load the data
#'  data("BD_Insitu", package = "RFplus")
#'  data("Cords_Insitu", package = "RFplus")
#'
#' # Load the covariates
#' Covariates <- list(
#'  MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "RFplus")),
#'  CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "RFplus")),
#'  DEM = terra::rast(system.file("extdata/DEM.nc", package = "RFplus"))
#'  )
#'
#'  # Apply the RFplus bias correction model
#' model = RFplus(BD_Insitu, Cords_Insitu, Covariates, n_round = 1, wet.day = 0.1,
#'         ntree = 2000, seed = 123, training = 1, Rain_threshold = 0.1, method = "RQUANT",
#'         ratio = 5, save_model = FALSE, name_save = NULL)
#' # Visualize the results
#' # Precipitation results within the study area
#' modelo_rainfall = model$Ensamble
#' # Validation statistic results
#' metrics = model$Validation
#' # Note: In the above example we used 80% of the data for training and 20% for # model validation.
#' }
#'
#' @return Returns a list containing two elements:
#'   \item{Ensamble}{A `SpatRaster` object containing the bias-corrected layers for each time step. The number of layers
#'   corresponds to the number of dates for which the correction is applied. This represents the corrected satellite data adjusted for bias.}
#'   \item{Validation}{A data frame or similar object containing the statistical results obtained from the validation process.
#'   These statistics assess the performance of the bias correction applied to the satellite data.}.
#' @author
#'  Jonnathan Augusto landi Bermeo, jonnathan.landi@outlook.com
#' @rdname RFplus
#' @export
RFplus <- function(BD_Insitu, Cords_Insitu, Covariates,...) {
  UseMethod("RFplus")
}

#' @rdname RFplus
#' @export
RFplus.default <- function(BD_Insitu, Cords_Insitu, Covariates, n_round = NULL, wet.day = FALSE,
                           ntree = 2000, seed = 123, training = 1, Rain_threshold = 0.1, method = c("RQUANT","QUANT","none"),
                           ratio = 15, save_model = FALSE, name_save = NULL, ...) {

  if (!inherits(BD_Insitu, "data.table")) stop("BD_Insitu must be a data.table.")

  RFplus.data.table(BD_Insitu = BD_Insitu, Cords_Insitu = Cords_Insitu, Covariates = Covariates,
                    n_round = n_round, wet.day = wet.day, ntree = ntree, seed = seed, training = training,
                    Rain_threshold = Rain_threshold,method = method, ratio = ratio, save_model = save_model,
                    name_save = name_save, ...)
}

#' @rdname RFplus
#' @export
RFplus.data.table <- function(BD_Insitu, Cords_Insitu, Covariates, n_round = NULL, wet.day = FALSE,
                              ntree = 2000, seed = 123, training = 1, Rain_threshold = 0.1, method = c("RQUANT","QUANT","none"),
                              ratio = 15, save_model = FALSE, name_save = NULL, ...) {

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

  # Check if there is a Date column
  date_column <- names(BD_Insitu)[which(sapply(BD_Insitu, function(x) inherits(x, c("Date", "IDate", "POSIXct"))))]
  if (length(date_column) == 0) stop("The Date column was not found in the observed data.")

  # Change the column name to match the full code
  if (date_column != "Date") setnames(BD_Insitu, date_column, "Date")

  # Check that the coordinate names appear in the observed data
  if (!all(Cords_Insitu$Cod %chin% setdiff(names(BD_Insitu), "Date"))) stop("The names of the coordinates do not appear in the observed data.")

  ##############################################################################
  #                Checking the input parameters for quantile mapping          #
  ##############################################################################
  # Verify the method
  method <- match.arg(method, choices = c("RQUANT", "QUANT", "none"))

  ##############################################################################
  #               Verify that there is a DEM and manage DEM layers.            #
  ##############################################################################
  # Check if there is a DEM layer
  nlyr_covs <- sapply(Covariates, function(x) terra::nlyr(x))
  index_dem <- which(nlyr_covs == 1)
  if (length(index_dem) == 0) stop("A single layer covariate was not found. Possibly the DEM was not entered.")

  # Replicating the DEM at covariate scale
  nlyrs_tots <- which(nlyr_covs != 1)
  nlyr_rep <- nlyr_covs[nlyrs_tots[1]]
  DEM <- Covariates[[index_dem]]
  DEM <- terra::rast(replicate(nlyr_rep, DEM))
  Covariates[[index_dem]] <- DEM

  # Verify the layers of the covariates.
  if (length(unique(sapply(Covariates, function(x) terra::nlyr(x)))) > 1) stop("The number of covariate layers does not match. Check the input data.")

  ##############################################################################
  #                          Verify if validation is to be done                #
  ##############################################################################
  if (training != 1) {
    message(paste("The training parameter has been introduced. The model will be trained with:", (training * 100), "%", "data and validated with:", (100 - (training * 100)), "%"))

    # Verify if the number entered in training is valid
    if (!(training %between% c(0, 1))) stop("The training parameter must be between 0 and 1.")

    # Exclude Date column and split remaining columns
    set.seed(seed)
    columns <- setdiff(names(BD_Insitu), "Date")

    #  Randomly select the % of the training columns
    train_columns <- sample(columns, size = floor(training * length(columns)))

    # Data train
    train_data <- BD_Insitu[, .SD, .SDcols = c("Date", train_columns)]

    # Data test
    test_data <- BD_Insitu[, .SD, .SDcols = setdiff(names(BD_Insitu), train_columns)]

    # Function used to validate data
    validate = function(data, Rain_threshold) {
      CC = round(hydroGOF::rPearson(data$Sim, data$Obs, na.rm = T),3)
      RMSE = round(hydroGOF::rmse(data$Sim, data$Obs, na.rm = T),3)
      KGE = round(hydroGOF::KGE(data$Sim, data$Obs, na.rm = T),3)
      PBIAS = round(hydroGOF::pbias(data$Sim, data$Obs, na.rm = T),3)

      # Calculations for other metrics
      data = na.omit(data)
      observed_fac = factor(ifelse(data$Obs > Rain_threshold, 1, 0), levels = c(1, 0))
      forecast_fac = factor(ifelse(data$Sim > Rain_threshold, 1, 0), levels = c(1, 0))
      contingency_table <- table(forecast_fac, observed_fac)

      # Extract components from the table
      H <- contingency_table["1", "1"]  # Hits (Forecast=1, Observed=1)
      FA <- contingency_table["1", "0"] # False alarms (Forecast=1, Observed=0)
      M <- contingency_table["0", "1"]  # Omissions (Forecast=0, Observed=1)
      CN <- contingency_table["0", "0"] # Negative corrections (Forecast=0, Observed=0)

      # functions for calculations
      calculate_pod <- function(H, M) {
        if (H + M == 0) return(NA)
        return(H / (H + M))
      }

      # Function for calculating False Alarm Rate (FAR)
      calculate_far <- function(H, FA) {
        if (H + FA == 0) return(NA)
        return(FA / (H + FA))
      }

      # Function for calculating the Critical Success Index (CSI)
      calculate_csi <- function(H, M, FA) {
        if (H + M + FA == 0) return(NA)
        return(H / (H + M + FA))
      }

      # Metrics
      pod <- calculate_pod(H, M)
      far <- calculate_far(H, FA)
      csi <- calculate_csi(H, M, FA)

      res_metrics = data.table(
        CC = CC,
        RMSE = RMSE,
        KGE = KGE,
        PBIAS = PBIAS,
        POD = pod,
        FAR = far,
        CSI = csi
      )
      return(res_metrics)
    }

    # Coordinates of the training data
    train_cords <- Cords_Insitu[Cod %chin% train_columns, ]

    # Coordinates of the test data
    test_cords <- Cords_Insitu[Cod %chin% setdiff(names(BD_Insitu), train_columns), ]

  } else {
    message("The training parameter was not entered. The model will be trained with all the data.")
    train_data <- BD_Insitu
    train_cords <- Cords_Insitu
  }

  ##############################################################################
  #                         Prepare data for training                          #
  ##############################################################################
  # Layer to sample
  Sample_lyrs <- DEM[[1]]

  # Data for training
  training_data <- melt(
    train_data,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  )[, ID := as.numeric(factor(Cod))]

  # Date of the data
  Dates_extracted <- unique(training_data[, Date])
  Points_Train <- merge(training_data, train_cords, by = "Cod")
  setDT(Points_Train)

  Points_Train <- unique(Points_Train, by = "Cod")[, .(ID, Cod, X, Y)]
  setorder(Points_Train, ID)

  Points_VectTrain <- terra::vect(Points_Train, geom = c("X", "Y"), crs = crs(Sample_lyrs))

  # Calculate the Distance Euclidean
  distance_ED <- setNames(lapply(1:nrow(Points_VectTrain), function(i) {
    terra::distance(Sample_lyrs, Points_VectTrain[i, ], rasterize = FALSE)
  }), Points_VectTrain$Cod)

  # Calculate altitude difference
  difference_altitude <- setNames(lapply(1:nrow(Points_VectTrain), function(i) {
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
      names(day_COV)[i] <- names(day_COV[[i]])
    }

    data_obs <- training_data[Date == as.Date(fecha), ]

    points_EstTrain = merge(
      data_obs[, .(ID, Cod)],
      Points_Train[, .(Cod, X, Y)],
      by = "Cod"
    )[order(ID)]

    points_EstTrain <- terra::vect(points_EstTrain, geom = c("X", "Y"), crs = crs(Sample_lyrs))

    # Covariates extras
    day_COV$dist_ED <- terra::rast(distance_ED[Points_Train$Cod])
    names(day_COV$dist_ED) = paste("dist_ED_", seq_along(Points_Train$Cod), sep = "")

    day_COV$diff_alt <- terra::rast(difference_altitude[Points_Train$Cod])
    names(day_COV$diff_alt) = paste("diff_alt_", seq_along(Points_Train$Cod), sep = "")

    #                   Training the model Random Forest (1)                   #
    data_cov <- Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE),
                       lapply(day_COV, terra::extract, y = points_EstTrain))

    dt.train <- merge(
      data_obs[, .(ID, Date, var)],
      data_cov,
      by = "ID"
    )

    dt.train <- dt.train[, setdiff(names(dt.train), "Date"), with = FALSE]

    # Train Model 1
    set.seed(seed)
    Model_P1 <- suppressWarnings(
      randomForest::randomForest(
        formula = var ~ .,
        data = dt.train[, setdiff(names(dt.train), "ID"), with = FALSE],
        ntree = ntree,
        na.action = stats::na.omit
      )
    )

    #                   Training the model Random Forest (2)                   #
    # Prediction of the Model 1
    val_P1 <- dt.train[, .(ID, Obs = var, sim = predict(Model_P1, .SD, na.rm = TRUE)), .SDcols = !"ID"]
    val_P1[, residuals := Obs - sim]

    dt.train_resi <- data.table(
      residuals = val_P1$residuals,
      dt.train[, setdiff(names(dt.train), c("ID", "var")), with = FALSE]
    )

    # Train Model 2
    Model_P2 <- suppressWarnings(randomForest::randomForest(
      formula = residuals ~ .,
      data = dt.train_resi,
      ntree = ntree,
      na.action = stats::na.omit
    )
    )
    # Create the corrected model 1
    cov_Sat <- terra::rast(day_COV)
    Ensamble <- predict(cov_Sat, Model_P1, na.rm = TRUE, fun = predict) +
      predict(cov_Sat, Model_P2, na.rm = TRUE, fun = predict)

    # Extra operations in case they have been established
    if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, function(x) round(x, n_round))
    if (wet.day != FALSE) Ensamble <- terra::app(Ensamble, function(x) ifelse(x < wet.day, 0, x))
    return(Ensamble)
  }

  pbapply::pboptions(type = "timer", use_lb = T, style = 1, char = "=")
  message("Analysis in progress: Stage 1 of 2. Please wait...")

  raster_Model <- pbapply::pblapply(Dates_extracted, function(fecha) {
    day_COV <- lapply(Covariates, function(x) x[[match(fecha, Dates_extracted)]])
    RF_Modelplus(day_COV, fecha)
  })

  Ensamble <- terra::rast(raster_Model)

  # Model of the QM or QDM correction ------------------------------------------
  if (method == "none") {
    message("Analysis completed, QUANT or RQUANT correction phase not applied.")
  } else if (method %in% c("RQUANT","QUANT")) {
    message(paste0("Analysis in progress: Stage 2 of 2. Correction by: ", method, ". Please wait..."))

    data_CM <- data.table(terra::extract(Ensamble, Points_VectTrain))
    data_CM[, ID := as.numeric(as.character(ID))]

    names_train = unique(training_data[, .(ID, Cod)])
    setkey(names_train, ID)

    data_CM[, ID := names_train[data_CM, on = "ID", Cod]]
    data_CM <- na.omit(data_CM)

    names <- as.character(data_CM$ID)

    data_CM <- data.table(t(data_CM[, -1]))
    colnames(data_CM) <- names
    data_CM <- data.table(training_data[,.(Date)], data_CM)

    common_columns <- setdiff(names(data_CM), c("ID", "Date"))

    res_interpolation <- lapply(common_columns, function(col) {
      dt <- merge(
        train_data[, .(Date, Obs = get(col))],
        data_CM[, .(Date, Sim = get(col))],
        by = "Date",
        all = FALSE
      )

      dt[, c("Obs", "Sim") := lapply(.SD, as.numeric), .SDcols = c("Obs", "Sim")]

      return(na.omit(dt))
    })

    names(res_interpolation) <- common_columns
    data_complete <- data.table(terra::as.data.frame(Ensamble, xy = TRUE))
    setnames(data_complete, new = c("x", "y", as.character(Dates_extracted)))

    points <- train_cords[Cod %chin% names, ]
    points <- terra::vect(points, geom = c("X", "Y"), crs = crs(Sample_lyrs))
    dat_final <- data.table()

    process_data <- function(method_fun, doQmap_fun) {
      cuantiles <- lapply(res_interpolation, function(x) method_fun(x$Obs, x$Sim, method = method, wet.day = wet.day))
      message("Applying correction method. This may take a while...")

      dat_final <- pbapply::pblapply(seq_len(nrow(data_complete)), function(i) {
        x <- data_complete[i, x]
        y <- data_complete[i, y]

        distances <- terra::distance(
          terra::vect(data.table(x, y), geom = c("x", "y"), crs = crs(Sample_lyrs)),
          points, unit = "km"
        )

        distances <- data.table(dist = as.vector(distances), Cod = points$Cod)

        if (any(distances$dist <= ratio)) {
          name <- distances[which.min(distances$dist), Cod]
          data <- data.table(Sim = t(data_complete[i, -c("x", "y")]))
          data_corregido <- doQmap_fun(data$Sim.V1, cuantiles[[name]])
          data_sat <- cbind(data_complete[i, c("x", "y")], t(data_corregido))
          colnames(data_sat) <- c("x", "y", as.character(Dates_extracted))

          return(data_sat)
        } else {
          return(data_complete[i])
        }
      })

      rbindlist(dat_final)
    }

    dat_final <- if (method == "QUANT") {
      process_data(fitQmapQUANT, doQmapQUANT)
    } else {
      process_data(method_fun = fitQmapRQUANT, doQmap_fun = doQmapRQUANT)
    }

    Ensamble <- terra::rast(dat_final, crs = crs(Sample_lyrs))
    if (!is.null(n_round)) Ensamble <- terra::app(Ensamble, \(x) round(x, n_round))
    message("Analysis completed.")
  }

  ##############################################################################
  #                           Perform validation if established                #
  ##############################################################################
  if (training != 1) {
    message("Validation process in progress. Please wait.")
    test_cords$ID <- seq_len(nrow(test_cords))
    Points_VectTest <- terra::vect(test_cords, geom = c("X", "Y"), crs = crs(Sample_lyrs))

    data_validation = data.table(terra::extract(Ensamble, Points_VectTest))
    setkey(data_validation, ID)

    testing_data <- melt(
      test_data,
      id.vars = "Date",
      variable.name = "Cod",
      value.name = "var"
    )[, ID := as.numeric(Cod)]

    names_test = unique(testing_data[, .(ID, Cod)])
    setkey(names_test, ID)

    data_validation$ID = names_test[data_validation, on = "ID", Cod]
    data_validation <- data.table(t(data_validation[, -1, with = FALSE]))
    setnames(data_validation, new = as.character(names_test$Cod))

    data_validation <- data.table(testing_data[, .(Date)], data_validation)

    # Store in a list
    common_columns = setdiff(names(data_validation), "Date")

    res_validate <- lapply(common_columns, function(col) {
      merge(test_data[, .(Date, Obs = get(col))], data_validation[, .(Date, Sim = get(col))], by = "Date", all = FALSE)
    })
    names(res_validate) <- common_columns
    validation_results <- lapply(res_validate, function(x) validate(x, Rain_threshold = Rain_threshold))

    final_results <- data.table::rbindlist(validation_results, idcol = "ID")
  } else {
    final_results <- NULL
  }

  ##############################################################################
  #                           Save the model if necessary                      #
  ##############################################################################
  if (save_model) {
    message("Saving model. Please wait.")
    if (is.null(name_save)) name_save = "Model_RFplus"
    name_saving <- paste0(name_save, ".nc")
    terra::writeCDF(Ensamble, filename = name_saving, overwrite=TRUE)
  }
  return(list(Ensamble = Ensamble, Validation = final_results))
}
