#' RFplus
#'
#' Descripción detallada de lo que hace la función. Explica cómo funciona,
#' qué parámetros acepta y qué resultados devuelve. Si es necesario, también
#' puedes poner ejemplos de cómo usarla.
#'
#' @param x Descripción del parámetro `x`. Puedes incluir el tipo de dato y lo
#' que representa en la función.
#' @return Qué valor devuelve la función y qué significa.
#' @examples
#' # Ejemplo de cómo usar la función:
#' mi_funcion(3)
#' mi_funcion(4)

RFplus = function(Covariates, BD.Insitu, Cords_Insitu, ntree = 2000) {
  ##############################################################################
  #                      Check the input data of the covariates                #
  ##############################################################################
  # Verify the class of covariates
  class_check = lapply(Covariates, class)
  if (!all(sapply(class_check, function(x) x == class_check[[1]]))) stop("The class of the covariates are different (all classes should be similar).")

  # Verify the extent of covariates
  ext_check = lapply(Covariates, function(x) terra::ext(x))
  if (!all(sapply(ext_check, function(x) x == ext_check[[1]]))) stop("The extension of the covariates are different (all extensions should be similar).")

  # Verify the crc of covariates
  crs_check = lapply(Covariates, function(x) terra::crs(x))
  if (!all(sapply(crs_check, function(x) x == crs_check[[1]]))) stop("The crs of the covariates are different (all crs should be similar).")

  ##############################################################################
  #                    Check input data from on-site stations                  #
  ##############################################################################
  # Verify the BD.Insitu is data.table
  if (!is.data.table(BD.Insitu)) stop("The data of the on-site stations should be a data.table.")

  # Verify the columns of Cords_Insitu
  if (!is.data.table(Cords_Insitu)) stop("The coordinate data of the on-site stations must be a data.table.")

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
  data_train = melt(
    BD.Insitu,
    id.vars = "Date",
    variable.name = "Cod",
    value.name = "var"
  ) %>%
    mutate(ID = as.numeric(factor(Cod)))

  # Extract Dates
  Dates_extracted = unique(data_train$Date)
  Points_Train = merge.data.table(data_train, Cords_Insitu, by = "Cod")
  Points_Train = unique(Points_Train[,.(ID, Cod, X, Y)])

  # Extraction of a sample to save the correction results
  Sample_lyrs = Covariates[[2]][[1]]
  Sample_lyrs = Sample_lyrs * 0

  # RandonForest function for bias correction
  RF_train = function(day_COV, fecha) {
    data_obs = data_train[Date == as.Date(fecha), ] # Filtering the observation 'i'
    data_obs = data_obs[!is.na(var)] # Method for training models independently of NAs

    points_EstTrain = merge.data.table(data_obs, Points_Train, by = c("ID", "Cod"))
    points_EstTrain = terra::vect(points_EstTrain, geom = c("X", "Y"), crs = terra::crs(Sample_lyrs))

    # Calculate the Euclidean distance of my training set
    distance_ED = terra::distance(Sample_lyrs, points_EstTrain, rasterize = FALSE)
    names(distance_ED) = "Distance_ED"

    # Adequacy of covariates for training.
    day_COV$distance_ED = distance_ED
    names_covs = names(day_COV)

    # Merge the covariates with the training data (BUG BUG BUG)
    data_cov = lapply(day_COV, function(x) terra::extract(x, points_EstTrain))
    data_cov = Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_cov)
    names(data_cov)[2:length(data_cov)] = names

    # Generate my dependent variable file
    obs = RF.data_Train[Date == as.Date(fecha), ] # Buscar metodologia luego para asignar fecha
    obs$ID = 1:nrow(obs)
    Fecha = unique(obs$Date)
    obs = obs[,.(ID, var)]
    dt.train = merge(obs, data_cov, by = "ID")


    dt.train = dt.train[, !("ID"), with = FALSE] # Usando data.table
    RF_model = randomForest::randomForest(var ~ ., data = dt.train, na.action = na.omit, ntree = 2000)

    # Primeras predicciones
    res_f = predict(RF_model, dt.train[,-1])
    res_f = data.table(var = dt.train$var, sim = res_f)
    res_f$residuals = res_f$var - res_f$sim
    res_f = res_f[,.(residuals)]

    dt.train_resi = cbind(res_f, dt.train)
    dt.train_resi = dt.train_resi[, !("var"), with = FALSE]
    RF_model_f = randomForest::randomForest(residuals ~ ., data = dt.train_resi, na.action = na.omit, ntree = 2000)


    cov.day = terra::rast(day_COV)

    # Realizar la predicción utilizando el modelo de Random Forest
    result_RF = terra::predict(cov.day, RF_model)
    result_RF_rsd = terra::predict(cov.day, RF_model_f)
    result_final = result_RF + result_RF_rsd

    #
    # plot(result)
    # vv = terra::extract(result_final, sample.points)
    # valido = merge(vv, obs, by = "ID")
    # valido$residuals = valido$lyr1 - valido$var
    # gof(valido$lyr1, valido$var)
    ## Version Prueba que Evalua en base a los residuos
    return(result_final)
  }


  # Creo mi archivo de puntos
  #   RF.dtrain.points = merge(RF.data_Train, cords_bqd, by = "Cod")
  #   sample.points = unique(RF.dtrain.points[,.(Cod, X, Y)])
  #
  # #  RF.dtrain.points = terra::vect(RF.dtrain.points, geom = c("lon", "lat"), crs = terra::crs(DEM)) # Cambiar luego por Dem dentro de la lista
  #   sample.points = terra::vect(sample.points, geom = c("X", "Y"), crs = terra::crs(Covariates[[1]]))
  ##############################################################################
  #                       Creo mi archivo de coordenadas                       #
  ##############################################################################


  ##############################################################################
  #                                  Entrenamiento diario                      #
  ##############################################################################
  # day_COV = sapply(Covariates, function(x) subset(x, n), simplify = FALSE) n = numero de capa
  #distance_ED = terra::distance(day_COV[[1]], sample.points, rasterize=FALSE)
  #  names(distance_ED) = "Distance_ED"
  # day_COV$distance_ED = distance_ED

  # day_predict=

  RF_train = function(day_COV, fecha) {

    names = names(day_COV)

    # Extract layer i from my covariables
    data_cov = lapply(day_COV, function(x) terra::extract(x, sample.points, method = "simple"))
    data_cov = Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_cov)
    names(data_cov)[2:length(data_cov)] = names

    # Generate my dependent variable file
    obs = RF.data_Train[Date == as.Date(fecha), ] # Buscar metodologia luego para asignar fecha
    obs$ID = 1:nrow(obs)
    Fecha = unique(obs$Date)
    obs = obs[,.(ID, var)]
    dt.train = merge(obs, data_cov, by = "ID")


    dt.train = dt.train[, !("ID"), with = FALSE] # Usando data.table
    RF_model = randomForest::randomForest(var ~ ., data = dt.train, na.action = na.omit, ntree = 2000)

    # Primeras predicciones
    res_f = predict(RF_model, dt.train[,-1])
    res_f = data.table(var = dt.train$var, sim = res_f)
    res_f$residuals = res_f$var - res_f$sim
    res_f = res_f[,.(residuals)]

    dt.train_resi = cbind(res_f, dt.train)
    dt.train_resi = dt.train_resi[, !("var"), with = FALSE]
    RF_model_f = randomForest::randomForest(residuals ~ ., data = dt.train_resi, na.action = na.omit, ntree = 2000)


    cov.day = terra::rast(day_COV)

    # Realizar la predicción utilizando el modelo de Random Forest
    result_RF = terra::predict(cov.day, RF_model)
    result_RF_rsd = terra::predict(cov.day, RF_model_f)
    result_final = result_RF + result_RF_rsd

    #
    # plot(result)
    # vv = terra::extract(result_final, sample.points)
    # valido = merge(vv, obs, by = "ID")
    # valido$residuals = valido$lyr1 - valido$var
    # gof(valido$lyr1, valido$var)
    ## Version Prueba que Evalua en base a los residuos
    return(result_final)
  }

  # Calculo de distancia euclidiana como Covariable.
  distance_ED = terra::distance(Covariates[[1]][[1]], sample.points, rasterize=FALSE)
  names(distance_ED) = "Distance_ED"
  res_prediction = list()
  for (i in 1:length(Dates_complete)) {
    message("Iteracion: ", i, " de ", length(Dates_complete))
    day_COV = sapply(Covariates, function(x) subset(x, i), simplify = FALSE)
    day_COV$distance_ED = distance_ED
    fecha = as.Date(Dates_complete[i])

    cap = RF_train(day_COV, fecha)
    names(cap) = fecha
    res_prediction[[i]] = cap
  }

  raster_complete = terra::rast(res_prediction)

  # # Metodo 3. Correccion con QM o DQM (por investigar los dos)
  # presc = terra::extract(raster_complete, sample.points, method = "simple")
  # dat_QM = terra::extract(raster_complete, sample.points, method = "simple")
  # nombres = names(dat_QM)
  # nombres = nombres[-1]
  # dat_QM = t(dat_QM)
  # dat_QM = data.table(dat_QM)
  # dat_QM = dat_QM[-1,]
  #
  # dat_QM = cbind(Date = nombres, dat_QM)
  # dat_QM$Date = as.IDate(dat_QM$Date)
  #
  # colnames(dat_QM) = names(bqd)
  #
  #
  # dat_QM = melt(
  #   dat_QM,
  #   id.vars = "Date", # Columna fija
  #   variable.name = "Cod", # Nombre para las variables originales
  #   value.name = "var" # Nombre para los valores
  # )
  #
  # merge_mf = merge(dat_QM, RF.data_Train, by = c("Date", "Cod"))
  # names(merge_mf) = c("Date", "Cod", "Sim", "Obs")
  #
  # # Cuantiles
  # cuantiles = fitQmapRQUANT(merge_mf$Obs, merge_mf$Sim, wet.day = FALSE, Qstep = 0.01,nboot = 100)
  # qfit = cuantiles$par$fitq
  #
  #
  # RF_cuantile = function(day_COV, fecha, n) {
  #
  #   names = names(day_COV)
  #
  #   # Extract layer i from my covariables
  #   data_cov = lapply(day_COV, function(x) terra::extract(x, sample.points, method = "simple"))
  #   data_cov = Reduce(function(x, y) merge(x, y, by = "ID", all = TRUE), data_cov)
  #   names(data_cov)[2:length(data_cov)] = names
  #
  #   # Generate my dependent variable file
  #   obs = rep(qfit[[n]], nrow(data_cov))
  #   obs = data.table(obs = obs)
  #
  #   obs$ID = 1:nrow(obs)
  #   Fecha = unique(obs$Date)
  #   obs = obs[,.(ID, obs)]
  #   dt.train = merge(obs, data_cov, by = "ID")
  #
  #
  #   dt.train = dt.train[, !("ID"), with = FALSE] # Usando data.table
  #   RF_model = randomForest::randomForest(obs ~ ., data = dt.train, na.action = na.omit, ntree = 2000)
  #
  #   cov.day = terra::rast(day_COV)
  #
  #   # Realizar la predicción utilizando el modelo de Random Forest
  #   result_RF = terra::predict(cov.day, RF_model)
  #
  #   #
  #   # plot(result)
  #   # vv = terra::extract(result_final, sample.points)
  #   # valido = merge(vv, obs, by = "ID")
  #   # valido$residuals = valido$lyr1 - valido$var
  #   # gof(valido$lyr1, valido$var)
  #   ## Version Prueba que Evalua en base a los residuos
  #   return(result_RF)
  # }
  #
  # cuantiles_list = list()
  # for (i in 1:length(Dates_complete)) {
  #   message("Iteracion: ", i, " de ", length(Dates_complete))
  #   day_COV = sapply(Covariates, function(x) subset(x, i), simplify = FALSE)
  #   day_COV$distance_ED = distance_ED
  #   fecha = as.Date(Dates_complete[i])
  #
  #   cap = RF_cuantile(day_COV, fecha, i)
  #   names(cap) = fecha
  #   cuantiles_list[[i]] = cap
  # }
  #
  # raster_cuantiles = terra::rast(cuantiles_list)

  #qm1 = doQmapRQUANT(merge_mf$Sim, cuantiles)



  # merge_f = cbind(merge_mf, cuantiles_obs)

  # Prediccion sobre cuantiles.


  return(raster_complete)
}

################################################################################
# Debug
library(data.table)
library(terra)
library(tidyr)
library(dplyr)
library(randomForest)
library(qmap)

BD.Insitu = fread("D:/ultimo_WINDOWS/Sequias_GIZ/RFmerge/Red_Actual/Archivos_Train/Est_trainPpts.csv")
BD.Insitu = BD.Insitu[c(1:100),]
Cords_Insitu = fread("C:/Users/Jonna/Desktop/Artículos/Repartición espacial/Scripts/Testing/PPGIs.csv")

CHIRPS = rast("C:/Users/Jonna/Desktop/Artículos/Repartición espacial/Scripts/Testing/CHIRPS.nc")
MSWEP = rast("C:/Users/Jonna/Desktop/Artículos/Repartición espacial/Scripts/Testing/MSWEP.nc")
DEM = rast("C:/Users/Jonna/Desktop/Artículos/Repartición espacial/Scripts/Testing/DEM.nc")

CHIRPS = CHIRPS[[1:100]]
MSWEP = MSWEP[[1:100]]

Covariates = list(DEM = DEM, CHIRPS = CHIRPS, MSWEP = MSWEP)
rm(CHIRPS, MSWEP, DEM)

probando = RFbias(Covariates, bqd, cords_bqd)
setwd("C:/Users/Jonna/Desktop/Artículos/Repartición espacial/Scripts/Testing/Results_1")
terra::writeCDF(modelo_1[[1]], filename= "CapaTest_1.nc", overwrite=TRUE)

# cords_sf = st_as_sf(cords_bqd, coords = c("lon", "lat"), crs = 4326)
# cords_utm_sf = st_transform(cords_sf, crs = 32717)
# cords_bqd[, X := st_coordinates(cords_utm_sf)[, 1]]
# cords_bqd[, Y := st_coordinates(cords_utm_sf)[, 2]]
# cords_bqd = cords_bqd[,.(Cod, X, Y)]
# write.csv(cords_bqd, "C:/Users/Jonna/Desktop/Artículos/Repartición espacial/Scripts/Testing/PPGIs.csv", row.names = FALSE)
# #
# n = length(colnames(bqd)) # -1 para quitar la columna timestam
# column_names <- c("TIMESTAMP", paste0("M", sprintf("%03d", 1:(n - 1))))
# names(Est_TrainPPts) = column_names
#
# Est_TrainPPts$TIMESTAMP <- as.Date(Est_TrainPPts$TIMESTAMP, format = "%Y-%m-%d", tz = "UTC")
# Est_TrainPPts = zoo::zoo(Est_TrainPPts[,-1], order.by = Est_TrainPPts$TIMESTAMP)
