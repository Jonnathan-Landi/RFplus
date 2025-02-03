test_that("RFplus works with different methods and included data", {

  Covariates = list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "RFplus")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "RFplus")),
    DEM = terra::rast(system.file("extdata/DEM.nc", package = "RFplus"))
  )

  BD_Insitu = data.table::fread(system.file("extdata/BD_Insitu.csv", package = "RFplus"))
  Cords_Insitu = data.table::fread(system.file("extdata/Cords_Insitu.csv", package = "RFplus"))

  # # Test with "QUANT" method
  result_quant = RFplus(BD_Insitu, Cords_Insitu, Covariates, n_round = 1, wet.day = 0.1,
                        ntree = 2000, seed = 123, method = "QUANT", ratio = 15,
                        save_model = FALSE, name_save = NULL)

  expect_true(inherits(result_quant, "SpatRaster"))
  expect_true(all(terra::values(result_quant, na.rm = T) >= 0))
})
