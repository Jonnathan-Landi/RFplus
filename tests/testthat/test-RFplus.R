test_that("RFplus works with included data", {

  Covariates = list(
    MSWEP = terra::rast(system.file("extdata/MSWEP.nc", package = "RFplus")),
    CHIRPS = terra::rast(system.file("extdata/CHIRPS.nc", package = "RFplus")),
    DEM = terra::rast(system.file("extdata/DEM.nc", package = "RFplus"))
  )

  BD_Insitu = data.table::fread(system.file("extdata/BD_Insitu.csv", package = "RFplus"))
  Cords_Insitu = data.table::fread(system.file("extdata/Cords_Insitu.csv", package = "RFplus"))

  # Call function RFplus
  result = RFplus(Covariates, BD_Insitu, Cords_Insitu, ntree = 2000, threshold = 0.1,
                  n_round = 1, save_model = F, name_save = NULL, seed = 123)

  # Test
  expect_s4_class(result, "SpatRaster")
  expect_true(all(terra::values(result) >= 0))
})
