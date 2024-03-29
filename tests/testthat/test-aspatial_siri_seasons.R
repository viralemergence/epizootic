test_that("Summer simulator works with valid inputs", {
  inputs <- list(
    populations = 6355,
    stages = 2,
    compartments = 4,
    abundance_threshold = 10,
    mortality = c(0.4, 0, 0.505, 0.105, 0.4, 0, 0.45, 0.05),
    mortality_unit = rep(1, 8),
    fecundity = 15,
    fecundity_unit = 1,
    fecundity_mask = c(0, 1, 0, 1, 0, 1, 0, 1),
    transmission = c(0.00002, 0.00001, 7.84e-06, 3.92e-06),
    transmission_unit = rep(0, 4),
    transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0),
    recovery = c(0.05714286, 0.05714286, 0.1, 0.1),
    recovery_unit = rep(0, 4),
    recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1),
    carrying_capacity = rep(150000, 6355),
    breeding_season_length = rep(100, 6355),
    segment_abundance = c(c(50000, 50000, 0, 1, 0, 0, 0, 0),
                          rep(c(50000, 50000, 0, 0, 0, 0, 0, 0), 6354)) |>
      matrix(nrow = 8),
    occupied_indices = c(1:1000)
  )
  expect_silent(siri_model_summer(inputs))
})

test_that("Winter simulator works with valid inputs", {
  inputs <- list(
    populations = 6355,
    stages = 2,
    compartments = 4,
    abundance_threshold = 10,
    mortality = c(0.4, 0, 0.505, 0.105, 0.4, 0, 0.45, 0.05),
    mortality_unit = rep(1, 8),
    transmission = c(0.00002, 0.00001, 7.84e-06, 3.92e-06),
    transmission_unit = rep(0, 4),
    transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0),
    recovery = c(0.05714286, 0.05714286, 0.1, 0.1),
    recovery_unit = rep(0, 4),
    recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1),
    carrying_capacity = rep(150000, 6355),
    breeding_season_length = rep(100, 6355),
    segment_abundance = c(c(50000, 50000, 0, 1, 0, 0, 0, 0),
                          rep(c(50000, 50000, 0, 0, 0, 0, 0, 0), 6354)) |>
      matrix(nrow = 8),
    occupied_indices = c(1:6355)
  )
  expect_silent(siri_model_winter(inputs))
})

test_that("aspatial_siri works with valid inputs", {
  expect_silent(aspatial_siri(
    initial_pop = c(50000, 50000, 0, 1, 0, 0, 0, 0),
    mortality = c(0.004, 0, 0.00505, 0.00105, 0.004, 0, 0.0045, 5e-04),
    fecundity = 15/182,
    transmission = c(0.00002, 0.00001, 7.84e-06, 3.92e-06),
    recovery = c(0.05714286, 0.05714286, 0.1, 0.1),
    carrying_capacity = 150000,
    abundance_threshold = 10,
    season = "breeding"
  ))
})
