test_that("Summer simulator works with valid inputs", {
  inputs <- list(
    populations = 63,
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
    carrying_capacity = rep(150000, 63),
    breeding_season_length = rep(100, 63),
    segment_abundance = c(c(50000, 50000, 0, 1, 0, 0, 0, 0),
                          rep(c(50000, 50000, 0, 0, 0, 0, 0, 0), 62)) |>
      matrix(nrow = 8),
    occupied_indices = c(1:63)
  )
  expect_silent(siri_model_summer(inputs))
})

test_that("Winter simulator works with valid inputs", {
  inputs <- list(
    populations = 63,
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
    carrying_capacity = rep(150000, 63),
    breeding_season_length = rep(100, 63),
    segment_abundance = c(c(50000, 50000, 0, 1, 0, 0, 0, 0),
                          rep(c(50000, 50000, 0, 0, 0, 0, 0, 0), 62)) |>
      matrix(nrow = 8),
    occupied_indices = c(1:63)
  )
  expect_silent(siri_model_winter(inputs))
})
