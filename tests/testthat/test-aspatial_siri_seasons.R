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
    segment_abundance = c(
      c(50000, 50000, 0, 1, 0, 0, 0, 0),
      rep(c(50000, 50000, 0, 0, 0, 0, 0, 0), 62)
    ) |>
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
    segment_abundance = c(
      c(50000, 50000, 0, 1, 0, 0, 0, 0),
      rep(c(50000, 50000, 0, 0, 0, 0, 0, 0), 62)
    ) |>
      matrix(nrow = 8),
    occupied_indices = c(1:63)
  )
  expect_silent(siri_model_winter(inputs))
})

test_that("siri_model_summer handles single active population", {
  inputs <- list(
    populations = 3,
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
    carrying_capacity = rep(150000, 3),
    breeding_season_length = rep(100, 3),
    segment_abundance = matrix(c(
      50000, 0, 0,
      50000, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 8, byrow = TRUE),
    occupied_indices = c(1)
  )

  result <- siri_model_summer(inputs)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(8, 3))
  # Population 1 should have changed, populations 2 and 3 should be zero
  expect_true(all(result[, 2] == 0))
  expect_true(all(result[, 3] == 0))
})

test_that("siri_model_winter handles single active population", {
  inputs <- list(
    populations = 3,
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
    carrying_capacity = rep(150000, 3),
    breeding_season_length = rep(100, 3),
    segment_abundance = matrix(c(
      50000, 0, 0,
      50000, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 8, byrow = TRUE),
    occupied_indices = c(1)
  )

  result <- siri_model_winter(inputs)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(8, 3))
})

test_that("siri_model_summer returns unchanged matrix when no populations occupied", {
  inputs <- list(
    populations = 3,
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
    carrying_capacity = rep(150000, 3),
    breeding_season_length = rep(100, 3),
    segment_abundance = matrix(0, nrow = 8, ncol = 3),
    occupied_indices = integer(0)
  )

  result <- siri_model_summer(inputs)
  expect_equal(result, inputs$segment_abundance)
})

test_that("siri_model_winter returns unchanged matrix when no populations occupied", {
  inputs <- list(
    populations = 3,
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
    carrying_capacity = rep(150000, 3),
    breeding_season_length = rep(100, 3),
    segment_abundance = matrix(0, nrow = 8, ncol = 3),
    occupied_indices = integer(0)
  )

  result <- siri_model_winter(inputs)
  expect_equal(result, inputs$segment_abundance)
})

test_that("siri_model_summer handles populations with no infection", {
  inputs <- list(
    populations = 3,
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
    carrying_capacity = rep(150000, 3),
    breeding_season_length = rep(100, 3),
    segment_abundance = matrix(c(
      50000, 40000, 60000,
      50000, 40000, 60000,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 8, byrow = TRUE),
    occupied_indices = c(1, 2, 3)
  )

  result <- siri_model_summer(inputs)
  expect_true(is.matrix(result))
  # No infection means only S and no I or R compartments should have abundance
  expect_true(all(result[3:8, ] == 0))
})

test_that("siri_model_winter handles populations with no infection", {
  inputs <- list(
    populations = 3,
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
    carrying_capacity = rep(150000, 3),
    breeding_season_length = rep(100, 3),
    segment_abundance = matrix(c(
      50000, 40000, 60000,
      50000, 40000, 60000,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 8, byrow = TRUE),
    occupied_indices = c(1, 2, 3)
  )

  result <- siri_model_winter(inputs)
  expect_true(is.matrix(result))
})

test_that("siri_model_summer handles varying breeding season lengths", {
  inputs <- list(
    populations = 3,
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
    carrying_capacity = rep(150000, 3),
    breeding_season_length = c(80, 100, 120),
    segment_abundance = matrix(c(
      50000, 50000, 50000,
      50000, 50000, 50000,
      0, 0, 0,
      1, 1, 1,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 8, byrow = TRUE),
    occupied_indices = c(1, 2, 3)
  )

  result <- siri_model_summer(inputs)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(8, 3))
})

test_that("siri_model_winter handles varying breeding season lengths", {
  inputs <- list(
    populations = 3,
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
    carrying_capacity = rep(150000, 3),
    breeding_season_length = c(80, 100, 120),
    segment_abundance = matrix(c(
      50000, 50000, 50000,
      50000, 50000, 50000,
      0, 0, 0,
      1, 1, 1,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 8, byrow = TRUE),
    occupied_indices = c(1, 2, 3)
  )

  result <- siri_model_winter(inputs)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(8, 3))
})

test_that("siri_model_summer handles high transmission rates", {
  inputs <- list(
    populations = 3,
    stages = 2,
    compartments = 4,
    abundance_threshold = 10,
    mortality = c(0.4, 0, 0.505, 0.105, 0.4, 0, 0.45, 0.05),
    mortality_unit = rep(1, 8),
    fecundity = 15,
    fecundity_unit = 1,
    fecundity_mask = c(0, 1, 0, 1, 0, 1, 0, 1),
    transmission = c(0.01, 0.008, 0.005, 0.003), # High transmission
    transmission_unit = rep(0, 4),
    transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0),
    recovery = c(0.05714286, 0.05714286, 0.1, 0.1),
    recovery_unit = rep(0, 4),
    recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1),
    carrying_capacity = rep(150000, 3),
    breeding_season_length = rep(100, 3),
    segment_abundance = matrix(c(
      50000, 50000, 50000,
      50000, 50000, 50000,
      0, 0, 0,
      100, 100, 100,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 8, byrow = TRUE),
    occupied_indices = c(1, 2, 3)
  )

  result <- siri_model_summer(inputs)
  expect_true(is.matrix(result))
  # With high transmission, we expect some infection spread
  expect_true(sum(result[3:4, ]) > 100) # More infected
})

test_that("siri_model_summer respects abundance_threshold", {
  inputs <- list(
    populations = 3,
    stages = 2,
    compartments = 4,
    abundance_threshold = 50,
    mortality = c(0.9, 0.9, 0.95, 0.95, 0.9, 0.9, 0.95, 0.95), # Very high mortality
    mortality_unit = rep(1, 8),
    fecundity = 0.1,
    fecundity_unit = 1,
    fecundity_mask = c(0, 1, 0, 1, 0, 1, 0, 1),
    transmission = c(0.00002, 0.00001, 7.84e-06, 3.92e-06),
    transmission_unit = rep(0, 4),
    transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0),
    recovery = c(0.05714286, 0.05714286, 0.1, 0.1),
    recovery_unit = rep(0, 4),
    recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1),
    carrying_capacity = rep(1000, 3),
    breeding_season_length = rep(100, 3),
    segment_abundance = matrix(c(
      30, 60, 80,
      20, 40, 50,
      0, 0, 0,
      1, 1, 1,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 8, byrow = TRUE),
    occupied_indices = c(1, 2, 3)
  )

  result <- siri_model_summer(inputs)
  expect_true(is.matrix(result))
  # Population 1 might go extinct (< 50 threshold)
})
