sir_inputs <- function(populations = 3,
                       segment_abundance = NULL,
                       occupied_indices = seq_len(populations),
                       ...) {
  if (is.null(segment_abundance)) {
    segment_abundance <- matrix(
      c(
        rep(50000, populations),
        rep(50000, populations),
        rep(0, populations),
        c(1, rep(0, populations - 1)),
        rep(0, populations),
        rep(0, populations)
      ),
      nrow = 6, byrow = TRUE
    )
  }
  defaults <- list(
    populations = populations,
    stages = 2,
    compartments = 3,
    abundance_threshold = 10,
    mortality = c(0.4, 0, 0.505, 0.105, 0.4, 0),
    mortality_unit = rep(1, 6),
    fecundity = 15,
    fecundity_unit = 1,
    fecundity_mask = c(0, 1, 0, 1, 0, 1),
    transmission = c(0.00002, 0.00001),
    transmission_unit = rep(0, 2),
    transmission_mask = c(1, 1, 0, 0, 0, 0),
    recovery = c(0.05714286, 0.05714286),
    recovery_unit = rep(0, 2),
    recovery_mask = c(0, 0, 1, 1, 0, 0),
    carrying_capacity = rep(150000, populations),
    breeding_season_length = rep(100, populations),
    segment_abundance = segment_abundance,
    occupied_indices = occupied_indices
  )
  utils::modifyList(defaults, list(...))
}

test_that("Summer SIR simulator works with valid inputs", {
  expect_silent(sir_model_summer(sir_inputs(populations = 63)))
})

test_that("Winter SIR simulator works with valid inputs", {
  inputs <- sir_inputs(populations = 63)
  inputs$fecundity <- NULL
  inputs$fecundity_unit <- NULL
  inputs$fecundity_mask <- NULL
  expect_silent(sir_model_winter(inputs))
})

test_that("sir_model_summer handles single active population", {
  inputs <- sir_inputs(
    segment_abundance = matrix(c(
      50000, 0, 0,
      50000, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE),
    occupied_indices = 1
  )

  result <- sir_model_summer(inputs)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(6, 3))
  # Population 1 should have changed, populations 2 and 3 should be zero
  expect_true(all(result[, 2] == 0))
  expect_true(all(result[, 3] == 0))
})

test_that("sir_model_winter handles single active population", {
  inputs <- sir_inputs(
    segment_abundance = matrix(c(
      50000, 0, 0,
      50000, 0, 0,
      0, 0, 0,
      1, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE),
    occupied_indices = 1
  )

  result <- sir_model_winter(inputs)
  expect_true(is.matrix(result))
  expect_equal(dim(result), c(6, 3))
})

test_that("SIR season functions return unchanged matrix when no populations occupied", {
  inputs <- sir_inputs(
    segment_abundance = matrix(0, nrow = 6, ncol = 3),
    occupied_indices = integer(0)
  )

  expect_equal(sir_model_summer(inputs), inputs$segment_abundance)
  expect_equal(sir_model_winter(inputs), inputs$segment_abundance)
})

test_that("sir_model_summer handles populations with no infection", {
  inputs <- sir_inputs(
    segment_abundance = matrix(c(
      50000, 40000, 60000,
      50000, 40000, 60000,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  )

  result <- sir_model_summer(inputs)
  expect_true(is.matrix(result))
  # No infection means the I and R compartments stay empty
  expect_true(all(result[3:6, ] == 0))
})

test_that("sir_model_winter handles populations with no infection", {
  inputs <- sir_inputs(
    segment_abundance = matrix(c(
      50000, 40000, 60000,
      50000, 40000, 60000,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  )

  result <- sir_model_winter(inputs)
  expect_true(is.matrix(result))
  expect_true(all(result[3:6, ] == 0))
})

test_that("SIR season functions handle varying breeding season lengths", {
  inputs <- sir_inputs(
    breeding_season_length = c(80, 100, 120),
    segment_abundance = matrix(c(
      50000, 50000, 50000,
      50000, 50000, 50000,
      0, 0, 0,
      1, 1, 1,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  )

  expect_equal(dim(sir_model_summer(inputs)), c(6, 3))
  expect_equal(dim(sir_model_winter(inputs)), c(6, 3))
})

test_that("sir_model_summer handles high transmission rates", {
  inputs <- sir_inputs(
    transmission = c(0.01, 0.008),
    segment_abundance = matrix(c(
      50000, 50000, 50000,
      50000, 50000, 50000,
      0, 0, 0,
      100, 100, 100,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  )

  result <- sir_model_summer(inputs)
  expect_true(is.matrix(result))
  # With high transmission, individuals should accumulate in the R compartment
  expect_true(sum(result[5:6, ]) > 0)
})

test_that("sir_model_summer respects abundance_threshold", {
  inputs <- sir_inputs(
    abundance_threshold = 50,
    mortality = c(0.9, 0.9, 0.95, 0.95, 0.9, 0.9),
    fecundity = 0.1,
    carrying_capacity = rep(1000, 3),
    segment_abundance = matrix(c(
      20, 60, 80,
      10, 40, 50,
      0, 0, 0,
      1, 1, 1,
      0, 0, 0,
      0, 0, 0
    ), nrow = 6, byrow = TRUE)
  )

  result <- sir_model_summer(inputs)
  # The first population starts below the threshold and must go extinct
  expect_true(all(result[, 1] == 0))
})

test_that("SIR season functions reject the wrong number of compartments", {
  inputs <- sir_inputs()
  inputs$compartments <- 4
  expect_error(sir_model_summer(inputs), "compartment")
  inputs <- sir_inputs()
  inputs$stages <- 3
  expect_error(sir_model_winter(inputs), "two-stage")
})

test_that("aspatial_sir conserves individuals in the absence of mortality", {
  result <- aspatial_sir(
    initial_pop = c(5000, 5000, 100, 100, 0, 0),
    season_length = 50,
    mortality = rep(0, 6),
    transmission = c(0.00002, 0.00001, 0, 0, 0, 0),
    recovery = c(0, 0, 0.05714286, 0.05714286, 0, 0),
    fecundity = rep(0, 6),
    abundance_threshold = 10,
    carrying_capacity = 150000,
    season = "non-breeding"
  )
  expect_length(result, 6)
  expect_equal(sum(result), 10200)
  # Recovery is permanent, so the recovered compartments must have filled
  expect_true(sum(result[5:6]) > 0)
})

test_that("aspatial_sir adds juveniles only in the breeding season", {
  args <- list(
    initial_pop = c(0, 5000, 0, 0, 0, 0),
    season_length = 30,
    mortality = rep(0, 6),
    transmission = rep(0, 6),
    recovery = rep(0, 6),
    fecundity = c(0, 0.1, 0, 0.1, 0, 0.1),
    abundance_threshold = 10,
    carrying_capacity = 150000
  )
  breeding <- do.call(aspatial_sir, c(args, season = "breeding"))
  non_breeding <- do.call(aspatial_sir, c(args, season = "non-breeding"))
  expect_true(breeding[1] > 0)
  expect_equal(non_breeding, c(0, 5000, 0, 0, 0, 0))
})
