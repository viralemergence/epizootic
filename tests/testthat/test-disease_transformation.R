test_that("User-defined transformation function behavior", {
  # Define a simple user-defined transformation function
  custom_transformation <- function(params) {
    params[["segment_abundance"]] <- params[["segment_abundance"]] * 2
    return(params)
  }

  params <- list(transformation = custom_transformation,
                 name = "test",
                 stages = 1)

  new_function <- disease_transformation(params)

  result <- new_function(
    carrying_capacity = rep(100, 4),
    segment_abundance = matrix(1:12, nrow = 3),
    occupied_indices = c(1, 2),
    breeding_season_length = rep(100, 4),
    r = 1,
    tm = 1
  )

  # Verify the transformation result
  expect_identical(result[["segment_abundance"]],
                   matrix(seq(2, 24, by = 2), nrow = 3))
  expect_identical(result[["carrying_capacity"]], rep(100, 4))
})

test_that("Error handling within user-defined transformation function", {
  # Define a user-defined transformation function that produces an error
  error_transformation <- function(params) {
    stop("Custom transformation function error")
  }

  params <- list(transformation = error_transformation,
                 name = "test",
                 stages = 1)

  new_function <- disease_transformation(params)

  # Expect an error when calling disease_transformation
  expect_error(
    new_function(
      carrying_capacity = rep(100, 4),
      segment_abundance = matrix(1:12, nrow = 3),
      occupied_indices = c(1, 2),
      r = 1,
      tm = 1
    )
  )
})

test_that_cli("Warnings for negative or non-finite values in transformed data",
              {
                # Define a user-defined transformation function that returns
                # non-finite values
                invalid_transformation <- function(params) {
                  params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                  params[["segment_abundance"]][1, 1] <- Inf
                  params[["carrying_capacity"]] <- 200
                  params[["carrying_capacity"]][2] <- -1
                  return(params)
                }

                params <- list(transformation = invalid_transformation,
                               name = "test",
                               stages = 1)

                new_function <- disease_transformation(params)

                expect_snapshot(
                  invisible(new_function(
                    carrying_capacity = rep(100, 4),
                    segment_abundance = matrix(1:12, nrow = 3),
                    occupied_indices = c(1, 2),
                    breeding_season_length = rep(100, 4),
                    r = 1,
                    tm = 1
                  ))
                )

                invalid_transformation <- function(params) {
                  params[["segment_abundance"]] <- matrix(1:12, nrow = 3)
                  params[["carrying_capacity"]] <- 200
                  params[["carrying_capacity"]][2] <- -1
                  return(params)
                }

                params <- list(transformation = invalid_transformation,
                               name = "test",
                               stages = 1)

                new_function <- disease_transformation(params)

                expect_snapshot(
                  invisible(new_function(
                    carrying_capacity = rep(100, 4),
                    segment_abundance = matrix(1:12, nrow = 3),
                    occupied_indices = c(1, 2),
                    breeding_season_length = rep(100, 4),
                    r = 1,
                    tm = 1
                  )
                ))
              })

test_that("The function works as expected with siri_model_summer", {

  params <- list(
    time_steps = 5,
    transformation = siri_model_summer,
    name = "Summer SIRI model",
    stages = 2,
    populations = 63,
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
    seasons = 2
  )

  new_function <- disease_transformation(params)

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
  set.seed(153)  # Set a fixed random seed

  # Call the first function
  result_from_new_function <- new_function(
    r = 1,
    tm = 1,
    carrying_capacity = rep(150000, 63),
    segment_abundance = c(c(50000, 50000, 0, 1, 0, 0, 0, 0),
                          rep(c(50000, 50000, 0, 0, 0, 0, 0, 0), 62)) |>
      matrix(nrow = 8),
    breeding_season_length = rep(100, 63),
    occupied_indices = c(1:63)
  )

  # Call the second function with the same random seed
  set.seed(153)
  result_from_siri_model_summer <- siri_model_summer(inputs)

  # Use expect_equal to compare the results
  expect_equal(result_from_new_function,
             list(segment_abundance = result_from_siri_model_summer),
             tolerance = 25)

})
