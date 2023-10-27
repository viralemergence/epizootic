test_that("User-defined transformation function behavior", {
  # Define a simple user-defined transformation function
  custom_transformation <- function(params) {
    params[["segment_abundance"]] <- params[["segment_abundance"]] * 2
    return(params)
  }

  new_function <- disease_transformation(
    transformation = custom_transformation,
    name = "test",
    density_stages = 1
  )

  result <- new_function(
    carrying_capacity = rep(100, 4),
    segment_abundance = matrix(1:12, nrow = 3),
    occupied_indices = c(1, 2),
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

  new_function <- disease_transformation(
    transformation = error_transformation,
    name = "test",
    density_stages = 1
  )

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

                new_function <- disease_transformation(
                  transformation = invalid_transformation,
                  name = "test",
                  density_stages = 1
                )

                expect_snapshot_warning(
                  new_function(
                    carrying_capacity = rep(100, 4),
                    segment_abundance = matrix(1:12, nrow = 3),
                    occupied_indices = c(1, 2),
                    r = 1,
                    tm = 1
                  )
                )

                invalid_transformation <- function(params) {
                  params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                  params[["carrying_capacity"]] <- 200
                  params[["carrying_capacity"]][2] <- -1
                  return(params)
                }

                expect_snapshot_warning(
                  new_function(
                    carrying_capacity = rep(100, 4),
                    segment_abundance = matrix(1:12, nrow = 3),
                    occupied_indices = c(1, 2),
                    r = 1,
                    tm = 1
                  )
                )
              })
