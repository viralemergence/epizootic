library(testthat)

test_that("disease_results returns a list", {
  result <- disease_results(
    replicates = 10,
    time_steps = 10,
    seasons = 4,
    coordinates = data.frame(x = c(1, 2, 3), y = c(1, 2, 3)),
    initial_abundance = matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3),
    results_selection = NULL,
    result_stages = 2,
    result_compartments = 2
  )
  expect_is(result, "list")
})

test_that("disease_results handles NULL inputs correctly", {
  result <- disease_results(
    replicates = 10,
    time_steps = 10,
    seasons = 4,
    coordinates = NULL,
    initial_abundance = matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3),
    results_selection = NULL,
    result_stages = NULL, 
    result_compartments = NULL
  )
  expect_is(result, "list")
})

test_that("disease_results handles non-NULL inputs correctly", {
  result <- disease_results(
    replicates = 10,
    time_steps = 10,
    seasons = 4,
    coordinates = data.frame(x = c(1, 2, 3), y = c(1, 2, 3)),
    initial_abundance = matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3),
    results_selection = c("abundance", "occupancy"),
    result_stages = 2,
    result_compartments = 2
  )
  expect_is(result, "list")
})

test_that("disease_results handles invalid inputs correctly", {
  expect_error(disease_results("invalid", "invalid", "invalid", "invalid"))
})