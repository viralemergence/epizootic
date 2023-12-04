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
  expect_type(result, "list")
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
  expect_type(result, "list")
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
  expect_type(result, "list")
  expect_true("initialize_attributes" %in% names(result))
  expect_true("initialize_replicate" %in% names(result))
  expect_true("calculate_at_season" %in% names(result))
  expect_true("calculate_at_replicate" %in% names(result))
  expect_true("finalize_attributes" %in% names(result))
})

test_that("disease_results handles invalid inputs correctly", {
  expect_error(disease_results("invalid", "invalid", "invalid", "invalid"))
})

test_that("initialize_attributes initializes result attributes correctly", {
  # Set up test data
  replicates <- 10
  time_steps <- 10
  seasons <- 4
  coordinates <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
  initial_abundance <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3)
  results_selection <- c("abundance", "occupancy")
  result_stages <- 2
  result_compartments <- 2

  # Call the disease_results function to generate the initialize_attributes function
  disease_results_functions <- disease_results(replicates, time_steps, seasons,
   coordinates, initial_abundance, results_selection, result_stages, 
   result_compartments)
  initialize_attributes <- disease_results_functions$initialize_attributes

  # Call the initialize_attributes function
  result <- initialize_attributes()

  # Verify the result attributes
  expect_type(result, "list")
  expect_true("abundance" %in% names(result))
  expect_true("all" %in% names(result))
  expect_true("occupancy" %in% names(result$all))
  expect_true("abundance" %in% names(result$all))
})

test_that("initialize_replicate initializes replicate result attributes correctly", {
  # Set up test data
  replicates <- 10
  time_steps <- 10
  seasons <- 4
  coordinates <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
  initial_abundance <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3)
  results_selection <- c("ema", "extinction_location")
  result_stages <- 2
  result_compartments <- 2

  # Call the disease_results function to generate the initialize_replicate function
  disease_results_functions <- disease_results(replicates, time_steps, seasons,
   coordinates, initial_abundance, results_selection, result_stages, 
   result_compartments)
  initialize_attributes <- disease_results_functions$initialize_attributes
  initialize_replicate <- disease_results_functions$initialize_replicate

  # Call the initialize_attributes function
  result <- initialize_attributes()

  # Call the initialize_replicate function
  rep_test <- initialize_replicate(result)

  # Verify the replicate result attributes
  expect_type(rep_test, "list")
  expect_true("abundance_count_min" %in% names(rep_test))
  expect_true("last_occupied_abundance_count" %in% names(rep_test))
})

test_that("calculate_at_season appends and calculates results correctly", {
  # Set up test data
  replicates <- 10
  time_steps <- 10
  seasons <- 4
  coordinates <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
  initial_abundance <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 3)
  results_selection <- c("abundance", "occupancy", "harvested")
  result_stages <- c(1, 0)
  result_compartments <- c(1, 0)

  # Call the disease_results function to generate the calculate_at_season function
  disease_results_functions <- disease_results(replicates, time_steps, seasons,
                                               coordinates, initial_abundance, results_selection, result_stages, 
                                               result_compartments)
  initialize_attributes <- disease_results_functions$initialize_attributes
  initialize_replicate <- disease_results_functions$initialize_replicate
  calculate_at_season <- disease_results_functions$calculate_at_season

  # Set up test inputs
  # Call the initialize_attributes function
  result <- initialize_attributes()
  # Call the initialize_replicate function
  rep_test <- initialize_replicate(result)
  r <- 1
  tm <- 1
  season <- 1
  segment_abundance <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 3)
  harvested <- matrix(c(1, 1, rep(0, 10)), nrow = 4, ncol = 3)

  # Call the calculate_at_season function
  result <- calculate_at_season(r, tm, season, segment_abundance, harvested, rep_test)

  # Verify the appended and calculated results
  expect_equal(result$abundance$mean, array(c(1, 1, 1, rep(0, 117)), dim = c(3, 10, 4)))
  expect_equal(result$abundance$sd, array(c(0, 0, 0, rep(0, 117)), dim = c(3, 10, 4)))
  expect_equal(result$abundance$min, array(c(1, 1, 1, rep(0, 117)), dim = c(3, 10, 4)))
  expect_equal(result$abundance$max, array(c(1, 1, 1, rep(0, 117)), dim = c(3, 10, 4)))
})