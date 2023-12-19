library(testthat)

test_that("disease_results returns a correct list", {
  result_functions <- disease_results(
    replicates = 10,
    time_steps = 10,
    seasons = 4,
    stages = 2,
    compartments = 2,
    coordinates = data.frame(x = c(1, 2, 3), y = c(1, 2, 3)),
    initial_abundance = matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3),
    results_selection = c("abundance", "ema", "extirpation", 
                          "extinction_location", "harvested", "occupancy"),
    results_breakdown = "pooled")
  expect_named(result_functions, c("initialize_attributes", 
                                   "initialize_replicate", 
                                   "calculate_at_season",
                                   "calculate_at_replicate", 
                                   "finalize_attributes"))
  expect_type(result_functions$initialize_attributes, "closure")
  expect_null(formals(result_functions$initialize_attributes))
  expect_type(result_functions$initialize_replicate, "closure")
  expect_named(formals(result_functions$initialize_replicate), c("results"))
  expect_type(result_functions$calculate_at_season, "closure")
  expect_named(formals(result_functions$calculate_at_season), 
               c("r", "tm", "season", "segment_abundance", "harvested",
                 "results"))
  expect_type(result_functions$calculate_at_replicate, "closure")
  expect_named(formals(result_functions$calculate_at_replicate), 
               c("r", "segment_abundance", "results"))
  expect_type(result_functions$finalize_attributes, "closure")
  expect_named(formals(result_functions$finalize_attributes), c("results"))
})

test_that("disease_results handles NULL inputs correctly", {
  result <- disease_results(
    replicates = 10,
    time_steps = 10,
    seasons = 4,
    stages = 2,
    compartments = 2,
    coordinates = NULL,
    initial_abundance = matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3),
    results_selection = NULL,
    results_breakdown = "pooled"
  )
  expect_type(result, "list")
})

test_that("disease_results handles non-NULL inputs correctly", {
  result <- disease_results(
    replicates = 10,
    time_steps = 10,
    seasons = 4,
    stages = 2,
    compartments = 2,
    coordinates = data.frame(x = c(1, 2, 3), y = c(1, 2, 3)),
    initial_abundance = matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3),
    results_selection = c("abundance", "occupancy"),
    results_breakdown = "pooled"
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
  stages <- 2
  compartments <- 2
  coordinates <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
  initial_abundance <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3)
  results_selection <- c("abundance", "occupancy")
  results_breakdown <- "pooled"

  # Call the disease_results function to generate the initialize_attributes
  # function
  disease_results_functions <- disease_results(
      replicates,
      time_steps,
      seasons,
      stages,
      compartments,
      coordinates,
      initial_abundance,
      results_selection,
      results_breakdown
    )
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

test_that("initialize_replicate initializes replicate result attributes
          correctly", {
            # Set up test data
            replicates <- 10
            time_steps <- 10
            seasons <- 4
            stages <- 2
            compartments <- 2
            coordinates <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
            initial_abundance <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 4, ncol = 3)
            results_selection <- c("ema", "extinction_location")
            results_breakdown <- "pooled"

            # Call the disease_results function to generate the
            # initialize_replicate function
            disease_results_functions <- disease_results(
              replicates,
              time_steps,
              seasons,
              stages,
              compartments,
              coordinates,
              initial_abundance,
              results_selection,
              results_breakdown
            )
            initialize_attributes <-
              disease_results_functions$initialize_attributes
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
  stages <- 2
  compartments <- 2
  coordinates <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
  initial_abundance <- matrix(c(1, 2, 3, 4), nrow = 4, ncol = 3)
  results_selection <- c("abundance", "occupancy", "harvested")
  results_breakdown <- "pooled"

  # Call the disease_results function to generate the calculate_at_season function
  disease_results_functions <- disease_results(
      replicates,
      time_steps,
      seasons,
      stages,
      compartments,
      coordinates,
      initial_abundance,
      results_selection,
      results_breakdown
    )
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
  result <-
    calculate_at_season(r, tm, season, segment_abundance, harvested, rep_test)

  # Verify the appended and calculated results
  # Abundance
  expect_equal(result$abundance$mean, array(c(10, 10, 10, rep(0, 117)), dim = c(3, 10, 4)))
  expect_equal(result$abundance$sd, array(c(0, 0, 0, rep(0, 117)), dim = c(3, 10, 4)))
  expect_equal(result$abundance$min, array(c(10, 10, 10, rep(0, 117)), dim = c(3, 10, 4)))
  expect_equal(result$abundance$max, array(c(10, 10, 10, rep(0, 117)), dim = c(3, 10, 4)))
  # Harvested
  expect_equal(result$harvested$mean, array(c(2, rep(0, 119)), dim = c(3, 10, 4)))
  expect_equal(result$harvested$sd, array(c(0, rep(0, 119)), dim = c(3, 10, 4)))
  expect_equal(result$harvested$min, array(c(2, rep(0, 119)), dim = c(3, 10, 4)))
  expect_equal(result$harvested$max, array(c(2, rep(0, 119)), dim = c(3, 10, 4)))
  # All
  expect_equal(result$all$abundance$mean, array(c(30, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$abundance$sd, array(c(0, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$abundance$min, array(c(30, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$abundance$max, array(c(30, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$harvested$mean, array(c(2, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$harvested$sd, array(c(0, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$harvested$min, array(c(2, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$harvested$max, array(c(2, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$occupancy$mean, array(c(3, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$occupancy$sd, array(c(0, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$occupancy$min, array(c(3, rep(0, 39)), dim = c(10, 4)))
  expect_equal(result$all$occupancy$max, array(c(3, rep(0, 39)), dim = c(10, 4)))
})

test_that("calculate_at_season combines stages properly", {
  # Set up test data
  replicates <- 10
  time_steps <- 10
  seasons <- 4
  stages <- 2
  compartments <- 4
  coordinates <- data.frame(x = c(1, 2, 3), y = c(1, 2, 3))
  initial_abundance <- matrix(c(1:8), nrow = 8, ncol = 3)
  results_selection <- c("abundance", "occupancy", "harvested")
  results_breakdown <- "stages"

  # Call the disease_results function to generate the calculate_at_season function
  disease_results_functions <- disease_results(
      replicates,
      time_steps,
      seasons,
      stages,
      compartments,
      coordinates,
      initial_abundance,
      results_selection,
      results_breakdown
    )
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
  segment_abundance <- matrix(c(1:8), nrow = 8, ncol = 3)
  harvested <- matrix(c(2:9), nrow = 8, ncol = 3)

  # Call the calculate_at_season function
  result <-
    calculate_at_season(r, tm, season, segment_abundance, harvested, rep_test)

  # Verify the appended and calculated results
  # All abundance
  abundance_comp1 <- matrix(c(108, rep(0, 39)), nrow = 10)
  abundance_comp2 <- matrix(rep(0, 40), nrow = 10)
  abundance_comp3 <- matrix(c(48, rep(0, 39)), nrow = 10)
  abundance_comp4 <- matrix(c(60, rep(0, 39)), nrow = 10)
  expect_equal(result$all$abundance$mean, abundance_comp1)
  expect_equal(result$all$abundance$sd, abundance_comp2)
  expect_equal(result$all$abundance$min, abundance_comp1)
  expect_equal(result$all$abundance$max, abundance_comp1)
  expect_equal(result$all$abundance_stages[[1]],
               list(mean = abundance_comp3, sd = abundance_comp2,
                    min = abundance_comp3, max = abundance_comp3))
  expect_equal(result$all$abundance_stages[[2]],
               list(mean = abundance_comp4, sd = abundance_comp2,
                    min = abundance_comp4, max = abundance_comp4))
  expect_equal(names(result$all$abundance_stages), c("stage_1", "stage_2"))
  # All harvested
  harvested_comp1 <- matrix(c(132, rep(0, 39)), nrow = 10)
  harvested_comp2 <- matrix(rep(0, 40), nrow = 10)
  harvested_comp3 <- matrix(c(60, rep(0, 39)), nrow = 10)
  harvested_comp4 <- matrix(c(72, rep(0, 39)), nrow = 10)
  expect_equal(result$all$harvested$mean, harvested_comp1)
  expect_equal(result$all$harvested$sd, harvested_comp2)
  expect_equal(result$all$harvested$min, harvested_comp1)
  expect_equal(result$all$harvested$max, harvested_comp1)
  expect_equal(result$all$harvested_stages[[1]],
               list(mean = harvested_comp3, sd = harvested_comp2,
                    min = harvested_comp3, max = harvested_comp3))
  expect_equal(result$all$harvested_stages[[2]],
                list(mean = harvested_comp4, sd = harvested_comp2,
                      min = harvested_comp4, max = harvested_comp4))
  expect_equal(names(result$all$harvested_stages), c("stage_1", "stage_2"))
  # All occupancy
  occupancy_comp1 <- matrix(c(3, rep(0, 39)), nrow = 10)
  occupancy_comp2 <- matrix(rep(0, 40), nrow = 10)
  expect_equal(result$all$occupancy$mean, occupancy_comp1)
  expect_equal(result$all$occupancy$sd, occupancy_comp2)
  expect_equal(result$all$occupancy$min, occupancy_comp1)
  expect_equal(result$all$occupancy$max, occupancy_comp1)
  # Abundance
  abundance_comp5 <- matrix(c(36, 36, 36, rep(0, 27)), nrow = 3)
  abundance_comp6 <- matrix(rep(0, 30), nrow = 3)
  abundance_comp7 <- matrix(c(16, 16, 16, rep(0, 27)), nrow = 3)
  abundance_comp8 <- matrix(c(20, 20, 20, rep(0, 27)), nrow = 3)
  expect_equal(result$abundance$mean[,,1], abundance_comp5)
  expect_equal(result$abundance$sd[,,1], abundance_comp6)
  expect_equal(result$abundance$min[,,1], abundance_comp5)
  expect_equal(result$abundance$max[,,1], abundance_comp5)
  expect_equal(result$abundance_stages[[1]] |> map(\(x) x[,,1]),
               list(mean = abundance_comp7, sd = abundance_comp6,
                    min = abundance_comp7, max = abundance_comp7))
  expect_equal(result$abundance_stages[[2]] |> map(\(x) x[,,1]),
                list(mean = abundance_comp8, sd = abundance_comp6,
                      min = abundance_comp8, max = abundance_comp8))
})
