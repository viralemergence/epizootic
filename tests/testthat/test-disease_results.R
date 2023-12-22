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

  results_breakdown <- "segments"

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
  result <- disease_results_functions$initialize_attributes()
  expect_equal(result$abundance_segments |> names(),
               c("stage_1_compartment_1", "stage_2_compartment_1",
                 "stage_1_compartment_2", "stage_2_compartment_2"))
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

test_that("calculate at season (replicate)", {
  coordinates <- array(c(1:4, 4:1), c(7, 2))
  initial_abundance <- matrix(c(7, 13, 0, 26, 0, 39, 47,
                                2,  0, 6,  8, 0, 12, 13,
                                0,  3, 4,  6, 0,  9, 10),
                              nrow = 3, ncol = 7, byrow = TRUE)
  harvested <- round(initial_abundance*0.3)
  results_selection <- c("abundance", "ema", "extirpation",
                         "extinction_location", "harvested", "occupancy")
  # Single replicate and included stages combined
  result_functions <- disease_results(replicates = 1, time_steps = 10,
                                      seasons = 2, stages = 1, compartments = 3,
                                      coordinates = coordinates,
                                      initial_abundance = initial_abundance,
                                      results_selection = results_selection,
                                      results_breakdown = "pooled")
  results <- result_functions$initialize_replicate(result_functions$initialize_attributes())
  expected_results <- results
  expected_results$all$abundance[2, 1] <- sum(initial_abundance)
  expected_results$all$ema[2, 1] <- sum(initial_abundance)
  expected_results$all$harvested[2, 1] <- sum(harvested)
  expected_results$all$occupancy[2, 1] <- sum(+(colSums(initial_abundance) > 0))
  expected_results$abundance[, 2, 1] <- colSums(initial_abundance)
  expected_results$harvested[, 2, 1] <- colSums(harvested)
  expect_equal(result_functions$calculate_at_season(r = 1, tm = 2, season = 1,
                                                    segment_abundance = initial_abundance,
                                                    harvested, results),
               expected_results)
  # Calculate abundance and harvested separately
  expected_results_abundance <- expected_results
  expected_results_abundance$all$harvested <- results$all$harvested
  expected_results_abundance$harvested <- results$harvested
  expect_equal(result_functions$calculate_at_season(r = 1, tm = 2, season = 1,
                                                    segment_abundance = initial_abundance,
                                                    harvested = NULL, results),
               expected_results_abundance)
  expected_results_harvested <- results
  expected_results_harvested$all$harvested <- expected_results$all$harvested
  expected_results_harvested$harvested <- expected_results$harvested
  expect_equal(result_functions$calculate_at_season(r = 1, tm = 2, season = 1,
                                                      segment_abundance = NULL,
                                                      harvested, results),
               expected_results_harvested)
  # Multiple replicates and separated stage combinations
  results_selection <- c("abundance", "ema", "extirpation",
                         "extinction_location", "harvested", "occupancy",
                         "replicate") # "summarize" (default) or "replicate"
  result_functions <- disease_results(replicates = 2, time_steps = 10,
                                         seasons = 2, stages = 1,
                                         compartments = 3,
                                         coordinates = coordinates,
                                         initial_abundance = initial_abundance,
                                         results_selection = results_selection,
                                         results_breakdown = "compartments")
  results <- result_functions$initialize_replicate(result_functions$initialize_attributes())
  results <- result_functions$calculate_at_season(r = 1, tm = 2, season = 2,
                                                  segment_abundance = initial_abundance,
                                                  harvested, results)
  results <- result_functions$initialize_replicate(results)
  expected_results <- results
  abundance_r2 <- initial_abundance
  abundance_r2[, 1:6] <- 0
  expected_results$all$abundance[3, 2, 2] <- sum(abundance_r2)
  expected_results$all$abundance_compartments$compartment_1[3, 2, 2] <- sum(abundance_r2[1,])
  expected_results$all$abundance_compartments$compartment_2[3, 2, 2] <- sum(abundance_r2[2,])
  expected_results$all$abundance_compartments$compartment_3[3, 2, 2] <- sum(abundance_r2[3,])
  expected_results$all$ema[3, 2] <- mean(c(0, sum(abundance_r2)))
  expected_results$all$harvested[3, 2, 2] <- sum(harvested)
  expected_results$all$harvested_compartments$compartment_1[3, 2, 2] <- sum(harvested[1,])
  expected_results$all$harvested_compartments$compartment_2[3, 2, 2] <- sum(harvested[2,])
  expected_results$all$harvested_compartments$compartment_3[3, 2, 2] <- sum(harvested[3,])
  expected_results$all$occupancy[3, 2, 2] <- sum(+(colSums(abundance_r2) > 0))
  expected_results$abundance[, 3, 2, 2] <- colSums(abundance_r2)
  expected_results$abundance_compartments$compartment_1[, 3, 2, 2] <- abundance_r2[1,]
  expected_results$abundance_compartments$compartment_2[, 3, 2, 2] <- abundance_r2[2,]
  expected_results$abundance_compartments$compartment_3[, 3, 2, 2] <- abundance_r2[3,]
  expected_results$extirpation[c(1:4, 6), 2] <- 3
  expected_results$last_occupied_abundance_count <- array(colSums(abundance_r2))
  expected_results$abundance_count_min <- sum(abundance_r2)
  expected_results$harvested[, 3, 2, 2] <- colSums(harvested)
  expected_results$harvested_compartments$compartment_1[, 3, 2, 2] <- harvested[1,]
  expected_results$harvested_compartments$compartment_2[, 3, 2, 2] <- harvested[2,]
  expected_results$harvested_compartments$compartment_3[, 3, 2, 2] <- harvested[3,]
  results_2 <- result_functions$calculate_at_season(r = 2, tm = 3, season = 2,
                                                    segment_abundance = abundance_r2,
                                                    harvested, results)
  expect_equal(results_2, expected_results)
  expected_results_2 <- expected_results
  abundance_t4 <- abundance_r2*0
  expected_results_2$all$extirpation[2] <- 3.5
  expected_results_2$extirpation[7, 2] <- 3.5
  expected_results_2$abundance_count_min <- 0
  expect_equal(result_functions$calculate_at_season(r = 2, tm = 4, season = 1,
                                                    segment_abundance = abundance_t4,
                                                    harvested*0, results_2),
               expected_results_2)
  # Calculate abundance and harvested separately
  expected_results_abundance <- expected_results
  expected_results_abundance$all$harvested <- results$all$harvested
  expected_results_abundance$all$harvested_compartments <- results$all$harvested_compartments
  expected_results_abundance$harvested <- results$harvested
  expected_results_abundance$harvested_compartments <- results$harvested_compartments
  results_2_a <- result_functions$calculate_at_season(r = 2, tm = 3, season = 2,
                                                      segment_abundance = abundance_r2,
                                                      NULL, results)
  expect_equal(results_2_a, expected_results_abundance)
  expected_results_abundance <- expected_results_2
  expected_results_abundance$all$harvested <- expected_results$all$harvested
  expected_results_abundance$all$harvested_compartments <- expected_results$all$harvested_compartments
  expected_results_abundance$harvested <- expected_results$harvested
  expected_results_abundance$harvested_compartments <- expected_results$harvested_compartments
  expect_equal(result_functions$calculate_at_season(r = 2, tm = 4, season = 1,
                                                    segment_abundance = abundance_t4,
                                                    NULL, results_2),
               expected_results_2)
  expected_results_harvested <- results
  expected_results_harvested$all$harvested <- expected_results$all$harvested
  expected_results_harvested$all$harvested_compartments <- expected_results$all$harvested_compartments
  expected_results_harvested$harvested <- expected_results$harvested
  expected_results_harvested$harvested_compartments <- expected_results$harvested_compartments
  expect_equal(result_functions$calculate_at_season(r = 2, tm = 3, season = 2,
                                                    segment_abundance = NULL,
                                                    harvested, results),
               expected_results_harvested)
})

test_that("calculate at replicate", {
  coordinates = array(c(1:4, 4:1), c(7, 2))
  initial_abundance <- matrix(c(0,  0, 6,  8, 0,  0,  0,
                                0,  0, 4,  6, 0,  0,  0), nrow = 2, ncol = 7,
                              byrow = TRUE)
  results_selection <- c("extinction_location")
  # Single replicate and included stages combined
  result_functions <- disease_results(replicates = 1, time_steps = 10,
                                      seasons = 2, stages = 1, compartments = 2,
                                      coordinates = coordinates,
                                      initial_abundance = initial_abundance,
                                      results_selection = results_selection,
                                      results_breakdown = "pooled")
  results <- result_functions$initialize_replicate(result_functions$initialize_attributes())
  results <- result_functions$initialize_replicate(results)
  expected_results <- results
  expected_results$all$extinction_location[1,] <- 10/24*coordinates[3,] + 14/24*coordinates[4,]
  stage_abundance <- matrix(0, nrow = 2, ncol = 7)
  expect_equal(result_functions$calculate_at_replicate(1, stage_abundance, results),
               expected_results)
  # Multiple replicates
  result_functions <- disease_results(replicates = 3, time_steps = 10,
                                      seasons = 2, stages = 1,
                                      compartments = 2,
                                      coordinates = coordinates,
                                      initial_abundance = initial_abundance,
                                      results_selection = results_selection,
                                      results_breakdown = "pooled")
  results <- result_functions$initialize_replicate(result_functions$initialize_attributes())
  results <- result_functions$initialize_replicate(results)
  expected_results <- results
  expected_results$all$extinction_location[1,] <- 10/24*coordinates[3,] + 14/24*coordinates[4,]
  expected_results$all$extinction_location[2,] <- coordinates[4,]
  expected_results$last_occupied_abundance_count[3] <- 0
  results <- result_functions$calculate_at_replicate(1, stage_abundance, results)
  results$last_occupied_abundance_count[3] <- 0
  expect_equal(result_functions$calculate_at_replicate(2, stage_abundance, results),
               expected_results)
})

test_that("finalize_attributes", {
  coordinates = array(c(1:4, 4:1), c(7, 2))
  initial_abundance <- matrix(c(7, 13, 0, 26, 0, 39, 47,
                                2,  0, 6,  8, 0, 12, 13,
                                0,  3, 4,  6, 0,  9, 10), nrow = 3, ncol = 7, byrow = TRUE)
  harvested <- round(initial_abundance*0.3)
  results_selection <- c("abundance", "ema", "extirpation",
                         "extinction_location", "harvested", "occupancy", "summarize")
  # Summarized replicates and separated stage combinations
  result_functions <- disease_results(replicates = 5, time_steps = 10,
                                         seasons = 2, stages = 3,
                                         compartments = 1,
                                         coordinates, initial_abundance,
                                         results_selection = results_selection,
                                         results_breakdown = "stages")
  results <- result_functions$initialize_replicate(result_functions$initialize_attributes())
  # Run 5 replicates of a single time step + season
  abundance <- list(initial_abundance, initial_abundance, initial_abundance,
                    initial_abundance, initial_abundance*0)
  abundance[[2]][, 1] <- 0
  abundance[[3]][, 1:2] <- 0
  abundance[[4]][, 1:3] <- 0
  for (i in 1:5) {
    results <- result_functions$initialize_replicate(results)
    results <- result_functions$calculate_at_season(r = i, tm = 3, season = 1,
                                                      segment_abundance = abundance[[i]],
                                                      harvested + i, results)
  }
  expected_results <- results
  expected_results$all$abundance$sd[3, 1] <- sqrt(expected_results$all$abundance$sd[3, 1]/(5 - 1))
  expected_results$all$abundance_stages$stage_1$sd[3, 1] <- sqrt(expected_results$all$abundance_stages$stage_1$sd[3, 1]/(5 - 1))
  expected_results$all$abundance_stages$stage_2$sd[3, 1] <- sqrt(expected_results$all$abundance_stages$stage_2$sd[3, 1]/(5 - 1))
  expected_results$all$abundance_stages$stage_3$sd[3, 1] <- sqrt(expected_results$all$abundance_stages$stage_3$sd[3, 1]/(5 - 1))
  expected_results$all$harvested$sd[3] <- sqrt(expected_results$all$harvested$sd[3, 1]/(5 - 1))
  expected_results$all$harvested_stages$stage_1$sd[3, 1] <- sqrt(expected_results$all$harvested_stages$stage_1$sd[3, 1]/(5 - 1))
  expected_results$all$harvested_stages$stage_2$sd[3, 1] <- sqrt(expected_results$all$harvested_stages$stage_2$sd[3, 1]/(5 - 1))
  expected_results$all$harvested_stages$stage_3$sd[3, 1] <- sqrt(expected_results$all$harvested_stages$stage_3$sd[3, 1]/(5 - 1))
  expected_results$all$occupancy$sd[3, 1] <- sqrt(expected_results$all$occupancy$sd[3, 1]/(5 - 1))
  expected_results$abundance$sd[, 3, 1] <- sqrt(expected_results$abundance$sd[, 3, 1]/(5 - 1))
  expected_results$abundance_stages$stage_1$sd[, 3, 1] <- sqrt(expected_results$abundance_stages$stage_1$sd[, 3, 1]/(5 - 1))
  expected_results$abundance_stages$stage_2$sd[, 3, 1] <- sqrt(expected_results$abundance_stages$stage_2$sd[, 3, 1]/(5 - 1))
  expected_results$abundance_stages$stage_3$sd[, 3, 1] <- sqrt(expected_results$abundance_stages$stage_3$sd[, 3, 1]/(5 - 1))
  expected_results$harvested$sd[, 3, 1] <- sqrt(expected_results$harvested$sd[, 3, 1]/(5 - 1))
  expected_results$harvested_stages$stage_1$sd[, 3, 1] <- sqrt(expected_results$harvested_stages$stage_1$sd[, 3, 1]/(5 - 1))
  expected_results$harvested_stages$stage_2$sd[, 3, 1] <- sqrt(expected_results$harvested_stages$stage_2$sd[, 3, 1]/(5 - 1))
  expected_results$harvested_stages$stage_3$sd[, 3, 1] <- sqrt(expected_results$harvested_stages$stage_3$sd[, 3, 1]/(5 - 1))
  expected_results$occupancy$sd[, 3, 1] <- sqrt(expected_results$occupancy$sd[, 3, 1]/(5 - 1))
  expected_results$extirpation[which(is.na(expected_results$extirpation))] <- Inf
  expected_results$extirpation <- apply(expected_results$extirpation, 1, stats::fivenum)
  expected_results$extirpation[which(is.infinite(expected_results$extirpation))] <- NA
  expected_results$extirpation <- list(min = expected_results$extirpation[1,],
                                       q1 = expected_results$extirpation[2,],
                                       median = expected_results$extirpation[3,],
                                       q3 = expected_results$extirpation[4,],
                                       max = expected_results$extirpation[5,])
  expected_results$abundance_count_min <- NULL
  expected_results$last_occupied_abundance_count <- NULL
  expect_equal(result_functions$finalize_attributes(results), expected_results)
  # Full extirpation
  results <- result_functions$initialize_replicate(result_functions$initialize_attributes())
  for (i in 1:5) { # Run 5 replicates of a single time step + season
    results <- result_functions$initialize_replicate(results)
    results <- result_functions$calculate_at_season(r = i, tm = 3, season = 1,
                                                    segment_abundance = abundance[[i]],
                                                    harvested + i, results)
    results <- result_functions$calculate_at_season(r = i, tm = 5, season = 1,
                                                    segment_abundance = abundance[[i]]*0,
                                                    harvested + i, results)
  }
  expected_extirpation <- results$extirpation
  expected_extirpation <- list(mean = apply(expected_extirpation, 1, mean),
                               sd = apply(expected_extirpation, 1, stats::sd),
                               min = apply(expected_extirpation, 1, min),
                               max = apply(expected_extirpation, 1, max))
  expect_equal(result_functions$finalize_attributes(results)[["extirpation"]],
               expected_extirpation)
})
