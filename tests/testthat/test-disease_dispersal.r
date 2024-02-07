test_that("disease_dispersal returns NULL when no dispersal is specified", {
  result <- disease_dispersal(replicates = 10,
                              time_steps = 100,
                              populations = c(100, 200, 300),
                              demographic_stochasticity = TRUE,
                              dispersal = NULL,
                              dispersal_type = "pooled",
                              stages = 3,
                              compartments = 4,
                              simulator = ModelSimulator$new())

  expect_null(result)
})

test_that("disease_dispersal returns NULL when dispersal is not a valid type", {
  result <- disease_dispersal(replicates = 10,
                              time_steps = 100,
                              populations = c(100, 200, 300),
                              demographic_stochasticity = TRUE,
                              dispersal = c(100, 200, 300),
                              dispersal_type = "stages",
                              stages = 3,
                              compartments = 4,
                              simulator = ModelSimulator$new())

  expect_null(result)
})

test_that("disease_dispersal returns NULL when no dispersal is specified", {
  result <- disease_dispersal(replicates = 10,
                              time_steps = 100,
                              populations = c(100, 200, 300),
                              demographic_stochasticity = TRUE,
                              dispersal = NULL,
                              dispersal_type = "pooled",
                              stages = 3,
                              compartments = 4,
                              simulator = ModelSimulator$new())

  expect_null(result)
})

test_that_cli("disease_dispersal throws an error for invalid dispersal length in 'pooled' dispersal type", {
  expect_snapshot(disease_dispersal(replicates = 10,
                                 time_steps = 100,
                                 populations = c(100, 200, 300),
                                 demographic_stochasticity = TRUE,
                                 dispersal = c(0.2, 0.3),
                                 dispersal_type = "pooled",
                                 stages = 3,
                                 compartments = 4,
                                 simulator = ModelSimulator$new()),
                  error = TRUE)
})

test_that_cli("disease_dispersal throws an error for invalid dispersal length in 'stages' dispersal type", {
  expect_snapshot(disease_dispersal(replicates = 10,
                                 time_steps = 100,
                                 populations = c(100, 200, 300),
                                 demographic_stochasticity = TRUE,
                                 dispersal = c(0.2, 0.3, 0.4),
                                 dispersal_type = "stages",
                                 stages = 4,
                                 compartments = 4,
                                 simulator = ModelSimulator$new()),
                  error = TRUE)
})

test_that_cli("disease_dispersal throws an error for invalid dispersal length in 'compartments' dispersal type", {
  expect_snapshot(disease_dispersal(replicates = 10,
                                 time_steps = 100,
                                 populations = c(100, 200, 300),
                                 demographic_stochasticity = TRUE,
                                 dispersal = c(0.2, 0.3, 0.4),
                                 dispersal_type = "compartments",
                                 stages = 3,
                                 compartments = 5,
                                 simulator = "simulator"),
                  error = TRUE)
})

test_that_cli("disease_dispersal throws an error for invalid dispersal length in 'segments' dispersal type", {
  expect_snapshot(disease_dispersal(replicates = 10,
                                 time_steps = 100,
                                 populations = c(100, 200, 300),
                                 demographic_stochasticity = TRUE,
                                 dispersal = c(0.2, 0.3, 0.4),
                                 dispersal_type = "segments",
                                 stages = 3,
                                 compartments = 4,
                                 simulator = ModelSimulator$new()),
                  error = TRUE)
})

test_that("user-defined dispersal function behaves as expected", {
  simulator <- SimulatorReference$new()
  # User-defined dispersal as a function
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(function(params)
      0.33),
    dispersal_type = "pooled",
    dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
    dispersal_target_k = 5,
    dispersal_target_n = list(threshold = 10, cutoff = 15),
    stages = 2,
    compartments = 4,
    simulator = simulator
  )
  expect_type(dispersal_function, "closure")
  expect_named(
    formals(dispersal_function),
    c(
      "r",
      "tm",
      "carrying_capacity",
      "segment_abundance",
      "occupied_indices"
    )
  )
  expect_named(
    environment(dispersal_function)[["params"]],
    c(
      "replicates",
      "time_steps",
      "populations",
      "stages",
      "compartments",
      "dispersal_type",
      "demographic_stochasticity",
      "dispersal_source_n_k",
      "dispersal_target_k",
      "dispersal_target_n",
      "dispersal_target_n_k",
      "simulator"
    )
  )
  expect_equal(
    environment(dispersal_function)[["params"]][c(
      "replicates",
      "time_steps",
      "populations",
      "dispersal_type",
      "demographic_stochasticity",
      "dispersal_source_n_k",
      "dispersal_target_k",
      "dispersal_target_n",
      "stages",
      "compartments"
    )],
    list(
      replicates = 4,
      time_steps = 10,
      populations = 7,
      dispersal_type = "pooled",
      demographic_stochasticity = TRUE,
      dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
      dispersal_target_k = 5,
      dispersal_target_n = list(threshold = 10, cutoff = 15),
      stages = 2,
      compartments = 4
    )
  )

  # User-defined dispersal as list of functions with length > 1
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(function(params)
      0.33, function(params)
      0.44),
    dispersal_type = "stages",
    dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
    dispersal_target_k = 5,
    dispersal_target_n = list(threshold = 10, cutoff = 15),
    stages = 2,
    compartments = 4,
    simulator = simulator
  )
  expect_type(dispersal_function, "closure")
})

test_that("user-defined dispersal calculations", {
  simulator <- SimulatorReference$new()
  test_function <- function(params) {
    # mock dispersal
    params$simulator$attached$params <- params # attach to reference object
    emigrants <- round(params$segment_abundance * 0.4)
    return(params$segment_abundance - emigrants + emigrants[, c(7, 1:6)])
  }
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(test_function),
    dispersal_type = "pooled",
    dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
    dispersal_target_k = 5,
    dispersal_target_n = list(threshold = 10, cutoff = 15),
    stages = 3,
    compartments = 1,
    simulator = simulator
  )
  carrying_capacity <- rep(10, 7)
  segment_abundance <- matrix(
    c(7, 13, 0, 26, 0, 39, 47,
      2,  0, 6,  8, 0, 12, 13,
      0,  3, 4,  6, 0,  9, 10),
    nrow = 3,
    ncol = 7,
    byrow = TRUE
  )
  occupied_indices <- (1:7)[-5]
  emigrants <- round(segment_abundance * 0.4)
  expected_segment_abundance <- segment_abundance - emigrants + emigrants[, c(7, 1:6)]
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 6,
      carrying_capacity,
      segment_abundance,
      occupied_indices
    ),
    expected_segment_abundance
  )
  expect_named(
    simulator$attached$params,
    c(
      "replicates",
      "time_steps",
      "populations",
      "stages",
      "compartments",
      "dispersal_type",
      "demographic_stochasticity",
      "dispersal_source_n_k",
      "dispersal_target_k",
      "dispersal_target_n",
      "dispersal_target_n_k",
      "simulator",
      "r",
      "tm",
      "carrying_capacity",
      "segment_abundance",
      "occupied_indices"
    )
  )
  expect_equal(
    simulator$attached$params[c("r",
                                "tm",
                                "carrying_capacity",
                                "segment_abundance",
                                "occupied_indices")],
    list(
      r = 2,
      tm = 6,
      carrying_capacity = carrying_capacity,
      segment_abundance = segment_abundance,
      occupied_indices = occupied_indices
    )
  )
  # Errors and warnings
  test_function = function (params)
    stop("test error")
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(test_function),
    dispersal_type = "pooled",
    dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
    dispersal_target_k = 5,
    dispersal_target_n = list(threshold = 10, cutoff = 15),
    stages = 3,
    compartments = 1,
    simulator = simulator
  )
  expect_error(
    dispersal_function(
      r = 2,
      tm = 6,
      carrying_capacity,
      segment_abundance,
      occupied_indices
    ),
    "Error produced within user-defined dispersal function"
  )
  test_function <- function(params)
    params$segment_abundance
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(test_function),
    dispersal_type = "pooled",
    dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
    dispersal_target_k = 5,
    dispersal_target_n = list(threshold = 10, cutoff = 15),
    stages = 3,
    compartments = 1,
    simulator = simulator
  )
  segment_abundance[1] <- NA
  expect_warning(
    dispersal_function(
      r = 2,
      tm = 6,
      carrying_capacity,
      segment_abundance,
      occupied_indices
    ),
    "Non-finite abundances returned by user-defined dispersal function"
  )
  segment_abundance[1] <- -1
  expect_warning(
    dispersal_function(
      r = 2,
      tm = 6,
      carrying_capacity,
      segment_abundance,
      occupied_indices
    ),
    "Negative abundances returned by user-defined dispersal function"
  )
  # Multiple functions

})
