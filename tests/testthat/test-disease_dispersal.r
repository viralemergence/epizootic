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
      "segment_abundance"
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
  emigrants <- round(segment_abundance * 0.4)
  expected_segment_abundance <- segment_abundance - emigrants + emigrants[, c(7, 1:6)]
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 6,
      carrying_capacity,
      segment_abundance
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
                                "segment_abundance")],
    list(
      r = 2,
      tm = 6,
      carrying_capacity = carrying_capacity,
      segment_abundance = segment_abundance
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
      segment_abundance
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
      segment_abundance
    ),
    "Non-finite abundances returned by user-defined dispersal function"
  )
  segment_abundance[1] <- -1
  expect_warning(
    dispersal_function(
      r = 2,
      tm = 6,
      carrying_capacity,
      segment_abundance
    ),
    "Negative abundances returned by user-defined dispersal function"
  )
  # Multiple functions
  test_function_list <- list(
    function1 = function(params)
    params$segment_abundance,
    function2 = function(params) {
      emigrants <- round(params$segment_abundance * 0.4)
      return(params$segment_abundance - emigrants + emigrants[c(7, 1:6)])
    },
    function3 = function(params) {
      params$segment_abundance[c(6:7, 1:5)]
    }
  )
  segment_abundance[1] <- 1
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = test_function_list,
    dispersal_type = "stages",
    dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
    dispersal_target_k = 5,
    dispersal_target_n = list(threshold = 10, cutoff = 15),
    stages = 3,
    compartments = 1,
    simulator = simulator
  )
  expect_type(dispersal_function, "closure")
  emigrants <- round(segment_abundance[2,] * 0.4)
  expected_segment_abundance <- matrix(
    c(segment_abundance[1,],
      segment_abundance[2,] - emigrants + emigrants[c(7, 1:6)],
      segment_abundance[3, c(6:7, 1:5)]
    ), byrow = T, nrow = 3
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 6,
      carrying_capacity,
      segment_abundance
    ),
    expected_segment_abundance
  )
})

test_that("setup default function", {
  simulator <- SimulatorReference$new()
  region <- Region$new(coordinates = array(c(1:4, 4:1), c(7, 2)))
  conductance_raster <- raster::stack(replicate(10,+(region$region_raster > 0)))
  conductance_raster[[2]][11] <- 0
  dispersal_friction = DispersalFriction$new(region = region,
                                             conductance = conductance_raster)
  dispersal_gen1 <- DispersalGenerator$new(
    dispersal_friction = dispersal_friction,
    dispersal_proportion = 0.6,
    dispersal_breadth = 110,
    dispersal_max_distance = 300,
    distance_scale = 1000,
    distance_classes = seq(100, 400, 20)
  )
  dispersal_gen2 <- DispersalGenerator$new(
    dispersal_friction = dispersal_friction,
    dispersal_proportion = 0.4,
    dispersal_breadth = 110,
    dispersal_max_distance = 400,
    distance_scale = 1000,
    distance_classes = seq(100, 400, 20)
  )
  dispersal_gen1$calculate_distance_data()
  dispersal_gen1$calculate_dispersals(type = "matrix")
  dispersal_gen1$calculate_dispersals(type = "data")
  dispersal_gen2$calculate_distance_data()
  dispersal_gen2$calculate_dispersals(type = "matrix")
  dispersal_gen2$calculate_dispersals(type = "data")
  dispersal_data1 <- dispersal_gen1$dispersal_data[[1]]
  dispersal_data2 <- dispersal_gen2$dispersal_data[[1]]
  expected_dispersal_data_changes1 <- dispersal_gen1$dispersal_data
  expected_dispersal_data_changes2 <- dispersal_gen2$dispersal_data
  expected_dispersal_data_changes1[[1]] <- dispersal_data1[NULL, ]
  expected_dispersal_data_changes2[[1]] <- dispersal_data2[NULL, ]
  expected_compact_rows1 <- max(dispersal_data1[c("emigrant_row", "immigrant_row")])
  expected_compact_rows2 <- max(dispersal_data2[c("emigrant_row", "immigrant_row")])
  expected_compact_matrix1 <- array(0, c(expected_compact_rows1, 7))
  expected_compact_matrix2 <- array(0, c(expected_compact_rows2, 7))
  expected_compact_matrix1[as.matrix(dispersal_data1[, c("emigrant_row", "source_pop")])] <-
    dispersal_data1$dispersal_rate
  expected_compact_matrix2[as.matrix(dispersal_data2[, c("emigrant_row", "source_pop")])] <-
    dispersal_data2$dispersal_rate
  expected_target_pop_map1 <- array(0, c(expected_compact_rows1, 7))
  expected_target_pop_map2 <- array(0, c(expected_compact_rows2, 7))
  expected_target_pop_map1[as.matrix(dispersal_data1[, c("emigrant_row", "source_pop")])] <-
    dispersal_data1$target_pop
  expected_target_pop_map2[as.matrix(dispersal_data2[, c("emigrant_row", "source_pop")])] <-
  dispersal_data2$target_pop
  compact_indices1 <- array(1:(expected_compact_rows1 * 7), c(expected_compact_rows1, 7))
  compact_indices2 <- array(1:(expected_compact_rows2 * 7), c(expected_compact_rows2, 7))

  expected_immigrant_map1 <- array(0, c(expected_compact_rows1, 7))
  expected_immigrant_map2 <- array(0, c(expected_compact_rows2, 7))

  expected_immigrant_map1[as.matrix(dispersal_data1[, c("emigrant_row", "source_pop")])] <-
    compact_indices1[as.matrix(dispersal_data1[, c("immigrant_row", "target_pop")])]

  expected_immigrant_map2[as.matrix(dispersal_data2[, c("emigrant_row", "source_pop")])] <-
    compact_indices2[as.matrix(dispersal_data2[, c("immigrant_row", "target_pop")])]
  # Via dispersal data with temporal changes
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(dispersal_gen1$dispersal_data,
                     dispersal_gen2$dispersal_data),
    dispersal_type = "stages",
    dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
    dispersal_target_k = 5,
    dispersal_target_n = list(threshold = 10, cutoff = 15),
    stages = 2,
    compartments = 1,
    simulator = simulator
  )
  expect_type(dispersal_function, "closure")
  expect_equal(
    environment(dispersal_function)[["dispersal_compact_rows_list"]],
    list(expected_compact_rows1, expected_compact_rows2)
  )
  expect_equal(
    environment(dispersal_function)[["dispersal_compact_matrix_list"]],
    list(expected_compact_matrix1, expected_compact_matrix2)
  )
  expect_true(environment(dispersal_function)[["dispersals_change_over_time"]])
  expect_true(environment(dispersal_function)[["dispersal_depends_on_source_pop_n_k"]])
  expect_true(environment(dispersal_function)[["dispersal_depends_on_target_pop_k"]])
  expect_true(environment(dispersal_function)[["dispersal_depends_on_target_pop_n"]])
  expect_equal(
    environment(dispersal_function)[["dispersal_target_pop_map_list"]],
    list(expected_target_pop_map1, expected_target_pop_map2)
  )
  expect_equal(
    environment(dispersal_function)[["dispersal_data_changes_list"]],
    list(expected_dispersal_data_changes1, expected_dispersal_data_changes2)
  )
  expect_equal(
    environment(dispersal_function)[["dispersal_immigrant_map_list"]],
    list(expected_immigrant_map1, expected_immigrant_map2)
  )
  # Via dispersal matrix (no temporal changes)
  dispersal_gen1$distance_data <- NULL
  dispersal_gen1$dispersal_friction = DispersalFriction$new(region = region)
  dispersal_gen1$calculate_distance_data()
  dispersal_gen1$calculate_dispersals(type = "matrix")
  dispersal_gen1$calculate_dispersals(type = "data")
  dispersal_data1 <- dispersal_gen1$dispersal_data[[1]]
  expected_compact_rows1 <- max(dispersal_data1[c("emigrant_row", "immigrant_row")])
  expected_compact_matrix1 <- array(0, c(expected_compact_rows1, 7))
  expected_compact_matrix1[as.matrix(dispersal_data1[, c("emigrant_row", "source_pop")])] <- dispersal_data1$dispersal_rate
  compact_indices1 <- array(1:(expected_compact_rows1 * 7), c(expected_compact_rows1, 7))
  expected_immigrant_map1 <- array(0, c(expected_compact_rows1, 7))
  expected_immigrant_map1[as.matrix(dispersal_data1[, c("emigrant_row", "source_pop")])] <- compact_indices1[as.matrix(dispersal_data1[, c("immigrant_row", "target_pop")])]

  dispersal_gen2$distance_data <- NULL
  dispersal_gen2$dispersal_friction = DispersalFriction$new(region = region)
  dispersal_gen2$calculate_distance_data()
  dispersal_gen2$calculate_dispersals(type = "matrix")
  dispersal_gen2$calculate_dispersals(type = "data")
  dispersal_data2 <- dispersal_gen2$dispersal_data[[1]]
  expected_compact_rows2 <- max(dispersal_data2[c("emigrant_row", "immigrant_row")])
  expected_compact_matrix2 <- array(0, c(expected_compact_rows2, 7))
  expected_compact_matrix2[as.matrix(dispersal_data2[, c("emigrant_row", "source_pop")])] <- dispersal_data2$dispersal_rate
  compact_indices2 <- array(1:(expected_compact_rows2 * 7), c(expected_compact_rows2, 7))
  expected_immigrant_map2 <- array(0, c(expected_compact_rows2, 7))
  expected_immigrant_map2[as.matrix(dispersal_data2[, c("emigrant_row", "source_pop")])] <- compact_indices2[as.matrix(dispersal_data2[, c("immigrant_row", "target_pop")])]

  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(dispersal_gen1$dispersal_matrix,
                     dispersal_gen2$dispersal_matrix),
    dispersal_type = "stages",
    stages = 2,
    compartments = 1,
    simulator = simulator
  )
  expect_type(dispersal_function, "closure")
  expect_equal(
    environment(dispersal_function)[["dispersal_compact_rows_list"]],
    list(expected_compact_rows1, expected_compact_rows2)
  )
  expect_equal(
    environment(dispersal_function)[["dispersal_compact_matrix_list"]],
    list(expected_compact_matrix1, expected_compact_matrix2)
  )
  expect_false(environment(dispersal_function)[["dispersals_change_over_time"]])
  expect_false(environment(dispersal_function)[["dispersal_depends_on_source_pop_n_k"]])
  expect_false(environment(dispersal_function)[["dispersal_depends_on_target_pop_k"]])
  expect_false(environment(dispersal_function)[["dispersal_depends_on_target_pop_n"]])
  expect_equal(
    environment(dispersal_function)[["dispersal_immigrant_map_list"]],
    list(expected_immigrant_map1, expected_immigrant_map2)
  )
  expect_warning(
    dispersal_function <- disease_dispersal(
      replicates = 4,
      time_steps = 10,
      populations = 7,
      demographic_stochasticity = TRUE,
      dispersal = list(dispersal_gen1$dispersal_matrix,
                      dispersal_gen2$dispersal_matrix),
      dispersal_type = "stages",
      stages = 2,
      compartments = 1,
      simulator = simulator,
      dispersal_source_n_k = list(cutoff = 1.5, threshold = -0.5),
      dispersal_target_k = NULL,
      dispersal_target_n = NULL
    ),
    "Dispersal density dependence for source N/K threshold must be greater than cutoff."
  )
  expect_false(environment(dispersal_function)[["dispersal_depends_on_source_pop_n_k"]])
})

test_that("default dispersal calculations", {
  simulator <- SimulatorReference$new()
  region <- Region$new(coordinates = array(c(1:4, 4:1), c(7, 2)))
  conductance_raster <- raster::stack(replicate(10,+(region$region_raster > 0)))
  conductance_raster[[2]][11] <- 0
  dispersal_friction = DispersalFriction$new(region = region,
                                              conductance = conductance_raster)
  dispersal_gen1 <- DispersalGenerator$new(
    dispersal_friction = dispersal_friction,
    dispersal_proportion = 0.6,
    dispersal_breadth = 110,
    dispersal_max_distance = 300,
    distance_scale = 1000,
    distance_classes = seq(100, 400, 20)
  )
  dispersal_gen2 <- DispersalGenerator$new(
    dispersal_friction = dispersal_friction,
    dispersal_proportion = 0.4,
    dispersal_breadth = 110,
    dispersal_max_distance = 400,
    distance_scale = 1000,
    distance_classes = seq(100, 400, 20)
  )
  dispersal_gen1$calculate_distance_data()
  dispersal_gen1$calculate_dispersals(type = "matrix")
  dispersal_gen1$calculate_dispersals(type = "data")
  dispersal_gen2$calculate_distance_data()
  dispersal_gen2$calculate_dispersals(type = "matrix")
  dispersal_gen2$calculate_dispersals(type = "data")

  dispersal_matrix_tm1 <- list()
  dispersal_matrix_tm1[[1]] <- dispersal_matrix_tm1[[2]] <- dispersal_gen1$dispersal_matrix
  dispersal_matrix_tm1[[3]] <- dispersal_matrix_tm1[[4]] <- dispersal_gen1$dispersal_matrix
  dispersal_matrix_tm1[[2]][as.matrix(dispersal_gen1$dispersal_data[[2]][c("target_pop", "source_pop")])] <- 0
  dispersal_matrix_tm1[[4]] <- dispersal_matrix_tm1[[4]] * matrix(
    1 / colSums(dispersal_matrix_tm1[[4]]),
    # high dispersal
    nrow = 7,
    ncol = 7,
    byrow = TRUE
  )
  dispersal_gen1$dispersal_data[[4]] <- dispersal_gen1$dispersal_data[[1]]
  dispersal_gen1$dispersal_data[[4]]$dispersal_rate <- dispersal_matrix_tm1[[4]][as.matrix(dispersal_gen1$dispersal_data[[4]][c("target_pop", "source_pop")])]

  dispersal_matrix_tm2 <- list()
  dispersal_matrix_tm2[[1]] <- dispersal_matrix_tm2[[2]] <- dispersal_gen2$dispersal_matrix
  dispersal_matrix_tm2[[3]] <- dispersal_matrix_tm2[[4]] <- dispersal_gen2$dispersal_matrix
  dispersal_matrix_tm2[[2]][as.matrix(dispersal_gen2$dispersal_data[[2]][c("target_pop", "source_pop")])] <- 0
  dispersal_matrix_tm2[[4]] <- dispersal_matrix_tm2[[4]] * matrix(
    1 / colSums(dispersal_matrix_tm2[[4]]),
    # high dispersal
    nrow = 7,
    ncol = 7,
    byrow = TRUE
  )
  dispersal_gen2$dispersal_data[[4]] <- dispersal_gen2$dispersal_data[[1]]
  dispersal_gen2$dispersal_data[[4]]$dispersal_rate <- dispersal_matrix_tm2[[4]][as.matrix(dispersal_gen2$dispersal_data[[4]][c("target_pop", "source_pop")])]
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(dispersal_gen1$dispersal_data,
                      dispersal_gen2$dispersal_data),
    dispersal_type = "stages",
    dispersal_target_k = 5,
    stages = 2,
    compartments = 1,
    simulator = simulator
  )
  carrying_capacity_tm_1 <- carrying_capacity_tm_2 <- carrying_capacity_tm_4 <- rep(10, 7)
  carrying_capacity_tm_3 <- c(10, 0, 3, 10, 10, 10, 10) # some target k < 5
  dispersal_matrix_tm1[[3]][2, ] <- 0
  dispersal_matrix_tm1[[3]][3, ] <- 3 / 5 * dispersal_matrix_tm1[[3]][3, ]
  dispersal_matrix_tm2[[3]][2, ] <- 0
  dispersal_matrix_tm2[[3]][3, ] <- 3 / 5 * dispersal_matrix_tm2[[3]][3, ]
  segment_abundance <- matrix(
    c(7, 13, 0, 26, 0, 39, 47,
      2,  0, 6,  8, 0, 12, 13),
    nrow = 2,
    ncol = 7,
    byrow = TRUE
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 1,
      carrying_capacity = carrying_capacity_tm_1,
      segment_abundance = segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 2,
      carrying_capacity = carrying_capacity_tm_2,
      segment_abundance = segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 3,
      carrying_capacity = carrying_capacity_tm_3,
      segment_abundance = segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 4,
      carrying_capacity = carrying_capacity_tm_4,
      segment_abundance = segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
})

test_that("density dependent dispersal", {
  simulator <- SimulatorReference$new()
  region <- Region$new(coordinates = array(c(1:4, 4:1), c(7, 2)))
  dispersal_gen1 <- DispersalGenerator$new(
    dispersal_proportion = 0.6,
    dispersal_breadth = 110,
    dispersal_max_distance = 300,
    distance_scale = 1000,
    distance_classes = seq(100, 400, 20),
    region = region
  )
  dispersal_gen2 <- DispersalGenerator$new(
    dispersal_proportion = 0.4,
    dispersal_breadth = 110,
    dispersal_max_distance = 400,
    distance_scale = 1000,
    distance_classes = seq(100, 400, 20),
    region = region
  )
  dispersal_gen1$calculate_distance_data()
  dispersal_gen1$calculate_dispersals(type = "matrix")
  dispersal_gen1$calculate_dispersals(type = "data")
  dispersal_gen2$calculate_distance_data()
  dispersal_gen2$calculate_dispersals(type = "matrix")
  dispersal_gen2$calculate_dispersals(type = "data")
  # Target abundance N threshold < cutoff (avoids overcrowded cells)
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(dispersal_gen1$dispersal_data,
                     dispersal_gen2$dispersal_data),
    dispersal_type = "compartments",
    dispersal_source_n_k = NULL,
    dispersal_target_k = c(5, 5, 5, 10, 15, 20, 25),
    dispersal_target_n = list(threshold = 10, cutoff = 15),
    stages = 1,
    compartments = 2,
    simulator = simulator
  )
  carrying_capacity <- c(10, 0, 3, 10, 10, 10, 10) # some target k < 5
  segment_abundance <- matrix(
    c(7, 13, 0, 26, 0, 39, 47,
      2,  0, 6,  8, 0, 12, 13),
    nrow = 2,
    ncol = 7,
    byrow = TRUE
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 1,
      carrying_capacity,
      segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
  # Set abundance N cutoff as array
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(dispersal_gen1$dispersal_data,
                     dispersal_gen2$dispersal_data),
    dispersal_type = "compartments",
    dispersal_source_n_k = NULL,
    dispersal_target_k = c(5, 5, 5, 10, 15, 20, 25),
    dispersal_target_n = list(
      threshold = 10,
      cutoff = c(15, 15, 15, 15, 20, 25, 30)
    ),
    stages = 1,
    compartments = 2,
    simulator = simulator
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 1,
      carrying_capacity,
      segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
  # Target abundance N threshold > cutoff (seeks company)
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(dispersal_gen1$dispersal_data,
                     dispersal_gen2$dispersal_data),
    dispersal_type = "compartments",
    dispersal_source_n_k = NULL,
    dispersal_target_k = c(5, 5, 5, 10, 15, 20, 25),
    dispersal_target_n = list(
      threshold = c(5, 10, 15, 20, 20, 20, 25),
      cutoff = c(15, 15, 15, 15, 20, 25, 30)
    ),
    stages = 1,
    compartments = 2,
    simulator = simulator
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 1,
      carrying_capacity,
      segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
  # Source abundance divide by carrying capacity N/K threshold > cutoff
  # (leaves over-exploited cells)
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(dispersal_gen1$dispersal_data,
                     dispersal_gen2$dispersal_data),
    dispersal_type = "compartments",
    dispersal_source_n_k = list(cutoff = -0.5, threshold = 1.5),
    stages = 1,
    compartments = 2,
    simulator = simulator
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 1,
      carrying_capacity,
      segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
  # Target abundance divide by carrying capacity N/K threshold < cutoff
  # (avoids over-exploited cells)
  dispersal_function <- disease_dispersal(
    replicates = 4,
    time_steps = 10,
    populations = 7,
    demographic_stochasticity = TRUE,
    dispersal = list(dispersal_gen1$dispersal_matrix,
                     dispersal_gen2$dispersal_matrix),
    dispersal_type = "compartments",
    dispersal_target_n_k = list(cutoff = 0.115, threshold = 0.679),
    stages = 1,
    compartments = 2,
    simulator = simulator
  )
  expect_equal(
    dispersal_function(
      r = 2,
      tm = 1,
      carrying_capacity,
      segment_abundance
    ) |> sum(),
    segment_abundance |> sum()
  )
})
