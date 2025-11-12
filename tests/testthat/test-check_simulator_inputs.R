test_that_cli("minimal inputs", {
  expect_snapshot_error(check_simulator_inputs(list()))
  inputs <- list(
    random_seed = 42,
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = rep(10, 5),
    carrying_capacity = rep(20, 5),
    mortality = 0.1,
    fecundity = 0.11,
    transmission = 0,
    simulation_order = c("transition")
  )
  expect_silent({check_simulator_inputs(inputs)})
})

test_that("Replicates, seasons, stages and compartments take on default values
when not provided in inputs", {
inputs <- list(time_steps = 100, populations = 5,
               initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
               mortality = 0.1, fecundity = 0.11, transmission = 0,
               simulation_order = c("transition"))
  result <- check_simulator_inputs(inputs)
  expect_identical(result[["replicates"]], 1)
  expect_identical(result[["seasons"]], 2)
  expect_identical(result[["stages"]], 1)
  expect_identical(result[["compartments"]], 1)
})

test_that_cli("check error handling for initial abundance", {
  inputs <- list(replicates = 5, time_steps = 100, populations = 5,
                 initial_abundance = "boop", carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"), stages = 2)
  expect_snapshot_error(check_simulator_inputs(inputs))
  inputs <- list(replicates = 5, time_steps = 100, populations = 5,
                 initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"), stages = 2)
  expect_snapshot_error(check_simulator_inputs(inputs))
  inputs <- list(replicates = 5, time_steps = 100, populations = 6,
                 initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"))
  expect_snapshot_error(check_simulator_inputs(inputs))
})

test_that("vector and raster input for initial abundance are converted", {
  inputs <- list(random_seed = 42, replicates = 5, time_steps = 100,
                 populations = 5, initial_abundance = rep(10, 5),
                 carrying_capacity = rep(20, 5), mortality = 0.1,
                 fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"))
  result <- check_simulator_inputs(inputs)
  expect_contains(class(result[["initial_abundance"]]), "matrix")
  expect_identical(dim(result[["initial_abundance"]]), c(1L, 5L))
})

test_that_cli("check error handling for breeding season length", {
  inputs <- list(replicates = 5, time_steps = 100, populations = 5,
                 initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"), seasons = 2,
                 breeding_season_length = c(50, 60, 70, 80, NA))
  expect_snapshot_error(check_simulator_inputs(inputs))
  inputs <- list(replicates = 5, time_steps = 100, populations = 5,
                 initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"), seasons = 2,
                 season_lengths = c(300, 100, 65))
  expect_snapshot_error(check_simulator_inputs(inputs))
})

test_that_cli("check error handling for mortality", {
  inputs <- list(time_steps = 100, populations = 5, stages = 2,
                 initial_abundance = matrix(rep(10, 5*2), nrow = 2),
                 carrying_capacity = rep(20, 5),
                 fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"), seasons = 2,
                 mortality = list(c(0.01, 0.02), c(0.01)))
  expect_snapshot_error(check_simulator_inputs(inputs))
  region <- Region$new(coordinates = array(c(1:4, 4:1), c(7, 2)))
  inputs <- list(time_steps = 100, populations = 5, stages = 2,
                 initial_abundance = matrix(rep(10, 5*2), nrow = 2),
                 mortality = region$region_raster,
                 carrying_capacity = rep(20, 5),
                 fecundity = rep(0.11, 2), transmission = 0,
                 simulation_order = c("transition"), seasons = 2)
  expect_snapshot_error(check_simulator_inputs(inputs))
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*2*3), nrow = 2*3),
    carrying_capacity = rep(20, 5),
    fecundity = rep(0.11, 6),
    transmission = 0,
    simulation_order = c("transition"),
    mortality = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)),
    mortality_unit = list(c(1, 1, 1), c(0, 0)),
    seasons = 2,
    stages = 3,
    compartments = 2
  )
  expect_snapshot_error(check_simulator_inputs(inputs))
})

test_that("Valid mortality list input", {
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3), nrow = 3),
    carrying_capacity = rep(20, 5),
    fecundity = rep(0.11, 3),
    transmission = rep(0, 3),
    simulation_order = c("transition"),
    mortality = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)),
    seasons = 2,
    stages = 3,
    compartments = 1
  )
  result <- check_simulator_inputs(inputs)
  expect_contains(class(result[["mortality"]]), "list")
  expect_length(result[["mortality"]], 2)
  expect_identical(lengths(result[["mortality"]]), c(3L, 3L))
})

test_that("Named lists get sorted appropriately", {
    inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3), nrow = 3),
    carrying_capacity = rep(20, 5),
    fecundity = rep(0.11, 3),
    transmission = rep(0, 3),
    simulation_order = c("transition"),
    mortality = list(b = c(c = 0.1, b = 0.2, a = 0.3), a = c(0.4, 0.5, 0.6)),
    seasons = 2,
    stages = 3,
    compartments = 1,
    verbose = F
  )
  expected_list <- list(a = c(0.4, 0.5, 0.6), b = c(a = 0.3, b = 0.2, c = 0.1))
  result <- check_simulator_inputs(inputs)
  expect_equal(result[["mortality"]], expected_list)
})

test_that("Valid mortality vector input", {
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*2*3), nrow = 2*3),
    carrying_capacity = rep(20, 5),
    fecundity = rep(0.11, 6),
    transmission = rep(0, 6),
    simulation_order = c("transition"),
    mortality = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    seasons = 2,
    stages = 3,
    compartments = 2
  )
  result <- check_simulator_inputs(inputs)
  expect_contains(class(result[["mortality"]]), "list")
  expect_length(result[["mortality"]], 2)
  expect_identical(lengths(result[["mortality"]]), c(6L, 6L))
})

test_that("Valid 'mortality_unit' input", {
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3), nrow = 3),
    carrying_capacity = rep(20, 5),
    fecundity = rep(0.11, 3),
    transmission = rep(0, 3),
    simulation_order = c("transition"),
    mortality = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)),
    mortality_unit = list(c(1, 1, 1), c(0, 0, 0)),
    seasons = 2,
    stages = 3,
    compartments = 1
  )
  result <- check_simulator_inputs(inputs)
  expect_contains(class(result[["mortality_unit"]]), "list")
  expect_length(result[["mortality_unit"]], 2)
  expect_identical(lengths(result[["mortality_unit"]]), c(3L, 3L))
})

test_that("Valid fecundity list input", {
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3), nrow = 3),
    carrying_capacity = rep(20, 5),
    transmission = rep(0, 3),
    simulation_order = c("transition"),
    mortality = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)),
    fecundity = list(c(0, 1, 0), c(1, 0, 1)),
    seasons = 2,
    stages = 3,
    compartments = 1
  )
  result <- check_simulator_inputs(inputs)
  expect_contains(class(result[["fecundity"]]), "list")
  expect_length(result[["fecundity"]], 2)
  expect_identical(lengths(result[["fecundity"]]), c(3L, 3L))
})

test_that("Valid fecundity vector input", {
  inputs <- list(
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3*2), nrow = 6),
    carrying_capacity = rep(20, 5),
    transmission = rep(0, 6),
    simulation_order = c("transition"),
    mortality = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    fecundity = c(0, 1, 0, 1, 0, 1),
    seasons = 2,
    stages = 3,
    compartments = 2
  )
  result <- check_simulator_inputs(inputs)
  expect_contains(class(result[["fecundity"]]), "list")
  expect_length(result[["fecundity"]], 2)
  expect_identical(lengths(result[["fecundity"]]), c(6L, 6L))
})

test_that_cli("check error handling for fecundity", {
  inputs <- list(time_steps = 100, populations = 5, stages = 2,
                 initial_abundance = matrix(rep(10, 5*2), nrow = 2),
                 carrying_capacity = rep(20, 5),
                 mortality = c(0.11, 0.22), transmission = 0,
                 simulation_order = c("transition"), seasons = 2,
                 fecundity = list(c(0.01, 0.02), c(0.01)))
  expect_snapshot_error(check_simulator_inputs(inputs))
  region <- Region$new(coordinates = array(c(1:4, 4:1), c(7, 2)))
  inputs <- list(time_steps = 100, populations = 5, stages = 2,
                 initial_abundance = matrix(rep(10, 5*2), nrow = 2),
                 fecundity = region$region_raster,
                 carrying_capacity = rep(20, 5),
                 mortality = rep(0.11, 2), transmission = 0,
                 simulation_order = c("transition"), seasons = 2)
  expect_snapshot_error(check_simulator_inputs(inputs))
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3), nrow = 3),
    carrying_capacity = rep(20, 5),
    mortality = rep(0.11, 3),
    transmission = 0,
    simulation_order = c("transition"),
    fecundity = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)),
    fecundity_unit = list(c(1, 1, 1), c(0, 0)),
    seasons = 2,
    stages = 3,
    compartments = 1
  )
  expect_snapshot_error(check_simulator_inputs(inputs))
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3), nrow = 3),
    carrying_capacity = rep(20, 5),
    mortality = rep(0.11, 3),
    transmission = 0,
    simulation_order = c("transition"),
    fecundity = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)),
    fecundity_unit = list(c(1, 1, 1), c(0, 0, 2)),
    seasons = 2,
    stages = 3,
    compartments = 1
  )
  expect_snapshot_error(check_simulator_inputs(inputs))
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3), nrow = 3),
    carrying_capacity = rep(20, 5),
    mortality = rep(0.11, 3),
    transmission = rep(0, 3),
    simulation_order = c("transition"),
    fecundity = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)),
    fecundity_unit = list(c(1, 1, 1), c(0, 0, 1)),
    fecundity_mask = c(0, 1, 1, 1),
    seasons = 2,
    stages = 3,
    compartments = 1
  )
  expect_snapshot_error(check_simulator_inputs(inputs))
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3), nrow = 3),
    carrying_capacity = rep(20, 5),
    mortality = rep(0.11, 3),
    transmission = 0,
    simulation_order = c("transition"),
    fecundity = list(c(0.1, 0.2, 0.3), c(0.4, 0.5, 0.6)),
    fecundity_unit = list(c(1, 1, 1), c(0, 0, 1)),
    fecundity_mask = c(0, 1, 2),
    seasons = 2,
    stages = 3,
    compartments = 1
  )
  expect_snapshot_error(check_simulator_inputs(inputs))
})

test_that("Default transmission_mask assignment when NULL", {
  inputs <- list(
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3*2), nrow = 6),
    carrying_capacity = rep(20, 5),
    simulation_order = c("transition"),
    mortality = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    fecundity = c(0, 1, 0, 1, 0, 1),
    transmission = list(c(0, 1, 0, 0, 0, 0), c(1, 0, 1, 1, 1, 1)),
    seasons = 2,
    stages = 3,
    compartments = 2
  )

  result <- check_simulator_inputs(inputs)

  # Check that transmission_mask is assigned with default values
  expect_contains(class(result[["transmission_mask"]]), "list")
  expect_length(result[["transmission_mask"]], 2)  # Number of seasons
  # Check the first season's transmission_mask
  expect_identical(result[["transmission_mask"]][[1]],
               c(1, 1, 1, 0, 0, 0))
  # Check the second season's transmission_mask
  expect_identical(result[["transmission_mask"]][[2]],
               c(1, 1, 1, 0, 0, 0))
})

test_that("Default recovery_mask assignment when NULL", {
  inputs <- list(
    time_steps = 100,
    populations = 5,
    initial_abundance = matrix(rep(10, 5*3*2), nrow = 6),
    carrying_capacity = rep(20, 5),
    simulation_order = c("transition"),
    mortality = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
    fecundity = c(0, 1, 0, 1, 0, 1),
    transmission = list(c(0, 1, 0, 0, 0, 0), c(1, 0, 1, 1, 1, 1)),
    recovery = list(
      c(0, 1, 0, 1, 0, 1),
      c(1, 0, 1, 0, 1, 0)
    ),
    seasons = 2,
    stages = 2,
    compartments = 3
  )

  result <- check_simulator_inputs(inputs)

  # Check that recovery_mask is assigned with default values
  expect_contains(class(result[["recovery_mask"]]), "list")
  expect_length(result[["recovery_mask"]], 2)  # Number of seasons
  # Check the first season's recovery_mask
  expect_identical(result[["recovery_mask"]][[1]],
               c(0, 0, 1, 1, 0, 0))
  # Check the second season's recovery_mask
  expect_identical(result[["recovery_mask"]][[2]],
               c(0, 0, 1, 1, 0, 0))
})

test_that("Default dispersal_source_n_k assignment when NULL", {
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = rep(10, 5),
    carrying_capacity = rep(20, 5),
    mortality = 0.1,
    fecundity = 0.11,
    transmission = 0,
    simulation_order = c("transition"),
    dispersal_source_n_k = NULL,
    dispersal_source_n_k_cutoff = 0.1,
    dispersal_source_n_k_threshold = 0.2
  )

  result <- check_simulator_inputs(inputs)

  expect_contains(class(result[["dispersal_source_n_k"]]), "list")
  expect_identical(result[["dispersal_source_n_k"]][["cutoff"]], 0.1)
  expect_identical(result[["dispersal_source_n_k"]][["threshold"]], 0.2)
})

test_that("Replicate simulation_order when a vector", {
  inputs <- list(
    replicates = 5,
    time_steps = 100,
    populations = 5,
    initial_abundance = rep(10, 5),
    carrying_capacity = rep(20, 5),
    mortality = 0.1,
    fecundity = 0.11,
    transmission = 0,
    simulation_order = c("season1", "season2"),
    seasons = 2
  )

  result <- check_simulator_inputs(inputs)

  expect_contains(class(result[["simulation_order"]]), "list")
  expect_length(result[["simulation_order"]], 2)  # Number of seasons
  expect_identical(result[["simulation_order"]][[1]],
                   c("season1", "season2", "results"))
  expect_identical(result[["simulation_order"]][[2]],
                   c("season1", "season2", "results"))
})
