test_that_cli("minimal inputs", {
  expect_snapshot_error(disease_simulator(list()))
  inputs <- list(random_seed = 42, replicates = 5, time_steps = 100, populations = 5,
                 initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"))
  expect_silent({disease_simulator(inputs)})
})

# test_that("Replicates, seasons, stages and compartments take on default values
# when not provided in inputs", {
# inputs <- list(time_steps = 100, populations = 5,
#                initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
#                mortality = 0.1, fecundity = 0.11, transmission = 0,
#                simulation_order = c("transition"))
#   result <- disease_simulator(inputs)
#   expect_equal(result[["replicates"]], 1)
#   expect_equal(result[["seasons"]], 2)
#   expect_equal(result[["stages"]], 1)
#   expect_equal(result[["compartments"]], 1)
# })

test_that_cli("check validity of initial abundance", {
  inputs <- list(replicates = 5, time_steps = 100, populations = 5,
                 initial_abundance = "boop", carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"), stages = 2)
  expect_snapshot_error(disease_simulator(inputs))
  inputs <- list(replicates = 5, time_steps = 100, populations = 5,
                 initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"), stages = 2)
  expect_snapshot_error(disease_simulator(inputs))
  inputs <- list(replicates = 5, time_steps = 100, populations = 6,
                 initial_abundance = rep(10, 5), carrying_capacity = rep(20, 5),
                 mortality = 0.1, fecundity = 0.11, transmission = 0,
                 simulation_order = c("transition"))
  expect_snapshot_error(disease_simulator(inputs))
})
