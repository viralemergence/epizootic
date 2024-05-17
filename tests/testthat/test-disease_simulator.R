test_that('disease_simulator works with valid inputs', {
  inputs <- list(
    time_steps = 5,
    seasons = 2,
    populations = 25,
    stages = 2,
    compartments = 4,
    coordinates = data.frame(x = rep(seq(177.01, 177.05, 0.01), 5),
                             y = rep(seq(-18.01, -18.05, -0.01), each = 5)),
    initial_abundance = c(c(5000, 5000, 0, 1, 0, 0, 0, 0),
                          rep(c(5000, 5000, 0, 0, 0, 0, 0, 0), 24)) |>
      matrix(nrow = 8),
    carrying_capacity = matrix(100000, nrow = 25, ncol = 5),
    breeding_season_length = rep(100, 25),
    mortality = c(0.4, 0, 0.505, 0.105, 0.4, 0, 0.45, 0.05),
    mortality_unit = 1,
    fecundity = 15,
    fecundity_unit = 1,
    fecundity_mask = c(0, 1, 0, 1, 0, 1, 0, 1),
    transmission = c(0.00002, 0.00001, 7.84e-06, 3.92e-06),
    transmission_unit = 0,
    transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0),
    recovery = c(0.05714286, 0.05714286, 0.1, 0.1),
    recovery_unit = rep(0, 8),
    recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1),
    season_functions = list(siri_model_summer, siri_model_winter),
    simulation_order = c("transition", "season_functions", "results")
  )
  expect_silent(disease_simulator(inputs))
})

test_that("disease_simulator works with minimal inputs", {
  expect_error(
    disease_simulator(list())
  )
  inputs <- list(
    time_steps = 10,
    populations = 1,
    initial_abundance = 10,
    carrying_capacity = 40,
    mortality = 0.1,
    fecundity = 0.1,
    transmission = 0.1,
    simulation_order = c("transition", "results")
  )
  expect_silent(disease_simulator(inputs))
})
