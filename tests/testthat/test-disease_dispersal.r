library(testthat)

# Test disease_dispersal function
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

test_that("disease_dispersal throws an error for invalid dispersal length in 'pooled' dispersal type", {
  expect_snapshot_error(disease_dispersal(replicates = 10,
                                 time_steps = 100,
                                 populations = c(100, 200, 300),
                                 demographic_stochasticity = TRUE,
                                 dispersal = c(0.2, 0.3),
                                 dispersal_type = "pooled",
                                 stages = 3,
                                 compartments = 4,
                                 simulator = ModelSimulator$new()))
})

test_that("disease_dispersal throws an error for invalid dispersal length in 'stages' dispersal type", {
  expect_snapshot_error(disease_dispersal(replicates = 10,
                                 time_steps = 100,
                                 populations = c(100, 200, 300),
                                 demographic_stochasticity = TRUE,
                                 dispersal = c(0.2, 0.3, 0.4),
                                 dispersal_type = "stages",
                                 stages = 4,
                                 compartments = 4,
                                 simulator = ModelSimulator$new()))
})

test_that("disease_dispersal throws an error for invalid dispersal length in 'compartments' dispersal type", {
  expect_snapshot_error(disease_dispersal(replicates = 10,
                                 time_steps = 100,
                                 populations = c(100, 200, 300),
                                 demographic_stochasticity = TRUE,
                                 dispersal = c(0.2, 0.3, 0.4),
                                 dispersal_type = "compartments",
                                 stages = 3,
                                 compartments = 5,
                                 simulator = "simulator"))
})

test_that("disease_dispersal throws an error for invalid dispersal length in 'segments' dispersal type", {
  expect_error(disease_dispersal(replicates = 10,
                                 time_steps = 100,
                                 populations = c(100, 200, 300),
                                 demographic_stochasticity = TRUE,
                                 dispersal = c(0.2, 0.3, 0.4),
                                 dispersal_type = "segments",
                                 stages = 3,
                                 compartments = 4,
                                 simulator = ModelSimulator$new()))
})

test_that()