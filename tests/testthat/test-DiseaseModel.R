test_that("DiseaseModel initializes correctly", {
  model <- DiseaseModel$new()
  expect_true(inherits(model, "DiseaseModel"))
  expect_true(inherits(model, "SimulationModel"))
  expect_equal(model$seasons, 1)
  expect_equal(model$results_breakdown, "segments")
  expect_equal(model$stages, 1)
  expect_equal(model$compartments, 1)
})

test_that("DiseaseModel sets and gets attributes correctly", {
  model <- DiseaseModel$new()
  model$populations <- 100
  expect_equal(model$populations, 100)
  
  model$initial_abundance <- matrix(1:100, nrow = 10)
  expect_equal(model$initial_abundance, matrix(1:100, nrow = 10))
  
  model$demographic_stochasticity <- FALSE
  expect_false(model$demographic_stochasticity)
  
  model$standard_deviation <- 0.1
  expect_equal(model$standard_deviation, 0.1)
  
  model$correlation <- list(type = "spatial", value = 0.5)
  expect_equal(model$correlation, list(type = "spatial", value = 0.5))
  
  model$stages <- 3
  expect_equal(model$stages, 3)
  
  model$compartments <- 2
  expect_equal(model$compartments, 2)
  
  model$results_breakdown <- "pooled"
  expect_equal(model$results_breakdown, "pooled")
  
  model$carrying_capacity <- matrix(1:100, nrow = 10)
  expect_equal(model$carrying_capacity, matrix(1:100, nrow = 10))
  
  model$density_dependence <- "logistic"
  expect_equal(model$density_dependence, "logistic")
  
  model$growth_rate_max <- 1.2
  expect_equal(model$growth_rate_max, 1.2)
  
  model$fecundity <- c(0.5, 0.6, 0.7)
  expect_equal(model$fecundity, c(0.5, 0.6, 0.7))
  
  model$density_stages <- c(1, 0.5, 0.2)
  expect_equal(model$density_stages, c(1, 0.5, 0.2))

  model$dispersal_source_n_k <- list(cutoff = 0.5, threshold = 0.2)
  expect_equal(model$dispersal_source_n_k, list(cutoff = 0.5, threshold = 0.2))
  
  model$dispersal_target_k <- 0.8
  expect_equal(model$dispersal_target_k, 0.8)
  
  model$dispersal_target_n <- list(threshold = 0.3, cutoff = 0.1)
  expect_equal(model$dispersal_target_n, list(threshold = 0.3, cutoff = 0.1))
  
  model$dispersal_target_n_k <- list(threshold = 0.4, cutoff = 0.2)
  expect_equal(model$dispersal_target_n_k, list(threshold = 0.4, cutoff = 0.2))
  
  model$abundance_threshold <- 10
  expect_equal(model$abundance_threshold, 10)
  
  model$seasons <- 4
  expect_equal(model$seasons, 4)
  
  model$results_selection <- c("abundance", "occupancy")
  expect_equal(model$results_selection, c("abundance", "occupancy"))
})

test_that("DiseaseModel handles attribute aliases correctly", {
  model <- DiseaseModel$new(attribute_aliases = list(dispersal_data = "dispersal"))
  model$set_attributes(dispersal_data = list(strategy = "random", rate = 0.1))
  expect_equal(model$dispersal, list(strategy = "random", rate = 0.1))
})

# Test set_sample_attributes method
test_that("DiseaseModel set_sample_attributes works correctly", {
  model <- DiseaseModel$new()
  
  model$set_sample_attributes(params = list(populations = 100, seasons = 2))
  expect_equal(model$populations, 100)
  expect_equal(model$seasons, 2)
})