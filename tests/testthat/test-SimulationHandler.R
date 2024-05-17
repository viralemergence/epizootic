test_that("Runs successfully with valid inputs", {
  model_template <- DiseaseModel$new(
    time_steps = 5,
    seasons = 2,
    populations = 25,
    stages = 2,
    compartments = 4,
    coordinates = data.frame(
      x = rep(seq(177.01, 177.05, 0.01), 5),
      y = rep(seq(-18.01, -18.05, -0.01), each = 5)
    ),
    initial_abundance = c(
      c(5000, 5000, 0, 1, 0, 0, 0, 0),
      rep(c(5000, 5000, 0, 0, 0, 0, 0, 0), 24)
    ) |> matrix(nrow = 8),
    breeding_season_length = rep(100, 25),
    mortality = c(0.4, 0, 0.505, 0.105, 0.4, 0, 0.45, 0.05),
    mortality_unit = 1,
    fecundity_mask = c(0, 1, 0, 1, 0, 1, 0, 1),
    transmission = c(0.00002, 0.00001, 7.84e-06, 3.92e-06),
    transmission_unit = 0,
    transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0),
    recovery = c(0.05714286, 0.05714286, 0.1, 0.1),
    recovery_unit = rep(0, 8),
    recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1),
    season_functions = list(siri_model_summer, siri_model_winter),
    simulation_order = c("transition", "season_functions", "results"),
    results_selection = "abundance",
    results_breakdown = "stages",
    dispersal_type = "stages",
    attribute_aliases = list(dispersal1 = "dispersal$a",
                             dispersal2 = "dispersal$b"),
    verbose = FALSE
  )

  model_simulator <- ModelSimulator$new(simulation_model = model_template,
                                        simulation_function = disease_simulator)

  results_dir <- tempdir()

  generator3 <- Generator$new(
    description = "Test generator",
    decimals = 0,
    inputs = "seed_number",
    outputs = "carrying_capacity"
  )

  generator3$add_generative_requirements(carrying_capacity = "function")

  generator3$add_function_template(
    "carrying_capacity",
    function_def = function(params) {
      matrix(params$seed_number, ncol = 5, nrow = 25)
    },
    call_params = c("seed_number")
  )

  distance_matrix <- geosphere::distm(model_template$coordinates,
                                       model_template$coordinates,
                                       fun = geosphere::distGeo)/1000
  dispersal_gen1 <- DispersalGenerator$new(coordinates = model_template$coordinates,
                                           distance_classes = seq(100, 600, 20))
  dispersal_gen1$calculate_distance_data(distance_matrix = distance_matrix)
  dispersal_gen1$set_attributes(proportion = 0.4, breadth = 110, max_distance = 300)
  dispersal_gen2 <- DispersalGenerator$new(coordinates = model_template$coordinates,
                                           distance_classes = seq(100, 600, 20))
  dispersal_gen2$calculate_distance_data(distance_matrix = distance_matrix)
  dispersal_gen2$set_attributes(proportion = 0.2, breadth = 110, max_distance = 500)

  sample_data <- data.frame(seed_number = 100000, fecundity = 15, fecundity_unit = 1)

  # Check with valid inputs
  expect_silent(SimulationHandler$new(sample_data = sample_data,
                model_template = model_template,
                model_simulator = model_simulator,
                generators = list(dispersal_gen1, dispersal_gen2, generator3),
                parallel_cores = 1, results_dir = test_path("test_results")))
  handler <- SimulationHandler$new(sample_data = sample_data,
                model_template = model_template,
                model_simulator = model_simulator,
                generators = list(dispersal_gen1, dispersal_gen2, generator3),
                parallel_cores = 1, results_dir = test_path("test_results"))
  run_output <- handler$run()
  expect_named(run_output, c("summary", "failed_indices", "warning_indices", "full_log"))
  expect_equal(run_output$summary, "1 of 1 sample models ran and saved results successfully")
  expect_length(run_output$failed_indices, 0)
  expect_equal(run_output$full_log, list(list(successful = TRUE,
                                              message = "Model sample 1 simulation ran successfully and the results were saved")))
  expect_true(all(c("sample_1_results.qs", "simulation_log.txt") %in% list.files(handler$results_dir)))

  # Check if it works when the generator input is in the model template
  sample_data <- data.frame(fecundity = 15, fecundity_unit = 1)
  model_template$set_attributes(seed_number = 100000)
  expect_silent(SimulationHandler$new(sample_data = sample_data,
                model_template = model_template,
                model_simulator = model_simulator,
                generators = list(dispersal_gen1, dispersal_gen2, generator3),
                parallel_cores = 1, results_dir = test_path("test_results")))
  handler <- SimulationHandler$new(sample_data = sample_data,
                model_template = model_template,
                model_simulator = model_simulator,
                generators = list(dispersal_gen1, dispersal_gen2, generator3),
                parallel_cores = 1, results_dir = test_path("test_results"))
  expect_silent(handler$run())
  run_output <- handler$run()
  expect_named(run_output, c("summary", "failed_indices", "warning_indices", "full_log"))
  expect_equal(run_output$summary, "1 of 1 sample models ran and saved results successfully")
  expect_length(run_output$failed_indices, 0)
  expect_equal(run_output$full_log, list(list(successful = TRUE,
                                              message = "Model sample 1 simulation ran successfully and the results were saved")))
  expect_true(all(c("sample_1_results.qs", "simulation_log.txt") %in% list.files(handler$results_dir)))
})

test_that("initialization and parameter setting", {
  # Default initialization
  sim_manager <- SimulationHandler$new()
  expect_true("ModelSimulator" %in% class(sim_manager$model_simulator))
  expect_null(sim_manager$model_simulator$simulation_function)
  sim_manager$model_template <- DiseaseModel$new()
  expect_true(is.function(sim_manager$model_simulator$simulation_function))
  sim_manager <- SimulationHandler$new(model_template = DiseaseModel$new())
  expect_true(is.function(sim_manager$model_simulator$simulation_function))
  # Invalid attributes
  expect_error(sim_manager$model_template <- "dummy")
  expect_error(sim_manager$nested_model <- "dummy")
  expect_error(sim_manager$model_simulator <- "dummy")
  expect_silent(sim_manager$model_template <- NULL)
  expect_silent(sim_manager$nested_model <- NULL)
  expect_silent(sim_manager$model_simulator <- NULL)
})

test_that("attempt run with incomplete attributes", {
  TEST_DIRECTORY <- test_path("test_results")
  sim_manager <- SimulationHandler$new()
  # No model template and samples
  sim_manager$model_simulator <- NULL
  expect_error(sim_manager$run())
  sim_manager$model_template <- SimulationModel$new(model_attributes = c("time_steps", "attr1", "attr2"))
  sim_manager$sample_data <- data.frame() # empty
  expect_error(sim_manager$run())
  sim_manager$sample_data <- data.frame(attr1 = 3:4, attr2 = 5:6)
  # No model simulator
  expect_error(sim_manager$run(), "The model simulator has not been set")
  sim_manager$model_simulator <- ModelSimulator$new()
  sim_manager$model_simulator$simulation_function <- NULL
  expect_error(sim_manager$run(), "The model simulator function has not been set")
  sim_manager$model_simulator$simulation_function <- "max"
  # No results output directory
  expect_error(sim_manager$run(), "No output directory set for results")
  sim_manager$results_dir <- "invalid_dir"
  expect_error(sim_manager$run(), "Could not find results directory invalid_dir")
  sim_manager$results_dir <- TEST_DIRECTORY
  # With incomplete model
  expect_null(sim_manager$nested_model)
  generator <- Generator$new(generative_requirements = list(attr3 = "function"),
                             inputs = c("attr2"), outputs = c("attr3"))
  generator$function_templates <- list(attr3 = list(function_def = function(params) return(params$a + 2),
                                                    call_params = c("attr2")))
  sim_manager$generators <- list(gen3 = generator)
  expect_error(sim_manager$run(), "Model attributes are incomplete: time_steps")
  expect_true("SimulationModel" %in% class(sim_manager$nested_model))
  expect_equal(sim_manager$nested_model$attached, list(sample_model_names = c("attr1", "attr2"),
                                                       sample_generative_names = list("attr3")))
  # Set model sample
  model_clone <- sim_manager$nested_model$clone()
  sim_manager$set_model_sample(model_clone, 1)
  expect_equal(model_clone$get_attributes(), list(attr1 = 3, attr2 = 5, sample_model_names = c("attr1", "attr2"),
                                                  sample_generative_names = list("attr3"), attr3 = 7))
})
