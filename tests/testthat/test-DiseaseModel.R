test_that("active get and set", {
  region = Region$new(coordinates = array(c(1:4, 4:1), c(7, 2)))
  disease_model <- DiseaseModel$new(time_steps = 10)
  expect_null(disease_model$populations)
  disease_model$initial_abundance <- seq(10, 60, by = 10)
  expect_equal(disease_model$populations, 6)
  disease_model$region <- region
  expect_equal(disease_model$populations, 7)
  # expect_equal(disease_model$get_attributes(
  #   c(
  #     "stages",
  #     "stage_matrix",
  #     "fecundity_mask",
  #     "density_affects",
  #     "density_stages",
  #     "dispersal_stages",
  #     "result_stages"
  #   )
  # ),
  # list(density_affects = "all"))
  # disease_model$stage_matrix <- array(c(0, 0.5, 0, 3, 0, 0.7, 4, 0, 0.8), c(3, 3))
  # expect_equal(disease_model$stages, 3)
  # expect_equal(disease_model$fecundity_mask, array(c(0, 0, 0, 1, 0, 0, 1, 1, 0), c(3, 3)))
  # expect_equal(disease_model$density_affects, array(c(0, 1, 0, 1, 0, 1, 1, 0, 1), c(3, 3)))
  # expect_equal(
  #   disease_model$get_attributes(c(
  #     "density_stages", "dispersal_stages", "result_stages"
  #   )),
  #   list(
  #     density_stages = array(1, 3),
  #     dispersal_stages = array(1, 3),
  #     result_stages = array(1, 3)
  #   )
  # )
  # expect_equal(
  #   disease_model$get_attribute_aliases(
  #     c(
  #       "dispersal_source_n_k",
  #       "dispersal_target_k",
  #       "dispersal_target_n"
  #     )
  #   ),
  #   c(
  #     "dispersal_source_n_k",
  #     "dispersal_target_k",
  #     "dispersal_target_n",
  #     "dispersal_n_k_cutoff",
  #     "dispersal_n_k_threshold",
  #     "dispersal_k_threshold",
  #     "dispersal_n_threshold",
  #     "dispersal_n_cutoff",
  #     "dispersal_target_n_k_threshold",
  #     "dispersal_target_n_k_cutoff"
  #   )
  # )
  # disease_model$set_attributes(list(
  #   dispersal_n_k_cutoff = -0.5,
  #   dispersal_n_k_threshold = 1.5
  # ))
  # expect_equal(disease_model$dispersal_source_n_k,
  #              list(cutoff = -0.5, threshold = 1.5))
  # disease_model$set_attributes(list(
  #   dispersal_n_threshold = 1,
  #   dispersal_n_cutoff = 2
  # ))
  # expect_equal(disease_model$dispersal_target_n,
  #              list(threshold = 1, cutoff = 2))
})