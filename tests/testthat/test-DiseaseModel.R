test_that("active get and set", {
  region = Region$new(coordinates = array(c(1:4, 4:1), c(7, 2)))
  disease_model <- DiseaseModel$new(time_steps = 10)
  expect_null(disease_model$populations)
  disease_model$initial_abundance <- seq(10, 60, by = 10)
  expect_equal(disease_model$populations, 6)
  disease_model$region <- region
  expect_equal(disease_model$populations, 7)
})
