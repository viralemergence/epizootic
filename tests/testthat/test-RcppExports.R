test_that("aspatial_siri works with valid inputs", {
  expect_silent(aspatial_siri(
    initial_pop = c(50000, 50000, 0, 1, 0, 0, 0, 0),
    season_length = 100,
    mortality = c(0.004, 0, 0.00505, 0.00105, 0.004, 0, 0.0045, 5e-04),
    fecundity = c(0, 15/182, 0, 15/182, 0, 15/182, 0, 15/182),
    transmission = c(0.00002, 0.00001, 0, 0, 7.84e-06, 3.92e-06, 0, 0),
    recovery = c(0, 0, 0.05714286, 0.05714286, 0, 0, 0.1, 0.1),
    carrying_capacity = 150000,
    abundance_threshold = 10,
    season = "breeding"
  ))
})