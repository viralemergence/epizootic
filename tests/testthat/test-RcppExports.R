test_that("aspatial_siri works with valid inputs", {
  expect_silent(aspatial_siri(
    pop_list = list(c(50000, 50000, 0, 1, 0, 0, 0, 0)),
    mortality_list = list(c(0.004, 0, 0.00505, 0.00105, 0.004, 0, 0.0045, 5e-04) |>
      (`*`)(100)),
    fecundity_list = list(0.08241758 |> (`*`)(100)),
    transmission_list = list(c(0.00002, 0.00001, 7.84e-06, 3.92e-06) |>
      (`*`)(100)),
    recovery_list = list(rep(1, 4)),
    carrying_capacity_list = 150000,
    abundance_threshold = 10,
    season = "breeding"
  ))
  expect_silent(aspatial_siri(
    pop_list = list(c(50000, 50000, 0, 1, 0, 0, 0, 0)),
    mortality_list = list(c(0.004, 0, 0.00505, 0.00105, 0.004, 0, 0.0045, 5e-04) |>
                            (`*`)(100)),
    fecundity_list = list(0.08241758 |> (`*`)(100)),
    transmission_list = list(c(0.00002, 0.00001, 7.84e-06, 3.92e-06) |>
                               (`*`)(100)),
    recovery_list = list(rep(1, 4)),
    carrying_capacity_list = 150000,
    abundance_threshold = 10,
    season = "non-breeding"
  ))
})
