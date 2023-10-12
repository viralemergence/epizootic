library(terra)
library(purrr)

test_that("raster interpolates", {
  m1 <- matrix(1:25, nrow=5, ncol=5) |> rast()
  m2 <- matrix(2:26, nrow=5, ncol=5) |> rast()
  m3 <- matrix(3:27, nrow=5, ncol=5) |> rast()
  r <- c(m1, m2, m3)
  names(r) <- c("y1", "y2", "y3")
  test_r <- interpolate_raster(c(m1, m3), c(1, 3), 1:3, time_label = "y")
  expect_equal(values(test_r), values(r))
  expect_equal(names(test_r), names(r))
})
