# Test disease_transitions function

test_that("disease_transitions works with single stage and compartment", {
  # Test the simplest case: 1 stage, 1 compartment
  transition_fn <- disease_transitions(stages = 1, compartments = 1)

  # Create a simple abundance matrix: 1 row (1 stage × 1 compartment), 3 populations
  segment_abundance <- matrix(c(10, 20, 30), nrow = 1)
  occupied_indices <- c(1, 2, 3)

  result <- transition_fn(segment_abundance, occupied_indices)

  # With only 1 stage, transitions should not change anything
  expect_equal(result, segment_abundance)
  expect_equal(dim(result), dim(segment_abundance))
})

test_that("disease_transitions advances juveniles to adults with 2 stages", {
  # Test with 2 stages (juvenile, adult), 1 compartment
  transition_fn <- disease_transitions(stages = 2, compartments = 1)

  # 2 rows: juvenile, adult; 3 populations
  segment_abundance <- matrix(c(
    100, 200, 150, # juveniles
    50, 100, 75 # adults
  ), nrow = 2, byrow = TRUE)
  occupied_indices <- c(1, 2, 3)

  result <- transition_fn(segment_abundance, occupied_indices)

  # Juveniles should become adults (row 1 -> row 2)
  # Row 1 (juveniles) should be 0 after transition
  expect_equal(result[1, ], c(0, 0, 0))

  # Row 2 (adults) should have original adults + transitioned juveniles
  expect_equal(result[2, ], c(150, 300, 225))
})

test_that("disease_transitions works with 2 stages and 2 compartments (SI model)", {
  # Test with 2 stages, 2 compartments (Susceptible, Infected)
  transition_fn <- disease_transitions(stages = 2, compartments = 2)

  # 4 rows: S_juv, S_adult, I_juv, I_adult; 2 populations
  segment_abundance <- matrix(c(
    100, 150, # S_juv
    50, 75, # S_adult
    10, 20, # I_juv
    5, 10 # I_adult
  ), nrow = 4, byrow = TRUE)
  occupied_indices <- c(1, 2)

  result <- transition_fn(segment_abundance, occupied_indices)

  # Juveniles should transition to adults within each compartment
  # S_juv (row 1) -> S_adult (row 2)
  expect_equal(result[1, ], c(0, 0))
  expect_equal(result[2, ], c(150, 225)) # 50 + 100, 75 + 150

  # I_juv (row 3) -> I_adult (row 4)
  expect_equal(result[3, ], c(0, 0))
  expect_equal(result[4, ], c(15, 30)) # 5 + 10, 10 + 20
})

test_that("disease_transitions works with 2 stages and 4 compartments (SIRI model)", {
  # Test with 2 stages, 4 compartments (S, I1, R, I2)
  transition_fn <- disease_transitions(stages = 2, compartments = 4)

  # 8 rows: S_juv, S_adult, I1_juv, I1_adult, R_juv, R_adult, I2_juv, I2_adult
  # 3 populations
  segment_abundance <- matrix(c(
    1000, 2000, 1500, # S_juv
    500, 1000, 750, # S_adult
    50, 100, 75, # I1_juv
    25, 50, 40, # I1_adult
    20, 40, 30, # R_juv
    10, 20, 15, # R_adult
    5, 10, 8, # I2_juv
    2, 5, 4 # I2_adult
  ), nrow = 8, byrow = TRUE)
  occupied_indices <- c(1, 2, 3)

  result <- transition_fn(segment_abundance, occupied_indices)

  # All juveniles (odd rows) should become 0
  expect_equal(result[1, ], c(0, 0, 0))
  expect_equal(result[3, ], c(0, 0, 0))
  expect_equal(result[5, ], c(0, 0, 0))
  expect_equal(result[7, ], c(0, 0, 0))

  # Adults should accumulate juveniles
  expect_equal(result[2, ], c(1500, 3000, 2250)) # S_adult
  expect_equal(result[4, ], c(75, 150, 115)) # I1_adult
  expect_equal(result[6, ], c(30, 60, 45)) # R_adult
  expect_equal(result[8, ], c(7, 15, 12)) # I2_adult
})

test_that("disease_transitions respects occupied_indices", {
  # Test that only occupied populations are transitioned
  transition_fn <- disease_transitions(stages = 2, compartments = 1)

  # 2 stages, 5 populations
  segment_abundance <- matrix(c(
    100, 200, 300, 400, 500, # juveniles
    50, 100, 150, 200, 250 # adults
  ), nrow = 2, byrow = TRUE)

  # Only populations 1, 3, and 5 are occupied
  occupied_indices <- c(1, 3, 5)

  result <- transition_fn(segment_abundance, occupied_indices)

  # Occupied populations should transition
  expect_equal(result[1, 1], 0)
  expect_equal(result[2, 1], 150) # 50 + 100

  expect_equal(result[1, 3], 0)
  expect_equal(result[2, 3], 450) # 150 + 300

  expect_equal(result[1, 5], 0)
  expect_equal(result[2, 5], 750) # 250 + 500

  # Unoccupied populations (2, 4) should become 0
  expect_equal(result[1, 2], 0)
  expect_equal(result[2, 2], 0)
  expect_equal(result[1, 4], 0)
  expect_equal(result[2, 4], 0)
})

test_that("disease_transitions works with 3 stages", {
  # Test with 3 stages (juvenile, yearling, adult), 1 compartment
  transition_fn <- disease_transitions(stages = 3, compartments = 1)

  # 3 rows, 2 populations
  segment_abundance <- matrix(c(
    100, 150, # juveniles
    50, 75, # yearlings
    200, 300 # adults
  ), nrow = 3, byrow = TRUE)
  occupied_indices <- c(1, 2)

  result <- transition_fn(segment_abundance, occupied_indices)

  # Juveniles -> yearlings
  expect_equal(result[1, ], c(0, 0))
  expect_equal(result[2, ], c(100, 150))

  # Yearlings -> adults, adults stay adults
  expect_equal(result[3, ], c(250, 375)) # 200 + 50, 300 + 75
})

test_that("disease_transitions works with 3 stages and 3 compartments", {
  # Test with 3 stages, 3 compartments (SIR model)
  transition_fn <- disease_transitions(stages = 3, compartments = 3)

  # 9 rows: S_juv, S_year, S_adult, I_juv, I_year, I_adult, R_juv, R_year, R_adult
  # 2 populations
  segment_abundance <- matrix(c(
    100, 150, # S_juv
    50, 75, # S_year
    200, 300, # S_adult
    10, 20, # I_juv
    5, 10, # I_year
    30, 40, # I_adult
    8, 12, # R_juv
    4, 6, # R_year
    20, 30 # R_adult
  ), nrow = 9, byrow = TRUE)
  occupied_indices <- c(1, 2)

  result <- transition_fn(segment_abundance, occupied_indices)

  # Check juveniles all become 0
  expect_equal(result[1, ], c(0, 0)) # S_juv
  expect_equal(result[4, ], c(0, 0)) # I_juv
  expect_equal(result[7, ], c(0, 0)) # R_juv

  # Check yearlings get juveniles
  expect_equal(result[2, ], c(100, 150)) # S_year
  expect_equal(result[5, ], c(10, 20)) # I_year
  expect_equal(result[8, ], c(8, 12)) # R_year

  # Check adults accumulate yearlings
  expect_equal(result[3, ], c(250, 375)) # S_adult: 200 + 50, 300 + 75
  expect_equal(result[6, ], c(35, 50)) # I_adult: 30 + 5, 40 + 10
  expect_equal(result[9, ], c(24, 36)) # R_adult: 20 + 4, 30 + 6
})

test_that("disease_transitions handles empty populations", {
  # Test with no occupied populations
  transition_fn <- disease_transitions(stages = 2, compartments = 2)

  segment_abundance <- matrix(c(
    100, 0, 50,
    50, 0, 25,
    10, 0, 5,
    5, 0, 2
  ), nrow = 4, byrow = TRUE)

  # Only population 1 is occupied
  occupied_indices <- c(1)

  result <- transition_fn(segment_abundance, occupied_indices)

  # Population 1 should transition
  expect_equal(result[1, 1], 0)
  expect_equal(result[2, 1], 150)
  expect_equal(result[3, 1], 0)
  expect_equal(result[4, 1], 15)

  # Populations 2 and 3 should be all zeros
  expect_equal(result[, 2], c(0, 0, 0, 0))
  expect_equal(result[, 3], c(0, 0, 0, 0))
})

test_that("disease_transitions preserves matrix dimensions", {
  # Test that output has same dimensions as input
  transition_fn <- disease_transitions(stages = 2, compartments = 3)

  segment_abundance <- matrix(1:30, nrow = 6, ncol = 5)
  occupied_indices <- c(1, 2, 3, 4, 5)

  result <- transition_fn(segment_abundance, occupied_indices)

  expect_equal(nrow(result), nrow(segment_abundance))
  expect_equal(ncol(result), ncol(segment_abundance))
})

test_that("disease_transitions handles decimal abundances", {
  # Test with non-integer abundances (can occur with certain demographic models)
  transition_fn <- disease_transitions(stages = 2, compartments = 1)

  segment_abundance <- matrix(c(
    10.5, 20.3, 15.7,
    5.2, 10.8, 7.9
  ), nrow = 2, byrow = TRUE)
  occupied_indices <- c(1, 2, 3)

  result <- transition_fn(segment_abundance, occupied_indices)

  expect_equal(result[1, ], c(0, 0, 0))
  expect_equal(result[2, ], c(15.7, 31.1, 23.6), tolerance = 1e-10)
})

test_that("disease_transitions works with large number of stages and compartments", {
  # Test scalability
  transition_fn <- disease_transitions(stages = 5, compartments = 5)

  # 25 rows (5 stages × 5 compartments), 10 populations
  # Set seed for reproducible test
  set.seed(123)
  segment_abundance <- matrix(runif(250, 0, 100), nrow = 25, ncol = 10)
  occupied_indices <- 1:10

  result <- transition_fn(segment_abundance, occupied_indices)

  # Check dimensions preserved
  expect_equal(dim(result), dim(segment_abundance))

  # In a 5-stage, 5-compartment system:
  # Rows are ordered: stage1_comp1, stage2_comp1, ..., stage5_comp1, stage1_comp2, ...
  # Adult indices (stage 5) are: 5, 10, 15, 20, 25
  adult_indices <- seq(5, 25, by = 5)

  # Stage 1 individuals (rows 1, 6, 11, 16, 21) should become 0 (they transition to stage 2)
  stage_1_indices <- seq(1, 21, by = 5)
  for (i in stage_1_indices) {
    expect_equal(result[i, ], rep(0, 10))
  }

  # Adults should have accumulated: original adults + stage 4 individuals
  for (row in adult_indices) {
    # Get the stage 4 row in this compartment (adult_row - 1)
    stage_4_row <- row - 1
    expected_min <- segment_abundance[row, ] + segment_abundance[stage_4_row, ]
    expect_true(all(result[row, ] >= segment_abundance[row, ]))
  }
})

test_that("disease_transitions returns a function", {
  result <- disease_transitions(stages = 2, compartments = 3)
  expect_true(is.function(result))
})

test_that("disease_transitions closure captures stages and compartments correctly", {
  # Test that the returned function properly captures the stage/compartment values
  fn1 <- disease_transitions(stages = 2, compartments = 2)
  fn2 <- disease_transitions(stages = 3, compartments = 1)

  # Different input matrices for each function
  segment_abundance1 <- matrix(c(
    100, 50,
    50, 25,
    30, 15,
    20, 10
  ), nrow = 4, byrow = TRUE)

  segment_abundance2 <- matrix(c(
    100, 50,
    50, 25,
    30, 15
  ), nrow = 3, byrow = TRUE)

  occupied_indices <- c(1, 2)

  result1 <- fn1(segment_abundance1, occupied_indices)
  result2 <- fn2(segment_abundance2, occupied_indices)

  # Results should be different because they use different stage/compartment structures
  expect_false(identical(result1, result2))

  # fn1 treats this as 2 stages × 2 compartments (4 segments)
  expect_equal(nrow(result1), 4)

  # fn2 treats this as 3 stages × 1 compartment (3 segments)
  expect_equal(nrow(result2), 3)
})
