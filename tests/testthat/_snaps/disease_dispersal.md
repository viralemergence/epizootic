# disease_dispersal throws an error for invalid dispersal length in 'pooled' dispersal type [plain]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3),
      dispersal_type = "pooled", stages = 3, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      Error in `disease_dispersal()`:
      ! Error: Dispersal length should be 1 for 'pooled' dispersal type.
      x Dispersal length is 2.

# disease_dispersal throws an error for invalid dispersal length in 'pooled' dispersal type [ansi]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3),
      dispersal_type = "pooled", stages = 3, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      [1m[33mError[39m in `disease_dispersal()`:[22m
      [1m[22m[33m![39m Error: Dispersal length should be 1 for 'pooled' dispersal type.
      [31mx[39m Dispersal length is 2.

# disease_dispersal throws an error for invalid dispersal length in 'pooled' dispersal type [unicode]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3),
      dispersal_type = "pooled", stages = 3, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      Error in `disease_dispersal()`:
      ! Error: Dispersal length should be 1 for 'pooled' dispersal type.
      ✖ Dispersal length is 2.

# disease_dispersal throws an error for invalid dispersal length in 'pooled' dispersal type [fancy]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3),
      dispersal_type = "pooled", stages = 3, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      [1m[33mError[39m in `disease_dispersal()`:[22m
      [1m[22m[33m![39m Error: Dispersal length should be 1 for 'pooled' dispersal type.
      [31m✖[39m Dispersal length is 2.

# disease_dispersal throws an error for invalid dispersal length in 'stages' dispersal type [plain]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "stages", stages = 4, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      Error in `disease_dispersal()`:
      ! Error: Dispersal length should be equal to the number of stages for 'stages' dispersal type.
      x There are 4 stages and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'stages' dispersal type [ansi]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "stages", stages = 4, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      [1m[33mError[39m in `disease_dispersal()`:[22m
      [1m[22m[33m![39m Error: Dispersal length should be equal to the number of stages for 'stages' dispersal type.
      [31mx[39m There are 4 stages and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'stages' dispersal type [unicode]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "stages", stages = 4, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      Error in `disease_dispersal()`:
      ! Error: Dispersal length should be equal to the number of stages for 'stages' dispersal type.
      ✖ There are 4 stages and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'stages' dispersal type [fancy]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "stages", stages = 4, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      [1m[33mError[39m in `disease_dispersal()`:[22m
      [1m[22m[33m![39m Error: Dispersal length should be equal to the number of stages for 'stages' dispersal type.
      [31m✖[39m There are 4 stages and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'compartments' dispersal type [plain]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "compartments", stages = 3, compartments = 5, simulator = "simulator")
    Condition
      Error in `disease_dispersal()`:
      ! Error: Dispersal length should be equal to the number of compartments for 'compartments' dispersal type.
      x There are 5 compartments and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'compartments' dispersal type [ansi]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "compartments", stages = 3, compartments = 5, simulator = "simulator")
    Condition
      [1m[33mError[39m in `disease_dispersal()`:[22m
      [1m[22m[33m![39m Error: Dispersal length should be equal to the number of compartments for 'compartments' dispersal type.
      [31mx[39m There are 5 compartments and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'compartments' dispersal type [unicode]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "compartments", stages = 3, compartments = 5, simulator = "simulator")
    Condition
      Error in `disease_dispersal()`:
      ! Error: Dispersal length should be equal to the number of compartments for 'compartments' dispersal type.
      ✖ There are 5 compartments and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'compartments' dispersal type [fancy]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "compartments", stages = 3, compartments = 5, simulator = "simulator")
    Condition
      [1m[33mError[39m in `disease_dispersal()`:[22m
      [1m[22m[33m![39m Error: Dispersal length should be equal to the number of compartments for 'compartments' dispersal type.
      [31m✖[39m There are 5 compartments and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'segments' dispersal type [plain]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "segments", stages = 3, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      Error in `disease_dispersal()`:
      ! Error: Dispersal length should be equal to the number of stages multiplied by the number of compartments for 'segments' dispersal type.
      x There are 3 stages and 4 compartments and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'segments' dispersal type [ansi]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "segments", stages = 3, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      [1m[33mError[39m in `disease_dispersal()`:[22m
      [1m[22m[33m![39m Error: Dispersal length should be equal to the number of stages multiplied by the number of compartments for 'segments' dispersal type.
      [31mx[39m There are 3 stages and 4 compartments and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'segments' dispersal type [unicode]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "segments", stages = 3, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      Error in `disease_dispersal()`:
      ! Error: Dispersal length should be equal to the number of stages multiplied by the number of compartments for 'segments' dispersal type.
      ✖ There are 3 stages and 4 compartments and dispersal length is 3.

# disease_dispersal throws an error for invalid dispersal length in 'segments' dispersal type [fancy]

    Code
      disease_dispersal(replicates = 10, time_steps = 100, populations = c(100, 200,
        300), demographic_stochasticity = TRUE, dispersal = c(0.2, 0.3, 0.4),
      dispersal_type = "segments", stages = 3, compartments = 4, simulator = ModelSimulator$
        new())
    Condition
      [1m[33mError[39m in `disease_dispersal()`:[22m
      [1m[22m[33m![39m Error: Dispersal length should be equal to the number of stages multiplied by the number of compartments for 'segments' dispersal type.
      [31m✖[39m There are 3 stages and 4 compartments and dispersal length is 3.

