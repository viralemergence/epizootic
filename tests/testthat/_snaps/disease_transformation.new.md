# Warnings for negative or non-finite values in transformed data [plain]

    Code
      new_function(carrying_capacity = rep(100, 4), segment_abundance = matrix(1:12,
      nrow = 3), occupied_indices = c(1, 2), r = 1, tm = 1)
    Condition
      Warning:
      Non-finite segment abundances returned by user-defined test function at indices 1.
      Warning:
      Negative carrying capacities returned by user-defined test function at indices 1 and 2
    Output
      $abundance_threshold
      
      
      $additional_attributes
      list()
      
      $demographic_stochasticity
      
      
      $density_stages
      [1] 1
      
      $fecundity
      
      
      $fecundity_mask
      
      
      $fecundity_unit
      
      
      $mortality
      
      
      $mortality_unit
      
      
      $name
      [1] "test"
      
      $populations
      
      
      $recovery
      
      
      $recovery_mask
      
      
      $recovery_unit
      
      
      $replicates
      
      
      $seasons
      
      
      $simulator
      
      
      $stages
      [1] 1
      
      $time_steps
      
      
      $transformation
      function(params) {
                        params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                        params[["segment_abundance"]][1, 1] <- Inf
                        params[["carrying_capacity"]] <- 200
                        params[["carrying_capacity"]][2] <- -1
                        return(params)
                      }
      <environment: 0x7ff4c33d2150>
      
      $transmission
      
      
      $transmission_mask
      
      
      $transmission_unit
      
      
      $r
      [1] 1
      
      $tm
      [1] 1
      
      $carrying_capacity
      [1] 200  -1
      
      $segment_abundance
           [,1] [,2] [,3]
      [1,]  Inf    4    7
      [2,]    2    5    8
      [3,]    3    6    9
      
      $occupied_indices
      [1] 1 2
      

---

    Code
      new_function(carrying_capacity = rep(100, 4), segment_abundance = matrix(1:12,
      nrow = 3), occupied_indices = c(1, 2), r = 1, tm = 1)
    Condition
      Warning:
      Non-finite segment abundances returned by user-defined test function at indices 1.
      Warning:
      Negative carrying capacities returned by user-defined test function at indices 1 and 2
    Output
      $abundance_threshold
      
      
      $additional_attributes
      list()
      
      $demographic_stochasticity
      
      
      $density_stages
      [1] 1
      
      $fecundity
      
      
      $fecundity_mask
      
      
      $fecundity_unit
      
      
      $mortality
      
      
      $mortality_unit
      
      
      $name
      [1] "test"
      
      $populations
      
      
      $recovery
      
      
      $recovery_mask
      
      
      $recovery_unit
      
      
      $replicates
      
      
      $seasons
      
      
      $simulator
      
      
      $stages
      [1] 1
      
      $time_steps
      
      
      $transformation
      function(params) {
                        params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                        params[["segment_abundance"]][1, 1] <- Inf
                        params[["carrying_capacity"]] <- 200
                        params[["carrying_capacity"]][2] <- -1
                        return(params)
                      }
      <environment: 0x7ff4c33d2150>
      
      $transmission
      
      
      $transmission_mask
      
      
      $transmission_unit
      
      
      $r
      [1] 1
      
      $tm
      [1] 1
      
      $carrying_capacity
      [1] 200  -1
      
      $segment_abundance
           [,1] [,2] [,3]
      [1,]  Inf    4    7
      [2,]    2    5    8
      [3,]    3    6    9
      
      $occupied_indices
      [1] 1 2
      

# Warnings for negative or non-finite values in transformed data [ansi]

    Code
      new_function(carrying_capacity = rep(100, 4), segment_abundance = matrix(1:12,
      nrow = 3), occupied_indices = c(1, 2), r = 1, tm = 1)
    Condition
      [1m[33mWarning[39m:[22m
      [1m[22mNon-finite segment abundances returned by user-defined test function at indices 1.
      [1m[33mWarning[39m:[22m
      [1m[22mNegative carrying capacities returned by user-defined test function at indices 1 and 2
    Output
      $abundance_threshold
      
      
      $additional_attributes
      list()
      
      $demographic_stochasticity
      
      
      $density_stages
      [1] 1
      
      $fecundity
      
      
      $fecundity_mask
      
      
      $fecundity_unit
      
      
      $mortality
      
      
      $mortality_unit
      
      
      $name
      [1] "test"
      
      $populations
      
      
      $recovery
      
      
      $recovery_mask
      
      
      $recovery_unit
      
      
      $replicates
      
      
      $seasons
      
      
      $simulator
      
      
      $stages
      [1] 1
      
      $time_steps
      
      
      $transformation
      function(params) {
                        params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                        params[["segment_abundance"]][1, 1] <- Inf
                        params[["carrying_capacity"]] <- 200
                        params[["carrying_capacity"]][2] <- -1
                        return(params)
                      }
      <environment: 0x7ff4a3cd6748>
      
      $transmission
      
      
      $transmission_mask
      
      
      $transmission_unit
      
      
      $r
      [1] 1
      
      $tm
      [1] 1
      
      $carrying_capacity
      [1] 200  -1
      
      $segment_abundance
           [,1] [,2] [,3]
      [1,]  Inf    4    7
      [2,]    2    5    8
      [3,]    3    6    9
      
      $occupied_indices
      [1] 1 2
      

---

    Code
      new_function(carrying_capacity = rep(100, 4), segment_abundance = matrix(1:12,
      nrow = 3), occupied_indices = c(1, 2), r = 1, tm = 1)
    Condition
      [1m[33mWarning[39m:[22m
      [1m[22mNon-finite segment abundances returned by user-defined test function at indices 1.
      [1m[33mWarning[39m:[22m
      [1m[22mNegative carrying capacities returned by user-defined test function at indices 1 and 2
    Output
      $abundance_threshold
      
      
      $additional_attributes
      list()
      
      $demographic_stochasticity
      
      
      $density_stages
      [1] 1
      
      $fecundity
      
      
      $fecundity_mask
      
      
      $fecundity_unit
      
      
      $mortality
      
      
      $mortality_unit
      
      
      $name
      [1] "test"
      
      $populations
      
      
      $recovery
      
      
      $recovery_mask
      
      
      $recovery_unit
      
      
      $replicates
      
      
      $seasons
      
      
      $simulator
      
      
      $stages
      [1] 1
      
      $time_steps
      
      
      $transformation
      function(params) {
                        params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                        params[["segment_abundance"]][1, 1] <- Inf
                        params[["carrying_capacity"]] <- 200
                        params[["carrying_capacity"]][2] <- -1
                        return(params)
                      }
      <environment: 0x7ff4a3cd6748>
      
      $transmission
      
      
      $transmission_mask
      
      
      $transmission_unit
      
      
      $r
      [1] 1
      
      $tm
      [1] 1
      
      $carrying_capacity
      [1] 200  -1
      
      $segment_abundance
           [,1] [,2] [,3]
      [1,]  Inf    4    7
      [2,]    2    5    8
      [3,]    3    6    9
      
      $occupied_indices
      [1] 1 2
      

# Warnings for negative or non-finite values in transformed data [unicode]

    Code
      new_function(carrying_capacity = rep(100, 4), segment_abundance = matrix(1:12,
      nrow = 3), occupied_indices = c(1, 2), r = 1, tm = 1)
    Condition
      Warning:
      Non-finite segment abundances returned by user-defined test function at indices 1.
      Warning:
      Negative carrying capacities returned by user-defined test function at indices 1 and 2
    Output
      $abundance_threshold
      
      
      $additional_attributes
      list()
      
      $demographic_stochasticity
      
      
      $density_stages
      [1] 1
      
      $fecundity
      
      
      $fecundity_mask
      
      
      $fecundity_unit
      
      
      $mortality
      
      
      $mortality_unit
      
      
      $name
      [1] "test"
      
      $populations
      
      
      $recovery
      
      
      $recovery_mask
      
      
      $recovery_unit
      
      
      $replicates
      
      
      $seasons
      
      
      $simulator
      
      
      $stages
      [1] 1
      
      $time_steps
      
      
      $transformation
      function(params) {
                        params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                        params[["segment_abundance"]][1, 1] <- Inf
                        params[["carrying_capacity"]] <- 200
                        params[["carrying_capacity"]][2] <- -1
                        return(params)
                      }
      <environment: 0x7ff4c30cc038>
      
      $transmission
      
      
      $transmission_mask
      
      
      $transmission_unit
      
      
      $r
      [1] 1
      
      $tm
      [1] 1
      
      $carrying_capacity
      [1] 200  -1
      
      $segment_abundance
           [,1] [,2] [,3]
      [1,]  Inf    4    7
      [2,]    2    5    8
      [3,]    3    6    9
      
      $occupied_indices
      [1] 1 2
      

---

    Code
      new_function(carrying_capacity = rep(100, 4), segment_abundance = matrix(1:12,
      nrow = 3), occupied_indices = c(1, 2), r = 1, tm = 1)
    Condition
      Warning:
      Non-finite segment abundances returned by user-defined test function at indices 1.
      Warning:
      Negative carrying capacities returned by user-defined test function at indices 1 and 2
    Output
      $abundance_threshold
      
      
      $additional_attributes
      list()
      
      $demographic_stochasticity
      
      
      $density_stages
      [1] 1
      
      $fecundity
      
      
      $fecundity_mask
      
      
      $fecundity_unit
      
      
      $mortality
      
      
      $mortality_unit
      
      
      $name
      [1] "test"
      
      $populations
      
      
      $recovery
      
      
      $recovery_mask
      
      
      $recovery_unit
      
      
      $replicates
      
      
      $seasons
      
      
      $simulator
      
      
      $stages
      [1] 1
      
      $time_steps
      
      
      $transformation
      function(params) {
                        params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                        params[["segment_abundance"]][1, 1] <- Inf
                        params[["carrying_capacity"]] <- 200
                        params[["carrying_capacity"]][2] <- -1
                        return(params)
                      }
      <environment: 0x7ff4c30cc038>
      
      $transmission
      
      
      $transmission_mask
      
      
      $transmission_unit
      
      
      $r
      [1] 1
      
      $tm
      [1] 1
      
      $carrying_capacity
      [1] 200  -1
      
      $segment_abundance
           [,1] [,2] [,3]
      [1,]  Inf    4    7
      [2,]    2    5    8
      [3,]    3    6    9
      
      $occupied_indices
      [1] 1 2
      

# Warnings for negative or non-finite values in transformed data [fancy]

    Code
      new_function(carrying_capacity = rep(100, 4), segment_abundance = matrix(1:12,
      nrow = 3), occupied_indices = c(1, 2), r = 1, tm = 1)
    Condition
      [1m[33mWarning[39m:[22m
      [1m[22mNon-finite segment abundances returned by user-defined test function at indices 1.
      [1m[33mWarning[39m:[22m
      [1m[22mNegative carrying capacities returned by user-defined test function at indices 1 and 2
    Output
      $abundance_threshold
      
      
      $additional_attributes
      list()
      
      $demographic_stochasticity
      
      
      $density_stages
      [1] 1
      
      $fecundity
      
      
      $fecundity_mask
      
      
      $fecundity_unit
      
      
      $mortality
      
      
      $mortality_unit
      
      
      $name
      [1] "test"
      
      $populations
      
      
      $recovery
      
      
      $recovery_mask
      
      
      $recovery_unit
      
      
      $replicates
      
      
      $seasons
      
      
      $simulator
      
      
      $stages
      [1] 1
      
      $time_steps
      
      
      $transformation
      function(params) {
                        params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                        params[["segment_abundance"]][1, 1] <- Inf
                        params[["carrying_capacity"]] <- 200
                        params[["carrying_capacity"]][2] <- -1
                        return(params)
                      }
      <environment: 0x7ff4a793c208>
      
      $transmission
      
      
      $transmission_mask
      
      
      $transmission_unit
      
      
      $r
      [1] 1
      
      $tm
      [1] 1
      
      $carrying_capacity
      [1] 200  -1
      
      $segment_abundance
           [,1] [,2] [,3]
      [1,]  Inf    4    7
      [2,]    2    5    8
      [3,]    3    6    9
      
      $occupied_indices
      [1] 1 2
      

---

    Code
      new_function(carrying_capacity = rep(100, 4), segment_abundance = matrix(1:12,
      nrow = 3), occupied_indices = c(1, 2), r = 1, tm = 1)
    Condition
      [1m[33mWarning[39m:[22m
      [1m[22mNon-finite segment abundances returned by user-defined test function at indices 1.
      [1m[33mWarning[39m:[22m
      [1m[22mNegative carrying capacities returned by user-defined test function at indices 1 and 2
    Output
      $abundance_threshold
      
      
      $additional_attributes
      list()
      
      $demographic_stochasticity
      
      
      $density_stages
      [1] 1
      
      $fecundity
      
      
      $fecundity_mask
      
      
      $fecundity_unit
      
      
      $mortality
      
      
      $mortality_unit
      
      
      $name
      [1] "test"
      
      $populations
      
      
      $recovery
      
      
      $recovery_mask
      
      
      $recovery_unit
      
      
      $replicates
      
      
      $seasons
      
      
      $simulator
      
      
      $stages
      [1] 1
      
      $time_steps
      
      
      $transformation
      function(params) {
                        params[["segment_abundance"]] <- matrix(1:9, nrow = 3)
                        params[["segment_abundance"]][1, 1] <- Inf
                        params[["carrying_capacity"]] <- 200
                        params[["carrying_capacity"]][2] <- -1
                        return(params)
                      }
      <environment: 0x7ff4a793c208>
      
      $transmission
      
      
      $transmission_mask
      
      
      $transmission_unit
      
      
      $r
      [1] 1
      
      $tm
      [1] 1
      
      $carrying_capacity
      [1] 200  -1
      
      $segment_abundance
           [,1] [,2] [,3]
      [1,]  Inf    4    7
      [2,]    2    5    8
      [3,]    3    6    9
      
      $occupied_indices
      [1] 1 2
      

