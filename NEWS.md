# v 0.2.1 (2 Apr 2024)

- Important bug fix: the number of occupied indices now updates more frequently within the 
`disease_simulator` function, which allows species range to expand properly.
- Documentation update: all references to features in development that are not yet
implemented have been removed.

# v 0.2.0 (28 Mar 2024)

- The aspatial_siri function now kills off populations at or below the abundance threshold (previously it killed off populations below the abundance threshold only)
- The disease_dispersal function can now handle the case when one or more segments of the population does not disperse.

# v 0.1.1

- Switched parallelization engine from future/furrr to doParallel/foreach
- Patched errors in dispersal (overcrowded cells were not re-dispersed properly)

# v 0.1.0

Initial release
