# v 0.2.0 (28 Mar 2024)

- The aspatial_siri function now kills off populations at or below the abundance threshold (previously it killed off populations below the abundance threshold only)
- The disease_dispersal function can now handle the case when one or more segments of the population does not disperse.

# v 0.1.1

- Switched parallelization engine from future/furrr to doParallel/foreach
- Patched errors in dispersal (overcrowded cells were not re-dispersed properly)

# v 0.1.0

Initial release