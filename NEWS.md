# v 1.0.0 (27 Sep 2024)

## Enhancements

- A vignette has been added using `epizootic` to model house finch conjunctivitis.
- Examples have been added to two of the functions.
- `aspatial_siri` now runs much faster.
- `check_aspatial_siri_inputs` now checks inputs much more thoroughly.
- There are now unit tests to ensure that no dispersers are created or destroyed.
- New tests have been added for `DiseaseModel`.

## Bug fixes
- Broken links to other packages have been fixed.
- `siri_model` functions can now handle one active population.
- `disease_simulator` no longer breaks when demographic rates are given as lists and masks are given as vectors.
- Non-breeding season length is now calculated correctly.
- `disease_dispersal` can now handle the case when dispersal does not change over time.
- `disease_dispersal` was not working properly when there was density dependent dispersal and matrix-formatted dispersal data. This is now fixed.

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
