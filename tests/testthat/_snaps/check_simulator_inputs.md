# minimal inputs [plain]

    Minimal inputs required to run simulation are missing.
    x Your input does not include time_steps, populations, initial_abundance, carrying_capacity, fecundity, mortality, transmission, and simulation_order.

# minimal inputs [ansi]

    [1m[22mMinimal inputs required to run simulation are missing.
    [31mx[39m Your input does not include time_steps, populations, initial_abundance, carrying_capacity, fecundity, mortality, transmission, and simulation_order.

# minimal inputs [unicode]

    Minimal inputs required to run simulation are missing.
    ✖ Your input does not include time_steps, populations, initial_abundance, carrying_capacity, fecundity, mortality, transmission, and simulation_order.

# minimal inputs [fancy]

    [1m[22mMinimal inputs required to run simulation are missing.
    [31m✖[39m Your input does not include time_steps, populations, initial_abundance, carrying_capacity, fecundity, mortality, transmission, and simulation_order.

# check error handling for initial abundance [plain]

    initial_abundance must be a numeric vector, raster, matrix, or array.
    x initial_abundance is a string.

---

    initial_abundance has 1 row.
    x There should be 2 rows (compartments x stages).

---

    initial_abundance has 5 columns.
    x There should be 6 columns.

# check error handling for initial abundance [ansi]

    [1m[22minitial_abundance must be a numeric vector, raster, matrix, or array.
    [31mx[39m initial_abundance is a string.

---

    [1m[22minitial_abundance has 1 row.
    [31mx[39m There should be 2 rows (compartments x stages).

---

    [1m[22minitial_abundance has 5 columns.
    [31mx[39m There should be 6 columns.

# check error handling for initial abundance [unicode]

    initial_abundance must be a numeric vector, raster, matrix, or array.
    ✖ initial_abundance is a string.

---

    initial_abundance has 1 row.
    ✖ There should be 2 rows (compartments x stages).

---

    initial_abundance has 5 columns.
    ✖ There should be 6 columns.

# check error handling for initial abundance [fancy]

    [1m[22minitial_abundance must be a numeric vector, raster, matrix, or array.
    [31m✖[39m initial_abundance is a string.

---

    [1m[22minitial_abundance has 1 row.
    [31m✖[39m There should be 2 rows (compartments x stages).

---

    [1m[22minitial_abundance has 5 columns.
    [31m✖[39m There should be 6 columns.

# check error handling for breeding season length [plain]

    There are 1 missing values in the breeding_season_length object.

---

    The length of `season_lengths` must equal `seasons`.
    i `seasons` = 2.
    x `season_lengths` = 300, 100, and 65.

# check error handling for breeding season length [ansi]

    [1m[22mThere are 1 missing values in the breeding_season_length object.

---

    [1m[22mThe length of `season_lengths` must equal `seasons`.
    [36mi[39m `seasons` = 2.
    [31mx[39m `season_lengths` = 300, 100, and 65.

# check error handling for breeding season length [unicode]

    There are 1 missing values in the breeding_season_length object.

---

    The length of `season_lengths` must equal `seasons`.
    ℹ `seasons` = 2.
    ✖ `season_lengths` = 300, 100, and 65.

# check error handling for breeding season length [fancy]

    [1m[22mThere are 1 missing values in the breeding_season_length object.

---

    [1m[22mThe length of `season_lengths` must equal `seasons`.
    [36mℹ[39m `seasons` = 2.
    [31m✖[39m `season_lengths` = 300, 100, and 65.

# check error handling for mortality [plain]

    Each vector inside the `mortality` list must have 2 elements, one for each combination of stage & compartment.

---

    `mortality` must be a vector or list, not a <RasterLayer> object.

---

    Each vector inside the `mortality` list must have 6 elements, one for each combination of stage & compartment.

# check error handling for mortality [ansi]

    [1m[22mEach vector inside the `mortality` list must have 2 elements, one for each combination of stage & compartment.

---

    [1m[22m`mortality` must be a vector or list, not a [34m<RasterLayer>[39m object.

---

    [1m[22mEach vector inside the `mortality` list must have 6 elements, one for each combination of stage & compartment.

# check error handling for mortality [unicode]

    Each vector inside the `mortality` list must have 2 elements, one for each combination of stage & compartment.

---

    `mortality` must be a vector or list, not a <RasterLayer> object.

---

    Each vector inside the `mortality` list must have 6 elements, one for each combination of stage & compartment.

# check error handling for mortality [fancy]

    [1m[22mEach vector inside the `mortality` list must have 2 elements, one for each combination of stage & compartment.

---

    [1m[22m`mortality` must be a vector or list, not a [34m<RasterLayer>[39m object.

---

    [1m[22mEach vector inside the `mortality` list must have 6 elements, one for each combination of stage & compartment.

# check error handling for fecundity [plain]

    Each vector inside the `fecundity` list must have 2 elements, one for each combination of stage and compartment. Or, provide a fecundity mask so that fecundity values may be assigned to the appropriate stages and compartments.

---

    `fecundity` must be a vector or list, not a <RasterLayer> object.

---

    vectors inside `fecundity` and `fecundity_unit` must be the same length.
    * `fecundity` vectors are lengths 3 and 3.
    * `fecundity_unit` vectors are lengths 3 and 2.

---

    fecundity_unit values must be 0 or 1. Unit values are 1, 0, and 2.

---

    Vectors inside `fecundity` and `fecundity_mask` must be the same length.
    * `fecundity` has vectors of length 3.
    * `fecundity_mask` has vectors of length 4.

---

    fecundity_mask values must be 0 or 1

# check error handling for fecundity [ansi]

    [1m[22mEach vector inside the `fecundity` list must have 2 elements, one for each combination of stage and compartment. Or, provide a fecundity mask so that fecundity values may be assigned to the appropriate stages and compartments.

---

    [1m[22m`fecundity` must be a vector or list, not a [34m<RasterLayer>[39m object.

---

    [1m[22mvectors inside `fecundity` and `fecundity_unit` must be the same length.
    [36m*[39m `fecundity` vectors are lengths 3 and 3.
    [36m*[39m `fecundity_unit` vectors are lengths 3 and 2.

---

    [1m[22mfecundity_unit values must be 0 or 1. Unit values are 1, 0, and 2.

---

    [1m[22mVectors inside `fecundity` and `fecundity_mask` must be the same length.
    [36m*[39m `fecundity` has vectors of length 3.
    [36m*[39m `fecundity_mask` has vectors of length 4.

---

    [1m[22mfecundity_mask values must be 0 or 1

# check error handling for fecundity [unicode]

    Each vector inside the `fecundity` list must have 2 elements, one for each combination of stage and compartment. Or, provide a fecundity mask so that fecundity values may be assigned to the appropriate stages and compartments.

---

    `fecundity` must be a vector or list, not a <RasterLayer> object.

---

    vectors inside `fecundity` and `fecundity_unit` must be the same length.
    • `fecundity` vectors are lengths 3 and 3.
    • `fecundity_unit` vectors are lengths 3 and 2.

---

    fecundity_unit values must be 0 or 1. Unit values are 1, 0, and 2.

---

    Vectors inside `fecundity` and `fecundity_mask` must be the same length.
    • `fecundity` has vectors of length 3.
    • `fecundity_mask` has vectors of length 4.

---

    fecundity_mask values must be 0 or 1

# check error handling for fecundity [fancy]

    [1m[22mEach vector inside the `fecundity` list must have 2 elements, one for each combination of stage and compartment. Or, provide a fecundity mask so that fecundity values may be assigned to the appropriate stages and compartments.

---

    [1m[22m`fecundity` must be a vector or list, not a [34m<RasterLayer>[39m object.

---

    [1m[22mvectors inside `fecundity` and `fecundity_unit` must be the same length.
    [36m•[39m `fecundity` vectors are lengths 3 and 3.
    [36m•[39m `fecundity_unit` vectors are lengths 3 and 2.

---

    [1m[22mfecundity_unit values must be 0 or 1. Unit values are 1, 0, and 2.

---

    [1m[22mVectors inside `fecundity` and `fecundity_mask` must be the same length.
    [36m•[39m `fecundity` has vectors of length 3.
    [36m•[39m `fecundity_mask` has vectors of length 4.

---

    [1m[22mfecundity_mask values must be 0 or 1

