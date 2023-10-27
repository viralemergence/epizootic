# minimal inputs [plain]

    Minimal inputs required to run simulation are missing.
    x Your input does not include time_steps, populations, initial_abundance, carrying_capacity, and fecundity.

# minimal inputs [ansi]

    [1m[22mMinimal inputs required to run simulation are missing.
    [31mx[39m Your input does not include time_steps, populations, initial_abundance, carrying_capacity, and fecundity.

# minimal inputs [unicode]

    Minimal inputs required to run simulation are missing.
    ✖ Your input does not include time_steps, populations, initial_abundance, carrying_capacity, and fecundity.

# minimal inputs [fancy]

    [1m[22mMinimal inputs required to run simulation are missing.
    [31m✖[39m Your input does not include time_steps, populations, initial_abundance, carrying_capacity, and fecundity.

# check validity of initial abundance [plain]

    `initial_abundance` must be a numeric vector, raster, matrix, or array.
    x `initial_abundance` is a a string.

---

    `initial_abundance` has 1 row.
    x There should be 2 rows (compartments x stages).

---

    `initial_abundance` has 5 columns.
    x There should be 6 columns (`populations`).

# check validity of initial abundance [ansi]

    [1m[22m`initial_abundance` must be a numeric vector, raster, matrix, or array.
    [31mx[39m `initial_abundance` is a a string.

---

    [1m[22m`initial_abundance` has 1 row.
    [31mx[39m There should be 2 rows (compartments x stages).

---

    [1m[22m`initial_abundance` has 5 columns.
    [31mx[39m There should be 6 columns (`populations`).

# check validity of initial abundance [unicode]

    `initial_abundance` must be a numeric vector, raster, matrix, or array.
    ✖ `initial_abundance` is a a string.

---

    `initial_abundance` has 1 row.
    ✖ There should be 2 rows (compartments x stages).

---

    `initial_abundance` has 5 columns.
    ✖ There should be 6 columns (`populations`).

# check validity of initial abundance [fancy]

    [1m[22m`initial_abundance` must be a numeric vector, raster, matrix, or array.
    [31m✖[39m `initial_abundance` is a a string.

---

    [1m[22m`initial_abundance` has 1 row.
    [31m✖[39m There should be 2 rows (compartments x stages).

---

    [1m[22m`initial_abundance` has 5 columns.
    [31m✖[39m There should be 6 columns (`populations`).

# check validity of breeding season length [plain]

    There are 1 missing values in the breeding_season_length object.

---

    The length of `season_lengths` must equal `seasons`.
    i `seasons` = 2.
    x `season_lengths` = 300, 100, and 65.

# check validity of breeding season length [ansi]

    [1m[22mThere are 1 missing values in the breeding_season_length object.

---

    [1m[22mThe length of `season_lengths` must equal `seasons`.
    [36mi[39m `seasons` = 2.
    [31mx[39m `season_lengths` = 300, 100, and 65.

# check validity of breeding season length [unicode]

    There are 1 missing values in the breeding_season_length object.

---

    The length of `season_lengths` must equal `seasons`.
    ℹ `seasons` = 2.
    ✖ `season_lengths` = 300, 100, and 65.

# check validity of breeding season length [fancy]

    [1m[22mThere are 1 missing values in the breeding_season_length object.

---

    [1m[22mThe length of `season_lengths` must equal `seasons`.
    [36mℹ[39m `seasons` = 2.
    [31m✖[39m `season_lengths` = 300, 100, and 65.

