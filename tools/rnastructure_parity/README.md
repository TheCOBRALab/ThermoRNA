# RNAstructure parity check

Proves that ThermoRNA's RNAstructure loader produces exactly the tables efn2 works with.

`rnastructure_dump` links RNAstructure's own `rna_library.cpp` and runs `datatable::opendat()`
(plus `ScaleToTemperature()` when a temperature is given), i.e. what efn2 does before it
evaluates a structure. `thermorna_dump` does the same through `AlphabetTables::load()`. Both
print every scalar (floats as raw bits), every special hairpin and every table entry that is not
`INFINITE_ENERGY`, so equal output means equal tables.

```bash
cmake -S . -B build -DTHERMORNA_RNASTRUCTURE=ON \
      -DTHERMORNA_RNASTRUCTURE_SOURCE=/path/to/RNAstructure
cmake --build build
ctest --test-dir build -R RNAstructureParity --output-on-failure
```

The default cases are every shipped alphabet at 37 °C, `rna` at eight temperatures between
273.15 K and 400 K, and alphabets without `.dh` files at other temperatures (both loaders must
refuse those). Pass `-DCASES="rna:310.15;m6A:310.15"` to `compare.cmake` for other cases.

Checked against RNAstructure 6.6 (April 2026): all cases identical.
