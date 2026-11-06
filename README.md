# ThermoRNA (DO NOT MAKE THIS REPO PRIVATE)

This repository is REQUIRED for Knotergy to run. Setting it to private will prevent people from compiling it.

## What is ThermoRNA
It includes:
- Lightweight package of ViennaRNA param loader + low level eval functions
- Lightweight package of RNAstructure's param loader: the tables efn2 uses, with any number of modified bases at once
- (TODO) RNAstructure low level eval functions

## How to use
In your cmake file you include fetch content, and can switch flags on or off for which subpackages you want
In the example below, only ViennaRNA files are included

```cmake
include(FetchContent)
#------------------ ThermoRNA -----------------
set(THERMORNA_VIENNARNA ON CACHE BOOL "" FORCE)
set(THERMORNA_RNASTRUCTURE ON CACHE BOOL "" FORCE)  # optional, see below

FetchContent_Declare(ThermoRNA
    GIT_REPOSITORY https://github.com/TheCOBRALab/ThermoRNA.git
    GIT_TAG v0.2.0 # Set desired version
)
FetchContent_MakeAvailable(ThermoRNA)
```

Include packages within your program

```c++
#include <ViennaRNA/model.hpp>
#include <ViennaRNA/params/basic.hpp>
#include <ViennaRNA/params/io.hpp>

namespace viennarna = thermorna::viennarna;

viennarna::vrna_md_t md{};     ///< ViennaRNA model details.
viennarna::vrna_param_t* p{};  ///< ViennaRNA parameters.
viennarna::vrna_params_load_DNA_Mathews2004();

viennarna::vrna_md_set_default(&md);
p = viennarna::vrna_params(&md);
```

## RNAstructure parameters (efn2)

`THERMORNA_RNASTRUCTURE` builds a loader for the parameter tables efn2 reads
(`<prefix>.<table>.dg/.dh` and `<prefix>.specification.dat`). It needs uni-algo (header
only), which is fetched unless the parent project already provides `uni-algo::uni-algo`.

Each parameter set has a folder of its own. ThermoRNA ships RNAstructure 6.6's sets this way
(sources and details in [params/RNAstructure](params/RNAstructure/README.md)):

```
params/RNAstructure/
├── turner2004/   rna.specification.dat, rna.stack.dg, rna.stack.dh, ...
├── m6A/          m6A.specification.dat, m6A.stack.dg, ...
├── m1Ψ/          1.specification.dat, 1.stack.dg, ...
├── dna/  dmem/  Ψ/  Ψ_full/  m1Ψ_full/
└── pseudoknot/   pseudconst.dat
```

A set of your own is a folder with `specification.dat`, `stack.dg`, ... (a file prefix such as
`mine.stack.dg` is optional).

```c++
#include <RNAstructure/rnastructure.hpp>

namespace rs = thermorna::rnastructure;

rs::LoadOptions options;
options.data_directory = "params/RNAstructure/turner2004";    // folder of the base set
options.temperature = 310.15;                                 // efn2 --temperature (K)
options.modifications = {{"params/RNAstructure/m6A"},         // optional, any number
                         {"params/RNAstructure/Ψ"},
                         {"params/RNAstructure/m1Ψ"}};

const rs::Datatable data = rs::Datatable::load(options);
const std::vector<rs::Code> s = data.encode("GG6AAPUCC1");    // 6 = m6A, P = Ψ, 1 = m1Ψ

rs::Energy e = data.stack(s[0], s[9], s[1], s[8]);             // data->stack[..][..][..][..]
rs::Energy t = data.tstkh(s[1], s[8], s[2], s[7]);             // data->tstkh[..]...
rs::Energy m = data.misc().efn2a;                               // data->efn2a
```

A folder with one set needs no file prefix: it is read from the specification file's name,
and the folder's name becomes the set's name in `origin()` and messages ("m1Ψ" rather than
"1"). RNAstructure's own `data_tables` folder, which holds every set, works too; there a set is
picked by its prefix, like efn2's `--alphabet`:

```c++
options.data_directory = "/path/to/RNAstructure/data_tables";  // empty: $DATAPATH
options.alphabet = "rna";                                       // the default there
options.modifications = {{"", "m6A"}, {"", "PSU"}, {"", "1"}};  // {folder, prefix}
```

Values are in RNAstructure's units (0.1 kcal/mol, `rs::kInfiniteEnergy` = 14000 for forbidden
entries), and every accessor takes the same arguments as the corresponding array of
RNAstructure's `datatable`, so efn2's energy functions port line by line. Lookups are `const`
and thread-safe.

| efn2 (`datatable* data`, `structure* ct`)  | `rs::Datatable data`, `s = data.encode(seq)` |
| ------------------------------------------ | -------------------------------------------- |
| `ct->numseq[i]` (1-based)                  | `s[i]` (0-based)                             |
| `data->stack[a][b][c][d]`, `tstkh`, `tstki`, `tstki23`, `tstki1n`, `tstack`, `tstkm`, `coax`, `tstackcoax`, `coaxstack`, `helix_ends` | `data.stack(a, b, c, d)`, `data.tstkh(...)`, ... |
| `data->dangle[a][b][c][d]`                 | `data.dangle(a, b, c, d)`                    |
| `data->iloop11[a]..[f]`, `iloop21`, `iloop22` | `data.iloop11(a, ..., f)`, ...            |
| `data->hairpin[n]`, `bulge`, `inter`       | `data.hairpin()[n]`, ...                     |
| `data->eparam[k]`, `poppen`, `maxpen`, `efn2a/b/c`, `prelog`, `gubonus`, ... | `data.misc().eparam[k]`, ... |
| `data->RT`                                 | `data.RT()`                                  |
| tetraloop search over `data->tloop`        | `data.tloop(&s[i])` (`triloop`, `hexaloop`)  |
| `data->pairing[a][b]`                      | `data.can_pair(a, b)`                        |
| `ct->IsNuc(i, 'G')`                        | `data.has_symbol(s[i], "G")`                 |
| `ct->numseq[i] == 2` (C)                   | `s[i] == data.code("C")`                     |

### What is loaded

The same files, in the same order, with the same post-processing as `datatable::opendat()`:
`loop`, `dangle`, `stack`, `helix_ends`, `tstackh`, `tstacki`, `tstacki23`, `tstacki1n`,
`coaxial`, `tstackcoax`, `coaxstack`, `tstack`, `tstackm`, `int11`, `int21`, `int22`, `tloop`,
`triloop`, `hexaloop`, `miscloop`. At any temperature other than 310.15 K (±0.01 K, as in
RNAstructure) the `.dh` tables are read and every value is recomputed with RNAstructure's float
arithmetic (`dG_T`/`Tscale`). `PseudoknotConstants::read(folder)` reads `pseudconst.dat`, the
penalties efn2 uses for pseudoknotted structures.

`tools/rnastructure_parity` compares every entry against RNAstructure's own loader. For
RNAstructure 6.6 all shipped alphabets at 37 °C and `rna` between 273.15 K and 400 K are
identical, and alphabets without `.dh` files are refused at other temperatures, like efn2 does.

RNAstructure reads some malformed input leniently, and those cases are reproduced (and reported
in `warnings()`) rather than fixed: a row with a missing value reads that entry as 0
(`PSU.coaxstack.dg` and `1.coaxstack.dg` have 15 such rows), and a short line in `loop`,
`miscloop` or the special-hairpin files reuses the previous value. Input that makes
RNAstructure read out of bounds (unknown bases, truncated blocks) raises `rs::ParameterError`.
Two deliberate differences: RNAstructure's `dG_T` only rescales codes 0-5, which matters only
for custom alphabets with more than six bases and `.dh` files (all codes are rescaled here), and
pre-6.5 `miscloop` files with terminal AU lines are accepted (RNAstructure 6.6 misreads them).

### Several modified bases at once

efn2 supports one modification at a time because RNAstructure represents it as a separate
alphabet (`--alphabet m6A`) with complete tables. `Datatable` instead merges the extra bases of
any number of modification alphabets into one code space on top of a base alphabet. Each lookup
is answered from the alphabet that covers all of its nucleotides:

1. only base nucleotides: the base alphabet (so canonical contexts are exactly `rna`, even
   though `m6A.stack.dg` lists other values for some canonical GU stacks);
2. nucleotides a modification covers: that modification's tables, exactly as `efn2 -a <it>`
   would read them; if several cover the context, the last one listed wins;
3. modified bases from different alphabets in one context: each is replaced by its unmodified
   fallback and the base tables are used, the rule Knotergy already applies to modified bases
   without parameters for a context.

With no modifications a `Datatable` is efn2's datatable; with the m6A set as the base it is
`efn2 -a m6A`.

- **Fallbacks** are inferred from the Pairing sections (m6A → A, pseudouridine → U, m1Ψ → U)
  or set explicitly: `{"params/RNAstructure/m6A", "", {{"6", "A"}}}`.
- **Symbols** come from the specification files. A symbol that two modifications use for
  different bases is rejected in sequences with an explanation (RNAstructure uses `M` for both
  m6A and m1Ψ, so write `6` or `1`). Knotergy's codes `6`, `P` and `1` are symbols of the m6A,
  PSU and 1 alphabets.
- **Graphemes**: symbols are grapheme clusters everywhere (specification, table headers,
  sequences), so alphabets can use `Ψ`, and `Modification::aliases` adds symbols
  (`{"params/RNAstructure/Ψ", "", {}, {{"Ψ", "P"}}}`). uni-algo is only called for the spans of text that contain
  non-ASCII bytes; ASCII is handled byte by byte, and the result is identical to segmenting
  everything with uni-algo.

Memory is small because internal-loop tables store one dense block per closing-pair combination
that the files contain: `rna` takes about 0.25 MB and each modification alphabet 0.4-0.8 MB
(RNAstructure allocates N^6 to N^8 entries per table). Loading `rna` takes about 0.5 ms and each
modification alphabet 1-2 ms.

### Tests

```bash
cmake -S . -B build -DTHERMORNA_RNASTRUCTURE=ON
cmake --build build -j
ctest --test-dir build
```

The tests use the sets in `params/RNAstructure`. `-DCMAKE_BUILD_TYPE=Debug` builds them with
AddressSanitizer and UBSan. With `-DTHERMORNA_RNASTRUCTURE_SOURCE=/path/to/RNAstructure` (its
source tree), `RNAstructureParity` also compares every table entry with RNAstructure's own
loader (see `tools/rnastructure_parity`).
