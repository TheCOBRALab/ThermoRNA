# ThermoRNA (DO NOT MAKE THIS REPO PRIVATE)

This repository is REQUIRED for Knotergy to run. Setting it to private will prevent people from compiling it.

## What is ThermoRNA

It includes:

- Lightweight package of ViennaRNA param loader + low level eval functions
- Lightweight package of RNAstructure's param loader: the tables efn2 uses, with any number of modified bases at once
- RNAstructure low level eval functions: efn2's energy of every kind of loop and of single stems, from encoded nucleotides

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
    GIT_TAG v0.2.2 # Set desired version
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

| efn2 (`datatable* data`, `structure* ct`)                                                                                             | `rs::Datatable data`, `s = data.encode(seq)`     |
| ------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------ |
| `ct->numseq[i]` (1-based)                                                                                                             | `s[i]` (0-based)                                 |
| `data->stack[a][b][c][d]`, `tstkh`, `tstki`, `tstki23`, `tstki1n`, `tstack`, `tstkm`, `coax`, `tstackcoax`, `coaxstack`, `helix_ends` | `data.stack(a, b, c, d)`, `data.tstkh(...)`, ... |
| `data->dangle[a][b][c][d]`                                                                                                            | `data.dangle(a, b, c, d)`                        |
| `data->iloop11[a]..[f]`, `iloop21`, `iloop22`                                                                                         | `data.iloop11(a, ..., f)`, ...                   |
| `data->hairpin[n]`, `bulge`, `inter`                                                                                                  | `data.hairpin()[n]`, ...                         |
| `data->eparam[k]`, `poppen`, `maxpen`, `efn2a/b/c`, `prelog`, `gubonus`, ...                                                          | `data.misc().eparam[k]`, ...                     |
| `data->RT`                                                                                                                            | `data.RT()`                                      |
| tetraloop search over `data->tloop`                                                                                                   | `data.tloop(&s[i])` (`triloop`, `hexaloop`)      |
| `data->pairing[a][b]`                                                                                                                 | `data.can_pair(a, b)`                            |
| `ct->IsNuc(i, 'G')`                                                                                                                   | `data.has_symbol(s[i], "G")`                     |
| `ct->numseq[i] == 2` (C)                                                                                                              | `s[i] == data.code("C")`                         |

### Energy functions (efn2)

`RNAstructure/eval` has efn2's energy functions. Like ViennaRNA's low level functions
(`vrna_E_hairpin(size, type, si1, sj1, sequence, P)`), they take the encoded nucleotides a loop
is made of, not a sequence and positions. Where the nucleotides come from is up to the caller:
a sequence that is encoded once, a part of one, two strands.

```c++
//                                           0         1         2
//                                           0123456789012345678901234567
const std::vector<rs::Code> s = data.encode("GCGGACAAAAGUCCGGACAAAAGUCCGC");
//                                           ((((((....))))((((....))))))
using Element = rs::LoopElement;

// The stack of pairs 0-27 and 1-26 is a table entry.                  efn2 -w: -3.4
int stack = data.stack(s[0], s[27], s[1], s[26]);

// The hairpin closed by 5-10: its size, the nucleotides from 5 to 10, and those at 4, 3
// and 11.                                                                        4.1
int hairpin = rs::hairpin_energy(4, &s[5], s[4], s[3], s[11], data);

// The multibranch loop closed by 1-26: the pair, the nucleotides at 0 and 27, and what is
// in the loop, here two helices that start with 2-13 and 14-25.                  8.2
int multiloop = rs::multibranch_energy(s[1], s[26], s[0], s[27],
                                       {Element::helix(s[2], s[13], s[3], s[12]),
                                        Element::helix(s[14], s[25], s[15], s[24])},
                                       data);

// The exterior loop: one helix, which starts with 0-27.                          0.0
int exterior = rs::exterior_energy({Element::helix(s[0], s[27], s[1], s[26])}, data);
```

| efn2 (`src/structure.cpp`)      | here                                                                                                                         |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| `erg1`: stack                   | `data.stack(i, j, p, q)`                                                                                                     |
| `erg2`: internal loop and bulge | `rs::internal_loop_energy(n1, n2, i, j, p, q, si1, sj1, sp1, sq1, beyond_i, beyond_j, beyond_p, beyond_q, bulge_run, data)` |
| `erg3`: hairpin                 | `rs::hairpin_energy(size, loop, beyond_i, beyond_i2, beyond_j, data)`                                                        |
| `erg4`: dangling end            | `data.dangle(a, b, x, 1)` on the 3' side, `data.dangle(a, b, x, 2)` on the 5' side                                           |
| `penalty`: helix end            | `data.helix_ends(a, b, after_a, before_b)`                                                                                   |
| `ergexterior`                   | `rs::exterior_energy(loop, data, coaxial)`                                                                                   |
| `ergmulti`                      | `rs::multibranch_energy(i, j, beyond_i, beyond_j, loop, data, rule, coaxial)`                                                |

- Nucleotides are passed as codes (`rs::Code`, from `data.code("A")` or `data.encode(...)`).
  `size`, `n1`, `n2` and `bulge_run` are numbers of nucleotides. The first ten arguments of
  `internal_loop_energy` are those of `vrna_E_internal`, with two codes for a pair instead of
  a pair type.
- The `beyond_...` arguments are nucleotides next to the loop that efn2 looks at too: the end
  of a helix depends on its next pair, and a hairpin closed by G-U gets a bonus after two Gs.
  For the pair that closes a loop, `rs::kNoNucleotide` stands for one that does not exist,
  such as the one before the first nucleotide of a strand.
- Exterior and multibranch loops are given as a whole, as a list of `rs::LoopElement` in 5' to
  3' order, because efn2 picks the best way for helices and nucleotides to stack over the whole
  loop: `unpaired(x)`, `helix(p, q, after_p, before_q)`, and `other(p, q)` for a branch that
  is not a helix, such as a pseudoknot. Such a branch has no helix end and nothing stacks on
  it, which is a rule of ThermoRNA, not of efn2.
- Energies are `int`s in tenths of kcal/mol, like the tables. The arguments are not checked.
- A loop RNAstructure does not allow (a pair it cannot form, a hairpin of fewer than three
  nucleotides, N where a table has no entry) gives `rs::kInfiniteEnergy`. efn2 itself adds such
  entries up like numbers and prints energies around 1400 kcal/mol.
- `rule` is `rs::MultibranchRule::Efn2` (the default) or `Simple`, which is `efn2 --simple`.
- `coaxial` is `rs::CoaxialStacking::On` (the default, as in efn2) or `Off`, which leaves out
  the helices that stack on each other (see "Without coaxial stacking" below).
- Not included: what efn2 adds for a structure of two strands, such as the intermolecular
  initiation (`data.misc().init`).

The sum over the loops of a structure without pseudoknots is efn2's total energy. `efn2()` in
`test/test_rnastructure_eval.cpp` does this for a dot-bracket structure, and shows how to
gather each function's arguments from an encoded sequence.

#### Exterior and multibranch loops, stem by stem

In RNAstructure's model the end of a helix does at most one of these: nothing, stack with the
nucleotide before it (5' dangle), with the one after it (3' dangle), with both (terminal
mismatch), or stack on the helix next to it (coaxial stack). A nucleotide stacks only once. So
these loops are not a plain sum over their stems, as they are in ViennaRNA with `-d2`: the
energy is the lowest sum over the choices. `exterior_energy` and `multibranch_energy` make the
choice for a whole loop. A folding algorithm makes it in its recursion, with the energy of one
stem and one choice at a time:

| `rs::`                                                                               | what                                                              |
| ------------------------------------------------------------------------------------ | ----------------------------------------------------------------- |
| `exterior_stem_energy(p, q, after_p, before_q, n5, n3, stacking, data)`              | helix end + dangle or mismatch, like `vrna_E_exterior_stem`       |
| `multibranch_stem_energy(p, q, after_p, before_q, n5, n3, stacking, data)`           | the same + the penalty per helix, like `vrna_E_multibranch_stem`  |
| `coaxial_stack_energy(p1, q1, p2, q2, data)`                                         | two helices next to each other                                    |
| `coaxial_stack_mismatch5_energy(p1, q1, p2, q2, before, between, data)`              | two helices one nucleotide apart, mismatch before the first       |
| `coaxial_stack_mismatch3_energy(p1, q1, p2, q2, between, after, data)`               | the same, mismatch after the second                               |
| `data.misc().efn2a`, `data.misc().efn2b`                                             | multibranch loop: offset, and penalty per unpaired nucleotide     |

- A helix is given as the loop sees it, going through the loop 5' to 3': the loop reaches it at
  `p` and leaves it at `q`. `after_p` and `before_q` are on the side of the helix (its next
  pair if it has one), `n5` and `n3` are the nucleotides before `p` and after `q` in the loop.
  For a branch, `p` is the 5' nucleotide of its first pair. For the pair `i`-`j` that closes a
  multibranch loop, `p` is `j` and `q` is `i`.
- `stacking` says which of `n5` and `n3` stack on the helix: `rs::Stacking::None`, `Dangle5`,
  `Dangle3` or `Mismatch`. Two helices that stack coaxially have `None`.
- The lowest sum of these is efn2's energy of an exterior loop. For a multibranch loop it is
  the energy by the simple rule, the one RNAstructure's folding algorithms use. efn2's default
  adds three terms that depend on the loop as a whole (a logarithm beyond 8 unpaired
  nucleotides, an asymmetry penalty and a strain penalty), so it only exists for a whole loop.

```c++
// The multibranch loop of the example above by the simple rule, from its stems. Seen from
// inside, the closing pair 1-26 is a helix that the loop reaches at 26 and leaves at 1.
using rs::Stacking;
int closing = rs::multibranch_stem_energy(s[26], s[1], s[27], s[0], s[25], s[2], Stacking::None, data);
int first   = rs::multibranch_stem_energy(s[2], s[13], s[3], s[12], s[1], s[14], Stacking::None, data);
int second  = rs::multibranch_stem_energy(s[14], s[25], s[15], s[24], s[13], s[26], Stacking::None, data);
int stacked = rs::coaxial_stack_energy(s[2], s[13], s[14], s[25], data);  // the two branches
int simple  = data.misc().efn2a + closing + first + second + stacked;     // efn2 -w -s: 5.1
// = rs::multibranch_energy(s[1], s[26], s[0], s[27], {...}, data, rs::MultibranchRule::Simple)

// Without coaxial stacking the two branches do not stack, and there is no unpaired
// nucleotide to stack on a helix instead:                                            7.5
int without = data.misc().efn2a + closing + first + second;
// = rs::multibranch_energy(s[1], s[26], s[0], s[27], {...}, data, rs::MultibranchRule::Simple,
//                          rs::CoaxialStacking::Off)
```

**Without coaxial stacking.** With `rs::CoaxialStacking::Off`, `exterior_energy` and
`multibranch_energy` leave the last of the choices out: helices do not stack on each other,
and each helix takes nothing, a dangling end or a terminal mismatch. That is the scheme of
ViennaRNA's dangle model 1 (`-d1`), with RNAstructure's parameters. With
`rs::MultibranchRule::Simple` it is the energy RNAstructure's own folding minimizes when its
coaxial stacking is disabled (`Fold --disablecoax`), and the energy a recursion finds that
uses the two stem functions and no coaxial stack. efn2 itself always considers coaxial stacks,
so it has no option that prints this energy.

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
   though the m6A, Ψ and m1Ψ files list other values for some GU stacks and helix ends);
2. nucleotides a modification covers: that modification's tables, exactly as `efn2 -a <it>`
   would read them; if several cover the context, the last one listed wins;
3. modified bases from different alphabets in one context: each is replaced by its unmodified
   fallback and the base tables are used, the rule Knotergy already applies to modified bases
   without parameters for a context.

With no modifications a `Datatable` is efn2's datatable. A modified-base set can be the base
too: with the m6A set as the base it is `efn2 -a m6A`, m6A is then one of the base bases (every
context, unmodified ones included, comes from m6A's files), and contexts that mix m6A with an
added modification keep m6A and only replace the added one.

- **Fallbacks** are inferred from the Pairing sections (m6A → A, pseudouridine → U, m1Ψ → U)
  or set explicitly: `{"params/RNAstructure/m6A", "", {{"6", "A"}}}`.
- **Symbols** come from the specification files. Bases are the same when their primary
  symbols are (or each is a symbol of the other, like U and T); other shared symbols don't
  count, since RNAstructure writes both m6A and m1Ψ as `M`. A symbol that two added sets use
  for different bases is rejected in sequences with an explanation (write `6` or `1`); one the
  base set already uses keeps the base set's meaning. Knotergy's codes `6`, `P` and `1` are
  symbols of the m6A, PSU and 1 alphabets.
- **Graphemes**: symbols are grapheme clusters everywhere (specification, table headers,
  sequences), so alphabets can use `Ψ`, and `Modification::aliases` adds symbols
  (`{"params/RNAstructure/Ψ", "", {}, {{"Ψ", "P"}}}`). uni-algo is only called for the spans of
  text that contain non-ASCII bytes; ASCII is handled byte by byte, and the result is identical
  to segmenting everything with uni-algo.

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

The tests use the sets in `params/RNAstructure`. The energy functions are compared with what
`efn2 -w` prints for each loop of about 170 structures. `-DCMAKE_BUILD_TYPE=Debug` builds them
with AddressSanitizer and UBSan. With `-DTHERMORNA_RNASTRUCTURE_SOURCE=/path/to/RNAstructure` (its
source tree), `RNAstructureParity` also compares every table entry with RNAstructure's own
loader (see `tools/rnastructure_parity`).
