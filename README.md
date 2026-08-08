# ThermoRNA (DO NOT MAKE THIS REPO PRIVATE)

This repository is REQUIRED for Knotergy to run. Setting it to private will prevent people from compiling it.

## What is ThermoRNA
It includes:
- Lightweight package of ViennaRNA param loader + low level eval functions
- Lightweight package of Knotergy pseudoknot param loader
- (TODO) Lightweight package of RNAstructure's param loader + low level eval functions

## How to use
In your cmake file you include fetch content, and can switch flags on or off for which subpackages you want
In the example below, only ViennaRNA files are included

```cmake
include(FetchContent)
#------------------ ThermoRNA -----------------
set(THERMORNA_KNOTERGY OFF CACHE BOOL "" FORCE)
set(THERMORNA_VIENNARNA ON CACHE BOOL "" FORCE)

FetchContent_Declare(ThermoRNA
    GIT_REPOSITORY https://github.com/TheCOBRALab/ThermoRNA.git
    GIT_TAG v0.1.1 # Set desired version
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