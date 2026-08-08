/*
 * WBL 24 Aug 2018 Add AVX512 based on sources_034_578/modular_decomposition_id3.c
 * WBL 22 Aug 2018 by hand d3c17fd3e04e2419c147a1e097d3c4d2c5a6f11d lines 1355-1357
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/params/default.hpp"
#include "ViennaRNA/eval/multibranch.hpp"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

namespace thermorna::viennarna {
PUBLIC int vrna_E_multibranch_stem(unsigned int  type,
                        int           si1,
                        int           sj1,
                        vrna_param_t  *P)
{
  int energy = INF;

  if (P) {
    energy = P->MLintern[type];

    if (si1 >= 0 && sj1 >= 0)
      energy += P->mismatchM[type][si1][sj1];
    else if (si1 >= 0)
      energy += P->dangle5[type][si1];
    else if (sj1 >= 0)
      energy += P->dangle3[type][sj1];

    if (type > 2)
      energy += P->TerminalAU;
  }

  return energy;
}

}
