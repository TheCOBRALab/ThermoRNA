#pragma once

#include "ViennaRNA/loading_params/params/basic.hpp"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

namespace thermorna::viennarna {
/**
 *  @brief  Evaluate the free energy contribution of a stem branching off a multibranch loop
 *
 *  This function yields the free energy contribution for the terminal base
 *  pairs of a stem branching off a multibranch loop. In essence, this consists of
 *  (i) a terminal mismatch or dangling end contribution, (ii) the score for a stem
 *  according to the affine multibranch loop model, and (iii) a terminal AU/GU penalty,
 *  if applicable.
 *
 *  @note By default, terminal mismatch energies are applied that correspond to the
 *        neighboring nucleotides provided by their encodings @p n5d and @p n3d. Whenever
 *        the encodings are negative, the implementation switches to usage of dangling
 *        end energies (for the non-negative base). If both encodings are negative, no
 *        terminal mismatch contributions are added.
 *
 *  @see vrna_exp_E_multibranch_stem(), vrna_E_exterior_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  si1   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  sj1   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  P     The pre-computed energy parameters
 *  @return       The energy contribution of the introduced mutlibranch loop stem
 */
int
vrna_E_multibranch_stem(unsigned int  type,
                        int           si1,
                        int           sj1,
                        vrna_param_t  *P);
} // namespace thermorna::viennarna