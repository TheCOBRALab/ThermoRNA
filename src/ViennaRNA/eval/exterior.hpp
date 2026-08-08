#pragma once
#include "ViennaRNA/params/basic.hpp"

namespace thermorna::viennarna {
/**
 *  @brief  Evaluate a stem branching off the exterior loop
 *
 *  Given a base pair @f$(i,j)@f$ encoded by @em type, compute the energy contribution
 *  including dangling-end/terminal-mismatch contributions. Instead of returning the
 *  energy contribution per-se, this function returns the corresponding Boltzmann factor.
 *  If either of the adjacent nucleotides @f$(i - 1)@f$ and @f$(j+1)@f$ must not
 *  contribute stacking energy, the corresponding encoding must be @f$-1@f$.
 *
 *  @note By default, terminal mismatch energies are applied that correspond to the
 *        neighboring nucleotides provided by their encodings @p n5d and @p n3d. Whenever
 *        the encodings are negative, the implementation switches to usage of dangling
 *        end energies (for the non-negative base). If both encodings are negative, no
 *        terminal mismatch contributions are added.
 *
 *  @see vrna_exp_E_exterior_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  n5d   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be
 * -1)
 *  @param  n3d   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be
 * -1)
 *  @param  p     The pre-computed energy parameters
 *  @return       The energy contribution of the introduced exterior-loop stem
 */
int vrna_E_exterior_stem(unsigned int type, int n5d, int n3d, vrna_param_t* p);

}  // namespace thermorna::viennarna