#pragma once

#include "ViennaRNA/loading_params/params/basic.hpp"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

namespace thermorna::viennarna {

/**
 *  @brief Retrieve the energy of a hairpin-loop
 *
 *  To evaluate the free energy of a hairpin-loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:<BR>
 *  <PRE>
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *
 *  @warning  This function @b only evaluates the free energy of a hairpin loop according to the current
 *            Turner energy parameter set! No additional hard- or soft constraints are applied. See vrna_eval_hp_loop()
 *            for a function that also takes into account any user-supplied constraints!
 *
 *  @note     Whenever one of the mismatch base encodings @p si1 or @p sj1 is negative, terminal mismatch energies are not applied!
 *
 *  @note     The parameter @p sequence is a 0-terminated string of size @p size + 2 that contain the nucleic acid
 *            sequence of the loop in upper-case letters. This parameter is only required for loops of @p size below
 *            7, since it is used for look-up of unusually stable tri-, tetra- and hexa-loops, such as GNRA tetra loops.
 *            Those may have additional sequence-dependent tabulated free energies available.
 *
 *  @see      #vrna_param_t, vrna_eval_hp_loop()
 *
 *  @param  size      The size of the loop (number of unpaired nucleotides)
 *  @param  type      The pair type of the base pair closing the hairpin
 *  @param  si1       The 5'-mismatching nucleotide
 *  @param  sj1       The 3'-mismatching nucleotide
 *  @param  sequence  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P         The datastructure containing scaled energy parameters
 *  @return           The Free energy of the Hairpin-loop in dcal/mol
 */
int
vrna_E_hairpin(unsigned int size,
               unsigned int type,
               int          si1,
               int          sj1,
               const char   *sequence,
               vrna_param_t *P);

}