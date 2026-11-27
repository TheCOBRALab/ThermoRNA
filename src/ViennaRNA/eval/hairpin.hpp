#pragma once

#include "ViennaRNA/params/basic.hpp"
#include "ViennaRNA/params/constants.hpp"
#include "ViennaRNA/params/salt.hpp"

#include <cstring>

#ifdef __GNUC__
#define INLINE inline
#else
#define INLINE
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
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as
 * there are six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is a1
 * while the 3' mismatch is a6. The nucleotide sequence of this loop is
 * &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *
 *  @warning  This function @b only evaluates the free energy of a hairpin loop according to the
 * current Turner energy parameter set! No additional hard- or soft constraints are applied. See
 * vrna_eval_hp_loop() for a function that also takes into account any user-supplied constraints!
 *
 *  @note     Whenever one of the mismatch base encodings @p si1 or @p sj1 is negative, terminal
 * mismatch energies are not applied!
 *
 *  @note     The parameter @p sequence is a 0-terminated string of size @p size + 2 that contain
 * the nucleic acid sequence of the loop in upper-case letters. This parameter is only required for
 * loops of @p size below 7, since it is used for look-up of unusually stable tri-, tetra- and
 * hexa-loops, such as GNRA tetra loops. Those may have additional sequence-dependent tabulated free
 * energies available.
 *
 *  @see      #vrna_param_t, vrna_eval_hp_loop()
 *
 *  @param  size      The size of the loop (number of unpaired nucleotides)
 *  @param  type      The pair type of the base pair closing the hairpin
 *  @param  si1       The 5'-mismatching nucleotide
 *  @param  sj1       The 3'-mismatching nucleotide
 *  @param  sequence  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size +
 * 2@f$ long)
 *  @param  P         The datastructure containing scaled energy parameters
 *  @return           The Free energy of the Hairpin-loop in dcal/mol
 */
// Defined here so that recursions can inline it into their innermost loops
// which makes the inlined code faster.
inline int vrna_E_hairpin(unsigned int size, unsigned int type, int si1, int sj1,
                          const char* sequence, vrna_param_t* P) {
    int energy, salt_correction;

    energy = INF;

    if (P) {
        salt_correction = 0;

        if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
            if (size <= MAXLOOP)
                salt_correction = P->SaltLoop[size + 1];
            else
                salt_correction =
                    vrna_salt_loop_int(size + 1, P->model_details.salt, P->temperature + K0,
                                       P->model_details.backbone_length);
        }

        if (size <= 30)
            energy = P->hairpin[size];
        else
            energy = P->hairpin[30] + (int) (P->lxc * log((size) / 30.));

        energy += salt_correction;

        if (size < 3) return energy; /* should only be the case when folding alignments */

        if ((sequence) && (P->model_details.special_hp)) {
            if (size == 4) {
                /* check for tetraloop bonus */
                char tl[7] = {0}, *ts;
                memcpy(tl, sequence, sizeof(char) * 6);
                tl[6] = '\0';
                if ((ts = strstr(P->Tetraloops, tl)))
                    return P->Tetraloop_E[(ts - P->Tetraloops) / 7] + salt_correction;
            } else if (size == 6) {
                char tl[9] = {0}, *ts;
                memcpy(tl, sequence, sizeof(char) * 8);
                tl[8] = '\0';
                if ((ts = strstr(P->Hexaloops, tl)))
                    return P->Hexaloop_E[(ts - P->Hexaloops) / 9] + salt_correction;
            } else if (size == 3) {
                char tl[6] = {0}, *ts;
                memcpy(tl, sequence, sizeof(char) * 5);
                tl[5] = '\0';
                if ((ts = strstr(P->Triloops, tl)))
                    return P->Triloop_E[(ts - P->Triloops) / 6] + salt_correction;

                return energy + (type > 2 ? P->TerminalAU : 0);
            }
        }

        if ((si1 >= 0) && (sj1 >= 0)) energy += P->mismatchH[type][si1][sj1];
    }

    return energy;
}

}  // namespace thermorna::viennarna