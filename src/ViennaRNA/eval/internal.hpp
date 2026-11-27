#pragma once

#include "ViennaRNA/params/basic.hpp"
#include "ViennaRNA/params/default.hpp"
#include "ViennaRNA/params/salt.hpp"

#include <algorithm>

#ifdef __GNUC__
#define INLINE inline
#else
#define INLINE
#endif
namespace thermorna::viennarna {
/**
 *  @brief  Compute the Energy of an internal loop
 *
 *  This function computes the free energy @f$ E @f$ of an internal-loop with the
 *  following structure: <BR>
 *  <PRE>
 *        3'  5'
 *        |   |
 *        U - V
 *    a_n       b_1
 *     .        .
 *     .        .
 *     .        .
 *    a_1       b_m
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  This general structure depicts an internal-loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop. In this example, the length of the internal-loop is @f$(n+m)@f$
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:<BR>
 *  5'-mismatch: a_1<BR>
 *  3'-mismatch: b_m<BR>
 *  and for the enclosed base pair (V,U):<BR>
 *  5'-mismatch: b_1<BR>
 *  3'-mismatch: a_n<BR>
 *
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
 *        must be 'turned arround' when evaluating the free energy of the internal-loop<br>
 *        This function is threadsafe
 *
 *  @see vrna_exp_E_internal()
 *
 *  @param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the internal loop
 *  @param  type_2  The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the internal loop in dcal/mol
 */
// Defined here so that recursions can inline it into their innermost loops
// which makes the inlined code faster.
inline int vrna_E_internal(unsigned int n1, unsigned int n2, unsigned int type, unsigned int type_2,
                           int si1, int sj1, int sp1, int sq1, vrna_param_t* P) {
    /* compute energy of degree 2 loop (stack bulge or internal) */
    unsigned int nl, ns, u, backbones;
    int energy;

    if (!P) return INF;

    if (n1 > n2) {
        nl = n1;
        ns = n2;
    } else {
        nl = n2;
        ns = n1;
    }

    if (nl == 0) /* stack */
        return P->stack[type][type_2] + P->SaltStack;

    energy = 0;

    switch (ns) {
        case 0:
            /* bulge */
            energy =
                (nl <= MAXLOOP) ? P->bulge[nl] : (P->bulge[30] + (int) (P->lxc * log(nl / 30.)));
            if (nl == 1) {
                /* add stacking energy for 1-bulges */
                energy += P->stack[type][type_2];
            } else {
                if (type > 2) energy += P->TerminalAU;

                if (type_2 > 2) energy += P->TerminalAU;
            }

            break;

        case 1:
            if (nl == 1) {
                /* 1x1 loop */
                energy = P->int11[type][type_2][si1][sj1];
            } else if (nl == 2) {
                /* 2x1 loop */
                if (n1 == 1)
                    energy = P->int21[type][type_2][si1][sq1][sj1];
                else
                    energy = P->int21[type_2][type][sq1][si1][sp1];
            } else {
                /* 1xn loop */
                energy = (nl + 1 <= MAXLOOP)
                             ? (P->internal_loop[nl + 1])
                             : (P->internal_loop[30] + (int) (P->lxc * log((nl + 1) / 30.)));
                energy += MIN2(MAX_NINIO, (int) (nl - ns) * P->ninio[2]);
                energy += P->mismatch1nI[type][si1][sj1] + P->mismatch1nI[type_2][sq1][sp1];
            }

            break;

        case 2:
            if (nl == 2) {
                /* 2x2 loop */
                energy = P->int22[type][type_2][si1][sp1][sq1][sj1];
                break;
            } else if (nl == 3) {
                /* 2x3 loop */
                energy = P->internal_loop[5] + P->ninio[2];
                energy += P->mismatch23I[type][si1][sj1] + P->mismatch23I[type_2][sq1][sp1];
                break;
            }

            /* fall through */

        default:
            /* generic internal loop (no else here!)*/
            u = nl + ns;
            energy += (u <= MAXLOOP) ? (P->internal_loop[u])
                                     : (P->internal_loop[30] + (int) (P->lxc * log((u) / 30.)));

            energy += std::min(MAX_NINIO, (int) (nl - ns) * P->ninio[2]);

            energy += P->mismatchI[type][si1][sj1] + P->mismatchI[type_2][sq1][sp1];

            break;
    }

    if ((P->model_details.noGUclosure) &&
        ((type_2 == 3) || (type_2 == 4) || (type == 3) || (type == 4)))
        return INF;

    if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
        /* salt correction for loop */
        backbones = nl + ns + 2;
        if (backbones <= MAXLOOP + 1)
            energy += P->SaltLoop[backbones];
        else
            energy += vrna_salt_loop_int(backbones, P->model_details.salt, P->temperature + K0,
                                         P->model_details.backbone_length);
    }

    return energy;
}
}  // namespace thermorna::viennarna