#pragma once

#include "ViennaRNA/params/basic.hpp"

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
int vrna_E_internal(unsigned int n1, unsigned int n2, unsigned int type, unsigned int type_2,
                    int si1, int sj1, int sp1, int sq1, vrna_param_t* P);
}  // namespace thermorna::viennarna