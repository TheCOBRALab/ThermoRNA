#pragma once

/**
 * @file RNAstructure/eval/multibranch.hpp
 * @brief efn2's multibranch loop energy, for a whole loop and for one of its helices.
 *
 * Like every function of RNAstructure/eval, these take encoded nucleotides (Datatable::code or
 * Datatable::encode) rather than a sequence and positions, return tenths of kcal/mol, and do
 * not check their arguments.
 */

#include "RNAstructure/constants.hpp"
#include "RNAstructure/eval/loop_element.hpp"
#include "RNAstructure/eval/stacking.hpp"
#include "RNAstructure/params/datatable.hpp"

#include <vector>

namespace thermorna::rnastructure {

/// Which of efn2's two energies a multibranch loop gets.
enum class MultibranchRule {
    /// efn2's default. It has terms that depend on the loop as a whole (a logarithm, an
    /// asymmetry and a strain penalty), so it is not a sum over the helices.
    Efn2,
    /// efn2 --simple, the rule of RNAstructure's folding algorithms: an offset, a penalty per
    /// unpaired nucleotide, and multibranch_stem_energy for every helix.
    Simple,
};

/**
 * @brief The energy of one helix of a multibranch loop with what stacks on its end, like
 * ViennaRNA's vrna_E_multibranch_stem.
 *
 * @code
 *              n5  p  after_p
 *     loop         |             helix
 *              n3  q  before_q
 * @endcode
 * The penalty per helix (misc().efn2c), the end of the helix (efn2's penalty), and the
 * dangling end or the terminal mismatch that @p stacking names (efn2's erg4, or its tstkm
 * table). In RNAstructure's model this is one choice among several (see
 * RNAstructure/eval/stacking.hpp). With MultibranchRule::Simple the energy of a multibranch
 * loop is
 * @code
 *   misc().efn2a + misc().efn2b * unpaired nucleotides
 *   + the lowest, over the ways the helices can stack, of the sum of multibranch_stem_energy
 *     over the helices and of the coaxial stacks between them
 * @endcode
 * where the pair that closes the loop is one of the helices. multibranch_energy does that for
 * a whole loop. RNAstructure's folding algorithms read the three penalties from
 * misc().eparam[5], [6] and [10], which have the values of efn2a, efn2b and efn2c in every set
 * shipped with RNAstructure 6.6.
 *
 * @param p, q The pair as the loop sees it: the loop reaches it at p and leaves it at q. For
 * a branch it is the first pair of the helix, and p is its 5' nucleotide. For the pair i-j
 * that closes the loop, p is j and q is i.
 * @param after_p The nucleotide after p, on the side of the helix: with before_q, the next
 * pair of the helix if it has one. The end of the helix depends on them. For the closing pair
 * they are the nucleotides after j and before i.
 * @param before_q The nucleotide before q.
 * @param n5 The nucleotide before p, in the loop, or kNoNucleotide.
 * @param n3 The nucleotide after q, in the loop, or kNoNucleotide.
 * @param stacking Which of n5 and n3 stack on the helix.
 * @param d The parameters.
 * @return The energy; kInfiniteEnergy when RNAstructure does not allow it: a pair it cannot
 * form, a forbidden table entry (N as a 5' dangling end, for example), or a nucleotide that is
 * to stack and is not there.
 *
 * When after_p or before_q does not exist (a closing pair at the end of the strand), efn2
 * looks at the pair from the other side: the end of the helix then depends on n3 and n5.
 */
[[nodiscard]] int multibranch_stem_energy(Code p, Code q, Code after_p, Code before_q, Code n5,
                                          Code n3, Stacking stacking, const Datatable& d) noexcept;

/**
 * @brief The energy of the multibranch loop closed by pair i-j (efn2's ergmulti).
 *
 * The loop is a circle of helices with unpaired nucleotides between them: the pair that closes
 * it, seen from inside, and its branches. Its energy is the sum of
 * - an offset, a penalty per helix and a penalty per unpaired nucleotide, which is logarithmic
 *   beyond 8 nucleotides;
 * - an asymmetry penalty, for how much the numbers of unpaired nucleotides on the two sides of
 *   a helix differ, averaged over the helices;
 * - a strain penalty when an odd number of helices has at most one unpaired nucleotide among
 *   them (a three-way junction without unpaired nucleotides, for example);
 * - the helix end of every helix;
 * - the lowest energy over the ways the helices and the nucleotides next to them can stack, as
 *   in an exterior loop (see exterior_energy), with the mismatches of multibranch loops.
 *
 * That is efn2's default, MultibranchRule::Efn2. MultibranchRule::Simple is efn2 --simple,
 * which leaves out the logarithm, the asymmetry and the strain like RNAstructure's folding
 * algorithms do: what is left is a sum over the helices (see multibranch_stem_energy).
 *
 * With CoaxialStacking::Off the helices do not stack on each other, whichever the rule. With
 * MultibranchRule::Simple that leaves the offset, the penalty per unpaired nucleotide and the
 * lowest sum of multibranch_stem_energy over the helices: the multibranch loop of a folding
 * algorithm without coaxial stacking, like RNAstructure's Fold --disablecoax, and the scheme
 * of ViennaRNA's dangle model 1. It is not an energy efn2 calculates.
 *
 * @param i, j The closing pair.
 * @param beyond_i The nucleotide before i, or kNoNucleotide. With beyond_j, it is what the
 * helix end of the closing pair depends on: the next pair of its helix if it has one. When one
 * of the two does not exist, efn2 uses the nucleotides after i and before j instead.
 * @param beyond_j The nucleotide after j, or kNoNucleotide.
 * @param loop What is between i and j, in 5' to 3' order: unpaired nucleotides and at least
 * one branch (see LoopElement).
 * @param d The parameters.
 * @param rule efn2's default energy, or the simple one.
 * @param coaxial Whether helices can stack on each other.
 * @return The energy; kInfiniteEnergy when RNAstructure does not allow the end of one of the
 * helices (a pair it cannot form, for example).
 */
[[nodiscard]] int multibranch_energy(Code i, Code j, Code beyond_i, Code beyond_j,
                                     const std::vector<LoopElement>& loop, const Datatable& d,
                                     MultibranchRule rule = MultibranchRule::Efn2,
                                     CoaxialStacking coaxial = CoaxialStacking::On);

}  // namespace thermorna::rnastructure
