#pragma once

/**
 * @file RNAstructure/eval/exterior.hpp
 * @brief efn2's exterior loop energy, for a whole loop and for one of its helices.
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

/**
 * @brief The energy of one helix of an exterior loop with what stacks on its end, like
 * ViennaRNA's vrna_E_exterior_stem.
 *
 * @code
 *              n5  p  after_p
 *     loop         |             helix
 *              n3  q  before_q
 * @endcode
 * The end of the helix (efn2's penalty) plus the dangling end or the terminal mismatch that
 * @p stacking names (efn2's erg4, or its tstack table). In RNAstructure's model this is one
 * choice among several, and the exterior loop has the lowest energy over the choices of its
 * helices (see RNAstructure/eval/stacking.hpp). exterior_energy does that for a whole loop.
 *
 * @param p, q The first pair of the helix: the loop reaches it at p, its 5' nucleotide, and
 * leaves it at q.
 * @param after_p The nucleotide after p, inside the helix: with before_q, the next pair of
 * the helix if it has one. The end of the helix depends on them.
 * @param before_q The nucleotide before q.
 * @param n5 The nucleotide before p, in the loop, or kNoNucleotide.
 * @param n3 The nucleotide after q, in the loop, or kNoNucleotide.
 * @param stacking Which of n5 and n3 stack on the helix.
 * @param d The parameters.
 * @return The energy; kInfiniteEnergy when RNAstructure does not allow it: a pair it cannot
 * form, a forbidden table entry (N as a 5' dangling end, for example), or a nucleotide that is
 * to stack and is not there.
 *
 * When after_p or before_q does not exist (the pair is the last of its strand), efn2 looks at
 * the pair from the other side: the end of the helix then depends on n3 and n5.
 */
[[nodiscard]] int exterior_stem_energy(Code p, Code q, Code after_p, Code before_q, Code n5,
                                       Code n3, Stacking stacking, const Datatable& d) noexcept;

/**
 * @brief The energy of an exterior loop (efn2's ergexterior).
 *
 * The helix end of every helix, plus the lowest energy over the ways the helices and the
 * nucleotides next to them can stack. Going through the loop 5' to 3', a helix does at most
 * one of these, and a nucleotide is used at most once:
 * - take the nucleotide before it as a 5' dangle, or the one after it as a 3' dangle;
 * - take both as a terminal mismatch;
 * - stack coaxially on the next helix when no nucleotide separates them;
 * - stack coaxially on the next helix across the single nucleotide between them, which forms
 *   a mismatch with the nucleotide before the first helix or after the second.
 * In other words: the lowest sum of exterior_stem_energy over the helices and of the coaxial
 * stacks between them (RNAstructure/eval/stacking.hpp).
 *
 * With CoaxialStacking::Off the helices do not stack on each other: only the first two ways
 * are left, the lowest sum of exterior_stem_energy over the helices. That is not efn2's
 * energy. It is the exterior loop of a folding algorithm without coaxial stacking, like
 * RNAstructure's Fold --disablecoax, and the scheme of ViennaRNA's dangle model 1.
 *
 * @param loop The loop in 5' to 3' order: its unpaired nucleotides and its branches (see
 * LoopElement).
 * @param d The parameters.
 * @param coaxial Whether helices can stack on each other.
 * @return The energy; kInfiniteEnergy when RNAstructure does not allow the end of a helix (a
 * pair it cannot form, for example).
 */
[[nodiscard]] int exterior_energy(const std::vector<LoopElement>& loop, const Datatable& d,
                                  CoaxialStacking coaxial = CoaxialStacking::On);

}  // namespace thermorna::rnastructure
