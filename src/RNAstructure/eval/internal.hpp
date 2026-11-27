#pragma once

/**
 * @file RNAstructure/eval/internal.hpp
 * @brief efn2's energy of two pairs with nothing but unpaired nucleotides between them: an
 * internal loop, a bulge or a stack.
 *
 * Like every function of RNAstructure/eval, it takes encoded nucleotides (Datatable::code or
 * Datatable::encode) rather than a sequence and positions, returns tenths of kcal/mol, and
 * does not check its arguments.
 */

#include "RNAstructure/constants.hpp"
#include "RNAstructure/params/datatable.hpp"

#include <cstddef>

namespace thermorna::rnastructure {

/**
 * @brief The energy of an internal loop or bulge (efn2's erg2).
 *
 * @code
 *              si1 .. sp1
 *     5'  i1  i            p  p1
 *             |            |
 *     3'  j1  j            q  q1
 *              sj1 .. sq1
 * @endcode
 * The pairs are i-j outside and p-q inside, with n1 unpaired nucleotides between i and p and
 * n2 between q and j. The first ten arguments are those of ViennaRNA's vrna_E_internal, with
 * two codes for a pair instead of a pair type. The others are what efn2 looks at in addition.
 *
 * - A single bulged nucleotide: the two pairs stack across it, and the bulge gets a bonus of
 *   RT ln(n) for the n positions it can take among identical neighbors, plus the bonus of a C
 *   bulge next to another C.
 * - Longer bulges: initiation by size and the helix ends of both pairs.
 * - 1x1, 2x1 and 2x2 internal loops: their own tables.
 * - Other internal loops: initiation by size, the mismatches next to both pairs (with separate
 *   tables for 1xn and 2x3 loops) and the asymmetry penalty.
 * Internal loops also get the helix ends of both pairs. Sizes above 30 are extrapolated.
 * Without unpaired nucleotides the two pairs are a stack: the stack table entry is returned.
 *
 * @param n1 Unpaired nucleotides on the 5' side, between i and p.
 * @param n2 Unpaired nucleotides on the 3' side, between q and j.
 * @param i, j The outer pair.
 * @param p, q The inner pair.
 * @param si1 The nucleotide after i. Like the next three, it is only looked at when its side
 * has unpaired nucleotides.
 * @param sj1 The nucleotide before j.
 * @param sp1 The nucleotide before p.
 * @param sq1 The nucleotide after q.
 * @param beyond_i The nucleotide before i (i1), or kNoNucleotide. With beyond_j, it is what
 * the helix end of the outer pair depends on: that helix's next pair if it has one. When one
 * of the two does not exist, efn2 uses the nucleotides after i and before j instead.
 * @param beyond_j The nucleotide after j (j1), or kNoNucleotide.
 * @param beyond_p The nucleotide after p (p1), for the helix end of the inner pair. Like
 * beyond_q, it has to be a nucleotide; within one strand both always exist.
 * @param beyond_q The nucleotide before q (q1).
 * @param bulge_run For a single bulged nucleotide: how many identical nucleotides in a row it
 * is one of, paired or not (1 when its neighbors differ from it). Not used otherwise.
 * @param d The parameters.
 * @return The energy; kInfiniteEnergy when RNAstructure does not allow the loop: a pair it
 * cannot form, or a forbidden table entry (N in a 1x1, 2x1 or 2x2 loop, or beside a helix of
 * one pair).
 */
[[nodiscard]] int internal_loop_energy(std::size_t n1, std::size_t n2, Code i, Code j, Code p,
                                       Code q, Code si1, Code sj1, Code sp1, Code sq1,
                                       Code beyond_i, Code beyond_j, Code beyond_p, Code beyond_q,
                                       std::size_t bulge_run, const Datatable& d) noexcept;

}  // namespace thermorna::rnastructure
