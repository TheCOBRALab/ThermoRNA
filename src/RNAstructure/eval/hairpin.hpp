#pragma once

/**
 * @file RNAstructure/eval/hairpin.hpp
 * @brief efn2's hairpin loop energy.
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
 * @brief The energy of a hairpin loop (efn2's erg3).
 *
 * @code
 *        a2 a3
 *      a1     a4
 *        i - j
 *       i1   j1
 *       i2
 *       5'   3'
 * @endcode
 * Covers the loop initiation (extrapolated beyond 30 nucleotides), the first mismatch, special
 * tri-, tetra- and hexaloops, the bonus of a G-U closure after two Gs, all-C loops, and the
 * helix end of the closing pair, which efn2 counts as part of the hairpin.
 *
 * @param size The number of unpaired nucleotides (4 in the drawing).
 * @param loop The loop from i to j: size + 2 codes (i a1 a2 a3 a4 j).
 * @param beyond_i The nucleotide before i (i1), or kNoNucleotide.
 * @param beyond_i2 The nucleotide before that one (i2), or kNoNucleotide.
 * @param beyond_j The nucleotide after j (j1), or kNoNucleotide.
 * @param d The parameters.
 * @return The energy; kInfiniteEnergy when RNAstructure does not allow the loop: fewer than 3
 * nucleotides, a closing pair it cannot form, or a forbidden table entry (N beside a closing
 * helix of one pair).
 *
 * The helix end depends on i1 and j1: the next pair of the closing helix if it has one. When
 * one of them does not exist, efn2 uses a1 and a4 instead.
 */
[[nodiscard]] int hairpin_energy(std::size_t size, const Code* loop, Code beyond_i, Code beyond_i2,
                                 Code beyond_j, const Datatable& d) noexcept;

}  // namespace thermorna::rnastructure
