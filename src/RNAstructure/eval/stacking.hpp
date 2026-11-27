#pragma once

/**
 * @file RNAstructure/eval/stacking.hpp
 * @brief What can stack on the end of a helix in an exterior or multibranch loop: the unpaired
 * nucleotides next to it, or the helix next to it.
 *
 * In RNAstructure's model the end of a helix does at most one of these, and an unpaired
 * nucleotide stacks at most once. So the energy of a loop is the lowest over the ways its
 * helices can stack, not a sum of independent terms. exterior_stem_energy,
 * multibranch_stem_energy and the functions here give the energy of one such choice, for a
 * program that makes the choice itself, like the recursion of a folding algorithm.
 * exterior_energy and multibranch_energy make it for a whole loop. With CoaxialStacking::Off
 * they leave out the helices that stack on each other, which is the model of a folding
 * algorithm without coaxial stacking.
 *
 * A helix is given as the loop sees it, going through the loop 5' to 3': by the pair p-q where
 * the loop reaches the helix at p and leaves it at q. For a branch of the loop, p is the 5'
 * nucleotide of its first pair. For the pair i-j that closes a multibranch loop, seen from
 * inside the loop, p is j and q is i.
 *
 * Like every function of RNAstructure/eval, these take encoded nucleotides (Datatable::code or
 * Datatable::encode), return tenths of kcal/mol, and do not check their arguments.
 */

#include "RNAstructure/constants.hpp"
#include "RNAstructure/params/datatable.hpp"

namespace thermorna::rnastructure {

/// Which unpaired nucleotides stack on the end of a helix.
enum class Stacking {
    None,      ///< none
    Dangle5,   ///< the nucleotide before the helix: a 5' dangling end
    Dangle3,   ///< the nucleotide after the helix: a 3' dangling end
    Mismatch,  ///< both: a terminal mismatch
};

/// Whether the helices of an exterior or multibranch loop can stack on each other, for
/// exterior_energy and multibranch_energy.
enum class CoaxialStacking {
    /// efn2's model: the end of a helix stacks on the nucleotides next to it or on the helix
    /// next to it (the coaxial_stack functions below).
    On,
    /// Only unpaired nucleotides stack on the helices, as dangling ends and terminal
    /// mismatches. What is left is the lowest sum of stem energies, each nucleotide stacking at
    /// most once: the scheme of ViennaRNA's dangle model 1, and with MultibranchRule::Simple
    /// what RNAstructure's folding minimizes when its coaxial stacking is disabled (Fold
    /// --disablecoax).
    Off,
};

/**
 * @brief Two helices of a loop that stack end to end with no nucleotide between them: a flush
 * coaxial stack (RNAstructure's ergcoaxflushbases).
 *
 * @code
 *   5'  ..  p1 == q1  p2 == q2  ..  3'
 * @endcode
 * This is the stacking alone. Each of the two helices also has its stem energy, with
 * Stacking::None.
 *
 * @param p1, q1 The first helix the loop meets, as the loop sees it.
 * @param p2, q2 The second one: p2 comes right after q1.
 * @param d The parameters.
 * @return The energy; kInfiniteEnergy where RNAstructure's table has none (a pair it cannot
 * form).
 */
[[nodiscard]] int coaxial_stack_energy(Code p1, Code q1, Code p2, Code q2,
                                       const Datatable& d) noexcept;

/**
 * @brief Two helices that stack across the single nucleotide between them, which forms a
 * mismatch with the nucleotide before the first helix (RNAstructure's ergcoaxinterbases1).
 *
 * @code
 *   5'  ..  before  p1 == q1  between  p2 == q2  ..  3'
 * @endcode
 * This is the stacking alone, as for coaxial_stack_energy. Both nucleotides are used up: they
 * cannot stack on anything else.
 *
 * The lookup is efn2's, in ergexterior and ergmulti. RNAstructure's ergcoaxinterbases1, which
 * its folding algorithms call, reads the coaxstack table the other way round
 * ([between][before][p2][q2] for [q2][p2][before][between]). That is another value for 16 of
 * the 1296 entries of the rna set, all of them with G and U.
 *
 * @param p1, q1 The first helix the loop meets, as the loop sees it.
 * @param p2, q2 The second one.
 * @param before The unpaired nucleotide before p1.
 * @param between The unpaired nucleotide between q1 and p2.
 * @param d The parameters.
 * @return The energy; kInfiniteEnergy where RNAstructure's tables have none.
 */
[[nodiscard]] int coaxial_stack_mismatch5_energy(Code p1, Code q1, Code p2, Code q2, Code before,
                                                 Code between, const Datatable& d) noexcept;

/**
 * @brief Two helices that stack across the single nucleotide between them, which forms a
 * mismatch with the nucleotide after the second helix (RNAstructure's ergcoaxinterbases2).
 *
 * @code
 *   5'  ..  p1 == q1  between  p2 == q2  after  ..  3'
 * @endcode
 * See coaxial_stack_mismatch5_energy.
 *
 * @param p1, q1 The first helix the loop meets, as the loop sees it.
 * @param p2, q2 The second one.
 * @param between The unpaired nucleotide between q1 and p2.
 * @param after The unpaired nucleotide after q2.
 * @param d The parameters.
 * @return The energy; kInfiniteEnergy where RNAstructure's tables have none.
 */
[[nodiscard]] int coaxial_stack_mismatch3_energy(Code p1, Code q1, Code p2, Code q2, Code between,
                                                 Code after, const Datatable& d) noexcept;

}  // namespace thermorna::rnastructure
