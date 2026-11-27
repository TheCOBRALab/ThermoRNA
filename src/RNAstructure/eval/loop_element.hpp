#pragma once

/**
 * @file RNAstructure/eval/loop_element.hpp
 * @brief How an exterior or multibranch loop is given to the energy functions: as the list of
 * what one meets going through it 5' to 3'.
 */

#include "RNAstructure/constants.hpp"

namespace thermorna::rnastructure {

/**
 * @brief One element of an exterior or multibranch loop: an unpaired nucleotide or a branch.
 *
 * efn2 scores these loops as a whole, because which nucleotide stacks on which helix, and
 * which helices stack on each other, is a choice over the whole loop. So the loop is given as
 * its elements in 5' to 3' order, each made of encoded values. A loop that goes through
 * nucleotide a, the helix that starts with pair p-q, nucleotide b and the helix that starts
 * with pair r-t is
 * @code
 *   {unpaired(a), helix(p, q, p1, q1), unpaired(b), helix(r, t, r1, t1)}
 * @endcode
 * where p is the 5' nucleotide of its pair, and p1 and q1 are the nucleotides after p and
 * before q, inside the helix: its next pair if it has one.
 */
struct LoopElement {
    enum class Type {
        Unpaired,  ///< an unpaired nucleotide
        Helix,     ///< a branch that starts with a base pair
        Other,     ///< a branch that is not a helix
    };

    Type type = Type::Unpaired;
    Code reached = 0;  ///< where the loop reaches the element: the nucleotide, or p of a branch
    Code left = 0;     ///< where the loop leaves it: the nucleotide, or q of a branch
    Code after_reached = kNoNucleotide;  ///< Helix: the nucleotide after p (inside the helix)
    Code before_left = kNoNucleotide;    ///< Helix: the nucleotide before q (inside the helix)

    /// An unpaired nucleotide of the loop.
    [[nodiscard]] static constexpr LoopElement unpaired(Code nucleotide) noexcept {
        return {Type::Unpaired, nucleotide, nucleotide, kNoNucleotide, kNoNucleotide};
    }

    /// A helix that leaves the loop: its first pair p-q (p on the 5' side) and the nucleotides
    /// after p and before q, which its helix-end term depends on. Both have to be nucleotides;
    /// within one strand they always exist.
    [[nodiscard]] static constexpr LoopElement helix(Code p, Code q, Code after_p,
                                                     Code before_q) noexcept {
        return {Type::Helix, p, q, after_p, before_q};
    }

    /// A branch from nucleotide p to nucleotide q that is not a helix, such as a pseudoknot.
    /// efn2 has no rule for those. Here it has no helix end and nothing stacks on it or across
    /// it, but a multibranch loop still counts it as one of its helices.
    [[nodiscard]] static constexpr LoopElement other(Code p, Code q) noexcept {
        return {Type::Other, p, q, kNoNucleotide, kNoNucleotide};
    }
};

}  // namespace thermorna::rnastructure
