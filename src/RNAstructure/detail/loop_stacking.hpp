#pragma once

/**
 * @file RNAstructure/detail/loop_stacking.hpp
 * @brief What efn2's loop energies have in common: the helix end, what stacks on it, and for
 * exterior and multibranch loops the lowest energy over the ways helices and nucleotides can
 * stack.
 */

#include "RNAstructure/constants.hpp"
#include "RNAstructure/eval/loop_element.hpp"
#include "RNAstructure/eval/stacking.hpp"
#include "RNAstructure/params/datatable.hpp"

#include <cstddef>
#include <vector>

namespace thermorna::rnastructure::detail {

/// The last index of RNAstructure's dangle table: which side of the helix the nucleotide is on.
inline constexpr std::size_t kDangle3 = 1;
inline constexpr std::size_t kDangle5 = 2;

/// efn2 adds a forbidden table entry (kInfiniteEnergy) up like any other term, so a forbidden
/// loop comes out near a multiple of it (efn2 prints 1399.9, 1400.5, 2800.5, ...). No allowed
/// loop is anywhere near half of it. Such a sum becomes kInfiniteEnergy itself.
[[nodiscard]] constexpr int forbidden_as_infinite(int energy) noexcept {
    return energy >= kInfiniteEnergy / 2 ? static_cast<int>(kInfiniteEnergy) : energy;
}

/// efn2's penalty(a, b): the end of a helix whose last pair a-b a loop meets at a first. It
/// depends on the nucleotides after a and before b, the helix's next pair if it has one. When
/// one of them does not exist, efn2 looks at the pair from the other side: b-a with the
/// nucleotides after b and before a. A pair with neither has no helix end.
[[nodiscard]] inline int helix_end(Code a, Code b, Code after_a, Code before_b, Code after_b,
                                   Code before_a, const Datatable& d) noexcept {
    if (after_a != kNoNucleotide && before_b != kNoNucleotide) {
        return d.helix_ends(a, b, after_a, before_b);
    }
    if (after_b != kNoNucleotide && before_a != kNoNucleotide) {
        return d.helix_ends(b, a, after_b, before_a);
    }
    return kInfiniteEnergy;
}

/// What stacks on the end of helix p-q, which a loop reaches at p and leaves at q: n5, the
/// nucleotide before p, and n3, the nucleotide after q, as @p stacking says. @p mismatch is
/// the terminal mismatch table of the loop (Table4D::Tstack or Table4D::Tstkm). The tables
/// have the helix with q first. kInfiniteEnergy when a nucleotide that is to stack is not
/// there.
[[nodiscard]] inline int stem_stacking(Code p, Code q, Code n5, Code n3, Stacking stacking,
                                       Table4D mismatch, const Datatable& d) noexcept {
    const bool has5 = n5 != kNoNucleotide;
    const bool has3 = n3 != kNoNucleotide;
    switch (stacking) {
        case Stacking::None:    return 0;
        case Stacking::Dangle5: return has5 ? d.dangle(q, p, n5, kDangle5) : kInfiniteEnergy;
        case Stacking::Dangle3: return has3 ? d.dangle(q, p, n3, kDangle3) : kInfiniteEnergy;
        case Stacking::Mismatch:
            return has5 && has3 ? d.table4(mismatch, q, p, n3, n5) : kInfiniteEnergy;
    }
    return 0;
}

/// The helix ends of the branches of a loop.
[[nodiscard]] int helix_ends(const std::vector<LoopElement>& elements, const Datatable& d) noexcept;

/// The minimisation of efn2's ergexterior and ergmulti: the lowest energy over the ways the
/// helices of @p elements and the nucleotides next to them can stack. @p mismatch is the
/// terminal mismatch table of the loop (Table4D::Tstack or Table4D::Tstkm). With
/// CoaxialStacking::Off the helices do not stack on each other.
[[nodiscard]] int lowest_stacking_energy(const std::vector<LoopElement>& elements, Table4D mismatch,
                                         const Datatable& d, CoaxialStacking coaxial);

}  // namespace thermorna::rnastructure::detail
