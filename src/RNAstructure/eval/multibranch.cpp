#include "RNAstructure/eval/multibranch.hpp"

#include "RNAstructure/detail/loop_stacking.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace thermorna::rnastructure {

int multibranch_stem_energy(Code p, Code q, Code after_p, Code before_q, Code n5, Code n3,
                            Stacking stacking, const Datatable& d) noexcept {
    return detail::forbidden_as_infinite(
        d.misc().efn2c + detail::helix_end(p, q, after_p, before_q, n3, n5, d) +
        detail::stem_stacking(p, q, n5, n3, stacking, Table4D::Tstkm, d));
}

int multibranch_energy(Code i, Code j, Code beyond_i, Code beyond_j,
                       const std::vector<LoopElement>& loop, const Datatable& d,
                       MultibranchRule rule, CoaxialStacking coaxial) {
    const MiscLoop& m = d.misc();

    // The helix ends: the closing pair seen from inside, and the branches
    int energy =
        detail::helix_end(j, i, beyond_j, beyond_i, loop.front().reached, loop.back().left, d) +
        detail::helix_ends(loop, d);

    // Unpaired nucleotides before each helix, going round the loop; the closing pair is last.
    std::vector<int> gaps;
    int gap = 0;
    for (const LoopElement& e : loop) {
        if (e.type == LoopElement::Type::Unpaired) {
            ++gap;
        } else {
            gaps.push_back(gap);
            gap = 0;
        }
    }
    gaps.push_back(gap);
    const int helices = static_cast<int>(gaps.size());
    const int unpaired = static_cast<int>(loop.size()) + 1 - helices;

    // Offset, penalty per helix and penalty per unpaired nucleotide
    energy += m.efn2a + helices * m.efn2c;
    if (rule == MultibranchRule::Simple) {
        energy += unpaired * m.efn2b;
    } else {
        // Beyond 8 nucleotides efn2 switches to a logarithm with a fixed coefficient of 1.1
        // kcal/mol.
        if (unpaired > 8) {
            energy += 8 * m.efn2b + static_cast<int>(11.0 * std::log(unpaired / 8.0) + 0.5);
        } else {
            energy += unpaired * m.efn2b;
        }

        // Strain: an odd number of helices with at most one unpaired nucleotide among them
        if (helices % 2 != 0 && unpaired <= 1) energy += m.strain;

        // Asymmetry: how much the gaps on the two sides of a helix differ, averaged over the
        // helices and capped at 2. In single precision like efn2, so that it rounds the same
        // way.
        int differences = 0;
        for (std::size_t k = 0; k < gaps.size(); ++k) {
            differences += std::abs(gaps[k] - gaps[(k + 1) % gaps.size()]);
        }
        const float asymmetry =
            std::min(2.0F, static_cast<float>(differences) / static_cast<float>(helices));
        energy += static_cast<int>(static_cast<float>(m.mlasym) * asymmetry + 0.5);
    }

    // Stacking, as in an exterior loop. The loop is a circle: the closing pair, which the loop
    // reaches at j and leaves at i, then what is inside it. It has to be opened somewhere to
    // go through it, and what stacks together covers up to 4 elements, so efn2 opens it in 4
    // places next to each other and keeps the best: one of them does not cut through a stack.
    std::vector<LoopElement> circle;
    circle.reserve(loop.size() + 1);
    circle.push_back(LoopElement::helix(j, i, kNoNucleotide, kNoNucleotide));
    circle.insert(circle.end(), loop.begin(), loop.end());
    int stacking = 0;
    for (int opening = 0; opening < 4; ++opening) {
        stacking =
            std::min(stacking, detail::lowest_stacking_energy(circle, Table4D::Tstkm, d, coaxial));
        std::rotate(circle.begin(), circle.begin() + 1, circle.end());
    }

    return detail::forbidden_as_infinite(energy + stacking);
}

}  // namespace thermorna::rnastructure
