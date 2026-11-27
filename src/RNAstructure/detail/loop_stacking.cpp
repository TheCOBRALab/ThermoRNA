#include "RNAstructure/detail/loop_stacking.hpp"

#include <algorithm>

namespace thermorna::rnastructure::detail {

int helix_ends(const std::vector<LoopElement>& elements, const Datatable& d) noexcept {
    int energy = 0;
    for (const LoopElement& e : elements) {
        if (e.type != LoopElement::Type::Helix) continue;
        energy += d.helix_ends(e.reached, e.left, e.after_reached, e.before_left);
    }
    return energy;
}

// best[n] is the lowest energy of the first n elements; each option uses up the elements it
// stacks. The terms are the ones the stem-level functions return (RNAstructure/eval/stacking.hpp,
// exterior_stem_energy, multibranch_stem_energy), without the helix ends.
int lowest_stacking_energy(const std::vector<LoopElement>& elements, Table4D mismatch,
                           const Datatable& d, CoaxialStacking coaxial) {
    const bool helices_stack = coaxial == CoaxialStacking::On;
    std::vector<int> best(elements.size() + 1, 0);
    for (std::size_t n = 1; n <= elements.size(); ++n) {
        const std::size_t k = n - 1;
        const LoopElement& e = elements[k];
        const auto helix = [&](std::size_t back) {
            return n > back && elements[k - back].type == LoopElement::Type::Helix;
        };
        const auto unpaired = [&](std::size_t back) {
            return n > back && elements[k - back].type == LoopElement::Type::Unpaired;
        };
        int lowest = best[n - 1];  // nothing stacks here
        const auto option = [&](std::size_t used, int energy) {
            lowest = std::min(lowest, best[n - used] + energy);
        };

        if (e.type == LoopElement::Type::Unpaired) {
            if (helix(1)) {
                const LoopElement& h = elements[k - 1];
                // 3' dangle on the helix before this nucleotide
                option(2, stem_stacking(h.reached, h.left, kNoNucleotide, e.left, Stacking::Dangle3,
                                        mismatch, d));
                // Terminal mismatch with the nucleotide before that helix
                if (unpaired(2)) {
                    option(3, stem_stacking(h.reached, h.left, elements[k - 2].left, e.left,
                                            Stacking::Mismatch, mismatch, d));
                }
                // The two helices before this nucleotide stack across the nucleotide between
                // them, which forms a mismatch with this one.
                if (helices_stack && unpaired(2) && helix(3)) {
                    const LoopElement& between = elements[k - 2];
                    const LoopElement& h2 = elements[k - 3];
                    option(4, coaxial_stack_mismatch3_energy(h2.reached, h2.left, h.reached, h.left,
                                                             between.left, e.left, d));
                }
            }
        } else if (e.type == LoopElement::Type::Helix) {
            // 5' dangle of the nucleotide before this helix
            if (unpaired(1)) {
                option(2, stem_stacking(e.reached, e.left, elements[k - 1].left, kNoNucleotide,
                                        Stacking::Dangle5, mismatch, d));
            }
            // Flush coaxial stack on the helix right before this one
            if (helices_stack && helix(1)) {
                const LoopElement& h = elements[k - 1];
                option(2, coaxial_stack_energy(h.reached, h.left, e.reached, e.left, d));
            }
            // This helix and the one before it stack across the nucleotide between them, which
            // forms a mismatch with the nucleotide before that helix.
            if (helices_stack && unpaired(1) && helix(2) && unpaired(3)) {
                const LoopElement& between = elements[k - 1];
                const LoopElement& h = elements[k - 2];
                const LoopElement& before = elements[k - 3];
                option(4, coaxial_stack_mismatch5_energy(h.reached, h.left, e.reached, e.left,
                                                         before.left, between.left, d));
            }
        }
        best[n] = lowest;
    }
    return best.back();
}

}  // namespace thermorna::rnastructure::detail
