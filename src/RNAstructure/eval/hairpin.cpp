#include "RNAstructure/eval/hairpin.hpp"

#include "RNAstructure/detail/loop_stacking.hpp"

#include <cmath>
#include <optional>

namespace thermorna::rnastructure {

int hairpin_energy(std::size_t size, const Code* loop, Code beyond_i, Code beyond_i2, Code beyond_j,
                   const Datatable& d) noexcept {
    if (size < 3) return kInfiniteEnergy;

    const MiscLoop& m = d.misc();
    const Code i = loop[0];
    const Code j = loop[size + 1];
    if (!d.can_pair(i, j)) return kInfiniteEnergy;

    // The end of the closing helix, seen from inside the loop
    const int end = detail::helix_end(j, i, beyond_j, beyond_i, loop[1], loop[size], d);

    // Special hairpins are looked up by the whole loop, closing pair included, and replace
    // every other term.
    if (size == 3 || size == 4 || size == 6) {
        const std::optional<int> special = size == 3   ? d.triloop(loop)
                                           : size == 4 ? d.tloop(loop)
                                                       : d.hexaloop(loop);
        if (special) return detail::forbidden_as_infinite(*special + end);
    }

    // Loop initiation; beyond 30 nucleotides it is extrapolated logarithmically.
    int energy = m.eparam[4];  // always 0 in RNAstructure
    if (size <= 30) {
        energy += d.hairpin()[size];
    } else {
        energy +=
            d.hairpin()[30] +
            static_cast<int>(std::round(m.prelog * std::log(static_cast<double>(size) / 30.0)));
    }

    // First mismatch; loops of 3 have none.
    if (size > 3) energy += d.tstkh(i, j, loop[1], loop[size]);

    // GU closure preceded by two Gs (5'-GGG...U-3')
    if (d.has_symbol(i, "G") && d.has_symbol(j, "U") && d.has_symbol(beyond_i, "G") &&
        d.has_symbol(beyond_i2, "G")) {
        energy += m.gubonus;
    }

    energy += end;

    // All-C loops
    const Code C = d.code("C");
    for (std::size_t k = 1; k <= size; ++k) {
        if (loop[k] != C) return detail::forbidden_as_infinite(energy);
    }
    energy += size == 3 ? m.c3 : m.cint + static_cast<int>(size) * m.cslope;
    return detail::forbidden_as_infinite(energy);
}

}  // namespace thermorna::rnastructure
