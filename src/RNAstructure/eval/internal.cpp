#include "RNAstructure/eval/internal.hpp"

#include "RNAstructure/detail/loop_stacking.hpp"

#include <algorithm>
#include <cmath>

namespace thermorna::rnastructure {

int internal_loop_energy(std::size_t n1, std::size_t n2, Code i, Code j, Code p, Code q, Code si1,
                         Code sj1, Code sp1, Code sq1, Code beyond_i, Code beyond_j, Code beyond_p,
                         Code beyond_q, std::size_t bulge_run, const Datatable& d) noexcept {
    const std::size_t size = n1 + n2;
    if (size == 0) return d.stack(i, j, p, q);

    const MiscLoop& m = d.misc();
    if (!d.can_pair(i, j) || !d.can_pair(p, q)) return kInfiniteEnergy;

    // Loop initiation by size; beyond 30 nucleotides it is extrapolated logarithmically.
    const auto initiation = [&](const std::vector<Energy>& by_size) {
        if (size <= 30) return static_cast<int>(by_size[size]);
        return by_size[30] +
               static_cast<int>(std::round(m.prelog * std::log(static_cast<double>(size) / 30.0)));
    };
    // The helix ends that face the loop: the outer pair from inside, the inner from outside.
    const auto helix_ends = [&] {
        return detail::helix_end(j, i, beyond_j, beyond_i, n1 > 0 ? si1 : p, n2 > 0 ? sj1 : q, d) +
               d.helix_ends(p, q, beyond_p, beyond_q);
    };

    if (n1 == 0 || n2 == 0) {
        // ---------------- Bulge ----------------
        if (size > 1) {
            return detail::forbidden_as_infinite(initiation(d.bulge()) + m.eparam[2] +
                                                 helix_ends());
        }

        // A single bulged nucleotide does not interrupt the helix: the two pairs stack. The
        // bulge can sit at any of the identical nucleotides next to it.
        int energy = d.stack(i, j, p, q) + d.bulge()[1] + m.eparam[2];
        const Code bulged = n1 == 1 ? si1 : sj1;
        if (d.has_symbol(bulged, "C") && bulge_run > 1) energy += m.singlecbulge;
        energy -=
            static_cast<int>(std::round(d.RT() * 10 * std::log(static_cast<double>(bulge_run))));
        return detail::forbidden_as_infinite(energy);
    }

    // ---------------- Internal loop ----------------
    int energy = 0;
    if (n1 == 2 && n2 == 2) {
        energy = d.iloop22(i, p, j, q, si1, sp1, sj1, sq1);
    } else if (n1 == 1 && n2 == 2) {
        energy = d.iloop21(i, j, si1, sj1, sq1, p, q);
    } else if (n1 == 2 && n2 == 1) {
        energy = d.iloop21(q, p, sq1, sp1, si1, j, i);
    } else if (size == 2) {
        energy = d.iloop11(i, si1, p, j, sj1, q);
    } else {
        // Mismatches next to both pairs; 1xn and 2x3 loops have tables of their own.
        Table4D mismatch = Table4D::Tstki;
        if (n1 == 1 || n2 == 1) {
            mismatch = Table4D::Tstki1n;
        } else if ((n1 == 2 && n2 == 3) || (n1 == 3 && n2 == 2)) {
            mismatch = Table4D::Tstki23;
        }
        energy = d.table4(mismatch, i, j, si1, sj1) + d.table4(mismatch, q, p, sq1, sp1) +
                 initiation(d.inter()) + m.eparam[3];

        // Asymmetry (Ninio) penalty
        const int lopsided = static_cast<int>(n1 > n2 ? n1 - n2 : n2 - n1);
        const std::size_t shorter = std::min<std::size_t>(2, std::min(n1, n2));
        energy += std::min<int>(m.maxpen, lopsided * m.poppen[shorter]);
    }
    return detail::forbidden_as_infinite(energy + helix_ends());
}

}  // namespace thermorna::rnastructure
