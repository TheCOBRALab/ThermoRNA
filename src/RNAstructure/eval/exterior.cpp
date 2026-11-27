#include "RNAstructure/eval/exterior.hpp"

#include "RNAstructure/detail/loop_stacking.hpp"

namespace thermorna::rnastructure {

int exterior_stem_energy(Code p, Code q, Code after_p, Code before_q, Code n5, Code n3,
                         Stacking stacking, const Datatable& d) noexcept {
    return detail::forbidden_as_infinite(
        detail::helix_end(p, q, after_p, before_q, n3, n5, d) +
        detail::stem_stacking(p, q, n5, n3, stacking, Table4D::Tstack, d));
}

int exterior_energy(const std::vector<LoopElement>& loop, const Datatable& d,
                    CoaxialStacking coaxial) {
    return detail::forbidden_as_infinite(
        detail::helix_ends(loop, d) +
        detail::lowest_stacking_energy(loop, Table4D::Tstack, d, coaxial));
}

}  // namespace thermorna::rnastructure
