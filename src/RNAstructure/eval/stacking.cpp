#include "RNAstructure/eval/stacking.hpp"

#include "RNAstructure/detail/loop_stacking.hpp"

// The tables have a helix as the loop sees it, but with the nucleotide where the loop leaves
// it first (its 3' neighbor is in the loop), then the one where the loop reaches it.

namespace thermorna::rnastructure {

int coaxial_stack_energy(Code p1, Code q1, Code p2, Code q2, const Datatable& d) noexcept {
    return d.coax(q1, p1, p2, q2);
}

int coaxial_stack_mismatch5_energy(Code p1, Code q1, Code p2, Code q2, Code before, Code between,
                                   const Datatable& d) noexcept {
    return detail::forbidden_as_infinite(d.coaxstack(q2, p2, before, between) +
                                         d.tstackcoax(q1, p1, between, before));
}

int coaxial_stack_mismatch3_energy(Code p1, Code q1, Code p2, Code q2, Code between, Code after,
                                   const Datatable& d) noexcept {
    return detail::forbidden_as_infinite(d.coaxstack(q1, p1, between, after) +
                                         d.tstackcoax(q2, p2, after, between));
}

}  // namespace thermorna::rnastructure
