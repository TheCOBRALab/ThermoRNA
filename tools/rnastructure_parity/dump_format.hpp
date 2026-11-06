#pragma once

// Shared text format for comparing parameter tables loaded by RNAstructure and by ThermoRNA.
// Both dumpers print every scalar, every special hairpin (as RNAstructure's seqtonum key) and
// every table entry that is not INFINITE_ENERGY, followed by per-table entry counts.

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>

namespace parity {

inline std::string float_bits(float value) {
    std::uint32_t bits = 0;
    std::memcpy(&bits, &value, sizeof bits);
    char buffer[32];
    std::snprintf(buffer, sizeof buffer, "%08x(%.9g)", bits, static_cast<double>(value));
    return buffer;
}

constexpr const char* kTable4Names[] = {"stack",     "helix_ends", "tstkh", "tstki",
                                        "tstki23",   "tstki1n",    "coax",  "tstackcoax",
                                        "coaxstack", "tstack",     "tstkm"};

}  // namespace parity
