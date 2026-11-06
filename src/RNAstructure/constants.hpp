#pragma once

/**
 * @file RNAstructure/constants.hpp
 * @brief Units and sentinels shared by the RNAstructure parameter loader.
 *
 * Every value mirrors a definition in RNAstructure's defines.h / rna_library.h so that tables
 * loaded here hold exactly the numbers efn2 works with. Names are prefixed with k (and never
 * spelled like the ViennaRNA macros INF, UNIT, ...), so this header can be included next to
 * the ViennaRNA headers of ThermoRNA.
 */

#include <cstddef>
#include <cstdint>

namespace thermorna::rnastructure {

/// Free energies and enthalpies in tenths of kcal/mol, as RNAstructure stores them (`short`).
using Energy = std::int16_t;

/// A nucleotide code: the index of a base in an alphabet (RNAstructure's `numseq` values).
using Code = std::uint8_t;

/// RNAstructure's INFINITE_ENERGY: the value of a "." entry and of every entry no file sets.
inline constexpr Energy kInfiniteEnergy = 14000;

/// RNAstructure's conversionfactor: file values in kcal/mol are multiplied by this and rounded.
inline constexpr int kConversionFactor = 10;

/// 37 degrees Celsius in Kelvin (RNAstructure's TEMP_37C and T37inK).
inline constexpr double kTemperature37C = 310.15;

/// Offset between Celsius and Kelvin.
inline constexpr double kCelsiusToKelvin = 273.15;

/// Temperatures closer than this to 310.15 K are not rescaled (RNA_class TOLERANCE).
inline constexpr double kTemperatureTolerance = 0.01;

/// RNAstructure's RT_37C (kcal/mol), stored as a float exactly like datatable::RT.
inline constexpr float kRT37C = 0.61633f;

/// Default pseudoknot parameters used by efn2 (DEFAULT_PSEUDOKNOT_P1 / _P2, kcal/mol).
inline constexpr double kDefaultPseudoknotP1 = 0.35;
inline constexpr double kDefaultPseudoknotP2 = 0.65;

/// Returned by symbol lookups for symbols that are not part of the alphabet.
inline constexpr Code kInvalidCode = 0xFF;

/// Largest number of distinct nucleotides one alphabet (or a merged Datatable) can hold.
inline constexpr std::size_t kMaxCodes = 0xFE;

/// True for the "infinite" sentinel.
[[nodiscard]] constexpr bool is_infinite(Energy e) noexcept { return e == kInfiniteEnergy; }

/// Converts a table value to kcal/mol.
[[nodiscard]] constexpr double to_kcal(Energy e) noexcept {
    return static_cast<double>(e) / kConversionFactor;
}

/// Converts a table value to ViennaRNA/Knotergy units (dcal/mol, i.e. 0.01 kcal/mol).
[[nodiscard]] constexpr int to_dcal(Energy e) noexcept { return static_cast<int>(e) * 10; }

}  // namespace thermorna::rnastructure
