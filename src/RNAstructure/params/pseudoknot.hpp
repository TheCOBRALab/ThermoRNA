#pragma once

/**
 * @file RNAstructure/params/pseudoknot.hpp
 * @brief The pseudoknot penalty constants efn2 reads for pseudoknotted structures.
 *
 * efn2 (CalculateFreeEnergywithPseudoknot) reads 16 constants from <DATAPATH>/pseudconst.dat,
 * each one written after an "-->" marker, and uses P1 = 0.35 and P2 = 0.65 kcal/mol unless told
 * otherwise. See ReadPseudoParam() in RNAstructure's PseudoParser.cpp.
 */

#include "RNAstructure/constants.hpp"

#include <array>
#include <string>
#include <string_view>

namespace thermorna::rnastructure {

struct PseudoknotConstants {
    /// Constants in file order: a, b, c2 ... c14, c15 (ShapeKnots penalty model).
    std::array<double, 16> values{};
    double p1 = kDefaultPseudoknotP1;
    double p2 = kDefaultPseudoknotP2;

    /// Parses the text of a pseudconst.dat file. @p source is used in error messages.
    [[nodiscard]] static PseudoknotConstants parse(std::string_view text,
                                                   const std::string& source);

    /// Reads <directory>/pseudconst.dat; an empty directory means the DATAPATH variable.
    [[nodiscard]] static PseudoknotConstants read(const std::string& directory = "");
};

}  // namespace thermorna::rnastructure
