#pragma once

/**
 * @file RNAstructure/rnastructure.hpp
 * @brief Everything of ThermoRNA's RNAstructure component: efn2's parameter tables and its
 * energy functions.
 *
 * @code
 * namespace rs = thermorna::rnastructure;
 *
 * rs::LoadOptions options;
 * options.data_directory = "params/RNAstructure/turner2004";  // a set's folder
 * options.temperature = 37.0 + rs::kCelsiusToKelvin;         // efn2 --temperature (Kelvin)
 * options.modifications = {{"params/RNAstructure/m6A"},      // any number of modified bases
 *                          {"params/RNAstructure/Ψ"}};
 *
 * const rs::Datatable data = rs::Datatable::load(options);
 * const std::vector<rs::Code> s = data.encode("GG6AAPUCC");  // '6' = m6A, 'P' = pseudouridine
 * rs::Energy e = data.stack(s[0], s[8], s[1], s[7]);          // data->stack[..][..][..][..]
 * int hairpin = rs::hairpin_energy(3, &s[2], s[1], s[0], s[7], data);  // efn2's erg3: the loop
 *                                                // 6AAPU, with G and G before it and C after
 * @endcode
 *
 * RNAstructure's data_tables folder, which holds every set, works as well: there the set is
 * chosen by file prefix (options.alphabet = "rna", modifications {{"", "m6A"}, {"", "PSU"}}).
 */

#include "RNAstructure/constants.hpp"
#include "RNAstructure/error.hpp"
#include "RNAstructure/eval/exterior.hpp"
#include "RNAstructure/eval/hairpin.hpp"
#include "RNAstructure/eval/internal.hpp"
#include "RNAstructure/eval/loop_element.hpp"
#include "RNAstructure/eval/multibranch.hpp"
#include "RNAstructure/eval/stacking.hpp"
#include "RNAstructure/params/alphabet.hpp"
#include "RNAstructure/params/alphabet_tables.hpp"
#include "RNAstructure/params/datatable.hpp"
#include "RNAstructure/params/parameter_set.hpp"
#include "RNAstructure/params/pseudoknot.hpp"
#include "RNAstructure/utils/graphemes.hpp"
