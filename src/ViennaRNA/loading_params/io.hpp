#pragma once

#include <vector>

namespace thermorna::viennarna {

/**
 * @brief Energy parameter file I/O.
 *
 * The class is stateless from the caller's perspective, so the API is exposed
 * through static member functions. Implementation helpers remain private to
 * io.cpp through an anonymous namespace.
 */
class ParameterIO {
public:
  static constexpr unsigned int PARAMETER_FORMAT_DEFAULT = 0;

  /**
   * @brief Load energy parameters from a file.
   *
   * @param fname   Path to the file containing energy parameters.
   * @param options File-format bit mask.
   * @return Non-zero on success, 0 on failure.
   */
  static int
  vrna_params_load(const char   fname[],
                   unsigned int options = PARAMETER_FORMAT_DEFAULT);

  /**
   * @brief Save energy parameters to a file.
   *
   * @param fname   Path to the output file.
   * @param options File-format bit mask.
   * @return Non-zero on success, 0 on failure.
   */
  static int
  vrna_params_save(const char   fname[],
                   unsigned int options = PARAMETER_FORMAT_DEFAULT);

  /**
   * @brief Load energy parameters from a string.
   *
   * @param string  Zero-terminated string containing energy parameters.
   * @param name    Optional name for the parameter set.
   * @param options File-format bit mask.
   * @return Non-zero on success, 0 on failure.
   */
  static int
  vrna_params_load_from_string(const char   *string,
                               const char   *name,
                               unsigned int options = PARAMETER_FORMAT_DEFAULT);

  /** @brief Load the default RNA parameter set (Turner 2004). */
  static int
  vrna_params_load_defaults();

  /** @brief Load the Turner 2004 RNA parameter set. */
  static int
  vrna_params_load_RNA_Turner2004();

  /** @brief Load the Turner 1999 RNA parameter set. */
  static int
  vrna_params_load_RNA_Turner1999();

  /** @brief Load the Andronescu 2007 RNA parameter set. */
  static int
  vrna_params_load_RNA_Andronescu2007();

  /** @brief Load the Langdon 2018 RNA parameter set. */
  static int
  vrna_params_load_RNA_Langdon2018();

  /** @brief Load the miscellaneous special-hairpin RNA parameter set. */
  static int
  vrna_params_load_RNA_misc_special_hairpins();

  /** @brief Load the Mathews 2004 DNA parameter set. */
  static int
  vrna_params_load_DNA_Mathews2004();

  /** @brief Load the Mathews 1999 DNA parameter set. */
  static int
  vrna_params_load_DNA_Mathews1999();
};

} // namespace thermorna::viennarna
