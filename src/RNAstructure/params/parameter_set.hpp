#pragma once

/**
 * @file RNAstructure/params/parameter_set.hpp
 * @brief Finding a parameter set's files on disk.
 *
 * RNAstructure keeps every parameter set in one folder (data_tables) and tells them apart by
 * file prefix: rna.stack.dg, m6A.stack.dg, 1.stack.dg, ... ThermoRNA also accepts a folder per
 * set, e.g.
 *
 *     params/RNAstructure/turner2004/rna.stack.dg, rna.specification.dat, ...
 *     params/RNAstructure/m1Ψ/1.stack.dg, 1.specification.dat, ...
 *     params/RNAstructure/mine/stack.dg, specification.dat, ...   (no prefix)
 *
 * In a folder that holds a single set the prefix is found from its specification file, and the
 * folder's name becomes the set's name in messages ("m1Ψ" rather than "1").
 */

#include <string>
#include <vector>

namespace thermorna::rnastructure {

/// Where one parameter set's files are and how they are named.
struct ParameterSetLocation {
    std::string directory{};  ///< folder holding the files
    std::string prefix{};     ///< file prefix ("rna", "m6A"); empty for unprefixed file names
    std::string name{};       ///< the folder's name if it holds only this set, else the prefix

    /// Path of one of the set's files, e.g. file("stack", ".dg").
    [[nodiscard]] std::string file(const std::string& table, const std::string& extension) const;
};

/// True when @p path is an existing directory.
[[nodiscard]] bool is_directory(const std::string& path);

/// File prefixes of the parameter sets in @p directory, sorted; "" stands for a set with
/// unprefixed file names (specification.dat). Throws ParameterError if the folder cannot be read.
[[nodiscard]] std::vector<std::string> list_parameter_sets(const std::string& directory);

/// Locates set @p prefix in @p directory. An empty prefix means the only set in the folder; a
/// folder holding several sets then needs the prefix named. Throws ParameterError when the
/// folder or the set does not exist.
[[nodiscard]] ParameterSetLocation locate_parameter_set(const std::string& directory,
                                                        const std::string& prefix = "");

/// The sets of @p root, a folder with one sub-folder per set (like Knotergy's
/// params/RNAstructure): every sub-folder that holds exactly one set, sorted by folder name.
/// Hidden folders and folders without a set (e.g. one with pseudconst.dat only) are skipped.
[[nodiscard]] std::vector<ParameterSetLocation> list_set_folders(const std::string& root);

/// Finds set @p name in @p root (see list_set_folders). @p name is a folder name
/// ("turner2004", "m1Ψ") or a file prefix, i.e. RNAstructure's alphabet name ("rna" finds
/// turner2004/ when it holds rna.*). Exact matches win over matches that ignore ASCII case
/// ("m6a"). Throws ParameterError, listing the sets of @p root, when none or several match.
[[nodiscard]] ParameterSetLocation find_set_folder(const std::string& root,
                                                   const std::string& name);

}  // namespace thermorna::rnastructure
