#include "RNAstructure/params/parameter_set.hpp"

#include "RNAstructure/error.hpp"

// POSIX directory access (Linux and macOS, the platforms Knotergy supports; std::filesystem is
// avoided for older macOS toolchains, as in Knotergy).
#include <sys/stat.h>

#include <dirent.h>

#include <algorithm>

namespace thermorna::rnastructure {

namespace {

constexpr const char* kSpecificationSuffix = ".specification.dat";
constexpr const char* kUnprefixedSpecification = "specification.dat";

bool ends_with(const std::string& text, const std::string& suffix) {
    return text.size() >= suffix.size() &&
           text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

// Last path component, or "" when there is none worth showing ("", ".", "..").
std::string folder_name(const std::string& directory) {
    std::string path = directory;
    while (path.size() > 1 && path.back() == '/') path.pop_back();
    const std::size_t slash = path.find_last_of('/');
    std::string name = slash == std::string::npos ? path : path.substr(slash + 1);
    if (name == "." || name == ".." || name == "/") name.clear();
    return name;
}

std::string join_path(const std::string& directory, const std::string& name) {
    if (directory.empty() || directory.back() == '/') return directory + name;
    return directory + "/" + name;
}

bool equal_ignoring_ascii_case(const std::string& a, const std::string& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t k = 0; k < a.size(); ++k) {
        char x = a[k];
        char y = b[k];
        if (x >= 'A' && x <= 'Z') x = static_cast<char>(x - 'A' + 'a');
        if (y >= 'A' && y <= 'Z') y = static_cast<char>(y - 'A' + 'a');
        if (x != y) return false;
    }
    return true;
}

// "turner2004 (rna)", "m6A": a set's folder name, with its file prefix when that differs.
std::string described(const ParameterSetLocation& set) {
    if (set.prefix.empty() || set.prefix == set.name) return "'" + set.name + "'";
    return "'" + set.name + "' (" + set.prefix + ")";
}

std::string described(const std::vector<const ParameterSetLocation*>& sets) {
    std::string out;
    for (const ParameterSetLocation* set : sets) {
        if (!out.empty()) out += ", ";
        out += described(*set);
    }
    return out;
}

std::string joined(const std::vector<std::string>& sets) {
    std::string out;
    for (const std::string& set : sets) {
        if (!out.empty()) out += ", ";
        out += set.empty() ? std::string("(unprefixed)") : "'" + set + "'";
    }
    return out;
}

}  // namespace

std::string ParameterSetLocation::file(const std::string& table,
                                       const std::string& extension) const {
    std::string path = directory;
    if (!path.empty() && path.back() != '/') path += '/';
    if (!prefix.empty()) path += prefix + ".";
    return path + table + extension;
}

bool is_directory(const std::string& path) {
    struct stat info {};
    return !path.empty() && ::stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
}

std::vector<std::string> list_parameter_sets(const std::string& directory) {
    DIR* dir = ::opendir(directory.c_str());
    if (dir == nullptr) {
        throw ParameterError(is_directory(directory)
                                 ? "Cannot read the parameter folder '" + directory + "'."
                                 : "The parameter folder '" + directory + "' does not exist.");
    }
    std::vector<std::string> sets;
    while (const dirent* entry = ::readdir(dir)) {
        const std::string file = entry->d_name;
        if (file == kUnprefixedSpecification) {
            sets.emplace_back();
        } else if (ends_with(file, kSpecificationSuffix) &&
                   file.size() > std::char_traits<char>::length(kSpecificationSuffix)) {
            sets.push_back(
                file.substr(0, file.size() - std::char_traits<char>::length(kSpecificationSuffix)));
        }
    }
    ::closedir(dir);
    std::sort(sets.begin(), sets.end());
    return sets;
}

ParameterSetLocation locate_parameter_set(const std::string& directory, const std::string& prefix) {
    const std::vector<std::string> sets = list_parameter_sets(directory);

    ParameterSetLocation location;
    location.directory = directory;
    location.prefix = prefix;
    if (prefix.empty()) {
        if (sets.empty()) {
            throw ParameterError("No RNAstructure parameter set in '" + directory +
                                 "' (expected specification.dat or <name>.specification.dat).");
        }
        if (sets.size() > 1) {
            throw ParameterError("The folder '" + directory + "' holds several parameter sets (" +
                                 joined(sets) + "); name the one to use.");
        }
        location.prefix = sets.front();
    } else if (!std::binary_search(sets.begin(), sets.end(), prefix)) {
        throw ParameterError(
            "No parameter set '" + prefix + "' in '" + directory + "': " + prefix +
            kSpecificationSuffix + " is missing" +
            (sets.empty() ? std::string(".") : " (the folder holds " + joined(sets) + ")."));
    }

    // A folder of its own names the set; in a shared folder (RNAstructure's data_tables) the
    // prefix does.
    const bool own_folder = sets.size() == 1 && sets.front() == location.prefix;
    location.name = own_folder ? folder_name(directory) : location.prefix;
    if (location.name.empty())
        location.name = location.prefix.empty() ? directory : location.prefix;
    return location;
}

std::vector<ParameterSetLocation> list_set_folders(const std::string& root) {
    DIR* dir = ::opendir(root.c_str());
    if (dir == nullptr) {
        throw ParameterError(is_directory(root)
                                 ? "Cannot read the parameter folder '" + root + "'."
                                 : "The parameter folder '" + root + "' does not exist.");
    }
    std::vector<std::string> names;
    while (const dirent* entry = ::readdir(dir)) {
        const std::string name = entry->d_name;
        if (!name.empty() && name[0] != '.') names.push_back(name);
    }
    ::closedir(dir);
    std::sort(names.begin(), names.end());

    std::vector<ParameterSetLocation> sets;
    for (const std::string& name : names) {
        const std::string folder = join_path(root, name);
        if (!is_directory(folder)) continue;
        std::vector<std::string> prefixes;
        try {
            prefixes = list_parameter_sets(folder);
        } catch (const ParameterError&) {
            continue;  // unreadable folder
        }
        if (prefixes.size() != 1) continue;
        ParameterSetLocation set;
        set.directory = folder;
        set.prefix = prefixes.front();
        set.name = name;
        sets.push_back(std::move(set));
    }
    return sets;
}

ParameterSetLocation find_set_folder(const std::string& root, const std::string& name) {
    if (name.empty()) throw ParameterError("No parameter set name was given.");
    const std::vector<ParameterSetLocation> sets = list_set_folders(root);

    // Folder name, then file prefix; exactly, then ignoring ASCII case.
    for (int pass = 0; pass < 4; ++pass) {
        std::vector<const ParameterSetLocation*> found;
        for (const ParameterSetLocation& set : sets) {
            const std::string& key = pass % 2 == 0 ? set.name : set.prefix;
            if (pass < 2 ? key == name : equal_ignoring_ascii_case(key, name)) {
                found.push_back(&set);
            }
        }
        if (found.size() == 1) return *found.front();
        if (found.size() > 1) {
            throw ParameterError("'" + name + "' matches several parameter sets in '" + root +
                                 "': " + described(found) + ".");
        }
    }
    std::vector<const ParameterSetLocation*> all;
    for (const ParameterSetLocation& set : sets) all.push_back(&set);
    throw ParameterError("No parameter set '" + name + "' in '" + root + "'" +
                         (all.empty() ? std::string(".") : " (it has " + described(all) + ")."));
}

}  // namespace thermorna::rnastructure
