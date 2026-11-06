#include "RNAstructure/params/pseudoknot.hpp"

#include "RNAstructure/detail/text.hpp"
#include "RNAstructure/error.hpp"

#include <cstdlib>
#include <cstring>

namespace thermorna::rnastructure {

PseudoknotConstants PseudoknotConstants::parse(std::string_view text, const std::string& source) {
    PseudoknotConstants constants;

    // ReadPseudoParam() reads whitespace-separated words, skips to each "-->" and reads the
    // number after it.
    detail::LineCursor words(text);
    std::string_view word;
    for (std::size_t i = 0; i < constants.values.size(); ++i) {
        do {
            if (!words.next_token(word)) {
                throw ParameterError(source +
                                     ": expected 16 constants after '-->' markers, found " +
                                     std::to_string(i));
            }
        } while (word != "-->");
        if (!words.next_token(word)) {
            throw ParameterError(source + ": missing value after the last '-->' marker");
        }
        const std::string value(word);
        char* end = nullptr;
        constants.values[i] = std::strtod(value.c_str(), &end);
        if (end == value.c_str()) {
            throw ParameterError(source + ": '" + value + "' after '-->' is not a number");
        }
    }
    return constants;
}

PseudoknotConstants PseudoknotConstants::read(const std::string& directory) {
    std::string where = directory;
    if (where.empty()) {
        const char* environment = std::getenv("DATAPATH");
        if (environment == nullptr || *environment == '\0') {
            throw ParameterError("No data_tables directory was given and DATAPATH is not set.");
        }
        where = environment;
    }
    const std::string path = where + "/pseudconst.dat";
    std::string text;
    if (!detail::read_file(path, text)) {
        throw ParameterError("Missing pseudoknot constants file: " + path);
    }
    return parse(text, path);
}

}  // namespace thermorna::rnastructure
