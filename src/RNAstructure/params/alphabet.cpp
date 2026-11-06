#include "RNAstructure/params/alphabet.hpp"

#include "RNAstructure/detail/text.hpp"
#include "RNAstructure/error.hpp"
#include "RNAstructure/params/parameter_set.hpp"
#include "RNAstructure/utils/graphemes.hpp"

#include <algorithm>

namespace thermorna::rnastructure {

namespace {

enum class Stage { None, Bases, Pairing, Single, NonInteracting, Linker };

struct PendingSymbol {
    std::string symbol;
    std::size_t line = 0;
};

struct PendingPair {
    std::string first;
    std::string second;
    std::size_t line = 0;
};

// RNAstructure removes spaces, '=' and '\r' from every specification line before reading it
// character by character. All C-locale whitespace is removed here.
std::string clean_line(std::string_view line) {
    std::string cleaned;
    cleaned.reserve(line.size());
    for (char c : line) {
        if (c == '=' || detail::is_space(c)) continue;
        cleaned.push_back(c);
    }
    return cleaned;
}

}  // namespace

Alphabet Alphabet::parse(std::string_view text, std::string name, const std::string& source) {
    Alphabet alphabet;
    alphabet.name_ = std::move(name);

    if (text.size() >= 3 && text.compare(0, 3, "\xEF\xBB\xBF") == 0) text.remove_prefix(3);

    Stage stage = Stage::None;
    std::vector<PendingPair> pairs;
    std::vector<PendingSymbol> non_interacting;
    std::vector<PendingSymbol> linker;

    std::size_t number = 0;
    std::size_t start = 0;
    while (start < text.size()) {
        std::size_t end = text.find('\n', start);
        if (end == std::string_view::npos) end = text.size();
        ++number;
        const std::string line = clean_line(text.substr(start, end - start));
        start = end + 1;

        if (line.empty() || line.front() == '#') continue;

        if (line == "Bases") {
            stage = Stage::Bases;
        } else if (line == "Pairing") {
            stage = Stage::Pairing;
        } else if (line == "Single") {
            stage = Stage::Single;
        } else if (line == "Non-interacting") {
            stage = Stage::NonInteracting;
        } else if (line == "Linker") {
            stage = Stage::Linker;
        } else {
            const std::vector<std::string_view> symbols = graphemes::split(line);
            switch (stage) {
                case Stage::Bases: {
                    std::vector<std::string> base;
                    base.reserve(symbols.size());
                    for (std::string_view s : symbols) base.emplace_back(s);
                    alphabet.bases_.push_back(std::move(base));
                    break;
                }
                case Stage::Pairing:
                    if (symbols.size() < 2) {
                        throw ParameterError(detail::where(source, number) +
                                             ": a Pairing line needs two bases, got '" + line +
                                             "'");
                    }
                    pairs.push_back({std::string(symbols[0]), std::string(symbols[1]), number});
                    break;
                case Stage::Single: alphabet.single_.emplace_back(symbols[0]); break;
                case Stage::NonInteracting:
                    non_interacting.push_back({std::string(symbols[0]), number});
                    break;
                case Stage::Linker: linker.push_back({std::string(symbols[0]), number}); break;
                case Stage::None:   break;  // RNAstructure ignores lines before a section
            }
        }
    }

    const std::size_t n = alphabet.bases_.size();
    if (n == 0) throw ParameterError(source + ": the specification lists no Bases");
    if (n > kMaxCodes) {
        throw ParameterError(source + ": the specification lists " + std::to_string(n) +
                             " bases; at most " + std::to_string(kMaxCodes) + " are supported");
    }

    // Symbol lookup. Like basetonum(), the first base that lists a symbol owns it.
    alphabet.ascii_lookup_.fill(kInvalidCode);
    for (std::size_t i = 0; i < n; ++i) {
        const Code code = static_cast<Code>(i);
        for (const std::string& symbol : alphabet.bases_[i]) {
            if (symbol.size() == 1) {
                Code& slot = alphabet.ascii_lookup_[static_cast<unsigned char>(symbol[0])];
                if (slot == kInvalidCode) slot = code;
            } else if (alphabet.code(symbol) == kInvalidCode) {
                alphabet.multibyte_.emplace_back(symbol, code);
                std::sort(alphabet.multibyte_.begin(), alphabet.multibyte_.end());
            }
        }
    }

    auto resolve = [&](const std::string& symbol, std::size_t line, const char* section) {
        const Code code = alphabet.code(symbol);
        if (code == kInvalidCode) {
            throw ParameterError(detail::where(source, line) + ": '" + symbol + "' in the " +
                                 section + " section is not one of the Bases");
        }
        return code;
    };

    alphabet.pairing_.assign(n * n, 0);
    for (const PendingPair& pair : pairs) {
        const std::size_t a = resolve(pair.first, pair.line, "Pairing");
        const std::size_t b = resolve(pair.second, pair.line, "Pairing");
        alphabet.pairing_[a * n + b] = 1;
        alphabet.pairing_[b * n + a] = 1;
    }

    for (const PendingSymbol& entry : non_interacting) {
        const Code code = resolve(entry.symbol, entry.line, "Non-interacting");
        alphabet.non_interacting_.push_back(entry.symbol);
        alphabet.non_interacting_list_.push_back(code);
        alphabet.non_interacting_codes_[code] = 1;
    }
    for (const PendingSymbol& entry : linker) {
        const Code code = resolve(entry.symbol, entry.line, "Linker");
        alphabet.linker_.push_back(entry.symbol);
        alphabet.linker_list_.push_back(code);
        alphabet.linker_codes_[code] = 1;
    }

    return alphabet;
}

Alphabet Alphabet::read(const std::string& directory, const std::string& prefix) {
    const ParameterSetLocation location = locate_parameter_set(directory, prefix);
    const std::string path = location.file("specification", ".dat");
    std::string text;
    if (!detail::read_file(path, text)) {
        throw ParameterError("Missing alphabet specification file: " + path);
    }
    return parse(text, location.name, path);
}

Code Alphabet::code(std::string_view symbol) const noexcept {
    if (symbol.size() == 1) return ascii_lookup_[static_cast<unsigned char>(symbol[0])];
    const auto it =
        std::lower_bound(multibyte_.begin(), multibyte_.end(), symbol,
                         [](const std::pair<std::string, Code>& entry, std::string_view key) {
                             return std::string_view(entry.first) < key;
                         });
    if (it != multibyte_.end() && it->first == symbol) return it->second;
    return kInvalidCode;
}

bool Alphabet::has_symbol(Code code, std::string_view symbol) const noexcept {
    if (code >= bases_.size()) return false;
    const std::vector<std::string>& symbols = bases_[code];
    return std::find(symbols.begin(), symbols.end(), symbol) != symbols.end();
}

bool Alphabet::is_single_symbol(std::string_view symbol) const noexcept {
    return std::find(single_.begin(), single_.end(), symbol) != single_.end();
}

std::vector<Code> Alphabet::partners(Code code) const {
    std::vector<Code> out;
    for (std::size_t other = 0; other < size(); ++other) {
        if (can_pair(code, static_cast<Code>(other))) out.push_back(static_cast<Code>(other));
    }
    return out;
}

}  // namespace thermorna::rnastructure
