#include "RNAstructure/params/datatable.hpp"

#include "RNAstructure/error.hpp"
#include "RNAstructure/params/parameter_set.hpp"
#include "RNAstructure/utils/graphemes.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>
#include <sstream>

namespace thermorna::rnastructure {

namespace {

// Marks single-byte symbols that several modifications use for different bases.
constexpr Code kAmbiguousCode = 0xFE;

std::string resolve_data_directory(const std::string& requested) {
    if (!requested.empty()) return requested;
    const char* environment = std::getenv("DATAPATH");
    if (environment != nullptr && *environment != '\0') return environment;
    throw ParameterError(
        "No RNAstructure data_tables directory was given and the DATAPATH environment variable "
        "is not set.");
}

// Error for a modification folder that does not exist, with a hint when @p folder is the name
// of a set in @p base_directory or of a folder next to it (params/RNAstructure/<folder>).
std::string missing_folder_message(const std::string& folder, const std::string& base_directory) {
    std::string message = "The parameter folder '" + folder + "' of a modification does not exist.";
    const std::vector<std::string> sets = list_parameter_sets(base_directory);
    if (std::find(sets.begin(), sets.end(), folder) != sets.end()) {
        return message + " '" + base_directory + "' has a set with file prefix '" + folder +
               "'; name it with Modification::alphabet ({\"\", \"" + folder + "\"}).";
    }
    std::string parent = base_directory;
    while (parent.size() > 1 && parent.back() == '/') parent.pop_back();
    const std::size_t slash = parent.find_last_of('/');
    if (folder.find('/') == std::string::npos && slash != std::string::npos) {
        const std::string sibling = parent.substr(0, slash + 1) + folder;
        if (is_directory(sibling)) message += " Did you mean '" + sibling + "'?";
    }
    return message;
}

std::string quoted_list(const std::vector<std::string>& symbols, std::size_t skip_first) {
    std::string out;
    for (std::size_t k = skip_first; k < symbols.size(); ++k) {
        if (!out.empty()) out += ", ";
        out += "'" + symbols[k] + "'";
    }
    return out;
}

bool contains(const std::vector<Code>& codes, Code code) {
    return std::find(codes.begin(), codes.end(), code) != codes.end();
}

}  // namespace

// ----------------------------------------------------------------------------
// Loading and merging
// ----------------------------------------------------------------------------

Datatable Datatable::load(const LoadOptions& options) {
    const std::string directory = resolve_data_directory(options.data_directory);

    // efn2 defaults to "rna" in RNAstructure's data_tables; a folder of its own needs no name.
    std::string alphabet = options.alphabet;
    if (alphabet.empty()) {
        const std::vector<std::string> sets = list_parameter_sets(directory);
        if (sets.size() > 1 && std::find(sets.begin(), sets.end(), "rna") != sets.end()) {
            alphabet = "rna";
        }
    }
    AlphabetTables base =
        AlphabetTables::load(directory, alphabet, options.temperature, options.missing_enthalpy);

    std::vector<std::pair<AlphabetTables, Modification>> modifications;
    modifications.reserve(options.modifications.size());
    for (const Modification& modification : options.modifications) {
        if (modification.directory.empty() && modification.alphabet.empty()) {
            throw ParameterError(
                "A modification needs a folder (Modification::directory) or a file prefix "
                "(Modification::alphabet).");
        }
        if (!modification.directory.empty() && !is_directory(modification.directory)) {
            throw ParameterError(missing_folder_message(modification.directory, directory));
        }
        const std::string& where =
            modification.directory.empty() ? directory : modification.directory;
        modifications.emplace_back(
            AlphabetTables::load(where, modification.alphabet, options.temperature,
                                 options.missing_enthalpy),
            modification);
    }
    return merge(std::move(base), std::move(modifications));
}

Datatable Datatable::merge(AlphabetTables base,
                           std::vector<std::pair<AlphabetTables, Modification>> modifications) {
    Datatable table;
    table.base_ = std::move(base);
    const Alphabet& alphabet = table.base_.alphabet();

    table.modified_.fill(1);
    table.fallback_.fill(0);
    for (std::size_t c = 0; c < alphabet.size(); ++c) {
        const Code code = static_cast<Code>(c);
        table.symbols_.push_back(alphabet.symbols(code));
        table.owner_.push_back(alphabet.name());
        table.modified_[c] = 0;
        table.fallback_[c] = code;
        table.linker_[c] = alphabet.is_linker(code) ? 1 : 0;
        table.non_interacting_[c] = alphabet.is_non_interacting(code) ? 1 : 0;
    }
    table.warnings_ = table.base_.warnings();

    std::vector<Modification> specs;
    specs.reserve(modifications.size());
    for (auto& entry : modifications) {
        table.add_modification(std::move(entry.first), entry.second);
        specs.push_back(std::move(entry.second));
    }
    table.build_symbol_lookup(specs);
    table.build_pairing();
    return table;
}

void Datatable::add_modification(AlphabetTables tables, const Modification& spec) {
    const Alphabet& base_alphabet = base_.alphabet();
    const Alphabet& local = tables.alphabet();
    const std::string& name = tables.name();
    const std::size_t base_count = base_alphabet.size();

    Overlay overlay;
    overlay.to_local.fill(kInvalidCode);
    std::vector<Code> merged(local.size(), kInvalidCode);
    std::vector<bool> is_new(local.size(), false);
    bool covers_added_bases = false;

    for (std::size_t l = 0; l < local.size(); ++l) {
        const Code lc = static_cast<Code>(l);
        const std::vector<std::string>& symbols = local.symbols(lc);

        // A base of the base alphabet when they share a symbol (X = N, A = a, U = T, ...).
        std::vector<Code> matches;
        for (const std::string& symbol : symbols) {
            const Code c = base_alphabet.code(symbol);
            if (c != kInvalidCode && !contains(matches, c)) matches.push_back(c);
        }
        if (matches.size() > 1) {
            throw ParameterError("Base '" + symbols.front() + "' of alphabet '" + name +
                                 "' shares symbols with several bases of alphabet '" +
                                 base_alphabet.name() + "' ('" +
                                 base_alphabet.primary_symbol(matches[0]) + "' and '" +
                                 base_alphabet.primary_symbol(matches[1]) + "').");
        }
        if (matches.size() == 1) {
            const Code c = matches.front();
            if (local.is_linker(lc) != base_alphabet.is_linker(c) ||
                local.is_non_interacting(lc) != base_alphabet.is_non_interacting(c)) {
                throw ParameterError(
                    "Base '" + symbols.front() + "' of alphabet '" + name + "' shares a symbol " +
                    "with '" + base_alphabet.primary_symbol(c) + "' of alphabet '" +
                    base_alphabet.name() +
                    "', but one of them is a linker or non-interacting base and the other is " +
                    "not. Give the modified base symbols that alphabet '" + base_alphabet.name() +
                    "' does not use.");
            }
            merged[l] = c;
            continue;
        }

        // A base an earlier modification added, recognised by its primary symbol; this
        // alphabet then also provides parameters for it (and wins where it covers a context).
        Code existing = kInvalidCode;
        for (std::size_t e = base_count; e < symbols_.size(); ++e) {
            if (symbols_[e].front() == symbols.front()) {
                existing = static_cast<Code>(e);
                break;
            }
        }
        if (existing != kInvalidCode) {
            merged[l] = existing;
            covers_added_bases = true;
            for (const std::string& symbol : symbols) {
                std::vector<std::string>& known = symbols_[existing];
                if (std::find(known.begin(), known.end(), symbol) == known.end()) {
                    known.push_back(symbol);
                }
            }
            continue;
        }

        if (symbols_.size() >= kMaxCodes) {
            throw ParameterError("Too many distinct bases: at most " + std::to_string(kMaxCodes) +
                                 " are supported.");
        }
        const Code added = static_cast<Code>(symbols_.size());
        symbols_.push_back(symbols);
        owner_.push_back(name);
        modified_[added] = 1;
        linker_[added] = local.is_linker(lc) ? 1 : 0;
        non_interacting_[added] = local.is_non_interacting(lc) ? 1 : 0;
        merged[l] = added;
        is_new[l] = true;
        overlay.added.push_back(added);
        covers_added_bases = true;
    }

    for (std::size_t l = 0; l < local.size(); ++l) {
        Code& slot = overlay.to_local[merged[l]];
        if (slot != kInvalidCode) {
            throw ParameterError("Bases '" + local.primary_symbol(slot) + "' and '" +
                                 local.primary_symbol(static_cast<Code>(l)) + "' of alphabet '" +
                                 name + "' both correspond to base '" +
                                 symbols_[merged[l]].front() + "'.");
        }
        slot = static_cast<Code>(l);
    }

    // Fallbacks given explicitly.
    std::vector<bool> has_fallback(local.size(), false);
    for (const auto& [symbol, parent] : spec.fallbacks) {
        const Code lc = local.code(symbol);
        if (lc == kInvalidCode) {
            throw ParameterError("Fallback given for '" + symbol +
                                 "', which is not a base of alphabet '" + name + "'.");
        }
        const Code target = merged[lc];
        if (static_cast<std::size_t>(target) < base_count) {
            throw ParameterError("Fallback given for '" + symbol + "', but it is base '" +
                                 base_alphabet.primary_symbol(target) + "' of alphabet '" +
                                 base_alphabet.name() + "', not a modified base.");
        }
        const Code to = base_alphabet.code(parent);
        if (to == kInvalidCode) {
            throw ParameterError("Fallback '" + parent + "' for '" + symbol +
                                 "' is not a base of alphabet '" + base_alphabet.name() + "'.");
        }
        if (!is_new[lc] && fallback_[target] != to) {
            throw ParameterError("Base '" + symbols_[target].front() + "' already falls back to '" +
                                 base_alphabet.primary_symbol(fallback_[target]) + "'; alphabet '" +
                                 name + "' asks for '" + parent + "'.");
        }
        fallback_[target] = to;
        has_fallback[lc] = true;
    }

    // Inferred fallbacks: the interacting base of the base alphabet that pairs with exactly the
    // same bases (m6A pairs with U like A; pseudouridine pairs with A and G like U).
    for (std::size_t l = 0; l < local.size(); ++l) {
        if (!is_new[l] || has_fallback[l]) continue;
        const Code lc = static_cast<Code>(l);
        const Code target = merged[l];

        auto first_with = [&](bool linker) {
            for (std::size_t c = 0; c < base_count; ++c) {
                const Code code = static_cast<Code>(c);
                if (linker ? base_alphabet.is_linker(code)
                           : base_alphabet.is_non_interacting(code)) {
                    return code;
                }
            }
            return kInvalidCode;
        };
        if (local.is_linker(lc) || local.is_non_interacting(lc)) {
            const Code to = first_with(local.is_linker(lc));
            if (to == kInvalidCode) {
                throw ParameterError("Base '" + symbols_[target].front() + "' of alphabet '" +
                                     name + "' has no counterpart in alphabet '" +
                                     base_alphabet.name() + "'.");
            }
            fallback_[target] = to;
            continue;
        }

        std::vector<Code> partners;
        for (const Code partner : local.partners(lc)) {
            const Code m = merged[partner];
            if (static_cast<std::size_t>(m) < base_count && !contains(partners, m)) {
                partners.push_back(m);
            }
        }
        std::sort(partners.begin(), partners.end());

        std::vector<Code> candidates;
        if (!partners.empty()) {
            for (std::size_t c = 0; c < base_count; ++c) {
                const Code code = static_cast<Code>(c);
                if (!base_alphabet.can_interact(code)) continue;
                if (base_alphabet.partners(code) == partners) candidates.push_back(code);
            }
        }
        if (candidates.size() != 1) {
            throw ParameterError(
                "Cannot infer which unmodified base '" + symbols_[target].front() +
                "' of alphabet '" + name + "' falls back to (" +
                (candidates.empty() ? std::string("no base of '") + base_alphabet.name() +
                                          "' has the same pairing partners"
                                    : std::string("several bases have the same pairing partners")) +
                "). Set Modification::fallbacks, e.g. {\"" + symbols_[target].front() +
                "\", \"A\"}.");
        }
        fallback_[target] = candidates.front();
    }

    if (!covers_added_bases) {
        warnings_.push_back("Alphabet '" + name + "' adds no bases to alphabet '" +
                            base_alphabet.name() + "', so its tables are never used.");
    }
    if (std::abs(tables.temperature() - base_.temperature()) >= kTemperatureTolerance) {
        std::ostringstream message;
        message << "Alphabet '" << name << "' is at " << tables.temperature()
                << " K while alphabet '" << base_alphabet.name() << "' is at "
                << base_.temperature() << " K.";
        warnings_.push_back(message.str());
    }
    for (const std::string& warning : tables.warnings()) warnings_.push_back(warning);

    overlay.tables = std::move(tables);
    overlays_.push_back(std::move(overlay));
}

void Datatable::build_symbol_lookup(const std::vector<Modification>& specs) {
    const std::size_t base_count = base_.size();

    // symbol -> code; kAmbiguousCode once two added bases claim it. Base symbols come first
    // and, like basetonum(), the first base listing a symbol owns it.
    std::map<std::string, Code> registry;
    std::map<std::string, std::vector<Code>> claims;
    for (std::size_t c = 0; c < symbols_.size(); ++c) {
        const Code code = static_cast<Code>(c);
        for (const std::string& symbol : symbols_[c]) {
            auto [it, inserted] = registry.emplace(symbol, code);
            if (c < base_count) continue;
            std::vector<Code>& claimants = claims[symbol];
            if (!contains(claimants, code)) claimants.push_back(code);
            const bool claimed_by_other_added_base =
                !inserted && it->second != code && it->second != kAmbiguousCode &&
                static_cast<std::size_t>(it->second) >= base_count;
            if (claimed_by_other_added_base) it->second = kAmbiguousCode;
        }
    }

    // Aliases requested in Modification::aliases override symbols taken from the files.
    std::map<std::string, Code> requested;
    for (std::size_t o = 0; o < specs.size(); ++o) {
        const Alphabet& local = overlays_[o].tables.alphabet();
        for (const auto& [alias, target] : specs[o].aliases) {
            const Code lc = local.code(target);
            if (lc == kInvalidCode) {
                throw ParameterError("Alias '" + alias + "' refers to '" + target +
                                     "', which is not a base of alphabet '" + local.name() + "'.");
            }
            if (graphemes::split(alias).size() != 1) {
                throw ParameterError("Alias '" + alias + "' must be a single grapheme.");
            }
            Code code = kInvalidCode;
            for (std::size_t m = 0; m < overlays_[o].to_local.size(); ++m) {
                if (overlays_[o].to_local[m] == lc) code = static_cast<Code>(m);
            }
            const auto existing = registry.find(alias);
            if (existing != registry.end() && existing->second != kAmbiguousCode &&
                static_cast<std::size_t>(existing->second) < base_count &&
                existing->second != code) {
                throw ParameterError("Alias '" + alias + "' is already a symbol of base '" +
                                     symbols_[existing->second].front() + "' of alphabet '" +
                                     base_.alphabet().name() + "'.");
            }
            const auto earlier = requested.find(alias);
            if (earlier != requested.end() && earlier->second != code) {
                throw ParameterError("Alias '" + alias +
                                     "' was requested for two different bases.");
            }
            requested[alias] = code;
            registry[alias] = code;
            std::vector<std::string>& known = symbols_[code];
            if (std::find(known.begin(), known.end(), alias) == known.end()) known.push_back(alias);
        }
    }

    ascii_.fill(kInvalidCode);
    multibyte_.clear();
    ambiguous_.clear();
    for (const auto& [symbol, code] : registry) {
        if (code == kAmbiguousCode) {
            std::string explanation;
            for (const Code claimant : claims[symbol]) {
                if (!explanation.empty()) explanation += " and ";
                explanation +=
                    "'" + symbols_[claimant].front() + "' of alphabet '" + owner_[claimant] + "'";
                const std::string others = quoted_list(symbols_[claimant], 1);
                if (!others.empty()) explanation += " (also written " + others + ")";
            }
            ambiguous_.emplace_back(symbol, explanation);
            warnings_.push_back("Symbol '" + symbol + "' stands for " + explanation +
                                "; it is rejected in sequences, use another symbol of the base "
                                "you mean.");
        }
        if (symbol.size() == 1) {
            ascii_[static_cast<unsigned char>(symbol[0])] = code;
        } else if (code != kAmbiguousCode) {
            multibyte_.emplace_back(symbol, code);
        }
    }
    std::sort(multibyte_.begin(), multibyte_.end());
    std::sort(ambiguous_.begin(), ambiguous_.end());

    // Symbols forced single-stranded: the base alphabet's, plus those of every modification
    // that name a base the modification added.
    single_symbols_ = base_.alphabet().single_symbols();
    for (const Overlay& overlay : overlays_) {
        const Alphabet& local = overlay.tables.alphabet();
        for (const std::string& symbol : local.single_symbols()) {
            const Code lc = local.code(symbol);
            bool added = false;
            for (const Code code : overlay.added) added = added || overlay.to_local[code] == lc;
            if (added && !is_single_symbol(symbol)) single_symbols_.push_back(symbol);
        }
    }
}

void Datatable::build_pairing() {
    const std::size_t n = size();
    pairing_.assign(n * n, 0);
    for (std::size_t a = 0; a < n; ++a) {
        for (std::size_t b = 0; b < n; ++b) {
            std::array<Code, 2> k{static_cast<Code>(a), static_cast<Code>(b)};
            const AlphabetTables& tables = resolve(k);
            pairing_[a * n + b] = tables.alphabet().can_pair(k[0], k[1]) ? 1 : 0;
        }
    }
}

// ----------------------------------------------------------------------------
// Symbols
// ----------------------------------------------------------------------------

Code Datatable::code(std::string_view symbol) const noexcept {
    if (symbol.size() == 1) {
        const Code c = ascii_[static_cast<unsigned char>(symbol[0])];
        return c == kAmbiguousCode ? kInvalidCode : c;
    }
    const auto it =
        std::lower_bound(multibyte_.begin(), multibyte_.end(), symbol,
                         [](const std::pair<std::string, Code>& entry, std::string_view key) {
                             return std::string_view(entry.first) < key;
                         });
    if (it != multibyte_.end() && it->first == symbol) return it->second;
    return kInvalidCode;
}

bool Datatable::is_ambiguous(std::string_view symbol) const noexcept {
    const auto it =
        std::lower_bound(ambiguous_.begin(), ambiguous_.end(), symbol,
                         [](const std::pair<std::string, std::string>& entry,
                            std::string_view key) { return std::string_view(entry.first) < key; });
    return it != ambiguous_.end() && it->first == symbol;
}

bool Datatable::is_single_symbol(std::string_view symbol) const noexcept {
    return std::find(single_symbols_.begin(), single_symbols_.end(), symbol) !=
           single_symbols_.end();
}

bool Datatable::has_symbol(Code code, std::string_view symbol) const noexcept {
    if (static_cast<std::size_t>(code) >= symbols_.size()) return false;
    const std::vector<std::string>& known = symbols_[code];
    return std::find(known.begin(), known.end(), symbol) != known.end();
}

void Datatable::throw_symbol_error(std::string_view symbol, std::size_t position) const {
    const std::string where = "at index " + std::to_string(position);
    for (const auto& [ambiguous, explanation] : ambiguous_) {
        if (ambiguous == symbol) {
            throw ParameterError("Base '" + ambiguous + "' " + where +
                                 " is ambiguous: it stands for " + explanation + ".");
        }
    }
    std::string alphabets = "'" + base_.alphabet().name() + "'";
    for (const Overlay& overlay : overlays_)
        alphabets += ", '" + overlay.tables.alphabet().name() + "'";
    throw ParameterError("Base '" + std::string(symbol) + "' " + where +
                         " is not a symbol of the loaded alphabets (" + alphabets + ").");
}

std::vector<Code> Datatable::encode(std::string_view sequence) const {
    if (graphemes::is_ascii(sequence)) {
        std::vector<Code> codes(sequence.size());
        for (std::size_t i = 0; i < sequence.size(); ++i) {
            const Code c = ascii_[static_cast<unsigned char>(sequence[i])];
            if (c >= kMaxCodes) throw_symbol_error(sequence.substr(i, 1), i);
            codes[i] = c;
        }
        return codes;
    }

    // Only the spans around non-ASCII bytes are segmented by uni-algo.
    std::vector<Code> codes;
    codes.reserve(sequence.size());
    graphemes::for_each(sequence, [&](std::string_view g) {
        Code c = g.size() == 1 ? ascii_[static_cast<unsigned char>(g[0])] : code(g);
        if (c >= kMaxCodes) throw_symbol_error(g, codes.size());
        codes.push_back(c);
    });
    return codes;
}

std::vector<Code> Datatable::encode(const std::vector<std::string_view>& graphemes) const {
    std::vector<Code> codes;
    codes.reserve(graphemes.size());
    for (std::size_t i = 0; i < graphemes.size(); ++i) {
        const Code c = code(graphemes[i]);
        if (c == kInvalidCode) throw_symbol_error(graphemes[i], i);
        codes.push_back(c);
    }
    return codes;
}

// ----------------------------------------------------------------------------
// Lookups involving modified bases
// ----------------------------------------------------------------------------

const AlphabetTables& Datatable::resolve_modified(Code* codes, std::size_t count,
                                                  LookupSource* source) const noexcept {
    for (std::size_t o = overlays_.size(); o-- > 0;) {
        const Overlay& overlay = overlays_[o];
        bool covered = true;
        for (std::size_t k = 0; k < count && covered; ++k) {
            covered = overlay.to_local[codes[k]] != kInvalidCode;
        }
        if (!covered) continue;
        for (std::size_t k = 0; k < count; ++k) codes[k] = overlay.to_local[codes[k]];
        if (source != nullptr) *source = {LookupSource::Kind::Modification, o};
        return overlay.tables;
    }
    for (std::size_t k = 0; k < count; ++k) codes[k] = fallback_[codes[k]];
    if (source != nullptr) *source = {LookupSource::Kind::Fallback, 0};
    return base_;
}

LookupSource Datatable::source(const Code* codes, std::size_t count) const noexcept {
    std::array<Code, 16> copy{};
    count = std::min(count, copy.size());
    unsigned char modified = 0;
    for (std::size_t k = 0; k < count; ++k) {
        copy[k] = codes[k];
        modified |= modified_[codes[k]];
    }
    LookupSource result;
    if (modified != 0) (void) resolve_modified(copy.data(), count, &result);
    return result;
}

}  // namespace thermorna::rnastructure
