#pragma once

/**
 * @file RNAstructure/params/alphabet.hpp
 * @brief An RNAstructure folding alphabet (<prefix>.specification.dat).
 *
 * The specification lists the bases (each with its equivalent symbols), the canonical pairs,
 * the symbols forced single-stranded, the non-interacting bases and the intermolecular linker.
 * Codes are assigned in file order, exactly like RNAstructure's basetonum(), so code 0 is X/N,
 * 1 is A, 2 is C, 3 is G, 4 is U in the shipped alphabets.
 *
 * Symbols are grapheme clusters rather than bytes, so a specification may use "Ψ" or other
 * multi-byte symbols. For ASCII files this is the same as RNAstructure's character handling.
 */

#include "RNAstructure/constants.hpp"

#include <array>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace thermorna::rnastructure {

class Alphabet {
   public:
    Alphabet() = default;

    /// Parses the text of a specification file. @p source is used in error messages.
    [[nodiscard]] static Alphabet parse(std::string_view text, std::string name,
                                        const std::string& source);

    /// Reads and parses the specification of set @p prefix in @p directory
    /// (<prefix>.specification.dat); an empty prefix means the only set in the folder.
    [[nodiscard]] static Alphabet read(const std::string& directory,
                                       const std::string& prefix = "");

    /// The set's name: its folder's name when it has a folder of its own ("turner2004", "m1Ψ"),
    /// else its file prefix ("rna", "m6A").
    [[nodiscard]] const std::string& name() const noexcept { return name_; }

    /// Number of bases (codes are 0 .. size() - 1).
    [[nodiscard]] std::size_t size() const noexcept { return bases_.size(); }

    /// All symbols of a base in file order; the first is its primary symbol.
    [[nodiscard]] const std::vector<std::string>& symbols(Code code) const {
        return bases_.at(code);
    }

    /// Primary symbol of a base.
    [[nodiscard]] const std::string& primary_symbol(Code code) const {
        return symbols(code).front();
    }

    /// Code of a symbol, or kInvalidCode. The first base listing the symbol wins (basetonum).
    [[nodiscard]] Code code(std::string_view symbol) const noexcept;

    /// True when @p symbol is one of the symbols of @p code (RNAstructure's structure::IsNuc).
    [[nodiscard]] bool has_symbol(Code code, std::string_view symbol) const noexcept;

    /// The Pairing section.
    [[nodiscard]] bool can_pair(Code a, Code b) const noexcept {
        const std::size_t n = size();
        if (a >= n || b >= n) return false;
        return pairing_[static_cast<std::size_t>(a) * n + static_cast<std::size_t>(b)] != 0;
    }

    /// Non-interacting bases (neither pair nor stack, e.g. X/N).
    [[nodiscard]] bool is_non_interacting(Code code) const noexcept {
        return non_interacting_codes_[code] != 0;
    }

    /// Intermolecular linker bases (I).
    [[nodiscard]] bool is_linker(Code code) const noexcept { return linker_codes_[code] != 0; }

    /// can_interact as RNAstructure's processDat computes it: neither non-interacting nor linker.
    [[nodiscard]] bool can_interact(Code code) const noexcept {
        return !is_non_interacting(code) && !is_linker(code);
    }

    /// Symbols of the Single section (forced unpaired, e.g. lower case letters).
    [[nodiscard]] const std::vector<std::string>& single_symbols() const noexcept {
        return single_;
    }
    [[nodiscard]] bool is_single_symbol(std::string_view symbol) const noexcept;

    /// Symbols as listed in the Non-interacting and Linker sections.
    [[nodiscard]] const std::vector<std::string>& non_interacting_symbols() const noexcept {
        return non_interacting_;
    }
    [[nodiscard]] const std::vector<std::string>& linker_symbols() const noexcept {
        return linker_;
    }

    /// Codes listed in the Non-interacting / Linker sections, in file order (with repeats).
    [[nodiscard]] const std::vector<Code>& non_interacting_list() const noexcept {
        return non_interacting_list_;
    }
    [[nodiscard]] const std::vector<Code>& linker_list() const noexcept { return linker_list_; }

    /// Codes paired with @p code according to the Pairing section.
    [[nodiscard]] std::vector<Code> partners(Code code) const;

   private:
    std::string name_;
    std::vector<std::vector<std::string>> bases_;
    std::vector<unsigned char> pairing_;  // size() x size()
    std::vector<std::string> single_;
    std::vector<std::string> non_interacting_;
    std::vector<std::string> linker_;
    std::vector<Code> non_interacting_list_;
    std::vector<Code> linker_list_;
    std::array<unsigned char, 256> non_interacting_codes_{};
    std::array<unsigned char, 256> linker_codes_{};
    std::array<Code, 256> ascii_lookup_{};                 // single-byte symbols
    std::vector<std::pair<std::string, Code>> multibyte_;  // sorted, first occurrence wins
};

}  // namespace thermorna::rnastructure
