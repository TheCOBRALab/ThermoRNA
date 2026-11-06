#pragma once

/**
 * @file RNAstructure/params/datatable.hpp
 * @brief efn2's parameters for one base alphabet plus any number of modified-base alphabets.
 *
 * RNAstructure handles a modified nucleotide by switching to a whole alphabet that contains it
 * (efn2 --alphabet m6A), so efn2 can only use one modification at a time. A Datatable keeps the
 * base alphabet's tables (e.g. "rna", exactly as efn2 loads them) and merges the extra bases of
 * any number of modification alphabets ("m6A", "PSU", "1", custom ones) into one code space.
 *
 * Every lookup is answered from the alphabet that covers all of its nucleotides:
 *   1. only base nucleotides               -> the base alphabet's tables
 *   2. nucleotides one modification covers -> that modification's tables (the last one loaded
 *      wins when several cover the context)
 *   3. modified bases from different files -> each modified base is replaced by its unmodified
 *      fallback and the base tables are used, the rule Knotergy applies to modified bases that
 *      have no parameters for a context.
 * A Datatable with no modifications is identical to efn2's datatable; one whose base alphabet
 * is "m6A" is identical to `efn2 --alphabet m6A`.
 *
 * Each set can sit in a folder of its own (params/RNAstructure/m6A, .../m1Ψ) or in
 * RNAstructure's shared data_tables folder, where it is picked by file prefix.
 *
 * Codes 0 .. base_size()-1 are the base alphabet's codes (the same numbers as RNAstructure's
 * numseq), so constants such as C = 2 used in efn2 code stay valid. Accessors take the same
 * arguments as the arrays of RNAstructure's datatable. All lookups are const and thread-safe.
 */

#include "RNAstructure/constants.hpp"
#include "RNAstructure/params/alphabet_tables.hpp"

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace thermorna::rnastructure {

/// A modified-base parameter set to merge into a Datatable.
///
/// One folder per set: {"params/RNAstructure/m1Ψ"}. RNAstructure's shared data_tables folder:
/// {"", "1"} (the set with file prefix "1" in LoadOptions::data_directory).
struct Modification {
    /// Folder with the set's files; empty means LoadOptions::data_directory.
    std::string directory{};

    /// File prefix of the set ("m6A", "PSU", "1"); empty means the only set in the folder.
    std::string alphabet{};

    /// Unmodified base each new base falls back to, as {symbol, base symbol}, e.g. {"M", "A"}.
    /// When omitted it is inferred from the Pairing section (the canonical base with the same
    /// pairing partners: m6A -> A, pseudouridine -> U, m1Ψ -> U).
    std::vector<std::pair<std::string, std::string>> fallbacks{};

    /// Extra symbols, as {new symbol, existing symbol of this alphabet}, e.g. {"Ψ", "P"}.
    std::vector<std::pair<std::string, std::string>> aliases{};
};

/// What to load. The defaults are efn2's defaults.
struct LoadOptions {
    /// Folder with the base set (e.g. params/RNAstructure/turner2004), or RNAstructure's
    /// data_tables folder. Empty means the DATAPATH environment variable.
    std::string data_directory{};

    /// File prefix of the base set, like efn2 --alphabet ("rna", "dna", "m6A", ...). Empty
    /// means the only set in data_directory, or "rna" (efn2's default) when it holds several.
    std::string alphabet{};

    /// Temperature in Kelvin (efn2 --temperature). Tables are rescaled with the *.dh files
    /// when it differs from 310.15 K by 0.01 K or more.
    double temperature = kTemperature37C;

    /// Modified-base alphabets merged on top of the base alphabet, in priority order (later
    /// ones win where several cover the same context).
    std::vector<Modification> modifications{};

    /// What to do at temperatures other than 37 degrees C when an alphabet has no *.dh files.
    MissingEnthalpy missing_enthalpy = MissingEnthalpy::Fail;
};

/// Which tables answered a lookup.
struct LookupSource {
    enum class Kind {
        Base,          ///< the base alphabet
        Modification,  ///< modifications()[modification]
        Fallback,      ///< base alphabet after replacing modified bases by their fallbacks
    };
    Kind kind = Kind::Base;
    std::size_t modification = 0;
};

class Datatable {
   public:
    Datatable() = default;

    /// Reads the base alphabet and every modification (efn2's opendat + ScaleToTemperature).
    [[nodiscard]] static Datatable load(const LoadOptions& options);

    /// Merges tables that were already read. Modification tables should be at the temperature
    /// of @p base; a warning is recorded when they are not.
    [[nodiscard]] static Datatable merge(
        AlphabetTables base, std::vector<std::pair<AlphabetTables, Modification>> modifications);

    // ------------------------------- Alphabet ----------------------------------
    /// Number of codes (base bases plus every added modified base).
    [[nodiscard]] std::size_t size() const noexcept { return symbols_.size(); }
    [[nodiscard]] std::size_t base_size() const noexcept { return base_.size(); }

    /// True for bases added by a modification alphabet.
    [[nodiscard]] bool is_modified(Code code) const noexcept {
        const std::size_t c = code;
        return c < size() && c >= base_size();
    }

    /// Base code used for @p code when no loaded alphabet covers a context (itself for base
    /// codes).
    [[nodiscard]] Code fallback(Code code) const noexcept { return fallback_[code]; }

    /// Symbols of a code, primary symbol first.
    [[nodiscard]] const std::vector<std::string>& symbols(Code code) const {
        return symbols_.at(code);
    }
    [[nodiscard]] const std::string& symbol(Code code) const { return symbols(code).front(); }

    /// Name of the alphabet that defined @p code ("rna" for base codes, "m6A", ...).
    [[nodiscard]] const std::string& origin(Code code) const { return owner_.at(code); }

    /// Code of one grapheme, or kInvalidCode when it is unknown or ambiguous.
    [[nodiscard]] Code code(std::string_view symbol) const noexcept;

    /// True when several modifications use @p symbol for different bases (e.g. "M" is m6A in
    /// "m6A" and m1Ψ in "1"); such symbols must be written with an unambiguous alternative.
    [[nodiscard]] bool is_ambiguous(std::string_view symbol) const noexcept;

    /// Symbol is listed in a Single section (RNAstructure forces those nucleotides unpaired).
    [[nodiscard]] bool is_single_symbol(std::string_view symbol) const noexcept;

    /// RNAstructure's structure::IsNuc(i, c): @p symbol is one of the symbols of @p code.
    [[nodiscard]] bool has_symbol(Code code, std::string_view symbol) const noexcept;

    /// Canonical pairing, resolved like the tables (base, covering modification or fallback).
    [[nodiscard]] bool can_pair(Code a, Code b) const noexcept {
        const std::size_t n = size();
        const std::size_t x = a;
        const std::size_t y = b;
        return x < n && y < n && pairing_[x * n + y] != 0;
    }

    [[nodiscard]] bool is_linker(Code code) const noexcept { return linker_[code] != 0; }
    [[nodiscard]] bool is_non_interacting(Code code) const noexcept {
        return non_interacting_[code] != 0;
    }

    /// Encodes a sequence. Pure-ASCII input is mapped byte by byte; grapheme segmentation is
    /// only used when the sequence contains non-ASCII bytes. Throws ParameterError (with the
    /// 0-based position) for unknown or ambiguous symbols.
    [[nodiscard]] std::vector<Code> encode(std::string_view sequence) const;

    /// Encodes a sequence that is already split into graphemes (e.g. Knotergy's
    /// ProcessedRNAEntry::get_modified_sequence()).
    [[nodiscard]] std::vector<Code> encode(const std::vector<std::string_view>& graphemes) const;

    // ------------------------------- Scalars (base alphabet) -------------------
    [[nodiscard]] const MiscLoop& misc() const noexcept { return base_.misc(); }
    [[nodiscard]] const std::vector<Energy>& inter() const noexcept { return base_.inter(); }
    [[nodiscard]] const std::vector<Energy>& bulge() const noexcept { return base_.bulge(); }
    [[nodiscard]] const std::vector<Energy>& hairpin() const noexcept { return base_.hairpin(); }
    [[nodiscard]] float RT() const noexcept { return base_.RT(); }
    [[nodiscard]] double temperature() const noexcept { return base_.temperature(); }

    // ------------------------------- Tables ------------------------------------
    [[nodiscard]] Energy table4(Table4D table, Code a, Code b, Code c, Code d) const noexcept {
        std::array<Code, 4> k{a, b, c, d};
        return resolve(k).table4(table, k[0], k[1], k[2], k[3]);
    }
    [[nodiscard]] Energy stack(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::Stack, a, b, c, d);
    }
    [[nodiscard]] Energy helix_ends(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::HelixEnds, a, b, c, d);
    }
    [[nodiscard]] Energy tstkh(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::Tstkh, a, b, c, d);
    }
    [[nodiscard]] Energy tstki(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::Tstki, a, b, c, d);
    }
    [[nodiscard]] Energy tstki23(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::Tstki23, a, b, c, d);
    }
    [[nodiscard]] Energy tstki1n(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::Tstki1n, a, b, c, d);
    }
    [[nodiscard]] Energy coax(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::Coax, a, b, c, d);
    }
    [[nodiscard]] Energy tstackcoax(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::TstackCoax, a, b, c, d);
    }
    [[nodiscard]] Energy coaxstack(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::CoaxStack, a, b, c, d);
    }
    [[nodiscard]] Energy tstack(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::Tstack, a, b, c, d);
    }
    [[nodiscard]] Energy tstkm(Code a, Code b, Code c, Code d) const noexcept {
        return table4(Table4D::Tstkm, a, b, c, d);
    }

    /// dangle[a][b][c][d]: d = 1 is a 3' dangle, d = 2 a 5' dangle.
    [[nodiscard]] Energy dangle(Code a, Code b, Code c, std::size_t d) const noexcept {
        std::array<Code, 3> k{a, b, c};
        return resolve(k).dangle(k[0], k[1], k[2], d);
    }

    /// iloop11[a][b][c][d][e][f]: 5'abc3' / 3'def5'.
    [[nodiscard]] Energy iloop11(Code a, Code b, Code c, Code d, Code e, Code f) const noexcept {
        std::array<Code, 6> k{a, b, c, d, e, f};
        return resolve(k).iloop11(k[0], k[1], k[2], k[3], k[4], k[5]);
    }
    /// iloop21[a][b][c][d][e][f][g]: 5'ac f3' / 3'bdeg5'.
    [[nodiscard]] Energy iloop21(Code a, Code b, Code c, Code d, Code e, Code f,
                                 Code g) const noexcept {
        std::array<Code, 7> k{a, b, c, d, e, f, g};
        return resolve(k).iloop21(k[0], k[1], k[2], k[3], k[4], k[5], k[6]);
    }
    /// iloop22[a][b][c][d][e][f][g][h]: 5'aefb3' / 3'cghd5'.
    [[nodiscard]] Energy iloop22(Code a, Code b, Code c, Code d, Code e, Code f, Code g,
                                 Code h) const noexcept {
        std::array<Code, 8> k{a, b, c, d, e, f, g, h};
        return resolve(k).iloop22(k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]);
    }

    /// Special hairpins; @p codes points at the closing base i of the loop i..j
    /// (6 codes for tloop, 5 for triloop, 8 for hexaloop).
    [[nodiscard]] std::optional<int> tloop(const Code* codes) const noexcept {
        std::array<Code, 6> k{};
        return resolve_copy(codes, k).find_tloop(k.data());
    }
    [[nodiscard]] std::optional<int> triloop(const Code* codes) const noexcept {
        std::array<Code, 5> k{};
        return resolve_copy(codes, k).find_triloop(k.data());
    }
    [[nodiscard]] std::optional<int> hexaloop(const Code* codes) const noexcept {
        std::array<Code, 8> k{};
        return resolve_copy(codes, k).find_hexaloop(k.data());
    }

    /// Which tables a lookup involving @p count (at most 16) codes starting at @p codes is
    /// answered from.
    [[nodiscard]] LookupSource source(const Code* codes, std::size_t count) const noexcept;

    // ------------------------------- Introspection -----------------------------
    [[nodiscard]] const AlphabetTables& base() const noexcept { return base_; }
    [[nodiscard]] std::size_t modification_count() const noexcept { return overlays_.size(); }
    [[nodiscard]] const AlphabetTables& modification(std::size_t index) const {
        return overlays_.at(index).tables;
    }
    [[nodiscard]] const std::string& modification_name(std::size_t index) const {
        return overlays_.at(index).tables.alphabet().name();
    }
    /// Codes a modification added to the code space (in the order it listed them).
    [[nodiscard]] const std::vector<Code>& added_codes(std::size_t index) const {
        return overlays_.at(index).added;
    }
    /// Messages about unusual input (short table rows, ambiguous symbols, ...).
    [[nodiscard]] const std::vector<std::string>& warnings() const noexcept { return warnings_; }

   private:
    struct Overlay {
        AlphabetTables tables;
        std::array<Code, 256> to_local{};  ///< merged code -> code in tables, or kInvalidCode
        std::vector<Code> added;           ///< merged codes this alphabet introduced
    };

    template <std::size_t K>
    [[nodiscard]] const AlphabetTables& resolve(std::array<Code, K>& codes) const noexcept {
        unsigned char modified = 0;
        for (const Code c : codes) modified |= modified_[c];
        if (modified == 0) return base_;
        return resolve_modified(codes.data(), K, nullptr);
    }

    template <std::size_t K>
    [[nodiscard]] const AlphabetTables& resolve_copy(const Code* codes,
                                                     std::array<Code, K>& out) const noexcept {
        for (std::size_t k = 0; k < K; ++k) out[k] = codes[k];
        return resolve(out);
    }

    // Rewrites @p codes into the codes of the tables that answer the lookup.
    [[nodiscard]] const AlphabetTables& resolve_modified(Code* codes, std::size_t count,
                                                         LookupSource* source) const noexcept;

    [[noreturn]] void throw_symbol_error(std::string_view symbol, std::size_t position) const;

    // Merging steps (datatable.cpp).
    void add_modification(AlphabetTables tables, const Modification& spec);
    void build_symbol_lookup(const std::vector<Modification>& specs);
    void build_pairing();

    AlphabetTables base_;
    std::vector<Overlay> overlays_;
    std::vector<std::vector<std::string>> symbols_;
    std::vector<std::string> owner_;  // alphabet that defined each code, for messages
    std::vector<unsigned char> pairing_;
    std::array<unsigned char, 256> modified_{};  // 1 for added codes and for invalid codes
    std::array<Code, 256> fallback_{};
    std::array<unsigned char, 256> linker_{};
    std::array<unsigned char, 256> non_interacting_{};
    std::array<Code, 256> ascii_{};
    std::vector<std::pair<std::string, Code>> multibyte_;         // sorted
    std::vector<std::pair<std::string, std::string>> ambiguous_;  // symbol -> explanation
    std::vector<std::string> single_symbols_;
    std::vector<std::string> warnings_;
};

}  // namespace thermorna::rnastructure
