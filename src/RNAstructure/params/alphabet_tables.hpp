#pragma once

/**
 * @file RNAstructure/params/alphabet_tables.hpp
 * @brief Every nearest-neighbor table efn2 loads for one alphabet (parameter set).
 *
 * AlphabetTables::read() is RNAstructure's datatable::opendat(): it reads the 21 files
 * <prefix>.<table>.dg (or .dh) and applies the same post-processing (processDat). The set can
 * sit in RNAstructure's shared data_tables folder or in a folder of its own, with or without
 * the file prefix (see parameter_set.hpp). load() adds ScaleToTemperature(), i.e. it reads the
 * .dh tables and recomputes every free energy with the same float arithmetic as RNAstructure's
 * dG_T()/Tscale(). The resulting values are identical to the ones in efn2's datatable,
 * including RNAstructure's handling of short rows and of "." entries.
 *
 * Accessors take the same arguments, in the same order, as the arrays of RNAstructure's
 * datatable (data->stack[a][b][c][d] becomes tables.stack(a, b, c, d)), so efn2 code ports
 * directly. Values are in tenths of kcal/mol; kInfiniteEnergy marks forbidden entries.
 *
 * Storage is compact: the 4-index tables are dense, the internal-loop tables keep one dense
 * block per closing-pair combination present in the files (efn2 allocates N^6..N^8 entries).
 */

#include "RNAstructure/constants.hpp"
#include "RNAstructure/params/alphabet.hpp"
#include "RNAstructure/params/parameter_set.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace thermorna::rnastructure {

/// Which file set a table set comes from.
enum class TableKind {
    FreeEnergy,  ///< *.dg, free energies at 37 degrees C
    Enthalpy,    ///< *.dh, enthalpies used to rescale to other temperatures
};

/// Policy for a temperature other than 37 degrees C when an alphabet has no *.dh files.
enum class MissingEnthalpy {
    Fail,            ///< throw ParameterError, as efn2 does
    KeepFreeEnergy,  ///< keep the 37 degrees C values and record a warning
};

/// The 4-index tables, named after RNAstructure's datatable members.
enum class Table4D : std::uint8_t {
    Stack,       ///< stack       (<alphabet>.stack)
    HelixEnds,   ///< helix_ends  (<alphabet>.helix_ends)
    Tstkh,       ///< tstkh       (<alphabet>.tstackh)
    Tstki,       ///< tstki       (<alphabet>.tstacki)
    Tstki23,     ///< tstki23     (<alphabet>.tstacki23)
    Tstki1n,     ///< tstki1n     (<alphabet>.tstacki1n)
    Coax,        ///< coax        (<alphabet>.coaxial)
    TstackCoax,  ///< tstackcoax  (<alphabet>.tstackcoax)
    CoaxStack,   ///< coaxstack   (<alphabet>.coaxstack)
    Tstack,      ///< tstack      (<alphabet>.tstack)
    Tstkm,       ///< tstkm       (<alphabet>.tstackm)
};
inline constexpr std::size_t kTable4DCount = 11;

/// File name (without alphabet prefix and extension) of a 4-index table.
[[nodiscard]] const char* table_file_name(Table4D table) noexcept;

/// Scalars from <alphabet>.miscloop, with RNAstructure's member names and units.
struct MiscLoop {
    float prelog = 0.0f;  ///< large-loop extrapolation, tenths of kcal/mol (kept as float)
    Energy maxpen = 0;    ///< Ninio maximum
    std::array<Energy, 5> poppen{};   ///< Ninio f(m); [1]..[4] are read
    std::array<Energy, 11> eparam{};  ///< [5] offset, [6] per nucleotide, [10] per helix
                                      ///< (multibranch); [7] = 30, [8] = 30, [9] = -500
    Energy efn2a = 0;                 ///< efn2 multibranch offset
    Energy efn2b = 0;                 ///< efn2 multibranch per unpaired nucleotide
    Energy efn2c = 0;                 ///< efn2 multibranch per helix
    Energy mlasym = 0;                ///< multibranch asymmetry
    Energy strain = 0;                ///< three-way junction strain
    Energy gubonus = 0;               ///< GGG hairpin bonus
    Energy cslope = 0;                ///< poly-C hairpin slope
    Energy cint = 0;                  ///< poly-C hairpin intercept
    Energy c3 = 0;                    ///< poly-C triloop penalty
    Energy init = 0;                  ///< intermolecular initiation
    Energy singlecbulge = 0;          ///< single C bulge bonus
};

/// One entry of the tloop/triloop/hexaloop tables, in file order.
struct SpecialHairpin {
    std::string sequence;     ///< as written in the file (closing pair included)
    std::vector<Code> codes;  ///< sequence encoded with the alphabet
    int energy = 0;           ///< tenths of kcal/mol (RNAstructure keeps these as int)
};

class AlphabetTables {
   public:
    AlphabetTables() = default;

    /// datatable::opendat(): reads all tables of one kind of the set with file prefix
    /// @p alphabet ("rna", "m6A") in @p directory. An empty prefix means the only set in the
    /// folder (e.g. params/RNAstructure/m6A).
    [[nodiscard]] static AlphabetTables read(const std::string& directory,
                                             const std::string& alphabet, TableKind kind);

    /// opendat() at 37 degrees C followed by ScaleToTemperature(@p kelvin) when needed.
    [[nodiscard]] static AlphabetTables load(const std::string& directory,
                                             const std::string& alphabet = "",
                                             double kelvin = kTemperature37C,
                                             MissingEnthalpy missing = MissingEnthalpy::Fail);

    /// dG_T(): recomputes the free energies for @p kelvin from these 37 degrees C values and
    /// @p enthalpy (tables of the same alphabet read with TableKind::Enthalpy).
    void scale_to_temperature(const AlphabetTables& enthalpy, double kelvin);

    // ------------------------------- Description -------------------------------
    /// The set's name: its folder's name when the folder holds only this set ("m1Ψ"),
    /// otherwise its file prefix ("1"). Also alphabet().name().
    [[nodiscard]] const std::string& name() const noexcept { return alphabet_.name(); }
    /// File prefix ("rna", "1"; empty for unprefixed file names).
    [[nodiscard]] const std::string& prefix() const noexcept { return prefix_; }
    [[nodiscard]] const Alphabet& alphabet() const noexcept { return alphabet_; }
    [[nodiscard]] std::size_t size() const noexcept { return n_; }
    [[nodiscard]] TableKind kind() const noexcept { return kind_; }
    [[nodiscard]] const std::string& directory() const noexcept { return directory_; }
    /// Temperature (K) the free energies correspond to.
    [[nodiscard]] double temperature() const noexcept { return temperature_; }
    /// datatable::RT, scaled with the temperature like RNAstructure does.
    [[nodiscard]] float RT() const noexcept { return RT_; }
    /// Notes about unusual input that was read the way RNAstructure reads it.
    [[nodiscard]] const std::vector<std::string>& warnings() const noexcept { return warnings_; }

    // ------------------------------- Scalars -----------------------------------
    [[nodiscard]] const MiscLoop& misc() const noexcept { return misc_; }
    /// Loop initiation by size; index 0 is kInfiniteEnergy and index k is line k of the file.
    [[nodiscard]] const std::vector<Energy>& inter() const noexcept { return inter_; }
    [[nodiscard]] const std::vector<Energy>& bulge() const noexcept { return bulge_; }
    [[nodiscard]] const std::vector<Energy>& hairpin() const noexcept { return hairpin_; }

    // ------------------------------- 4-index tables ----------------------------
    [[nodiscard]] Energy table4(Table4D table, Code a, Code b, Code c, Code d) const noexcept {
        return table4_[table_offset(table) + index4(a, b, c, d)];
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

    /// dangle[a][b][c][d]: d = 1 is a 3' dangle and d = 2 a 5' dangle (d = 0 is unused).
    [[nodiscard]] Energy dangle(Code a, Code b, Code c, std::size_t d) const noexcept {
        return dangle_[index3(a, b, c) * 3 + d];
    }

    // ------------------------------- Internal loops -----------------------------
    /// iloop11[a][b][c][d][e][f]: 5'abc3' / 3'def5', pairs a-d and c-f.
    [[nodiscard]] Energy iloop11(Code a, Code b, Code c, Code d, Code e, Code f) const noexcept {
        return lookup(iloop11_, index4(a, c, d, f), index2(b, e));
    }
    /// iloop21[a][b][c][d][e][f][g]: 5'ac f3' / 3'bdeg5', pairs a-b and f-g.
    [[nodiscard]] Energy iloop21(Code a, Code b, Code c, Code d, Code e, Code f,
                                 Code g) const noexcept {
        return lookup(iloop21_, index4(a, b, f, g), index3(c, d, e));
    }
    /// iloop22[a][b][c][d][e][f][g][h]: 5'aefb3' / 3'cghd5', pairs a-c and b-d.
    [[nodiscard]] Energy iloop22(Code a, Code b, Code c, Code d, Code e, Code f, Code g,
                                 Code h) const noexcept {
        return lookup(iloop22_, index4(a, b, c, d), index4(e, f, g, h));
    }

    // ------------------------------- Special hairpins --------------------------
    [[nodiscard]] const std::vector<SpecialHairpin>& tloop() const noexcept {
        return tloop_.entries;
    }
    [[nodiscard]] const std::vector<SpecialHairpin>& triloop() const noexcept {
        return triloop_.entries;
    }
    [[nodiscard]] const std::vector<SpecialHairpin>& hexaloop() const noexcept {
        return hexaloop_.entries;
    }

    /// Energy of the special hairpin whose @p length codes start at @p codes (closing pair
    /// included: 6 for tloop, 5 for triloop, 8 for hexaloop). The first matching file entry
    /// wins, as in efn2.
    [[nodiscard]] std::optional<int> find_tloop(const Code* codes) const noexcept {
        return find(tloop_, codes, 6);
    }
    [[nodiscard]] std::optional<int> find_triloop(const Code* codes) const noexcept {
        return find(triloop_, codes, 5);
    }
    [[nodiscard]] std::optional<int> find_hexaloop(const Code* codes) const noexcept {
        return find(hexaloop_, codes, 8);
    }

   private:
    friend class TableReader;

    // Internal-loop table: one dense block per closing-pair combination present in the file.
    struct BlockTable {
        std::vector<std::int32_t> block_of;  // size^4 closing-pair combinations, -1 = none
        std::vector<Energy> values;          // block_size values per present block
        std::size_t block_size = 0;
    };

    struct SpecialHairpinTable {
        std::vector<SpecialHairpin> entries;
        std::vector<std::pair<std::uint64_t, std::uint32_t>> index;  // key -> entry, sorted
    };

    // Packs up to 8 codes into a lookup key (0 when @p length is larger).
    [[nodiscard]] static std::uint64_t pack_codes(const Code* codes, std::size_t length) noexcept;

    [[nodiscard]] std::size_t index2(Code a, Code b) const noexcept {
        return static_cast<std::size_t>(a) * n_ + static_cast<std::size_t>(b);
    }
    [[nodiscard]] std::size_t index3(Code a, Code b, Code c) const noexcept {
        return index2(a, b) * n_ + static_cast<std::size_t>(c);
    }
    [[nodiscard]] std::size_t index4(Code a, Code b, Code c, Code d) const noexcept {
        return index3(a, b, c) * n_ + static_cast<std::size_t>(d);
    }
    [[nodiscard]] std::size_t table_offset(Table4D table) const noexcept {
        return static_cast<std::size_t>(table) * n4_;
    }
    [[nodiscard]] static Energy lookup(const BlockTable& table, std::size_t block,
                                       std::size_t inner) noexcept {
        const std::int32_t b = table.block_of[block];
        if (b < 0) return kInfiniteEnergy;
        return table.values[static_cast<std::size_t>(b) * table.block_size + inner];
    }
    [[nodiscard]] static std::optional<int> find(const SpecialHairpinTable& table,
                                                 const Code* codes, std::size_t length) noexcept;

    void process_dat();

    Alphabet alphabet_;
    TableKind kind_ = TableKind::FreeEnergy;
    std::string directory_;
    std::string prefix_;
    double temperature_ = kTemperature37C;
    float RT_ = kRT37C;
    std::size_t n_ = 0;
    std::size_t n4_ = 0;

    MiscLoop misc_;
    std::vector<Energy> inter_;
    std::vector<Energy> bulge_;
    std::vector<Energy> hairpin_;
    std::vector<Energy> table4_;  // kTable4DCount x n^4
    std::vector<Energy> dangle_;  // n^3 x 3
    BlockTable iloop11_;
    BlockTable iloop21_;
    BlockTable iloop22_;
    SpecialHairpinTable tloop_;
    SpecialHairpinTable triloop_;
    SpecialHairpinTable hexaloop_;
    std::vector<std::string> warnings_;
};

}  // namespace thermorna::rnastructure
