// Unit tests of ThermoRNA's RNAstructure parameter loader.
//
// They read the parameter sets shipped in params/RNAstructure (one folder per set) and a copy
// of them in RNAstructure's own layout (one folder, sets told apart by file prefix), so they
// need no setup. tools/rnastructure_parity compares every table entry against RNAstructure's
// own loader.

#include <RNAstructure/detail/text.hpp>
#include <RNAstructure/rnastructure.hpp>
#include <gtest/gtest.h>
#include <sys/stat.h>

#include <dirent.h>
#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

namespace rs = thermorna::rnastructure;

namespace {

// ThermoRNA's params/RNAstructure (set by CMake).
const std::string kParams = THERMORNA_TEST_PARAMS;

// Folder of a shipped set, by folder name or by RNAstructure's name ("rna", "PSU", "1").
std::string set_dir(const std::string& name) {
    return rs::find_set_folder(kParams, name).directory;
}

std::string read_text(const std::string& path) {
    std::string text;
    EXPECT_TRUE(rs::detail::read_file(path, text)) << path;
    return text;
}

void write_text(const std::string& path, const std::string& text) {
    std::FILE* file = std::fopen(path.c_str(), "wb");
    ASSERT_NE(file, nullptr) << path;
    std::fwrite(text.data(), 1, text.size(), file);
    std::fclose(file);
}

void remove_tree(const std::string& path) {
    if (DIR* d = opendir(path.c_str())) {
        while (const dirent* e = readdir(d)) {
            const std::string name = e->d_name;
            if (name == "." || name == "..") continue;
            const std::string child = path + "/" + name;
            struct stat info {};
            if (lstat(child.c_str(), &info) == 0 && S_ISDIR(info.st_mode)) {
                remove_tree(child);
            } else {
                std::remove(child.c_str());
            }
        }
        closedir(d);
    }
    rmdir(path.c_str());
}

// Temporary folders of this test program, removed when it exits.
struct TempDirs {
    std::vector<std::string> paths;
    TempDirs() = default;
    TempDirs(const TempDirs&) = delete;
    TempDirs& operator=(const TempDirs&) = delete;
    ~TempDirs() {
        for (const std::string& path : paths) remove_tree(path);
    }
};

std::string make_temp_dir() {
    static TempDirs dirs;
    const char* tmp = std::getenv("TMPDIR");
    std::string pattern = std::string(tmp != nullptr && *tmp != '\0' ? tmp : "/tmp") +
                          "/thermorna_rnastructure_XXXXXX";
    const char* dir = mkdtemp(pattern.data());
    EXPECT_NE(dir, nullptr) << pattern;
    if (dir == nullptr) return std::string();
    dirs.paths.emplace_back(dir);
    return dir;
}

std::vector<std::string> files_with_prefix(const std::string& dir, const std::string& prefix) {
    std::vector<std::string> names;
    if (DIR* d = opendir(dir.c_str())) {
        while (dirent* e = readdir(d)) {
            const std::string name = e->d_name;
            if (name == "." || name == "..") continue;
            if (name.compare(0, prefix.size(), prefix) == 0) names.push_back(name);
        }
        closedir(d);
    }
    std::sort(names.begin(), names.end());
    return names;
}

// RNAstructure's layout (its data_tables): every shipped set and pseudconst.dat in one folder.
const std::string& flat_dir() {
    static const std::string dir = [] {
        const std::string flat = make_temp_dir();
        for (const rs::ParameterSetLocation& set : rs::list_set_folders(kParams)) {
            for (const std::string& name : files_with_prefix(set.directory, "")) {
                write_text(flat + "/" + name, read_text(set.directory + "/" + name));
            }
        }
        write_text(flat + "/pseudconst.dat", read_text(kParams + "/pseudoknot/pseudconst.dat"));
        return flat;
    }();
    return dir;
}

// Copies <prefix>.* from @p from to @p to as <new_prefix>.*, passing every line that is not a
// comment through @p edit.
template <typename Edit>
void copy_alphabet(const std::string& from, const std::string& prefix, const std::string& to,
                   const std::string& new_prefix, Edit edit) {
    for (const std::string& name : files_with_prefix(from, prefix + ".")) {
        const std::string text = read_text(from + "/" + name);
        std::string out;
        std::size_t start = 0;
        while (start < text.size()) {
            std::size_t end = text.find('\n', start);
            if (end == std::string::npos) end = text.size();
            std::string line = text.substr(start, end - start);
            const std::size_t first = line.find_first_not_of(" \t\r");
            if (first != std::string::npos && line[first] != '#') line = edit(name, line);
            out += line + "\n";
            start = end + 1;
        }
        write_text(to + "/" + new_prefix + name.substr(prefix.size()), out);
    }
}

std::string replace_all(std::string text, const std::string& from, const std::string& to) {
    for (std::size_t pos = 0; (pos = text.find(from, pos)) != std::string::npos; pos += to.size()) {
        text.replace(pos, from.size(), to);
    }
    return text;
}

// The shipped base set plus shipped modified-base sets, named as in RNAstructure.
rs::Datatable load(const std::vector<std::string>& modifications, double kelvin = 310.15,
                   const std::string& base = "rna") {
    rs::LoadOptions options;
    options.data_directory = set_dir(base);
    options.temperature = kelvin;
    for (const std::string& m : modifications) options.modifications.push_back({set_dir(m)});
    return rs::Datatable::load(options);
}

// The same from RNAstructure's layout, where the file prefix picks the set.
rs::Datatable load_flat(const std::vector<std::string>& modifications, double kelvin = 310.15) {
    rs::LoadOptions options;
    options.data_directory = flat_dir();  // the base set is "rna", efn2's default
    options.temperature = kelvin;
    for (const std::string& m : modifications) options.modifications.push_back({"", m});
    return rs::Datatable::load(options);
}

// Every 4-tuple over the codes of @p count.
template <typename F>
void for_each4(std::size_t count, F f) {
    for (std::size_t a = 0; a < count; ++a)
        for (std::size_t b = 0; b < count; ++b)
            for (std::size_t c = 0; c < count; ++c)
                for (std::size_t d = 0; d < count; ++d)
                    f(rs::Code(a), rs::Code(b), rs::Code(c), rs::Code(d));
}

constexpr rs::Table4D kAllTables[] = {
    rs::Table4D::Stack,     rs::Table4D::HelixEnds, rs::Table4D::Tstkh, rs::Table4D::Tstki,
    rs::Table4D::Tstki23,   rs::Table4D::Tstki1n,   rs::Table4D::Coax,  rs::Table4D::TstackCoax,
    rs::Table4D::CoaxStack, rs::Table4D::Tstack,    rs::Table4D::Tstkm};

}  // namespace

// ============================================================================
// Graphemes
// ============================================================================

TEST(Graphemes, AsciiIsSplitPerByte) {
    const std::vector<std::string_view> g = rs::graphemes::split("ACGU-x");
    ASSERT_EQ(g.size(), 6u);
    EXPECT_EQ(g[0], "A");
    EXPECT_EQ(g[5], "x");
    EXPECT_TRUE(rs::graphemes::is_ascii("ACGU"));
    EXPECT_FALSE(rs::graphemes::is_ascii("ACΨU"));
    EXPECT_EQ(rs::graphemes::first_length("\r\nA"), 2u);
}

TEST(Graphemes, MultiByteSymbols) {
    const std::vector<std::string_view> g =
        rs::graphemes::split("GΨA\xCC\x81\xF0\x9F\x87\xA8\xF0\x9F\x87\xA6U");
    ASSERT_EQ(g.size(), 5u);
    EXPECT_EQ(g[1], "Ψ");
    EXPECT_EQ(g[2], "A\xCC\x81");                         // A + combining acute
    EXPECT_EQ(g[3], "\xF0\x9F\x87\xA8\xF0\x9F\x87\xA6");  // flag (two regional indicators)
    EXPECT_EQ(rs::graphemes::first_length("ΨA"), 2u);
}

// The ASCII shortcut must give the same clusters as segmenting everything with uni-algo.
TEST(Graphemes, MatchesUniAlgoOnRandomText) {
    const std::vector<std::string> pieces = {
        "A",
        "C",
        "G",
        "U",
        "x",
        " ",
        "\t",
        "\r",
        "\n",
        "\r\n",
        "1",
        "#",
        "Ψ",
        "ψ",
        "\xCC\x81",                                      // combining acute accent
        "e\xCC\x88",                                     // e + diaeresis
        "\xF0\x9F\x91\x8D\xF0\x9F\x8F\xBD",              // thumbs up + skin tone
        "\xF0\x9F\x91\xA8\xE2\x80\x8D\xF0\x9F\x91\xA9",  // man ZWJ woman
        "\xE2\x80\x8D",                                  // lone ZWJ
        "\xF0\x9F\x87\xA8",
        "\xF0\x9F\x87\xA6",  // regional indicators
        "\xD8\x80",          // U+0600, a Prepend character
        "\xE1\x84\x80",
        "\xE1\x85\xA1",
        "\xE1\x86\xA8",              // Hangul L, V, T jamo
        "\xEA\xB0\x80",              // Hangul syllable
        "\xEF\xB8\x8F\xE2\x83\xA3",  // VS16 + combining keycap
        "\xCE",
        "\xFF",
        "\x80",                      // ill-formed UTF-8
        "\xE0\xA4\x95\xE0\xA5\x8D",  // Devanagari ka + virama
    };
    std::mt19937 rng(12345);
    std::uniform_int_distribution<std::size_t> piece(0, pieces.size() - 1);
    std::uniform_int_distribution<int> length(0, 24);
    for (int trial = 0; trial < 20000; ++trial) {
        std::string text;
        for (int k = length(rng); k > 0; --k) text += pieces[piece(rng)];
        std::vector<std::string_view> fast;
        std::vector<std::string_view> reference;
        rs::graphemes::split(text, fast);
        rs::graphemes::split_with_uni_algo(text, reference);
        ASSERT_EQ(fast, reference) << "text: " << text;
        if (!text.empty()) ASSERT_EQ(rs::graphemes::first_length(text), reference.front().size());
    }
}

// ============================================================================
// Number parsing (atof / floor_entry_to_short)
// ============================================================================

TEST(Text, ParseDoubleMatchesStrtod) {
    std::vector<std::string> tokens = {".",
                                       "-",
                                       "+",
                                       "+.",
                                       "-.",
                                       "1e",
                                       "1e+",
                                       "1e-3",
                                       "2E2",
                                       "0x10",
                                       "inf",
                                       "-nan",
                                       "abc",
                                       "1.5abc",
                                       ".5",
                                       "5.",
                                       "-0",
                                       "-.0",
                                       "+0.7",
                                       "000012.3400",
                                       "123456789012345678901234",
                                       "0.1000000000000000055511151231257827",
                                       "1e400",
                                       "1.25"};
    std::mt19937 rng(7);
    std::uniform_int_distribution<int> digit(0, 9);
    std::uniform_int_distribution<int> count(0, 6);
    std::uniform_int_distribution<int> coin(0, 3);
    for (int trial = 0; trial < 50000; ++trial) {
        std::string token;
        if (coin(rng) == 0) token += coin(rng) % 2 ? "-" : "+";
        for (int k = count(rng); k > 0; --k) token += char('0' + digit(rng));
        if (coin(rng) != 0) token += '.';
        for (int k = count(rng); k > 0; --k) token += char('0' + digit(rng));
        if (coin(rng) == 0) token += "e" + std::to_string(digit(rng) - 5);
        tokens.push_back(token);
    }
    for (const std::string& token : tokens) {
        const double expected = std::strtod(token.c_str(), nullptr);
        const double got = rs::detail::parse_double(token);
        if (std::isnan(expected)) {
            EXPECT_TRUE(std::isnan(got)) << token;
        } else {
            EXPECT_EQ(std::memcmp(&expected, &got, sizeof got), 0) << token << " " << expected;
        }
    }
}

TEST(Text, EntryToEnergyMatchesRNAstructure) {
    // RNAstructure: (short)floor(atof(x) * conversionfactor + .5), "." is INFINITE_ENERGY.
    auto reference = [](const std::string& token) -> rs::Energy {
        if (token == ".") return rs::kInfiniteEnergy;
        return static_cast<rs::Energy>(std::floor(std::atof(token.c_str()) * 10 + .5));
    };
    std::vector<std::string> tokens = {".",    "-.0",   "0.45",    "0.15",  "-1.25",
                                       "3.42", "-0.05", "1.079",   "+0.7",  "2.",
                                       "-",    "",      "12345.6", "9999.9"};
    std::mt19937 rng(11);
    std::uniform_int_distribution<int> value(-50000, 50000);
    std::uniform_int_distribution<int> decimals(0, 3);
    for (int trial = 0; trial < 50000; ++trial) {
        const int places = decimals(rng);
        char buffer[64];
        std::snprintf(buffer, sizeof buffer, "%.*f", places, value(rng) / std::pow(10.0, places));
        tokens.emplace_back(buffer);
    }
    for (const std::string& token : tokens) {
        if (std::fabs(std::atof(token.c_str())) > 3000) continue;  // outside short in tenths
        EXPECT_EQ(rs::detail::entry_to_energy(token), reference(token)) << token;
    }
}

TEST(Text, DataLinesFollowReadDatatableLines) {
    const auto lines = rs::detail::data_lines("\xEF\xBB\xBF# c\n  \n\t A  X \r\n#x\nlast");
    ASSERT_EQ(lines.size(), 2u);
    EXPECT_EQ(lines[0].text, "A  X");
    EXPECT_EQ(lines[0].number, 3u);
    EXPECT_EQ(lines[1].text, "last");
}

// ============================================================================
// Alphabet
// ============================================================================

TEST(Alphabet, ParsesSpecificationWithGraphemes) {
    const rs::Alphabet a = rs::Alphabet::parse(
        "# comment\nBases\nX = N = x\nA = a\nC\nG\nU = T\nΨ = ψ = P\nI\n"
        "Pairing\nA U\nG C\nG U\nA Ψ\nSingle\na\nψ\nNon-interacting\nX\nLinker\nI\n",
        "test", "test.specification.dat");
    ASSERT_EQ(a.size(), 7u);
    EXPECT_EQ(a.code("N"), 0);
    EXPECT_EQ(a.code("T"), 4);
    EXPECT_EQ(a.code("Ψ"), 5);
    EXPECT_EQ(a.code("ψ"), 5);
    EXPECT_EQ(a.code("P"), 5);
    EXPECT_EQ(a.code("Q"), rs::kInvalidCode);
    EXPECT_TRUE(a.can_pair(a.code("Ψ"), a.code("A")));
    EXPECT_FALSE(a.can_pair(a.code("Ψ"), a.code("G")));
    EXPECT_TRUE(a.is_single_symbol("ψ"));
    EXPECT_TRUE(a.is_non_interacting(0));
    EXPECT_TRUE(a.is_linker(6));
    EXPECT_FALSE(a.can_interact(6));
    EXPECT_EQ(a.primary_symbol(5), "Ψ");
}

TEST(Alphabet, RejectsUnknownPairingSymbols) {
    EXPECT_THROW((void) rs::Alphabet::parse("Bases\nA\nU\nPairing\nA Q\n", "bad", "bad.spec"),
                 rs::ParameterError);
}

// ============================================================================
// Pseudoknot constants
// ============================================================================

TEST(Pseudoknot, ParsesConstantsAfterArrows) {
    std::string text = "header words\n";
    for (int k = 0; k < 16; ++k)
        text += "c" + std::to_string(k) + " is --> " + std::to_string(k * 1.5) + "\n";
    const rs::PseudoknotConstants c = rs::PseudoknotConstants::parse(text, "pseudconst.dat");
    EXPECT_DOUBLE_EQ(c.values[0], 0.0);
    EXPECT_DOUBLE_EQ(c.values[15], 22.5);
    EXPECT_DOUBLE_EQ(c.p1, 0.35);
    EXPECT_THROW((void) rs::PseudoknotConstants::parse("--> 1 --> 2", "short"), rs::ParameterError);
}

TEST(Pseudoknot, ReadsShippedFile) {
    const rs::PseudoknotConstants c = rs::PseudoknotConstants::read(kParams + "/pseudoknot");
    EXPECT_DOUBLE_EQ(c.values[0], 6.5);
    EXPECT_DOUBLE_EQ(c.values[1], 15);
    EXPECT_DOUBLE_EQ(c.values[15], 30000);
}

// ============================================================================
// AlphabetTables (one alphabet, as efn2 loads it)
// ============================================================================

TEST(AlphabetTables, RnaAt37) {
    const rs::AlphabetTables t = rs::AlphabetTables::load(set_dir("rna"));
    ASSERT_EQ(t.size(), 6u);
    const rs::MiscLoop& m = t.misc();
    EXPECT_EQ(m.prelog, static_cast<float>(1.079 * 10));
    EXPECT_EQ(m.maxpen, 30);
    EXPECT_EQ(m.poppen[1], 6);
    EXPECT_EQ(m.eparam[5], 93);
    EXPECT_EQ(m.eparam[6], 0);
    EXPECT_EQ(m.eparam[7], 30);
    EXPECT_EQ(m.eparam[9], -500);
    EXPECT_EQ(m.eparam[10], -6);
    EXPECT_EQ(m.efn2a, 93);
    EXPECT_EQ(m.efn2c, -6);
    EXPECT_EQ(m.mlasym, 9);
    EXPECT_EQ(m.strain, 31);
    EXPECT_EQ(m.gubonus, -22);
    EXPECT_EQ(m.cslope, 3);
    EXPECT_EQ(m.cint, 16);
    EXPECT_EQ(m.c3, 15);
    EXPECT_EQ(m.init, 41);
    EXPECT_EQ(m.singlecbulge, -9);
    ASSERT_EQ(t.inter().size(), 31u);
    EXPECT_EQ(t.inter()[0], rs::kInfiniteEnergy);
    EXPECT_EQ(t.inter()[4], 11);
    EXPECT_EQ(t.bulge()[1], 38);
    EXPECT_EQ(t.hairpin()[3], 54);

    // X=0 A=1 C=2 G=3 U=4 I=5
    EXPECT_EQ(t.stack(1, 4, 1, 4), -9);
    EXPECT_EQ(t.stack(1, 4, 1, 1), rs::kInfiniteEnergy);
    EXPECT_EQ(t.helix_ends(1, 4, 1, 1), 5);
    EXPECT_EQ(t.helix_ends(2, 3, 1, 1), 0);
    EXPECT_EQ(t.dangle(1, 4, 1, 1), -8);                    // 3' dangling A on A-U
    EXPECT_EQ(t.dangle(1, 4, 0, 1), 0);                     // processDat: X dangles are free
    EXPECT_EQ(t.tstack(1, 4, 5, 1), t.dangle(1, 4, 1, 2));  // linker: a mismatch is a dangle

    const rs::Code cuacgg[] = {2, 4, 1, 2, 3, 3};
    ASSERT_TRUE(t.find_tloop(cuacgg).has_value());
    EXPECT_EQ(*t.find_tloop(cuacgg), 28);
    const rs::Code aaaaaa[] = {1, 1, 1, 1, 1, 1};
    EXPECT_FALSE(t.find_tloop(aaaaaa).has_value());
    EXPECT_TRUE(t.warnings().empty());
}

// Reference values from RNAstructure's own datatable after ScaleToTemperature(300).
TEST(AlphabetTables, RnaScaledTo300K) {
    const rs::AlphabetTables t = rs::AlphabetTables::load(set_dir("rna"), "", 300.0);
    std::uint32_t bits = 0;
    const float rt = t.RT();
    std::memcpy(&bits, &rt, sizeof bits);
    EXPECT_EQ(bits, 0x3f189df0u);
    const float prelog = t.misc().prelog;
    std::memcpy(&bits, &prelog, sizeof bits);
    EXPECT_EQ(bits, 0x4126fd7cu);
    EXPECT_EQ(t.misc().maxpen, 29);
    EXPECT_EQ(t.misc().efn2a, 100);
    EXPECT_EQ(t.misc().efn2c, -7);
    EXPECT_EQ(t.misc().strain, 39);
    EXPECT_EQ(t.misc().gubonus, -26);
    EXPECT_EQ(t.misc().poppen[1], 7);
    EXPECT_EQ(t.misc().eparam[7], 30);  // not rescaled
    EXPECT_EQ(t.inter()[5], 17);
    EXPECT_EQ(t.stack(1, 4, 1, 4), -11);
    EXPECT_EQ(t.stack(3, 2, 3, 2), -36);
    EXPECT_EQ(t.helix_ends(1, 4, 1, 1), 6);
    EXPECT_EQ(t.iloop22(1, 1, 4, 4, 1, 1, 1, 1), 20);
    EXPECT_EQ(t.tloop().front().energy, 24);
    EXPECT_DOUBLE_EQ(t.temperature(), 300.0);
}

TEST(AlphabetTables, MissingEnthalpyTablesFailLikeEfn2) {
    EXPECT_THROW((void) rs::AlphabetTables::load(set_dir("m6A"), "", 300.0), rs::ParameterError);
    const rs::AlphabetTables kept =
        rs::AlphabetTables::load(set_dir("m6A"), "", 300.0, rs::MissingEnthalpy::KeepFreeEnergy);
    EXPECT_DOUBLE_EQ(kept.temperature(), rs::kTemperature37C);
    ASSERT_FALSE(kept.warnings().empty());
    EXPECT_NE(kept.warnings().back().find("m6A.loop.dh"), std::string::npos);
    // 0.005 K away from 37 degrees C is within RNAstructure's tolerance: nothing to rescale.
    EXPECT_NO_THROW((void) rs::AlphabetTables::load(set_dir("m6A"), "", 310.155));
}

// PSU.coaxstack.dg has rows with a value missing; RNAstructure reads those entries as 0.
TEST(AlphabetTables, ShortRowsAreReadAsZeroLikeRNAstructure) {
    const rs::AlphabetTables t = rs::AlphabetTables::load(set_dir("PSU"));
    const rs::Code A = t.alphabet().code("A");
    const rs::Code C = t.alphabet().code("C");
    const rs::Code P = t.alphabet().code("P");
    EXPECT_EQ(t.coaxstack(A, A, C, P), 0);
    bool reported = false;
    for (const std::string& w : t.warnings())
        reported |= w.find("PSU.coaxstack.dg") != std::string::npos;
    EXPECT_TRUE(reported);
}

// dna.miscloop.dh still has the pre-6.5 terminal AU lines (RNAstructure 6.6 misreads it).
TEST(AlphabetTables, ReadsPre65Miscloop) {
    const std::string dir = make_temp_dir();
    copy_alphabet(set_dir("dna"), "dna", dir, "dna",
                  [](const std::string&, const std::string& l) { return l; });
    write_text(dir + "/dna.helix_ends.dh", read_text(set_dir("rna") + "/rna.helix_ends.dh"));
    const rs::AlphabetTables dh = rs::AlphabetTables::read(dir, "dna", rs::TableKind::Enthalpy);
    EXPECT_EQ(dh.misc().gubonus, 0);
    EXPECT_EQ(dh.misc().init, -72);
    bool reported = false;
    for (const std::string& w : dh.warnings()) reported |= w.find("pre-6.5") != std::string::npos;
    EXPECT_TRUE(reported);
    EXPECT_NO_THROW((void) rs::AlphabetTables::load(dir, "dna", 300.0));
}

TEST(AlphabetTables, MissingFilesAreReported) {
    try {
        (void) rs::AlphabetTables::read("/nonexistent", "rna", rs::TableKind::FreeEnergy);
        FAIL();
    } catch (const rs::ParameterError& e) {
        EXPECT_NE(std::string(e.what()).find("'/nonexistent' does not exist"), std::string::npos);
    }
    const std::string dir = make_temp_dir();
    write_text(dir + "/rna.stack.dg", "");
    try {
        (void) rs::AlphabetTables::read(dir, "rna", rs::TableKind::FreeEnergy);
        FAIL();
    } catch (const rs::ParameterError& e) {
        EXPECT_NE(std::string(e.what()).find("rna.specification.dat"), std::string::npos);
    }
}

// ============================================================================
// Datatable (base alphabet + modifications)
// ============================================================================

TEST(Datatable, WithoutModificationsEqualsEfn2Tables) {
    for (double kelvin : {310.15, 300.0}) {
        const rs::Datatable d = load({}, kelvin);
        const rs::AlphabetTables t = rs::AlphabetTables::load(set_dir("rna"), "", kelvin);
        ASSERT_EQ(d.size(), t.size());
        for (rs::Table4D table : kAllTables) {
            for_each4(d.size(), [&](rs::Code a, rs::Code b, rs::Code c, rs::Code e) {
                ASSERT_EQ(d.table4(table, a, b, c, e), t.table4(table, a, b, c, e));
            });
        }
        for_each4(d.size(), [&](rs::Code a, rs::Code b, rs::Code c, rs::Code e) {
            for (rs::Code f = 0; f < d.size(); ++f)
                for (rs::Code g = 0; g < d.size(); ++g) {
                    ASSERT_EQ(d.iloop11(a, b, c, e, f, g), t.iloop11(a, b, c, e, f, g));
                    ASSERT_EQ(d.iloop22(a, b, c, e, f, g, a, b), t.iloop22(a, b, c, e, f, g, a, b));
                }
        });
        EXPECT_EQ(d.misc().efn2a, t.misc().efn2a);
        EXPECT_EQ(d.RT(), t.RT());
    }
}

TEST(Datatable, MergesSeveralModifications) {
    const rs::Datatable d = load({"m6A", "PSU", "1"});
    ASSERT_EQ(d.base_size(), 6u);
    ASSERT_EQ(d.size(), 9u);

    const rs::Code m6a = d.code("6");
    const rs::Code psi = d.code("P");
    const rs::Code m1psi = d.code("1");
    EXPECT_EQ(m6a, 6);
    EXPECT_EQ(psi, 7);
    EXPECT_EQ(m1psi, 8);
    EXPECT_EQ(d.code("p"), psi);
    EXPECT_EQ(d.origin(m6a), "m6A");
    EXPECT_EQ(d.origin(psi), "Ψ");  // the folder's name (its files are PSU.*)
    EXPECT_TRUE(d.is_modified(m6a));
    EXPECT_FALSE(d.is_modified(d.code("A")));

    // Fallbacks inferred from the Pairing sections.
    EXPECT_EQ(d.fallback(m6a), d.code("A"));
    EXPECT_EQ(d.fallback(psi), d.code("U"));
    EXPECT_EQ(d.fallback(m1psi), d.code("U"));

    // "M" and "m" mean m6A in m6A.specification.dat and m1Ψ in 1.specification.dat.
    EXPECT_TRUE(d.is_ambiguous("M"));
    EXPECT_TRUE(d.is_ambiguous("m"));
    EXPECT_EQ(d.code("M"), rs::kInvalidCode);
    EXPECT_THROW((void) d.encode("GGMCC"), rs::ParameterError);
    EXPECT_THROW((void) d.encode("GGZCC"), rs::ParameterError);

    EXPECT_TRUE(d.can_pair(m6a, d.code("U")));
    EXPECT_FALSE(d.can_pair(m6a, d.code("G")));
    EXPECT_TRUE(d.can_pair(psi, d.code("G")));
    EXPECT_TRUE(d.can_pair(m6a, psi));  // no alphabet has both: decided as A-U
    EXPECT_TRUE(d.is_single_symbol("p"));
    EXPECT_TRUE(d.has_symbol(d.code("U"), "T"));

    const std::vector<rs::Code> codes = d.encode("GG6AAPUCC1");
    const std::vector<rs::Code> expected = {3, 3, m6a, 1, 1, psi, 4, 2, 2, m1psi};
    EXPECT_EQ(codes, expected);
}

// Contexts covered by one modification give exactly that alphabet's values (efn2 -a <it>).
TEST(Datatable, ModificationContextsUseTheirAlphabet) {
    const rs::Datatable d = load({"m6A", "PSU", "1"});
    for (std::size_t o = 0; o < d.modification_count(); ++o) {
        const rs::AlphabetTables& own = d.modification(o);
        const rs::Alphabet& local = own.alphabet();
        std::vector<rs::Code> to_merged(local.size());
        for (std::size_t l = 0; l < local.size(); ++l) {
            to_merged[l] = d.code(local.primary_symbol(rs::Code(l)));
            if (to_merged[l] == rs::kInvalidCode)
                to_merged[l] = d.code(local.symbols(rs::Code(l)).back());
            ASSERT_NE(to_merged[l], rs::kInvalidCode) << local.primary_symbol(rs::Code(l));
        }
        auto modified = [&](std::initializer_list<rs::Code> ls) {
            for (rs::Code l : ls)
                if (d.is_modified(to_merged[l])) return true;
            return false;
        };
        std::size_t checked = 0;
        for (rs::Table4D table : kAllTables) {
            for_each4(local.size(), [&](rs::Code a, rs::Code b, rs::Code c, rs::Code e) {
                if (!modified({a, b, c, e})) return;
                ASSERT_EQ(d.table4(table, to_merged[a], to_merged[b], to_merged[c], to_merged[e]),
                          own.table4(table, a, b, c, e))
                    << local.name() << " table " << int(table);
                ++checked;
            });
        }
        for_each4(local.size(), [&](rs::Code a, rs::Code b, rs::Code c, rs::Code e) {
            for (rs::Code dd = 1; dd <= 2; ++dd) {
                if (!modified({a, b, c})) continue;
                ASSERT_EQ(d.dangle(to_merged[a], to_merged[b], to_merged[c], dd),
                          own.dangle(a, b, c, dd));
            }
            for (rs::Code f = 0; f < local.size(); ++f)
                for (rs::Code g = 0; g < local.size(); ++g) {
                    if (!modified({a, b, c, e, f, g})) continue;
                    ASSERT_EQ(d.iloop11(to_merged[a], to_merged[b], to_merged[c], to_merged[e],
                                        to_merged[f], to_merged[g]),
                              own.iloop11(a, b, c, e, f, g));
                }
        });
        std::mt19937 rng(static_cast<unsigned>(o));
        std::uniform_int_distribution<int> pick(0, static_cast<int>(local.size()) - 1);
        for (int trial = 0; trial < 200000; ++trial) {
            rs::Code l[8];
            for (rs::Code& x : l) x = rs::Code(pick(rng));
            rs::Code m[8];
            for (int k = 0; k < 8; ++k) m[k] = to_merged[l[k]];
            if (modified({l[0], l[1], l[2], l[3], l[4], l[5], l[6]})) {
                ASSERT_EQ(d.iloop21(m[0], m[1], m[2], m[3], m[4], m[5], m[6]),
                          own.iloop21(l[0], l[1], l[2], l[3], l[4], l[5], l[6]));
            }
            if (!modified({l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7]})) continue;
            ASSERT_EQ(d.iloop22(m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7]),
                      own.iloop22(l[0], l[1], l[2], l[3], l[4], l[5], l[6], l[7]));
            const rs::LookupSource source = d.source(m, 8);
            ASSERT_EQ(source.kind, rs::LookupSource::Kind::Modification);
            ASSERT_EQ(source.modification, o);
        }
        EXPECT_GT(checked, 0u);
    }
}

// Canonical-only contexts come from the base alphabet, even though m6A.stack.dg lists other
// values for some canonical GU stacks.
TEST(Datatable, CanonicalContextsUseTheBaseAlphabet) {
    const rs::Datatable d = load({"m6A"});
    const rs::AlphabetTables rna = rs::AlphabetTables::load(set_dir("rna"));
    const rs::AlphabetTables m6a = rs::AlphabetTables::load(set_dir("m6A"));
    const rs::Code A = 1, U = 4, G = 3;
    EXPECT_NE(rna.stack(A, U, G, U), m6a.stack(A, U, G, U));
    EXPECT_EQ(d.stack(A, U, G, U), rna.stack(A, U, G, U));
    const rs::Code ctx[] = {A, U, G, U};
    EXPECT_EQ(d.source(ctx, 4).kind, rs::LookupSource::Kind::Base);
}

// Modified bases from different alphabets in one context fall back to their unmodified bases.
TEST(Datatable, MixedModificationsFallBack) {
    const rs::Datatable d = load({"m6A", "PSU"});
    const rs::Code m6a = d.code("6"), psi = d.code("P");
    const rs::Code A = 1, C = 2, G = 3, U = 4;
    for (rs::Table4D table : kAllTables) {
        EXPECT_EQ(d.table4(table, m6a, psi, C, G), d.table4(table, A, U, C, G));
        EXPECT_EQ(d.table4(table, psi, m6a, G, A), d.table4(table, U, A, G, A));
    }
    EXPECT_EQ(d.iloop11(m6a, psi, G, U, A, C), d.iloop11(A, U, G, U, A, C));
    const rs::Code ctx[] = {m6a, psi, C, G};
    EXPECT_EQ(d.source(ctx, 4).kind, rs::LookupSource::Kind::Fallback);

    // A special hairpin mixing both falls back too; one with Ψ only uses PSU's table.
    const rs::Code mixed[] = {C, psi, m6a, C, G, G};
    const rs::Code plain[] = {C, U, A, C, G, G};
    EXPECT_EQ(d.tloop(mixed), d.tloop(plain));
    const rs::Code cpucgg[] = {C, psi, U, C, G, G};  // PSU.tloop.dg: CPUCGG -0.2
    ASSERT_TRUE(d.tloop(cpucgg).has_value());
    EXPECT_EQ(*d.tloop(cpucgg), -2);
}

// A later modification that reuses a base's primary symbol overrides earlier parameters.
TEST(Datatable, LaterModificationsWin) {
    const std::string dir = make_temp_dir();
    copy_alphabet(set_dir("m6A"), "m6A", dir, "m6Av2",
                  [](const std::string& file, const std::string& l) {
                      // M-U stacked on U-M becomes -9.9 kcal/mol.
                      return file == "m6A.stack.dg" && l.rfind("U\t-0.7", 0) == 0
                                 ? std::string("U\t-0.7\t.\t-0.3\t.\t-9.9")
                                 : l;
                  });
    rs::LoadOptions options;
    options.data_directory = set_dir("rna");
    options.modifications = {{set_dir("m6A")}, {dir, "m6Av2"}};
    const rs::Datatable d = rs::Datatable::load(options);
    ASSERT_EQ(d.size(), 7u);  // m6Av2's M is the same base as m6A's M
    const rs::Code M = d.code("6"), U = 4;
    EXPECT_EQ(d.stack(U, M, U, M), -99);
    const rs::Code ctx[] = {U, M, U, M};
    EXPECT_EQ(d.source(ctx, 4).modification, 1u);
}

// A modification written with multi-byte symbols behaves exactly like the ASCII one.
TEST(Datatable, GraphemeSymbolsInParameterFiles) {
    const std::string dir = make_temp_dir();
    copy_alphabet(set_dir("PSU"), "PSU", dir, "PSUg",
                  [](const std::string& file, const std::string& l) {
                      if (file == "PSU.specification.dat") {
                          if (l == "P = p") return std::string("Ψ = ψ");
                          if (l == "A P") return std::string("A Ψ");
                          if (l == "G P") return std::string("G Ψ");
                          if (l == "p") return std::string("ψ");
                          return l;
                      }
                      return replace_all(replace_all(l, "P", "Ψ"), "p", "ψ");
                  });
    const rs::Datatable ascii = load({"PSU"});
    rs::LoadOptions options;
    options.data_directory = set_dir("rna");
    options.modifications = {{dir, "PSUg"}};
    const rs::Datatable unicode = rs::Datatable::load(options);
    ASSERT_EQ(unicode.size(), ascii.size());
    EXPECT_EQ(unicode.code("Ψ"), ascii.code("P"));
    EXPECT_EQ(unicode.code("P"), rs::kInvalidCode);
    for (rs::Table4D table : kAllTables) {
        for_each4(ascii.size(), [&](rs::Code a, rs::Code b, rs::Code c, rs::Code e) {
            ASSERT_EQ(unicode.table4(table, a, b, c, e), ascii.table4(table, a, b, c, e));
        });
    }
    for_each4(ascii.size(), [&](rs::Code a, rs::Code b, rs::Code c, rs::Code e) {
        for (rs::Code f = 0; f < ascii.size(); ++f) {
            ASSERT_EQ(unicode.iloop22(a, b, c, e, f, a, b, c),
                      ascii.iloop22(a, b, c, e, f, a, b, c));
            ASSERT_EQ(unicode.iloop21(a, b, c, e, f, a, b), ascii.iloop21(a, b, c, e, f, a, b));
        }
    });
    EXPECT_EQ(unicode.encode("GGΨΨCC"), ascii.encode("GGPPCC"));
    EXPECT_TRUE(unicode.is_single_symbol("ψ"));
}

TEST(Datatable, AliasesAndExplicitFallbacks) {
    rs::LoadOptions options;
    options.data_directory = set_dir("rna");
    options.modifications = {{set_dir("PSU"), "", {}, {{"Ψ", "P"}}},
                             {set_dir("m6A"), "", {{"6", "G"}}}};
    const rs::Datatable d = rs::Datatable::load(options);
    EXPECT_EQ(d.code("Ψ"), d.code("P"));
    EXPECT_EQ(d.fallback(d.code("6")), d.code("G"));
    EXPECT_EQ(d.encode("AΨ6U"), d.encode("AP6U"));

    // Knotergy's ProcessedRNAEntry hands over grapheme views.
    const std::string raw = "GGΨA6ΨCC";
    const std::vector<std::string_view> views = rs::graphemes::split(raw);
    EXPECT_EQ(d.encode(views), d.encode(raw));

    rs::LoadOptions clash = options;
    clash.modifications = {{set_dir("PSU"), "", {}, {{"A", "P"}}}};
    EXPECT_THROW((void) rs::Datatable::load(clash), rs::ParameterError);
}

TEST(Datatable, TemperatureAppliesToEveryAlphabet) {
    EXPECT_THROW((void) load({"m6A"}, 300.0), rs::ParameterError);
    rs::LoadOptions options;
    options.data_directory = set_dir("rna");
    options.temperature = 300.0;
    options.missing_enthalpy = rs::MissingEnthalpy::KeepFreeEnergy;
    options.modifications = {{set_dir("m6A")}};
    const rs::Datatable d = rs::Datatable::load(options);
    EXPECT_DOUBLE_EQ(d.temperature(), 300.0);
    EXPECT_DOUBLE_EQ(d.modification(0).temperature(), rs::kTemperature37C);
    EXPECT_EQ(d.stack(1, 4, 1, 4), -11);  // base contexts are rescaled
}

// ============================================================================
// One folder per parameter set (params/RNAstructure/<set>/)
// ============================================================================

namespace {

// Copies the shipped set's <prefix>.* files into @p folder, renamed to <new_prefix>.* (or to
// unprefixed names when @p new_prefix is empty).
void copy_set(const std::string& prefix, const std::string& folder, const std::string& new_prefix) {
    const std::string from = set_dir(prefix);
    for (const std::string& name : files_with_prefix(from, prefix + ".")) {
        const std::string rest = name.substr(prefix.size() + 1);
        write_text(folder + "/" + (new_prefix.empty() ? rest : new_prefix + "." + rest),
                   read_text(from + "/" + name));
    }
}

// A params/RNAstructure-like tree: turner2004/rna.*, m6A/m6A.*, m1Ψ/1.*, and Ψ/ with PSU's
// files under unprefixed names.
std::string make_set_tree() {
    const std::string root = make_temp_dir();
    for (const char* folder : {"/turner2004", "/m6A", "/m1Ψ", "/Ψ"}) {
        EXPECT_EQ(mkdir((root + folder).c_str(), 0755), 0) << folder;
    }
    copy_set("rna", root + "/turner2004", "rna");
    copy_set("m6A", root + "/m6A", "m6A");
    copy_set("1", root + "/m1Ψ", "1");
    copy_set("PSU", root + "/Ψ", "");
    return root;
}

template <typename F>
void expect_message(F f, const std::string& part) {
    try {
        f();
        ADD_FAILURE() << "no ParameterError (expected one mentioning \"" << part << "\")";
    } catch (const rs::ParameterError& e) {
        EXPECT_NE(std::string(e.what()).find(part), std::string::npos) << e.what();
    }
}

}  // namespace

TEST(ParameterSets, FolderPerSetLoadsLikeDataTables) {
    const std::string root = make_set_tree();
    rs::LoadOptions options;
    options.data_directory = root + "/turner2004";
    options.modifications = {{root + "/m6A"}, {root + "/Ψ/"}, {root + "/m1Ψ"}};
    const rs::Datatable folders = rs::Datatable::load(options);
    const rs::Datatable flat = load_flat({"m6A", "PSU", "1"});
    EXPECT_EQ(flat.base().name(), "rna");
    EXPECT_EQ(flat.modification_name(1), "PSU");

    // Names come from the folders; file prefixes are found in them.
    EXPECT_EQ(folders.base().name(), "turner2004");
    EXPECT_EQ(folders.base().prefix(), "rna");
    EXPECT_EQ(folders.modification_name(0), "m6A");
    EXPECT_EQ(folders.modification_name(1), "Ψ");
    EXPECT_EQ(folders.modification(1).prefix(), "");
    EXPECT_EQ(folders.modification_name(2), "m1Ψ");
    EXPECT_EQ(folders.modification(2).prefix(), "1");
    EXPECT_EQ(folders.origin(folders.code("1")), "m1Ψ");
    EXPECT_EQ(folders.origin(folders.code("A")), "turner2004");

    ASSERT_EQ(folders.size(), flat.size());
    for (rs::Code c = 0; c < flat.size(); ++c) EXPECT_EQ(folders.symbols(c), flat.symbols(c));
    for (rs::Table4D table : kAllTables) {
        for_each4(flat.size(), [&](rs::Code a, rs::Code b, rs::Code c, rs::Code e) {
            ASSERT_EQ(folders.table4(table, a, b, c, e), flat.table4(table, a, b, c, e));
        });
    }
    std::mt19937 rng(7);
    std::uniform_int_distribution<int> pick(0, static_cast<int>(flat.size()) - 1);
    for (int trial = 0; trial < 100000; ++trial) {
        rs::Code k[8];
        for (rs::Code& x : k) x = rs::Code(pick(rng));
        ASSERT_EQ(folders.dangle(k[0], k[1], k[2], 1 + std::size_t(k[3] % 2)),
                  flat.dangle(k[0], k[1], k[2], 1 + std::size_t(k[3] % 2)));
        ASSERT_EQ(folders.iloop11(k[0], k[1], k[2], k[3], k[4], k[5]),
                  flat.iloop11(k[0], k[1], k[2], k[3], k[4], k[5]));
        ASSERT_EQ(folders.iloop21(k[0], k[1], k[2], k[3], k[4], k[5], k[6]),
                  flat.iloop21(k[0], k[1], k[2], k[3], k[4], k[5], k[6]));
        ASSERT_EQ(folders.iloop22(k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]),
                  flat.iloop22(k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]));
        ASSERT_EQ(folders.tloop(k), flat.tloop(k));
        ASSERT_EQ(folders.hexaloop(k), flat.hexaloop(k));
    }
    EXPECT_EQ(folders.misc().efn2a, flat.misc().efn2a);
    EXPECT_EQ(folders.hairpin(), flat.hairpin());
    EXPECT_EQ(folders.encode("GG6APUCC1"), flat.encode("GG6APUCC1"));
}

TEST(ParameterSets, FolderPerSetAtOtherTemperatures) {
    const std::string root = make_set_tree();
    rs::LoadOptions options;
    options.data_directory = root + "/turner2004";
    options.temperature = 300.0;
    const rs::Datatable folders = rs::Datatable::load(options);
    const rs::Datatable flat = load_flat({}, 300.0);
    for (rs::Table4D table : kAllTables) {
        for_each4(flat.size(), [&](rs::Code a, rs::Code b, rs::Code c, rs::Code e) {
            ASSERT_EQ(folders.table4(table, a, b, c, e), flat.table4(table, a, b, c, e));
        });
    }
    EXPECT_EQ(folders.misc().efn2a, flat.misc().efn2a);
    EXPECT_EQ(folders.RT(), flat.RT());

    // m1Ψ has no enthalpy tables, so it cannot be rescaled (like efn2 -a 1 -t 300).
    options.modifications = {{root + "/m1Ψ"}};
    expect_message([&] { (void) rs::Datatable::load(options); }, "'m1Ψ' cannot be rescaled");
}

TEST(ParameterSets, FindsSetsAndReportsMistakes) {
    const std::vector<std::string> sets = rs::list_parameter_sets(flat_dir());
    for (const char* set : {"rna", "dna", "m6A", "PSU", "1"}) {
        EXPECT_TRUE(std::find(sets.begin(), sets.end(), set) != sets.end()) << set;
    }
    // RNAstructure's shared folder: the prefix names the set; efn2's default is "rna".
    expect_message([] { (void) rs::locate_parameter_set(flat_dir()); }, "name the one to use");
    expect_message([] { (void) rs::locate_parameter_set(flat_dir(), "rna2"); },
                   "rna2.specification.dat");
    rs::LoadOptions options;
    options.data_directory = flat_dir();
    EXPECT_EQ(rs::Datatable::load(options).base().name(), "rna");

    const std::string root = make_set_tree();
    const rs::ParameterSetLocation m1psi = rs::locate_parameter_set(root + "/m1Ψ");
    EXPECT_EQ(m1psi.prefix, "1");
    EXPECT_EQ(m1psi.name, "m1Ψ");
    EXPECT_EQ(m1psi.file("stack", ".dg"), root + "/m1Ψ/1.stack.dg");
    EXPECT_EQ(rs::locate_parameter_set(root + "/Ψ").file("stack", ".dg"), root + "/Ψ/stack.dg");
    expect_message([&] { (void) rs::locate_parameter_set(root); }, "No RNAstructure parameter set");

    // Sets of a params/RNAstructure-like folder, found by folder name or by file prefix.
    const std::vector<rs::ParameterSetLocation> folders = rs::list_set_folders(root);
    ASSERT_EQ(folders.size(), 4u);
    EXPECT_EQ(folders[0].name, "m1Ψ");
    EXPECT_EQ(folders[1].name, "m6A");
    EXPECT_EQ(folders[2].name, "turner2004");
    EXPECT_EQ(folders[3].name, "Ψ");
    EXPECT_EQ(rs::find_set_folder(root, "m1Ψ").directory, root + "/m1Ψ");
    EXPECT_EQ(rs::find_set_folder(root, "rna").name, "turner2004");
    EXPECT_EQ(rs::find_set_folder(root, "1").name, "m1Ψ");
    EXPECT_EQ(rs::find_set_folder(root + "/", "m6a").directory, root + "/m6A");
    EXPECT_EQ(rs::find_set_folder(root, "RNA").name, "turner2004");
    expect_message([&] { (void) rs::find_set_folder(root, "PSU"); }, "'turner2004' (rna)");
    expect_message([&] { (void) rs::find_set_folder(root, ""); }, "No parameter set name");
    // Folder names win over file prefixes, exact matches over matches ignoring case.
    EXPECT_EQ(mkdir((root + "/m6A_v2").c_str(), 0755), 0);
    copy_set("m6A", root + "/m6A_v2", "m6A");
    EXPECT_EQ(rs::find_set_folder(root, "m6A").directory, root + "/m6A");
    EXPECT_EQ(rs::find_set_folder(root, "m6A_V2").directory, root + "/m6A_v2");
    EXPECT_EQ(mkdir((root + "/other").c_str(), 0755), 0);
    copy_set("m6A", root + "/other", "M6A");
    EXPECT_EQ(rs::find_set_folder(root, "M6A").name, "other");

    // Two sets with file prefix "1": the prefix no longer names one of them.
    EXPECT_EQ(mkdir((root + "/m1Ψ_copy").c_str(), 0755), 0);
    copy_set("1", root + "/m1Ψ_copy", "1");
    expect_message([&] { (void) rs::find_set_folder(root, "1"); }, "several parameter sets");
    expect_message([&] { (void) rs::find_set_folder(root, "1"); }, "'m1Ψ_copy' (1)");
    EXPECT_EQ(rs::find_set_folder(root, "m1Ψ").name, "m1Ψ");

    // A folder name where an alphabet name was meant, and a modification naming nothing.
    options.modifications = {{"m6A"}};
    expect_message([&] { (void) rs::Datatable::load(options); }, "Modification::alphabet");
    options.data_directory = root + "/turner2004";
    expect_message([&] { (void) rs::Datatable::load(options); }, root + "/m6A");
    options.modifications = {{}};
    expect_message([&] { (void) rs::Datatable::load(options); }, "needs a folder");
}

// The sets shipped in params/RNAstructure, under their folder names.
TEST(ParameterSets, ShippedSets) {
    const std::vector<rs::ParameterSetLocation> sets = rs::list_set_folders(kParams);
    std::vector<std::string> names;
    std::vector<std::string> prefixes;
    for (const rs::ParameterSetLocation& set : sets) {
        names.push_back(set.name);
        prefixes.push_back(set.prefix);
        const rs::AlphabetTables tables = rs::AlphabetTables::load(set.directory);
        EXPECT_EQ(tables.name(), set.name);
        EXPECT_GE(tables.size(), 6u) << set.name;
    }
    const std::vector<std::string> expected_names = {"dmem", "dna",        "m1Ψ", "m1Ψ_full",
                                                     "m6A",  "turner2004", "Ψ",   "Ψ_full"};
    const std::vector<std::string> expected_prefixes = {"dmem", "dna", "1",   "Full1",
                                                        "m6A",  "rna", "PSU", "FullP"};
    EXPECT_EQ(names, expected_names);
    EXPECT_EQ(prefixes, expected_prefixes);

    // Only turner2004 has every enthalpy table.
    EXPECT_NO_THROW((void) rs::AlphabetTables::load(set_dir("turner2004"), "", 300.0));
    EXPECT_THROW((void) rs::AlphabetTables::load(set_dir("dna"), "", 300.0), rs::ParameterError);

    // Every U is m1Ψ in m1Ψ_full, which is meant as a base set.
    rs::LoadOptions options;
    options.data_directory = set_dir("m1Ψ_full");
    const rs::Datatable full = rs::Datatable::load(options);
    EXPECT_EQ(full.encode("GG1UC"), full.encode("GGUUC"));
}
