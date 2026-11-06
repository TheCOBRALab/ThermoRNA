#include "RNAstructure/params/alphabet_tables.hpp"

#include "RNAstructure/detail/text.hpp"
#include "RNAstructure/error.hpp"
#include "RNAstructure/utils/graphemes.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

namespace thermorna::rnastructure {

const char* table_file_name(Table4D table) noexcept {
    switch (table) {
        case Table4D::Stack:      return "stack";
        case Table4D::HelixEnds:  return "helix_ends";
        case Table4D::Tstkh:      return "tstackh";
        case Table4D::Tstki:      return "tstacki";
        case Table4D::Tstki23:    return "tstacki23";
        case Table4D::Tstki1n:    return "tstacki1n";
        case Table4D::Coax:       return "coaxial";
        case Table4D::TstackCoax: return "tstackcoax";
        case Table4D::CoaxStack:  return "coaxstack";
        case Table4D::Tstack:     return "tstack";
        case Table4D::Tstkm:      return "tstackm";
    }
    return "";
}

namespace {

// RNAstructure's Tscale(): dH - (short)floor((float)(dH - dG) * T / T37inK + 0.5).
[[nodiscard]] Energy tscale(float kelvin, Energy dG, Energy dH) noexcept {
    if (dG == kInfiniteEnergy) return kInfiniteEnergy;
    const float product = static_cast<float>(dH - dG) * kelvin;
    const Energy shift =
        detail::narrow_to_energy(std::floor(static_cast<double>(product) / kTemperature37C + 0.5));
    return static_cast<Energy>(dH - shift);
}

// Collects rows that RNAstructure reads in a lenient way, so they can be reported once per file.
class Diagnostics {
   public:
    explicit Diagnostics(const std::string& file) : file_(file) {}

    void missing_as_zero(std::size_t line) { add(zero_, line); }
    void missing_as_previous(std::size_t line) { add(previous_, line); }
    void extra_values(std::size_t line) { add(extra_, line); }

    void flush(std::vector<std::string>& warnings) const {
        report(warnings, zero_,
               "had fewer values than columns; the missing entries were read as 0, as "
               "RNAstructure does");
        report(warnings, previous_,
               "had fewer values than expected; the missing entries reuse the previous value, "
               "as RNAstructure does");
        report(warnings, extra_, "had more values than columns; the extra values were ignored");
    }

   private:
    struct Lines {
        std::size_t count = 0;
        std::vector<std::size_t> first;
    };

    static void add(Lines& lines, std::size_t line) {
        if (!lines.first.empty() && lines.first.back() == line) return;
        ++lines.count;
        if (lines.first.size() < 5) lines.first.push_back(line);
    }

    void report(std::vector<std::string>& warnings, const Lines& lines, const char* what) const {
        if (lines.count == 0) return;
        std::ostringstream message;
        message << file_ << ": " << lines.count << (lines.count == 1 ? " row " : " rows ") << what
                << " (line";
        if (lines.first.size() > 1) message << 's';
        for (std::size_t k = 0; k < lines.first.size(); ++k) {
            message << (k == 0 ? " " : ", ") << lines.first[k];
        }
        if (lines.count > lines.first.size()) message << ", ...";
        message << ").";
        warnings.push_back(message.str());
    }

    const std::string& file_;
    Lines zero_;
    Lines previous_;
    Lines extra_;
};

}  // namespace

// Reads the files of one alphabet into an AlphabetTables, mirroring the readers in
// RNAstructure's rna_library.cpp (read_loop, read_dangle, read_4D_table, ...).
class TableReader {
   public:
    TableReader(AlphabetTables& tables, const ParameterSetLocation& location, TableKind kind)
        : t_(tables),
          location_(location),
          name_(location.name),
          extension_(kind == TableKind::FreeEnergy ? ".dg" : ".dh") {}

    // Same files and order as datatable::opendat().
    void read_all() {
        const std::string specification = location_.file("specification", ".dat");
        std::string text;
        if (!detail::read_file(specification, text)) {
            throw ParameterError("Missing alphabet specification file: " + specification);
        }
        t_.alphabet_ = Alphabet::parse(text, name_, specification);
        const std::size_t n = t_.alphabet_.size();
        t_.n_ = n;
        t_.n4_ = n * n * n * n;
        t_.table4_.assign(kTable4DCount * t_.n4_, kInfiniteEnergy);
        t_.dangle_.assign(n * n * n * 3, kInfiniteEnergy);
        init_blocks(t_.iloop11_, t_.n4_, n * n);
        init_blocks(t_.iloop21_, t_.n4_, n * n * n);
        init_blocks(t_.iloop22_, t_.n4_, n * n * n * n);

        read_loop();
        read_dangle();
        read_table4(Table4D::Stack);
        read_table4(Table4D::HelixEnds);
        read_table4(Table4D::Tstkh);
        read_table4(Table4D::Tstki);
        read_table4(Table4D::Tstki23);
        read_table4(Table4D::Tstki1n);
        read_table4(Table4D::Coax);
        read_table4(Table4D::TstackCoax);
        read_table4(Table4D::CoaxStack);
        read_table4(Table4D::Tstack);
        read_table4(Table4D::Tstkm);
        read_int11();
        read_int21();
        read_int22();
        read_special_hairpins(t_.tloop_, "tloop", 6);
        read_special_hairpins(t_.triloop_, "triloop", 5);
        read_special_hairpins(t_.hexaloop_, "hexaloop", 8);
        read_miscloop();
    }

   private:
    using Lines = std::vector<detail::Line>;

    [[nodiscard]] std::string path(const char* table) const {
        return location_.file(table, extension_);
    }

    // Reads one file. The returned views stay valid until the next call.
    [[nodiscard]] Lines lines_of(const std::string& file) {
        if (!detail::read_file(file, content_)) {
            throw ParameterError("Missing data file: " + file);
        }
        return detail::data_lines(content_);
    }

    // Blocks are keyed by the four closing-pair bases (n^4 keys) and hold block_size values.
    static void init_blocks(AlphabetTables::BlockTable& table, std::size_t keys,
                            std::size_t block_size) {
        table.block_size = block_size;
        table.values.clear();
        table.block_of.assign(keys, -1);
    }

    [[nodiscard]] static std::size_t block(AlphabetTables::BlockTable& table, std::size_t key) {
        std::int32_t& slot = table.block_of[key];
        if (slot < 0) {
            slot = static_cast<std::int32_t>(table.values.size() / table.block_size);
            table.values.resize(table.values.size() + table.block_size, kInfiniteEnergy);
        }
        return static_cast<std::size_t>(slot);
    }

    [[noreturn]] void fail(const std::string& file, std::size_t line,
                           const std::string& message) const {
        throw ParameterError(detail::where(file, line) + ": " + message);
    }

    // Throws when a block needs a line past the end of the file (RNAstructure reads out of
    // bounds here).
    void need(const Lines& lines, std::size_t i, const std::string& file, const char* what) const {
        if (i >= lines.size()) {
            const std::size_t last = lines.empty() ? 0 : lines.back().number;
            fail(file, last, std::string("unexpected end of file, expected ") + what);
        }
    }

    [[nodiscard]] Code code_of(std::string_view symbol, const std::string& file,
                               std::size_t line) const {
        const Code code = t_.alphabet_.code(symbol);
        if (code == kInvalidCode) {
            fail(file, line,
                 "'" + std::string(symbol) + "' is not a base of alphabet '" + name_ + "'");
        }
        return code;
    }

    [[nodiscard]] Code first_code(const detail::Line& line, const std::string& file) const {
        const std::size_t length = graphemes::first_length(line.text);
        return code_of(line.text.substr(0, length), file, line.number);
    }

    // `stream >> char` repeated until the line is exhausted.
    [[nodiscard]] std::vector<Code> codes_on_line(const detail::Line& line,
                                                  const std::string& file) const {
        std::vector<Code> codes;
        detail::LineCursor cursor(line.text);
        std::string_view symbol;
        while (cursor.next_grapheme(symbol)) codes.push_back(code_of(symbol, file, line.number));
        return codes;
    }

    // Two `stream >> char` extractions.
    void two_codes(const detail::Line& line, const std::string& file, Code& first,
                   Code& second) const {
        detail::LineCursor cursor(line.text);
        std::string_view a;
        std::string_view b;
        if (!cursor.next_grapheme(a) || !cursor.next_grapheme(b)) {
            fail(file, line.number, "expected two bases, got '" + std::string(line.text) + "'");
        }
        first = code_of(a, file, line.number);
        second = code_of(b, file, line.number);
    }

    // A value inside a table row. RNAstructure declares a fresh std::string for every value,
    // so a value missing at the end of a short row is atof("") = 0.
    [[nodiscard]] static Energy row_value(detail::LineCursor& cursor, Diagnostics& diagnostics,
                                          std::size_t line) {
        std::string_view token;
        if (cursor.next_token(token)) return detail::entry_to_energy(token);
        diagnostics.missing_as_zero(line);
        return 0;
    }

    static void check_extra(detail::LineCursor& cursor, Diagnostics& diagnostics,
                            std::size_t line) {
        std::string_view token;
        if (cursor.next_token(token)) diagnostics.extra_values(line);
    }

    // ----------------------------- loop (read_loop) -----------------------------
    void read_loop() {
        const std::string file = path("loop");
        const Lines lines = lines_of(file);
        Diagnostics diagnostics(file);

        t_.inter_.assign(1, kInfiniteEnergy);
        t_.bulge_.assign(1, kInfiniteEnergy);
        t_.hairpin_.assign(1, kInfiniteEnergy);

        // temp0..temp3 are declared once in RNAstructure and keep their value on a short line.
        std::string_view temp[4] = {"", "", "", ""};
        for (const detail::Line& line : lines) {
            detail::LineCursor cursor(line.text);
            bool complete = true;
            for (std::string_view& value : temp) {
                std::string_view token;
                if (cursor.next_token(token)) {
                    value = token;
                } else {
                    complete = false;
                }
            }
            if (!complete) diagnostics.missing_as_previous(line.number);
            t_.inter_.push_back(detail::entry_to_energy(temp[1]));
            t_.bulge_.push_back(detail::entry_to_energy(temp[2]));
            t_.hairpin_.push_back(detail::entry_to_energy(temp[3]));
        }
        diagnostics.flush(t_.warnings_);
    }

    // ----------------------------- dangle (read_dangle) -------------------------
    void read_dangle() {
        const std::string file = path("dangle");
        const Lines lines = lines_of(file);
        Diagnostics diagnostics(file);

        for (std::size_t i = 0; i < lines.size();) {
            // "AX" is a 3' dangle (d = 1); anything else, e.g. "A", is a 5' dangle (d = 2).
            const std::vector<std::string_view> head = graphemes::split(lines[i].text);
            const Code a = code_of(head.front(), file, lines[i].number);
            const std::size_t d = (head.size() == 2 && head[1] == "X") ? 1 : 2;
            ++i;
            need(lines, i, file, "the second line of a dangle header");
            const Code b = first_code(lines[i], file);
            ++i;
            need(lines, i, file, "the dangling base line");
            const std::vector<Code> columns = codes_on_line(lines[i], file);
            ++i;
            need(lines, i, file, "a dangle data line");
            detail::LineCursor cursor(lines[i].text);
            for (const Code c : columns) {
                t_.dangle_[t_.index3(a, b, c) * 3 + d] =
                    row_value(cursor, diagnostics, lines[i].number);
            }
            check_extra(cursor, diagnostics, lines[i].number);
            ++i;
        }
        diagnostics.flush(t_.warnings_);
    }

    // ----------------------------- 4x4 tables (read_4D_table) -------------------
    void read_table4(Table4D table) {
        const std::string file = path(table_file_name(table));
        const Lines lines = lines_of(file);
        Diagnostics diagnostics(file);
        const std::size_t offset = t_.table_offset(table);

        for (std::size_t i = 0; i < lines.size();) {
            const Code a = first_code(lines[i], file);
            ++i;
            need(lines, i, file, "the second line of a block header");
            const Code b = first_code(lines[i], file);
            ++i;
            need(lines, i, file, "the column header line");
            const std::vector<Code> columns = codes_on_line(lines[i], file);
            ++i;
            for (std::size_t x = 0; x < columns.size(); ++x) {
                need(lines, i, file, "a table row");
                detail::LineCursor cursor(lines[i].text);
                std::string_view symbol;
                if (!cursor.next_grapheme(symbol)) fail(file, lines[i].number, "empty row");
                const Code c = code_of(symbol, file, lines[i].number);
                for (const Code d : columns) {
                    t_.table4_[offset + t_.index4(a, b, c, d)] =
                        row_value(cursor, diagnostics, lines[i].number);
                }
                check_extra(cursor, diagnostics, lines[i].number);
                ++i;
            }
        }
        diagnostics.flush(t_.warnings_);
    }

    // ----------------------------- 1x1 loops (read_6D_table) --------------------
    void read_int11() {
        const std::string file = path("int11");
        const Lines lines = lines_of(file);
        Diagnostics diagnostics(file);
        AlphabetTables::BlockTable& table = t_.iloop11_;

        for (std::size_t i = 0; i < lines.size();) {
            ++i;  // "X"
            need(lines, i, file, "the first closing-pair line");
            Code a = 0;
            Code c = 0;
            two_codes(lines[i], file, a, c);
            ++i;
            need(lines, i, file, "the second closing-pair line");
            Code d = 0;
            Code f = 0;
            two_codes(lines[i], file, d, f);
            i += 2;  // "Y"
            need(lines, i, file, "the column header line");
            const std::vector<Code> columns = codes_on_line(lines[i], file);
            ++i;
            const std::size_t base = block(table, t_.index4(a, c, d, f)) * table.block_size;
            for (std::size_t x = 0; x < columns.size(); ++x) {
                need(lines, i, file, "a table row");
                detail::LineCursor cursor(lines[i].text);
                std::string_view symbol;
                if (!cursor.next_grapheme(symbol)) fail(file, lines[i].number, "empty row");
                const Code b = code_of(symbol, file, lines[i].number);
                for (const Code e : columns) {
                    table.values[base + t_.index2(b, e)] =
                        row_value(cursor, diagnostics, lines[i].number);
                }
                check_extra(cursor, diagnostics, lines[i].number);
                ++i;
            }
        }
        diagnostics.flush(t_.warnings_);
    }

    // ----------------------------- 2x1 loops (read_7D_table) --------------------
    void read_int21() {
        const std::string file = path("int21");
        const Lines lines = lines_of(file);
        Diagnostics diagnostics(file);
        AlphabetTables::BlockTable& table = t_.iloop21_;

        for (std::size_t i = 0; i < lines.size();) {
            ++i;  // "X"
            need(lines, i, file, "the first closing-pair line");
            Code a = 0;
            Code f = 0;
            two_codes(lines[i], file, a, f);
            ++i;
            need(lines, i, file, "the second closing-pair line");
            Code b = 0;
            Code g = 0;
            two_codes(lines[i], file, b, g);
            ++i;
            need(lines, i, file, "the \"YA\" line");
            const std::vector<std::string_view> ya = detail::graphemes_without_space(lines[i].text);
            if (ya.size() < 2) fail(file, lines[i].number, "expected a line such as \"YA\"");
            const Code e = code_of(ya[1], file, lines[i].number);
            ++i;
            need(lines, i, file, "the column header line");
            const std::vector<Code> columns = codes_on_line(lines[i], file);
            ++i;
            const std::size_t base = block(table, t_.index4(a, b, f, g)) * table.block_size;
            for (std::size_t x = 0; x < columns.size(); ++x) {
                need(lines, i, file, "a table row");
                detail::LineCursor cursor(lines[i].text);
                std::string_view symbol;
                if (!cursor.next_grapheme(symbol)) fail(file, lines[i].number, "empty row");
                const Code c = code_of(symbol, file, lines[i].number);
                for (const Code d : columns) {
                    table.values[base + t_.index3(c, d, e)] =
                        row_value(cursor, diagnostics, lines[i].number);
                }
                check_extra(cursor, diagnostics, lines[i].number);
                ++i;
            }
        }
        diagnostics.flush(t_.warnings_);
    }

    // ----------------------------- 2x2 loops (read_8D_table) --------------------
    void read_int22() {
        const std::string file = path("int22");
        const Lines lines = lines_of(file);
        Diagnostics diagnostics(file);
        AlphabetTables::BlockTable& table = t_.iloop22_;

        auto pair_header = [&](const detail::Line& line, Code& first, Code& last) {
            std::vector<std::string_view> tokens;
            detail::LineCursor cursor(line.text);
            std::string_view token;
            while (cursor.next_token(token)) tokens.push_back(token);
            if (tokens.size() < 4) {
                fail(file, line.number,
                     "expected a header such as \"A X1 Y1 A\", got '" + std::string(line.text) +
                         "'");
            }
            first =
                code_of(tokens[0].substr(0, graphemes::first_length(tokens[0])), file, line.number);
            last =
                code_of(tokens[3].substr(0, graphemes::first_length(tokens[3])), file, line.number);
        };

        for (std::size_t i = 0; i < lines.size();) {
            Code a = 0;
            Code b = 0;
            pair_header(lines[i], a, b);
            ++i;
            need(lines, i, file, "the second closing-pair line");
            Code c = 0;
            Code d = 0;
            pair_header(lines[i], c, d);
            ++i;
            need(lines, i, file, "the first column header line");
            const std::vector<Code> y1 = codes_on_line(lines[i], file);
            ++i;
            need(lines, i, file, "the second column header line");
            const std::vector<Code> y2 = codes_on_line(lines[i], file);
            if (y2.size() < y1.size()) {
                fail(file, lines[i].number, "the two column header lines differ in length");
            }
            ++i;
            const std::size_t base = block(table, t_.index4(a, b, c, d)) * table.block_size;
            for (std::size_t x = 0; x < y1.size(); ++x) {
                need(lines, i, file, "a table row");
                detail::LineCursor cursor(lines[i].text);
                std::string_view label;
                if (!cursor.next_token(label)) fail(file, lines[i].number, "empty row");
                const std::vector<std::string_view> pair = graphemes::split(label);
                if (pair.size() < 2) {
                    fail(file, lines[i].number,
                         "a row must start with two bases, got '" + std::string(label) + "'");
                }
                const Code e = code_of(pair[0], file, lines[i].number);
                const Code g = code_of(pair[1], file, lines[i].number);
                for (std::size_t y = 0; y < y1.size(); ++y) {
                    table.values[base + t_.index4(e, y1[y], g, y2[y])] =
                        row_value(cursor, diagnostics, lines[i].number);
                }
                check_extra(cursor, diagnostics, lines[i].number);
                ++i;
            }
        }
        diagnostics.flush(t_.warnings_);
    }

    // ----------------------------- tloop, triloop, hexaloop (read_xloop) --------
    void read_special_hairpins(AlphabetTables::SpecialHairpinTable& table, const char* name,
                               std::size_t expected_length) {
        const std::string file = path(name);
        const Lines lines = lines_of(file);
        Diagnostics diagnostics(file);

        table.entries.clear();
        table.entries.reserve(lines.size());
        std::size_t wrong_length = 0;

        // seq and temp are declared once in RNAstructure and keep their value on a short line.
        std::string_view sequence = "";
        std::string_view energy = "";
        for (const detail::Line& line : lines) {
            detail::LineCursor cursor(line.text);
            std::string_view token;
            bool complete = true;
            if (cursor.next_token(token)) {
                sequence = token;
            } else {
                complete = false;
            }
            if (cursor.next_token(token)) {
                energy = token;
            } else {
                complete = false;
            }
            if (!complete) diagnostics.missing_as_previous(line.number);

            SpecialHairpin entry;
            entry.sequence = std::string(sequence);
            for (std::string_view symbol : graphemes::split(sequence)) {
                entry.codes.push_back(code_of(symbol, file, line.number));
            }
            entry.energy = detail::entry_to_energy(energy);
            if (entry.codes.size() != expected_length) ++wrong_length;
            table.entries.push_back(std::move(entry));
        }

        // Index by sequence; the first entry of a sequence wins, like efn2's linear search.
        table.index.clear();
        for (std::size_t k = 0; k < table.entries.size(); ++k) {
            const std::vector<Code>& codes = table.entries[k].codes;
            const std::uint64_t key = AlphabetTables::pack_codes(codes.data(), codes.size());
            if (key != 0) table.index.emplace_back(key, static_cast<std::uint32_t>(k));
        }
        std::sort(table.index.begin(), table.index.end());
        table.index.erase(
            std::unique(table.index.begin(), table.index.end(),
                        [](const auto& x, const auto& y) { return x.first == y.first; }),
            table.index.end());

        diagnostics.flush(t_.warnings_);
        if (wrong_length != 0) {
            t_.warnings_.push_back(file + ": " + std::to_string(wrong_length) +
                                   " entries are not " + std::to_string(expected_length) +
                                   " nucleotides long and can never match a loop.");
        }
    }

    // ----------------------------- miscloop (read_miscloop) ---------------------
    void read_miscloop() {
        const std::string file = path("miscloop");
        const Lines lines = lines_of(file);
        Diagnostics diagnostics(file);
        MiscLoop& misc = t_.misc_;

        // RNAstructure 6.5+ files have 13 data lines. Older ones have two more lines (terminal
        // AU penalty and whether it applies to GU) that helix_ends now replaces.
        const bool legacy = lines.size() == 15 && (lines[8].text.substr(0, 4) == "true" ||
                                                   lines[8].text.substr(0, 4) == "fals");
        if (lines.size() != 13 && !legacy) {
            throw ParameterError(file +
                                 ": expected 13 data lines (or 15 in the pre-6.5 format), "
                                 "found " +
                                 std::to_string(lines.size()));
        }

        // RNAstructure reads every value through one std::string `temp`, so a value missing
        // from a short line reuses whatever was read last.
        std::string_view temp = "";
        std::size_t i = 0;
        auto whole_line = [&]() {
            temp = lines[i].text;
            ++i;
            return temp;
        };
        auto value = [&](detail::LineCursor& cursor, std::size_t line) {
            std::string_view token;
            if (cursor.next_token(token)) {
                temp = token;
            } else {
                diagnostics.missing_as_previous(line);
            }
            return detail::entry_to_energy(temp);
        };

        misc.prelog = static_cast<float>(detail::parse_double(whole_line()) * kConversionFactor);
        misc.maxpen = detail::entry_to_energy(whole_line());
        {
            detail::LineCursor cursor(lines[i].text);
            const std::size_t number = lines[i].number;
            ++i;
            for (std::size_t k = 1; k <= 4; ++k) misc.poppen[k] = value(cursor, number);
        }
        {
            detail::LineCursor cursor(lines[i].text);
            const std::size_t number = lines[i].number;
            ++i;
            misc.eparam[1] = 0;
            misc.eparam[2] = 0;
            misc.eparam[3] = 0;
            misc.eparam[4] = 0;
            misc.eparam[5] = value(cursor, number);
            misc.eparam[6] = value(cursor, number);
            misc.eparam[7] = 30;
            misc.eparam[8] = 30;
            misc.eparam[9] = -500;
            misc.eparam[10] = value(cursor, number);
        }
        {
            detail::LineCursor cursor(lines[i].text);
            const std::size_t number = lines[i].number;
            ++i;
            misc.efn2a = value(cursor, number);
            misc.efn2b = value(cursor, number);
            misc.efn2c = value(cursor, number);
        }
        misc.mlasym = detail::entry_to_energy(whole_line());
        misc.strain = detail::entry_to_energy(whole_line());
        if (legacy) {
            const std::size_t number = lines[i].number;
            i += 2;
            t_.warnings_.push_back(detail::where(file, number) +
                                   ": terminal AU penalty lines of the pre-6.5 format were "
                                   "ignored; terminal penalties come from helix_ends.");
        }
        misc.gubonus = detail::entry_to_energy(whole_line());
        misc.cslope = detail::entry_to_energy(whole_line());
        misc.cint = detail::entry_to_energy(whole_line());
        misc.c3 = detail::entry_to_energy(whole_line());
        misc.init = detail::entry_to_energy(whole_line());
        misc.singlecbulge = detail::entry_to_energy(whole_line());
        diagnostics.flush(t_.warnings_);
    }

    AlphabetTables& t_;
    const ParameterSetLocation& location_;
    const std::string& name_;
    const char* extension_;
    std::string content_;
};

// ----------------------------------------------------------------------------
// AlphabetTables
// ----------------------------------------------------------------------------

AlphabetTables AlphabetTables::read(const std::string& directory, const std::string& alphabet,
                                    TableKind kind) {
    const ParameterSetLocation location = locate_parameter_set(directory, alphabet);
    AlphabetTables tables;
    tables.kind_ = kind;
    tables.directory_ = directory;
    tables.prefix_ = location.prefix;
    TableReader(tables, location, kind).read_all();
    tables.process_dat();
    return tables;
}

AlphabetTables AlphabetTables::load(const std::string& directory, const std::string& alphabet,
                                    double kelvin, MissingEnthalpy missing) {
    AlphabetTables tables = read(directory, alphabet, TableKind::FreeEnergy);
    if (std::abs(kelvin - kTemperature37C) < kTemperatureTolerance) return tables;

    AlphabetTables enthalpy;
    try {
        enthalpy = read(directory, tables.prefix(), TableKind::Enthalpy);
    } catch (const ParameterError& error) {
        std::ostringstream message;
        message << "Alphabet '" << tables.name() << "' cannot be rescaled to " << kelvin
                << " K because its enthalpy tables could not be read (" << error.what() << ")";
        if (missing == MissingEnthalpy::Fail) throw ParameterError(message.str());
        message << "; its free energies stay at 37 degrees C.";
        tables.warnings_.push_back(message.str());
        return tables;
    }

    for (const std::string& warning : enthalpy.warnings_) tables.warnings_.push_back(warning);
    tables.scale_to_temperature(enthalpy, kelvin);
    return tables;
}

void AlphabetTables::process_dat() {
    const std::size_t n = n_;
    const Alphabet& spec = alphabet_;

    auto set_terminal_stacks = [&](std::size_t index, Energy value) {
        table4_[table_offset(Table4D::Tstkh) + index] = value;
        table4_[table_offset(Table4D::Tstki) + index] = value;
        table4_[table_offset(Table4D::Tstki23) + index] = value;
        table4_[table_offset(Table4D::Tstki1n) + index] = value;
    };

    // Non-interacting bases (X): zero dangles and terminal mismatches.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const Code ci = static_cast<Code>(i);
            const Code cj = static_cast<Code>(j);
            if (!spec.can_interact(ci) || !spec.can_interact(cj)) continue;
            for (const Code k : spec.non_interacting_list()) {
                dangle_[index3(ci, cj, k) * 3 + 0] = 0;
                dangle_[index3(ci, cj, k) * 3 + 1] = 0;
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = 0; k < n; ++k) {
                for (std::size_t l = 0; l < n; ++l) {
                    const Code ci = static_cast<Code>(i);
                    const Code cj = static_cast<Code>(j);
                    const Code ck = static_cast<Code>(k);
                    const Code cl = static_cast<Code>(l);
                    if (spec.can_interact(ck) && spec.can_interact(cl)) continue;
                    if (!spec.can_interact(ci) || !spec.can_interact(cj)) continue;
                    const std::size_t index = index4(ci, cj, ck, cl);
                    set_terminal_stacks(index, 0);
                    table4_[table_offset(Table4D::Tstack) + index] = 0;
                    table4_[table_offset(Table4D::Tstkm) + index] = 0;
                }
            }
        }
    }

    // Intermolecular linker (I): a mismatch next to the linker is really a dangle.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            const Code ci = static_cast<Code>(i);
            const Code cj = static_cast<Code>(j);
            if (!spec.can_interact(ci) || !spec.can_interact(cj)) continue;
            for (const Code k : spec.linker_list()) {
                dangle_[index3(ci, cj, k) * 3 + 0] = 0;
                dangle_[index3(ci, cj, k) * 3 + 1] = 0;
            }
        }
    }
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            for (std::size_t k = 0; k < n; ++k) {
                for (std::size_t l = 0; l < n; ++l) {
                    const Code ci = static_cast<Code>(i);
                    const Code cj = static_cast<Code>(j);
                    const Code ck = static_cast<Code>(k);
                    const Code cl = static_cast<Code>(l);
                    if (!spec.can_interact(ci) || !spec.can_interact(cj)) continue;
                    const bool linker_k = spec.is_linker(ck);
                    const bool linker_l = spec.is_linker(cl);
                    if (!linker_k && !linker_l) continue;

                    const std::size_t index = index4(ci, cj, ck, cl);
                    set_terminal_stacks(index, 0);
                    Energy value = 0;
                    if (linker_k && !linker_l) {
                        value = dangle_[index3(ci, cj, cl) * 3 + 2];
                    } else if (!linker_k) {
                        value = dangle_[index3(ci, cj, ck) * 3 + 1];
                    }
                    table4_[table_offset(Table4D::Tstack) + index] = value;
                    table4_[table_offset(Table4D::Tstkm) + index] = value;
                }
            }
        }
    }
}

void AlphabetTables::scale_to_temperature(const AlphabetTables& enthalpy, double kelvin) {
    if (enthalpy.n_ != n_) {
        throw ParameterError("Enthalpy tables of alphabet '" + enthalpy.alphabet_.name() +
                             "' have a different alphabet size than its free energy tables");
    }
    const AlphabetTables& dh = enthalpy;
    const float T = static_cast<float>(kelvin);
    const float t37 = static_cast<float>(kTemperature37C);

    misc_.prelog = misc_.prelog * T / t37;
    RT_ = RT_ * T / t37;

    auto scale_hairpins = [&](SpecialHairpinTable& table, const SpecialHairpinTable& other,
                              const char* name) {
        if (other.entries.size() < table.entries.size()) {
            throw ParameterError(std::string("The enthalpy ") + name + " table of alphabet '" +
                                 alphabet_.name() + "' has fewer entries than the free energy one");
        }
        for (std::size_t k = 0; k < table.entries.size(); ++k) {
            if (table.entries[k].codes != other.entries[k].codes) {
                throw ParameterError(std::string("The ") + name + " tables of alphabet '" +
                                     alphabet_.name() + "' list different sequences in entry " +
                                     std::to_string(k + 1) + " ('" + table.entries[k].sequence +
                                     "' vs '" + other.entries[k].sequence + "')");
            }
            table.entries[k].energy = tscale(T, static_cast<Energy>(table.entries[k].energy),
                                             static_cast<Energy>(other.entries[k].energy));
        }
    };
    scale_hairpins(tloop_, dh.tloop_, "tloop");
    scale_hairpins(triloop_, dh.triloop_, "triloop");
    scale_hairpins(hexaloop_, dh.hexaloop_, "hexaloop");

    for (std::size_t a = 1; a < 5; ++a)
        misc_.poppen[a] = tscale(T, misc_.poppen[a], dh.misc_.poppen[a]);
    misc_.maxpen = tscale(T, misc_.maxpen, dh.misc_.maxpen);
    for (std::size_t a = 1; a < 7; ++a)
        misc_.eparam[a] = tscale(T, misc_.eparam[a], dh.misc_.eparam[a]);
    misc_.eparam[10] = tscale(T, misc_.eparam[10], dh.misc_.eparam[10]);

    for (std::size_t a = 1; a < 31; ++a) {
        if (a < inter_.size() && a < dh.inter_.size())
            inter_[a] = tscale(T, inter_[a], dh.inter_[a]);
        if (a < bulge_.size() && a < dh.bulge_.size())
            bulge_[a] = tscale(T, bulge_[a], dh.bulge_[a]);
        if (a < hairpin_.size() && a < dh.hairpin_.size()) {
            hairpin_[a] = tscale(T, hairpin_[a], dh.hairpin_[a]);
        }
    }

    misc_.gubonus = tscale(T, misc_.gubonus, dh.misc_.gubonus);
    misc_.cint = tscale(T, misc_.cint, dh.misc_.cint);
    misc_.cslope = tscale(T, misc_.cslope, dh.misc_.cslope);
    misc_.c3 = tscale(T, misc_.c3, dh.misc_.c3);
    misc_.efn2a = tscale(T, misc_.efn2a, dh.misc_.efn2a);
    misc_.efn2b = tscale(T, misc_.efn2b, dh.misc_.efn2b);
    misc_.efn2c = tscale(T, misc_.efn2c, dh.misc_.efn2c);
    misc_.init = tscale(T, misc_.init, dh.misc_.init);
    misc_.mlasym = tscale(T, misc_.mlasym, dh.misc_.mlasym);
    misc_.strain = tscale(T, misc_.strain, dh.misc_.strain);
    misc_.singlecbulge = tscale(T, misc_.singlecbulge, dh.misc_.singlecbulge);

    for (std::size_t k = 0; k < table4_.size(); ++k)
        table4_[k] = tscale(T, table4_[k], dh.table4_[k]);
    for (std::size_t k = 0; k < dangle_.size(); k += 3) {
        dangle_[k + 1] = tscale(T, dangle_[k + 1], dh.dangle_[k + 1]);
        dangle_[k + 2] = tscale(T, dangle_[k + 2], dh.dangle_[k + 2]);
    }

    auto scale_blocks = [&](BlockTable& table, const BlockTable& other) {
        for (std::size_t key = 0; key < table.block_of.size(); ++key) {
            if (table.block_of[key] < 0) continue;
            const std::size_t start =
                static_cast<std::size_t>(table.block_of[key]) * table.block_size;
            for (std::size_t inner = 0; inner < table.block_size; ++inner) {
                table.values[start + inner] =
                    tscale(T, table.values[start + inner], lookup(other, key, inner));
            }
        }
    };
    scale_blocks(iloop11_, dh.iloop11_);
    scale_blocks(iloop21_, dh.iloop21_);
    scale_blocks(iloop22_, dh.iloop22_);

    temperature_ = kelvin;
}

std::uint64_t AlphabetTables::pack_codes(const Code* codes, std::size_t length) noexcept {
    if (length == 0 || length > 8) return 0;
    std::uint64_t key = 0;
    for (std::size_t k = 0; k < length; ++k) {
        if (codes[k] >= kMaxCodes) return 0;
        key = (key << 8) | (static_cast<std::uint64_t>(codes[k]) + 1);
    }
    return key;
}

std::optional<int> AlphabetTables::find(const SpecialHairpinTable& table, const Code* codes,
                                        std::size_t length) noexcept {
    const std::uint64_t key = pack_codes(codes, length);
    if (key == 0) return std::nullopt;
    const auto it = std::lower_bound(table.index.begin(), table.index.end(), key,
                                     [](const std::pair<std::uint64_t, std::uint32_t>& entry,
                                        std::uint64_t k) { return entry.first < k; });
    if (it == table.index.end() || it->first != key) return std::nullopt;
    return table.entries[it->second].energy;
}

}  // namespace thermorna::rnastructure
