#pragma once

/**
 * @file RNAstructure/detail/text.hpp
 * @brief Low-level text handling that reproduces how RNAstructure reads its data tables.
 *
 * RNAstructure reads tables with std::getline, trimLeft/trimRight (C-locale isspace),
 * std::istringstream extraction and atof. The helpers here behave the same way, so a file is
 * read into exactly the values efn2 sees, but they work on string_views of one buffer instead
 * of allocating a string per line and per token.
 */

#include "RNAstructure/constants.hpp"
#include "RNAstructure/utils/graphemes.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace thermorna::rnastructure::detail {

/// Byte classes used by the tokenizer.
struct ByteClasses {
    bool space[256] = {};

    constexpr ByteClasses() {
        space[static_cast<unsigned char>(' ')] = true;
        space[static_cast<unsigned char>('\t')] = true;
        space[static_cast<unsigned char>('\n')] = true;
        space[static_cast<unsigned char>('\v')] = true;
        space[static_cast<unsigned char>('\f')] = true;
        space[static_cast<unsigned char>('\r')] = true;
    }
};
inline constexpr ByteClasses kByteClasses{};

/// C-locale isspace: space, \t, \n, \v, \f, \r.
[[nodiscard]] constexpr bool is_space(char c) noexcept {
    return kByteClasses.space[static_cast<unsigned char>(c)];
}

/// Reads a whole file into @p out. Returns false when the file cannot be opened.
[[nodiscard]] bool read_file(const std::string& path, std::string& out);

/// A data line of a table file, with its 1-based line number in the file.
struct Line {
    std::string_view text;
    std::size_t number = 0;
};

/**
 * @brief The data lines of a table file, as RNAstructure's read_datatable_lines returns them.
 *
 * Lines are trimmed on the left, dropped when empty or starting with '#', then trimmed on the
 * right. A leading UTF-8 byte order mark is ignored.
 */
[[nodiscard]] std::vector<Line> data_lines(std::string_view content);

/// Same result as atof()/strtod() in the C locale: parses the longest valid prefix, 0.0 if none.
[[nodiscard]] double parse_double(std::string_view token);

/// floor(atof(token) * 10 + 0.5) through doubles, exactly as RNAstructure computes it.
[[nodiscard]] Energy entry_to_energy_exact(std::string_view token);

/// RNAstructure's floor_entry_to_short(): "." is infinite, otherwise floor(atof(x) * 10 + 0.5).
///
/// Values with at most one decimal (nearly every table entry) are converted without doubles:
/// for them floor(atof(x) * 10 + 0.5) is exactly 10 * x, because the binary rounding error of x
/// stays far below 0.5 after the multiplication. Everything else takes the exact double path.
[[nodiscard]] inline Energy entry_to_energy(std::string_view token) {
    const std::size_t n = token.size();
    std::size_t i = 0;
    bool negative = false;
    if (n == 1 && token[0] == '.') return kInfiniteEnergy;
    if (i < n && (token[i] == '+' || token[i] == '-')) {
        negative = token[i] == '-';
        ++i;
    }
    int value = 0;
    int digits = 0;
    for (; i < n && token[i] >= '0' && token[i] <= '9'; ++i) {
        value = value * 10 + (token[i] - '0');
        if (++digits > 4) return entry_to_energy_exact(token);
    }
    int fraction = 0;
    if (i < n && token[i] == '.') {
        ++i;
        if (i < n && token[i] >= '0' && token[i] <= '9') {
            fraction = token[i] - '0';
            ++digits;
            ++i;
        }
    }
    if (digits == 0 || i != n) return entry_to_energy_exact(token);
    const int tenths = value * 10 + fraction;
    return static_cast<Energy>(negative ? -tenths : tenths);
}

/// Narrows a rounded value to Energy the way a C cast to short behaves on common targets.
[[nodiscard]] Energy narrow_to_energy(double rounded) noexcept;

/**
 * @brief Reads a line the way std::istringstream extraction does.
 *
 * next_token() behaves like `stream >> std::string` and next_grapheme() like `stream >> char`,
 * except that it extracts one grapheme cluster, so multi-byte symbols such as "Ψ" are read as one
 * base. As with a stream, once an extraction fails every later extraction fails as well.
 */
class LineCursor {
   public:
    explicit LineCursor(std::string_view line) noexcept : line_(line) {}

    [[nodiscard]] bool next_token(std::string_view& token) noexcept {
        if (failed_ || !skip_space()) return false;
        const std::size_t start = pos_;
        while (pos_ < line_.size() && !is_space(line_[pos_])) ++pos_;
        token = line_.substr(start, pos_ - start);
        return true;
    }

    [[nodiscard]] bool next_grapheme(std::string_view& grapheme) {
        if (failed_ || !skip_space()) return false;
        std::size_t length = 1;
        const auto next = pos_ + 1 < line_.size() ? static_cast<unsigned char>(line_[pos_ + 1]) : 0;
        if (static_cast<unsigned char>(line_[pos_]) >= 0x80 || next >= 0x80) {
            length = graphemes::first_length(line_.substr(pos_));
        }
        grapheme = line_.substr(pos_, length);
        pos_ += length;
        return true;
    }

    /// Graphemes left on the line (ignoring whitespace), without consuming them.
    [[nodiscard]] std::vector<std::string_view> remaining_graphemes() const;

   private:
    // Skips whitespace; returns false (and fails the cursor) at the end of the line.
    bool skip_space() noexcept {
        while (pos_ < line_.size() && is_space(line_[pos_])) ++pos_;
        if (pos_ < line_.size()) return true;
        failed_ = true;
        return false;
    }

    std::string_view line_;
    std::size_t pos_ = 0;
    bool failed_ = false;
};

/// The grapheme clusters of @p text with all whitespace removed.
[[nodiscard]] std::vector<std::string_view> graphemes_without_space(std::string_view text);

/// Formats "file:line" for messages.
[[nodiscard]] std::string where(const std::string& file, std::size_t line);

}  // namespace thermorna::rnastructure::detail
