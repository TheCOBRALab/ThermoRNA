#include "RNAstructure/detail/text.hpp"

#include "RNAstructure/utils/graphemes.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>

namespace thermorna::rnastructure::detail {

namespace {

[[nodiscard]] constexpr bool is_digit(char c) noexcept { return c >= '0' && c <= '9'; }

// Powers of ten that are exactly representable as doubles.
constexpr double kExactPowersOfTen[] = {1e0,  1e1,  1e2,  1e3,  1e4,  1e5,  1e6,  1e7,
                                        1e8,  1e9,  1e10, 1e11, 1e12, 1e13, 1e14, 1e15,
                                        1e16, 1e17, 1e18, 1e19, 1e20, 1e21, 1e22};
constexpr int kMaxExactPower = 22;
constexpr std::uint64_t kMaxExactMantissa = std::uint64_t{1} << 53;
constexpr std::uint64_t kMantissaLimit = 100000000000000000ULL;  // 1e17, leaves room for a digit

double strtod_fallback(std::string_view token) {
    const std::string copy(token);
    return std::strtod(copy.c_str(), nullptr);
}

}  // namespace

bool read_file(const std::string& path, std::string& out) {
    std::FILE* file = std::fopen(path.c_str(), "rb");
    if (file == nullptr) return false;

    out.clear();
    // Read in one go when the size is known; fall back to chunks otherwise (pipes, ...).
    if (std::fseek(file, 0, SEEK_END) == 0) {
        const long length = std::ftell(file);
        if (length > 0 && std::fseek(file, 0, SEEK_SET) == 0) {
            out.resize(static_cast<std::size_t>(length));
            const std::size_t got = std::fread(&out[0], 1, out.size(), file);
            out.resize(got);
        } else {
            std::rewind(file);
        }
    }
    char buffer[1 << 14];
    std::size_t got = 0;
    while ((got = std::fread(buffer, 1, sizeof(buffer), file)) > 0) out.append(buffer, got);

    const bool ok = std::ferror(file) == 0;
    std::fclose(file);
    return ok;
}

std::vector<Line> data_lines(std::string_view content) {
    std::vector<Line> lines;
    lines.reserve(content.size() / 16 + 1);

    if (content.size() >= 3 && content.compare(0, 3, "\xEF\xBB\xBF") == 0) {
        content.remove_prefix(3);
    }

    std::size_t number = 0;
    std::size_t start = 0;
    while (start < content.size()) {
        std::size_t end = content.find('\n', start);
        if (end == std::string_view::npos) end = content.size();
        ++number;

        std::string_view text = content.substr(start, end - start);
        start = end + 1;

        std::size_t left = 0;
        while (left < text.size() && is_space(text[left])) ++left;
        text.remove_prefix(left);
        if (text.empty() || text.front() == '#') continue;

        std::size_t right = text.size();
        while (right > 0 && is_space(text[right - 1])) --right;
        lines.push_back(Line{text.substr(0, right), number});
    }
    return lines;
}

double parse_double(std::string_view token) {
    std::size_t i = 0;
    const std::size_t n = token.size();
    while (i < n && is_space(token[i])) ++i;

    bool negative = false;
    if (i < n && (token[i] == '+' || token[i] == '-')) {
        negative = token[i] == '-';
        ++i;
    }

    // Hexadecimal floats, inf, nan and other unusual spellings are left to strtod itself.
    if (i < n && !is_digit(token[i]) && token[i] != '.') return strtod_fallback(token);
    if (i + 1 < n && token[i] == '0' && (token[i + 1] == 'x' || token[i + 1] == 'X')) {
        return strtod_fallback(token);
    }

    std::uint64_t mantissa = 0;
    int exponent = 0;
    bool any_digit = false;

    for (; i < n && is_digit(token[i]); ++i) {
        if (mantissa >= kMantissaLimit) return strtod_fallback(token);
        mantissa = mantissa * 10 + static_cast<std::uint64_t>(token[i] - '0');
        any_digit = true;
    }
    if (i < n && token[i] == '.') {
        ++i;
        for (; i < n && is_digit(token[i]); ++i) {
            if (mantissa >= kMantissaLimit) return strtod_fallback(token);
            mantissa = mantissa * 10 + static_cast<std::uint64_t>(token[i] - '0');
            --exponent;
            any_digit = true;
        }
    }
    if (!any_digit) return 0.0;  // strtod performs no conversion

    if (i < n && (token[i] == 'e' || token[i] == 'E')) {
        std::size_t j = i + 1;
        bool exponent_negative = false;
        if (j < n && (token[j] == '+' || token[j] == '-')) {
            exponent_negative = token[j] == '-';
            ++j;
        }
        if (j < n && is_digit(token[j])) {
            int value = 0;
            for (; j < n && is_digit(token[j]); ++j) {
                if (value > 10000) return strtod_fallback(token);
                value = value * 10 + (token[j] - '0');
            }
            exponent += exponent_negative ? -value : value;
        }
    }

    // Both operands are exact, so one IEEE operation gives the correctly rounded value strtod
    // returns. Anything outside that range is handed to strtod.
    if (mantissa > kMaxExactMantissa || exponent < -kMaxExactPower || exponent > kMaxExactPower) {
        return strtod_fallback(token);
    }
    const double magnitude = static_cast<double>(mantissa);
    const double value = exponent < 0 ? magnitude / kExactPowersOfTen[-exponent]
                                      : magnitude * kExactPowersOfTen[exponent];
    return negative ? -value : value;
}

Energy narrow_to_energy(double rounded) noexcept {
    // (short)floor(x) in RNAstructure. Within the range of int this is what every target
    // produces; out-of-range values (never present in real tables) are pinned like x86 does.
    std::int32_t wide = std::numeric_limits<std::int32_t>::min();
    if (rounded >= static_cast<double>(std::numeric_limits<std::int32_t>::min()) &&
        rounded <= static_cast<double>(std::numeric_limits<std::int32_t>::max())) {
        wide = static_cast<std::int32_t>(rounded);
    }
    return static_cast<Energy>(static_cast<std::uint16_t>(static_cast<std::uint32_t>(wide)));
}

Energy entry_to_energy_exact(std::string_view token) {
    if (token == ".") return kInfiniteEnergy;
    return narrow_to_energy(std::floor(parse_double(token) * kConversionFactor + .5));
}

std::vector<std::string_view> LineCursor::remaining_graphemes() const {
    if (failed_) return {};
    return graphemes_without_space(line_.substr(pos_));
}

std::vector<std::string_view> graphemes_without_space(std::string_view text) {
    std::vector<std::string_view> all = graphemes::split(text);
    std::vector<std::string_view> out;
    out.reserve(all.size());
    for (std::string_view g : all) {
        bool space = true;
        for (char c : g) space = space && is_space(c);
        if (!space) out.push_back(g);
    }
    return out;
}

std::string where(const std::string& file, std::size_t line) {
    return file + ":" + std::to_string(line);
}

}  // namespace thermorna::rnastructure::detail
