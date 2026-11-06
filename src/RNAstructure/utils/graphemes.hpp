#pragma once

/**
 * @file RNAstructure/utils/graphemes.hpp
 * @brief Grapheme splitting that only pays for Unicode segmentation where it is needed.
 *
 * Between two ASCII characters there is always a grapheme boundary (the only exception is
 * CR LF, Unicode rule GB3). The splitter therefore emits ASCII bytes one by one and hands only
 * the short spans that contain non-ASCII bytes to uni-algo. A pure-ASCII string never reaches
 * uni-algo, and a long sequence with a handful of modified bases such as "Ψ" segments only the
 * few bytes around them. The result is identical to segmenting the whole string with uni-algo.
 *
 * Views point into the input, which must outlive them.
 */

#include <cstddef>
#include <string_view>
#include <vector>

namespace thermorna::rnastructure::graphemes {

namespace detail {

[[nodiscard]] constexpr bool is_ascii_byte(unsigned char c) noexcept { return c < 0x80; }

[[nodiscard]] constexpr bool is_crlf(const unsigned char* p, std::size_t k) noexcept {
    return p[k - 1] == '\r' && p[k] == '\n';
}

// Length of a grapheme made of the ASCII byte at i, or 0 when the next byte is not ASCII
// (a combining mark or similar may extend the cluster, so uni-algo has to decide).
[[nodiscard]] constexpr std::size_t ascii_grapheme_length(const unsigned char* p, std::size_t i,
                                                          std::size_t n) noexcept {
    if (!is_ascii_byte(p[i])) return 0;
    if (i + 1 == n) return 1;
    if (!is_ascii_byte(p[i + 1])) return 0;
    return is_crlf(p, i + 1) ? 2 : 1;
}

// End of the span starting at the boundary i that must be segmented by uni-algo: the next
// position with ASCII bytes on both sides that are not CR LF, where none of the Unicode rules
// that join characters (GB3, GB9, GB9a, GB9b, GB9c, GB11, GB12/13) can apply.
[[nodiscard]] constexpr std::size_t unicode_span_end(const unsigned char* p, std::size_t i,
                                                     std::size_t n) noexcept {
    std::size_t k = i + 1;
    while (k < n && !(is_ascii_byte(p[k - 1]) && is_ascii_byte(p[k]) && !is_crlf(p, k))) ++k;
    return k;
}

}  // namespace detail

/// True when every byte is ASCII (< 0x80).
[[nodiscard]] bool is_ascii(std::string_view text) noexcept;

/// Segments @p text with uni-algo only, without the ASCII shortcut.
void split_with_uni_algo(std::string_view text, std::vector<std::string_view>& out);

/// Calls @p f with every grapheme cluster of @p text, in order.
template <typename F>
void for_each(std::string_view text, F&& f) {
    const auto* p = reinterpret_cast<const unsigned char*>(text.data());
    const std::size_t n = text.size();
    std::vector<std::string_view> span;
    std::size_t i = 0;
    while (i < n) {
        const std::size_t ascii_length = detail::ascii_grapheme_length(p, i, n);
        if (ascii_length != 0) {
            f(text.substr(i, ascii_length));
            i += ascii_length;
            continue;
        }
        const std::size_t end = detail::unicode_span_end(p, i, n);
        span.clear();
        split_with_uni_algo(text.substr(i, end - i), span);
        for (std::string_view g : span) f(g);
        i = end;
    }
}

/// Appends the grapheme clusters of @p text to @p out.
void split(std::string_view text, std::vector<std::string_view>& out);

/// Returns the grapheme clusters of @p text.
[[nodiscard]] std::vector<std::string_view> split(std::string_view text);

/// Byte length of the first grapheme cluster of @p text (0 for an empty string).
[[nodiscard]] std::size_t first_length(std::string_view text);

}  // namespace thermorna::rnastructure::graphemes
