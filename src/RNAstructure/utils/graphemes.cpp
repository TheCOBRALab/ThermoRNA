#include "RNAstructure/utils/graphemes.hpp"

// uni-algo is only ever included here, so its large tables are compiled into one object file.
#include <uni_algo/ranges_grapheme.h>

namespace thermorna::rnastructure::graphemes {

bool is_ascii(std::string_view text) noexcept {
    for (char c : text) {
        if (!detail::is_ascii_byte(static_cast<unsigned char>(c))) return false;
    }
    return true;
}

void split_with_uni_algo(std::string_view text, std::vector<std::string_view>& out) {
    for (std::string_view g : una::views::grapheme::utf8(text)) {
        out.push_back(g);
    }
}

void split(std::string_view text, std::vector<std::string_view>& out) {
    for_each(text, [&out](std::string_view g) { out.push_back(g); });
}

std::vector<std::string_view> split(std::string_view text) {
    std::vector<std::string_view> out;
    out.reserve(text.size());
    split(text, out);
    return out;
}

std::size_t first_length(std::string_view text) {
    if (text.empty()) return 0;

    const auto* p = reinterpret_cast<const unsigned char*>(text.data());
    const std::size_t n = text.size();

    const std::size_t ascii_length = detail::ascii_grapheme_length(p, 0, n);
    if (ascii_length != 0) return ascii_length;

    const std::string_view span = text.substr(0, detail::unicode_span_end(p, 0, n));
    for (std::string_view g : una::views::grapheme::utf8(span)) {
        return g.size();
    }
    return span.size();
}

}  // namespace thermorna::rnastructure::graphemes
