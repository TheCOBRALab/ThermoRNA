// Dumps the tables ThermoRNA's RNAstructure loader produces, in the format of dump_format.hpp.
// Usage: thermorna_dump <folder> <file prefix, or "" for a folder with one set> [kelvin]

#include "dump_format.hpp"

#include <RNAstructure/params/alphabet_tables.hpp>

#include <cstdio>
#include <cstdlib>
#include <exception>

namespace rs = thermorna::rnastructure;

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <folder> <prefix or \"\"> [kelvin]\n", argv[0]);
        return 2;
    }
    const double kelvin = argc > 3 ? std::atof(argv[3]) : rs::kTemperature37C;

    rs::AlphabetTables data;
    try {
        data = rs::AlphabetTables::load(argv[1], argv[2], kelvin);
    } catch (const std::exception& error) {
        std::printf("ERROR %s\n", error.what());
        return 1;
    }
    for (const std::string& warning : data.warnings())
        std::fprintf(stderr, "warning: %s\n", warning.c_str());

    const int n = static_cast<int>(data.size());
    const rs::MiscLoop& misc = data.misc();
    std::printf("alphabet %s size %d\n", argv[2], n);
    std::printf("temperature %.17g\n", data.temperature());
    std::printf("RT %s\n", parity::float_bits(data.RT()).c_str());
    std::printf("prelog %s\n", parity::float_bits(misc.prelog).c_str());
    std::printf(
        "misc maxpen %d efn2a %d efn2b %d efn2c %d mlasym %d strain %d gubonus %d "
        "cslope %d cint %d c3 %d init %d singlecbulge %d\n",
        misc.maxpen, misc.efn2a, misc.efn2b, misc.efn2c, misc.mlasym, misc.strain, misc.gubonus,
        misc.cslope, misc.cint, misc.c3, misc.init, misc.singlecbulge);
    std::printf("poppen");
    for (int k = 0; k < 5; ++k) std::printf(" %d", misc.poppen[static_cast<std::size_t>(k)]);
    std::printf("\neparam");
    for (int k = 0; k < 11; ++k) std::printf(" %d", misc.eparam[static_cast<std::size_t>(k)]);
    std::printf("\n");

    auto loop = [](const char* name, const std::vector<rs::Energy>& v) {
        std::printf("%s %zu", name, v.size());
        for (rs::Energy x : v) std::printf(" %d", x);
        std::printf("\n");
    };
    loop("inter", data.inter());
    loop("bulge", data.bulge());
    loop("hairpin", data.hairpin());

    auto xloop = [n](const char* name, const std::vector<rs::SpecialHairpin>& v) {
        std::printf("%s %zu\n", name, v.size());
        for (const rs::SpecialHairpin& entry : v) {
            int key = 0;
            int factor = 1;
            for (rs::Code c : entry.codes) {
                key += c * factor;
                factor *= n;
            }
            std::printf("%s %d %d\n", name, key, entry.energy);
        }
    };
    xloop("tloop", data.tloop());
    xloop("triloop", data.triloop());
    xloop("hexaloop", data.hexaloop());

    using C = rs::Code;
    long finite = 0;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            for (int c = 0; c < n; ++c)
                for (int d = 0; d < 3; ++d) {
                    const rs::Energy v = data.dangle(C(a), C(b), C(c), static_cast<std::size_t>(d));
                    if (v == rs::kInfiniteEnergy) continue;
                    std::printf("dangle %d %d %d %d %d\n", a, b, c, d, v);
                    ++finite;
                }
    std::printf("count dangle %ld\n", finite);

    for (int t = 0; t < 11; ++t) {
        finite = 0;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                for (int c = 0; c < n; ++c)
                    for (int d = 0; d < n; ++d) {
                        const rs::Energy v =
                            data.table4(static_cast<rs::Table4D>(t), C(a), C(b), C(c), C(d));
                        if (v == rs::kInfiniteEnergy) continue;
                        std::printf("%s %d %d %d %d %d\n", parity::kTable4Names[t], a, b, c, d, v);
                        ++finite;
                    }
        std::printf("count %s %ld\n", parity::kTable4Names[t], finite);
    }

    finite = 0;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            for (int c = 0; c < n; ++c)
                for (int d = 0; d < n; ++d)
                    for (int e = 0; e < n; ++e)
                        for (int f = 0; f < n; ++f) {
                            const rs::Energy v = data.iloop11(C(a), C(b), C(c), C(d), C(e), C(f));
                            if (v == rs::kInfiniteEnergy) continue;
                            std::printf("iloop11 %d %d %d %d %d %d %d\n", a, b, c, d, e, f, v);
                            ++finite;
                        }
    std::printf("count iloop11 %ld\n", finite);

    finite = 0;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            for (int c = 0; c < n; ++c)
                for (int d = 0; d < n; ++d)
                    for (int e = 0; e < n; ++e)
                        for (int f = 0; f < n; ++f)
                            for (int g = 0; g < n; ++g) {
                                const rs::Energy v =
                                    data.iloop21(C(a), C(b), C(c), C(d), C(e), C(f), C(g));
                                if (v == rs::kInfiniteEnergy) continue;
                                std::printf("iloop21 %d %d %d %d %d %d %d %d\n", a, b, c, d, e, f,
                                            g, v);
                                ++finite;
                            }
    std::printf("count iloop21 %ld\n", finite);

    finite = 0;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            for (int c = 0; c < n; ++c)
                for (int d = 0; d < n; ++d)
                    for (int e = 0; e < n; ++e)
                        for (int f = 0; f < n; ++f)
                            for (int g = 0; g < n; ++g)
                                for (int h = 0; h < n; ++h) {
                                    const rs::Energy v = data.iloop22(C(a), C(b), C(c), C(d), C(e),
                                                                      C(f), C(g), C(h));
                                    if (v == rs::kInfiniteEnergy) continue;
                                    std::printf("iloop22 %d %d %d %d %d %d %d %d %d\n", a, b, c, d,
                                                e, f, g, h, v);
                                    ++finite;
                                }
    std::printf("count iloop22 %ld\n", finite);
    return 0;
}
