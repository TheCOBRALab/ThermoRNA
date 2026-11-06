// Dumps the tables RNAstructure's own loader produces (datatable::opendat + ScaleToTemperature,
// exactly what efn2 uses) in the format of dump_format.hpp.
//
// Build against an RNAstructure source tree (no COUNTING):
//   c++ -O2 -std=c++17 -I<RNAstructure>/src rnastructure_dump.cpp \
//       <RNAstructure>/src/rna_library.cpp <RNAstructure>/src/common_utils.cpp -o rnastructure_dump
// Usage: rnastructure_dump <data_tables dir> <alphabet> [temperature in K]

#include "dump_format.hpp"
#include "rna_library.h"

#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <data_tables> <alphabet> [kelvin]\n", argv[0]);
        return 2;
    }
    const double kelvin = argc > 3 ? std::atof(argv[3]) : 310.15;

    datatable data;
    if (data.opendat(argv[1], argv[2]) == 0) {
        std::printf("ERROR opendat\n");
        return 1;
    }
    if (std::abs(kelvin - 310.15) >= 0.01) {
        const int error = data.ScaleToTemperature(kelvin);
        if (error != 0) {
            std::printf("ERROR scale %d\n", error);
            return 1;
        }
    }

    const int n = static_cast<int>(data.alphabet.size());
    std::printf("alphabet %s size %d\n", argv[2], n);
    std::printf("temperature %.17g\n", data.temperature);
    std::printf("RT %s\n", parity::float_bits(data.RT).c_str());
    std::printf("prelog %s\n", parity::float_bits(data.prelog).c_str());
    std::printf(
        "misc maxpen %d efn2a %d efn2b %d efn2c %d mlasym %d strain %d gubonus %d "
        "cslope %d cint %d c3 %d init %d singlecbulge %d\n",
        data.maxpen, data.efn2a, data.efn2b, data.efn2c, data.mlasym, data.strain, data.gubonus,
        data.cslope, data.cint, data.c3, data.init, data.singlecbulge);
    std::printf("poppen");
    for (int k = 0; k < 5; ++k) std::printf(" %d", data.poppen[k]);
    std::printf("\neparam");
    for (int k = 0; k < 11; ++k) std::printf(" %d", data.eparam[k]);
    std::printf("\n");

    auto loop = [](const char* name, const std::vector<short>& v) {
        std::printf("%s %zu", name, v.size());
        for (short x : v) std::printf(" %d", x);
        std::printf("\n");
    };
    loop("inter", data.inter);
    loop("bulge", data.bulge);
    loop("hairpin", data.hairpin);

    auto xloop = [](const char* name, const datatable::intVector2D& v) {
        std::printf("%s %zu\n", name, v.size());
        for (const auto& entry : v) std::printf("%s %d %d\n", name, entry[0], entry[1]);
    };
    xloop("tloop", data.tloop);
    xloop("triloop", data.triloop);
    xloop("hexaloop", data.hexaloop);

    long finite = 0;
    for (int a = 0; a < n; ++a)
        for (int b = 0; b < n; ++b)
            for (int c = 0; c < n; ++c)
                for (int d = 0; d < 3; ++d)
                    if (data.dangle[a][b][c][d] != INFINITE_ENERGY) {
                        std::printf("dangle %d %d %d %d %d\n", a, b, c, d, data.dangle[a][b][c][d]);
                        ++finite;
                    }
    std::printf("count dangle %ld\n", finite);

    const datatable::shortVector4D* tables[] = {
        &data.stack, &data.helix_ends, &data.tstkh,     &data.tstki,  &data.tstki23, &data.tstki1n,
        &data.coax,  &data.tstackcoax, &data.coaxstack, &data.tstack, &data.tstkm};
    for (int t = 0; t < 11; ++t) {
        const auto& m = *tables[t];
        finite = 0;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < n; ++b)
                for (int c = 0; c < n; ++c)
                    for (int d = 0; d < n; ++d)
                        if (m[a][b][c][d] != INFINITE_ENERGY) {
                            std::printf("%s %d %d %d %d %d\n", parity::kTable4Names[t], a, b, c, d,
                                        m[a][b][c][d]);
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
                            const short v = data.iloop11[a][b][c][d][e][f];
                            if (v == INFINITE_ENERGY) continue;
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
                                const short v = data.iloop21[a][b][c][d][e][f][g];
                                if (v == INFINITE_ENERGY) continue;
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
                                    const short v = data.iloop22[a][b][c][d][e][f][g][h];
                                    if (v == INFINITE_ENERGY) continue;
                                    std::printf("iloop22 %d %d %d %d %d %d %d %d %d\n", a, b, c, d,
                                                e, f, g, h, v);
                                    ++finite;
                                }
    std::printf("count iloop22 %ld\n", finite);
    return 0;
}
