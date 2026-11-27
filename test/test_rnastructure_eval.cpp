// Unit tests of ThermoRNA's RNAstructure energy functions (RNAstructure/eval).
//
// The expected values are what RNAstructure 6.6's `efn2 -w` (thermodynamic details) prints for
// each structure, in the units of the tables: tenths of kcal/mol. The tests use the parameter
// sets shipped in params/RNAstructure, so they need no setup.

#include <RNAstructure/rnastructure.hpp>
#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace rs = thermorna::rnastructure;

namespace {

// ThermoRNA's params/RNAstructure (set by CMake).
const std::string kParams = THERMORNA_TEST_PARAMS;

// A shipped set, by folder name or by RNAstructure's name ("rna", "PSU", "1"), with modified
// bases from other sets merged in.
rs::Datatable load(const std::string& set, const std::vector<std::string>& modifications = {},
                   double celsius = 37.0) {
    rs::LoadOptions options;
    options.data_directory = rs::find_set_folder(kParams, set).directory;
    for (const std::string& name : modifications) {
        options.modifications.push_back({rs::find_set_folder(kParams, name).directory});
    }
    options.temperature = celsius + rs::kCelsiusToKelvin;
    return rs::Datatable::load(options);
}

// The position each position of a dot-bracket structure is paired with (itself if unpaired).
std::vector<std::size_t> pair_table(const std::string& structure) {
    std::vector<std::size_t> pair(structure.size());
    std::vector<std::size_t> open;
    for (std::size_t k = 0; k < structure.size(); ++k) {
        pair[k] = k;
        if (structure[k] == '(') open.push_back(k);
        if (structure[k] == ')') {
            pair[k] = open.back();
            pair[open.back()] = k;
            open.pop_back();
        }
    }
    return pair;
}

// The functions take encoded nucleotides, not a sequence. For a structure of one sequence,
// these gather what each of them needs from the encoded sequence s.

// The nucleotide a number of positions before or after position k, if there is one.
rs::Code before(const std::vector<rs::Code>& s, std::size_t k, std::size_t steps = 1) {
    return k >= steps ? s[k - steps] : rs::kNoNucleotide;
}

rs::Code after(const std::vector<rs::Code>& s, std::size_t k) {
    return k + 1 < s.size() ? s[k + 1] : rs::kNoNucleotide;
}

// The hairpin loop closed by i-j.
int hairpin(const std::vector<rs::Code>& s, std::size_t i, std::size_t j, const rs::Datatable& d) {
    return rs::hairpin_energy(j - i - 1, &s[i], before(s, i), before(s, i, 2), after(s, j), d);
}

// How many identical nucleotides in a row the one at position k is part of.
std::size_t run_at(const std::vector<rs::Code>& s, std::size_t k) {
    std::size_t run = 1;
    for (std::size_t x = k; x > 0 && s[x - 1] == s[k]; --x) ++run;
    for (std::size_t x = k + 1; x < s.size() && s[x] == s[k]; ++x) ++run;
    return run;
}

// The internal loop, bulge or stack between i-j and p-q.
int internal(const std::vector<rs::Code>& s, std::size_t i, std::size_t j, std::size_t p,
             std::size_t q, const rs::Datatable& d) {
    const std::size_t n1 = p - i - 1;
    const std::size_t n2 = j - q - 1;
    const std::size_t run = n1 + n2 == 1 ? run_at(s, n1 == 1 ? i + 1 : j - 1) : 1;
    return rs::internal_loop_energy(n1, n2, s[i], s[j], s[p], s[q], s[i + 1], s[j - 1], s[p - 1],
                                    s[q + 1], before(s, i), after(s, j), s[p + 1], s[q - 1], run,
                                    d);
}

// The elements of the loop between two positions: its unpaired nucleotides and its helices.
std::vector<rs::LoopElement> elements(const std::vector<rs::Code>& s,
                                      const std::vector<std::size_t>& pair, std::size_t from,
                                      std::size_t to) {
    std::vector<rs::LoopElement> loop;
    for (std::size_t k = from; k < to; ++k) {
        if (pair[k] > k) {
            const std::size_t q = pair[k];
            loop.push_back(rs::LoopElement::helix(s[k], s[q], s[k + 1], s[q - 1]));
            k = q;
        } else {
            loop.push_back(rs::LoopElement::unpaired(s[k]));
        }
    }
    return loop;
}

// efn2 adds forbidden table entries up like numbers (it prints 1400.5, 2799.1, ...). The eval
// functions return kInfiniteEnergy for such a loop, and efn2() below for its structure.
constexpr int kInf = rs::kInfiniteEnergy;

// What efn2 -w prints for a structure: "Exterior loop", then the energy of the loop each pair
// closes ("Stack", "Hairpin", "Internal/bulge" or "Multiloop") in 5' to 3' order of the pairs,
// then "Total energy". The sequence is encoded once. With MultibranchRule::Simple: efn2 -w -s.
// With CoaxialStacking::Off the exterior and multibranch loops are no longer efn2's.
std::vector<int> efn2(const std::string& sequence, const std::string& structure,
                      const rs::Datatable& d, rs::MultibranchRule rule = rs::MultibranchRule::Efn2,
                      rs::CoaxialStacking coaxial = rs::CoaxialStacking::On) {
    const std::vector<rs::Code> s = d.encode(sequence);
    const std::vector<std::size_t> pair = pair_table(structure);

    std::vector<int> energies{rs::exterior_energy(elements(s, pair, 0, s.size()), d, coaxial)};
    for (std::size_t i = 0; i < s.size(); ++i) {
        const std::size_t j = pair[i];
        if (j <= i) continue;
        const std::vector<rs::LoopElement> inside = elements(s, pair, i + 1, j);
        const auto is_helix = [](const rs::LoopElement& e) {
            return e.type == rs::LoopElement::Type::Helix;
        };
        const auto helices = std::count_if(inside.begin(), inside.end(), is_helix);
        if (helices == 0) {
            energies.push_back(hairpin(s, i, j, d));
        } else if (helices == 1) {
            const std::size_t p =
                i + 1 +
                static_cast<std::size_t>(std::find_if(inside.begin(), inside.end(), is_helix) -
                                         inside.begin());
            energies.push_back(internal(s, i, j, p, pair[p], d));
        } else {
            energies.push_back(rs::multibranch_energy(s[i], s[j], before(s, i), after(s, j), inside,
                                                      d, rule, coaxial));
        }
    }

    int total = 0;
    for (const int energy : energies) {
        if (energy == kInf) total = kInf;
        if (total != kInf) total += energy;
    }
    energies.push_back(total);
    return energies;
}

struct Efn2Case {
    const char* what;
    std::string sequence;
    std::string structure;
    std::vector<int> energies;  // see efn2()
};

void expect_efn2(const std::vector<Efn2Case>& cases, const rs::Datatable& d,
                 rs::MultibranchRule rule = rs::MultibranchRule::Efn2) {
    for (const Efn2Case& c : cases) {
        EXPECT_EQ(efn2(c.sequence, c.structure, d, rule), c.energies)
            << c.what << ": " << c.sequence << " " << c.structure;
    }
}

// An exterior or multibranch loop put together from the stem-level functions alone, the way
// the recursion of a folding algorithm does it: every helix gets its stem energy with one way
// of stacking, or stacks coaxially on the next helix, and the lowest sum is kept. `loop` is
// the loop in the order it is gone through. A multibranch loop is a circle, opened somewhere:
// the nucleotides next to a helix then wrap around, but what stacks together does not.
// Without `coaxial` the helices do not stack on each other, as with CoaxialStacking::Off.
class StemSum {
   public:
    StemSum(const std::vector<rs::LoopElement>& loop, bool multibranch, const rs::Datatable& d,
            bool coaxial)
        : loop_(loop),
          multibranch_(multibranch),
          d_(d),
          coaxial_(coaxial),
          memo_(2 * loop.size()) {}

    // The lowest energy of the elements from k on. `free`: the element before k is an unpaired
    // nucleotide that nothing has stacked on.
    int lowest(std::size_t k = 0, bool free = false) {
        if (k >= loop_.size()) return 0;
        std::optional<int>& known = memo_[2 * k + (free ? 1 : 0)];
        if (known) return *known;

        int best = 0;
        if (unpaired(k)) {
            best = lowest(k + 1, true);
        } else if (!helix(k)) {  // a branch that is not a helix: only the penalty per helix
            best = (multibranch_ ? int{d_.misc().efn2c} : 0) + lowest(k + 1, false);
        } else {
            const rs::LoopElement& h = loop_[k];
            const int alone = stem(k, rs::Stacking::None);
            best = alone + lowest(k + 1, false);
            if (free) {
                best = std::min(best, stem(k, rs::Stacking::Dangle5) + lowest(k + 1, false));
            }
            if (unpaired(k + 1)) {
                best = std::min(best, stem(k, rs::Stacking::Dangle3) + lowest(k + 2, false));
            }
            if (free && unpaired(k + 1)) {
                best = std::min(best, stem(k, rs::Stacking::Mismatch) + lowest(k + 2, false));
            }
            // Coaxial stacks on the next helix: both helices then have nothing else on them.
            if (coaxial_ && helix(k + 1)) {
                const rs::LoopElement& next = loop_[k + 1];
                best = std::min(best, alone + stem(k + 1, rs::Stacking::None) +
                                          rs::coaxial_stack_energy(h.reached, h.left, next.reached,
                                                                   next.left, d_) +
                                          lowest(k + 2, false));
            }
            if (coaxial_ && unpaired(k + 1) && helix(k + 2)) {
                const rs::LoopElement& next = loop_[k + 2];
                const int both = alone + stem(k + 2, rs::Stacking::None);
                if (unpaired(k + 3)) {
                    best = std::min(best, both +
                                              rs::coaxial_stack_mismatch3_energy(
                                                  h.reached, h.left, next.reached, next.left,
                                                  loop_[k + 1].left, loop_[k + 3].left, d_) +
                                              lowest(k + 4, false));
                }
                if (free) {
                    best = std::min(best, both +
                                              rs::coaxial_stack_mismatch5_energy(
                                                  h.reached, h.left, next.reached, next.left,
                                                  loop_[k - 1].left, loop_[k + 1].left, d_) +
                                              lowest(k + 3, false));
                }
            }
        }
        known = best;
        return best;
    }

   private:
    bool helix(std::size_t k) const {
        return k < loop_.size() && loop_[k].type == rs::LoopElement::Type::Helix;
    }
    bool unpaired(std::size_t k) const {
        return k < loop_.size() && loop_[k].type == rs::LoopElement::Type::Unpaired;
    }
    // The stem energy of the helix at k, with the nucleotides before and after it in the loop
    int stem(std::size_t k, rs::Stacking stacking) const {
        const rs::LoopElement& h = loop_[k];
        const rs::Code none = rs::kNoNucleotide;
        const rs::Code n5 = k > 0 ? loop_[k - 1].left : (multibranch_ ? loop_.back().left : none);
        const rs::Code n3 = k + 1 < loop_.size() ? loop_[k + 1].reached
                                                 : (multibranch_ ? loop_.front().reached : none);
        return multibranch_ ? rs::multibranch_stem_energy(h.reached, h.left, h.after_reached,
                                                          h.before_left, n5, n3, stacking, d_)
                            : rs::exterior_stem_energy(h.reached, h.left, h.after_reached,
                                                       h.before_left, n5, n3, stacking, d_);
    }

    const std::vector<rs::LoopElement>& loop_;
    bool multibranch_;
    const rs::Datatable& d_;
    bool coaxial_;
    std::vector<std::optional<int>> memo_;
};

// A sum with a forbidden term in it, as the eval functions return it
int as_infinite(int energy) { return energy >= kInf / 2 ? kInf : energy; }

// The exterior loop as a sum of stem energies
int exterior_from_stems(const std::vector<rs::LoopElement>& loop, const rs::Datatable& d,
                        bool coaxial) {
    return as_infinite(StemSum(loop, false, d, coaxial).lowest());
}

// The multibranch loop closed by i-j, by the simple rule, as a sum of stem energies. `loop` is
// what is between i and j. Seen from inside, the closing pair is a helix that the loop reaches
// at j and leaves at i, and whose helix side is outside the pair.
int multibranch_from_stems(rs::Code i, rs::Code j, rs::Code beyond_i, rs::Code beyond_j,
                           const std::vector<rs::LoopElement>& loop, const rs::Datatable& d,
                           bool coaxial) {
    std::vector<rs::LoopElement> circle{rs::LoopElement::helix(j, i, beyond_j, beyond_i)};
    circle.insert(circle.end(), loop.begin(), loop.end());
    const auto unpaired = std::count_if(loop.begin(), loop.end(), [](const rs::LoopElement& e) {
        return e.type == rs::LoopElement::Type::Unpaired;
    });
    // Open the circle in every place: one of them does not cut through what stacks together.
    int lowest = kInf;
    for (std::size_t opening = 0; opening < circle.size(); ++opening) {
        lowest = std::min(lowest, as_infinite(StemSum(circle, true, d, coaxial).lowest()));
        std::rotate(circle.begin(), circle.begin() + 1, circle.end());
    }
    return as_infinite(d.misc().efn2a + d.misc().efn2b * static_cast<int>(unpaired) + lowest);
}

// Every exterior and multibranch loop of the structures: the sum of stem energies is what
// exterior_energy and multibranch_energy with MultibranchRule::Simple return, with the coaxial
// stacks and without them.
void expect_stem_sums(const std::vector<Efn2Case>& cases, const rs::Datatable& d) {
    for (const bool coaxial : {true, false}) {
        const rs::CoaxialStacking setting =
            coaxial ? rs::CoaxialStacking::On : rs::CoaxialStacking::Off;
        for (const Efn2Case& c : cases) {
            const std::string where = std::string(c.what) + ": " + c.sequence + " " + c.structure +
                                      (coaxial ? "" : ", without coaxial stacking");
            const std::vector<rs::Code> s = d.encode(c.sequence);
            const std::vector<std::size_t> pair = pair_table(c.structure);
            const std::vector<rs::LoopElement> outside = elements(s, pair, 0, s.size());
            EXPECT_EQ(exterior_from_stems(outside, d, coaxial),
                      rs::exterior_energy(outside, d, setting))
                << where;
            for (std::size_t i = 0; i < s.size(); ++i) {
                const std::size_t j = pair[i];
                if (j <= i) continue;
                const std::vector<rs::LoopElement> inside = elements(s, pair, i + 1, j);
                const auto helices =
                    std::count_if(inside.begin(), inside.end(), [](const rs::LoopElement& e) {
                        return e.type == rs::LoopElement::Type::Helix;
                    });
                if (helices < 2) continue;
                EXPECT_EQ(multibranch_from_stems(s[i], s[j], before(s, i), after(s, j), inside, d,
                                                 coaxial),
                          rs::multibranch_energy(s[i], s[j], before(s, i), after(s, j), inside, d,
                                                 rs::MultibranchRule::Simple, setting))
                    << where << ", loop closed at " << i;
            }
        }
    }
}

// The total energy RNAstructure's own folding recursion (the dynamic programming behind Fold)
// gives a structure when its coaxial stacking is disabled, as with Fold --disablecoax, and the
// structure is the only one it is allowed to find.
struct FoldCase {
    const char* what;  // the Efn2Case with the structure
    int energy;
};

// `expected` has the structures of `cases` in their order, without those RNAstructure's
// recursion has no energy for. With MultibranchRule::Simple and CoaxialStacking::Off the
// functions give the same total.
void expect_fold_without_coaxial_stacking(const std::vector<Efn2Case>& cases,
                                          const std::vector<FoldCase>& expected,
                                          const rs::Datatable& d) {
    auto next = expected.begin();
    for (const Efn2Case& c : cases) {
        if (next == expected.end() || std::string(c.what) != next->what) continue;
        EXPECT_EQ(
            efn2(c.sequence, c.structure, d, rs::MultibranchRule::Simple, rs::CoaxialStacking::Off)
                .back(),
            next->energy)
            << c.what << ": " << c.sequence << " " << c.structure;
        ++next;
    }
    EXPECT_EQ(next, expected.end()) << "structures of the table are missing or out of order";
}

// ------------------------------------------------------------------------------------------
// efn2 -w on hand-picked structures (generated: do not edit the numbers by hand)
// ------------------------------------------------------------------------------------------

const std::vector<Efn2Case> kHairpins = {
    {"initiation and first mismatch", "GGACAAAAGUCC", "((((....))))", {0, -33, -24, -22, 41, -38}},
    {"mismatch with a bonus", "GGACGAAAGUCC", "((((....))))", {0, -33, -24, -22, 33, -46}},
    {"five nucleotides", "GGACAGUAAGUCC", "((((.....))))", {0, -33, -24, -22, 42, -37}},
    {"six nucleotides", "GGACACAGUGGUCC", "((((......))))", {0, -33, -24, -22, 40, -39}},
    {"three nucleotides: no mismatch", "GGACGAAGUCC", "((((...))))", {0, -33, -24, -22, 54, -25}},
    {"special tetraloop", "GGACUUCGGUCC", "((((....))))", {0, -33, -24, -22, 37, -42}},
    {"special triloop", "GGACAACGUCC", "((((...))))", {0, -33, -24, -22, 68, -11}},
    {"special hexaloop", "GGCACAGUGCUGCC", "((((......))))", {0, -33, -34, -21, 29, -59}},
    {"all-C triloop", "GGAGCCCCUCC", "((((...))))", {0, -33, -24, -21, 69, -9}},
    {"all-C loop", "GGAGCCCCCCUCC", "((((.....))))", {0, -33, -24, -21, 81, 3}},
    {"GGG-U closure", "GGGAAAAUCC", "(((....)))", {0, -33, -15, 36, -12}},
    {"GGG-U closure after an unpaired A",
     "AGGGACAAGUCCU",
     ".(((.....))).",
     {-12, -33, -15, 32, -28}},
    {"35 nucleotides",
     "GGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGUCC",
     "((((...................................))))",
     {0, -33, -24, -22, 64, -15}},
    {"closing pair at both sequence ends", "CAAAAG", "(....)", {0, 41, 41}},
};

const std::vector<Efn2Case> kInternalLoops = {
    {"1x1",
     "GGACAGCCAAAAGGCAGUCC",
     "((((.(((....))).))))",
     {0, -33, -24, -22, 9, -34, -33, 41, -96}},
    {"1x1 between A-U pairs",
     "GGAUAAGCAAAAGCUAAUCC",
     "((((.(((....))).))))",
     {0, -33, -24, -11, 19, -21, -34, 41, -63}},
    {"1x2",
     "GGACAGCCAAAAGGCAAGUCC",
     "((((.(((....)))..))))",
     {0, -33, -24, -22, 23, -34, -33, 41, -82}},
    {"2x1",
     "GGACAAGCCAAAAGGCAGUCC",
     "((((..(((....))).))))",
     {0, -33, -24, -22, 23, -34, -33, 41, -82}},
    {"2x2",
     "GGACAAGCCAAAAGGCAAGUCC",
     "((((..(((....)))..))))",
     {0, -33, -24, -22, 12, -34, -33, 41, -93}},
    {"2x2 with tandem G-A mismatches",
     "GGACGAGCCAAAAGGCGAGUCC",
     "((((..(((....)))..))))",
     {0, -33, -24, -22, -7, -34, -33, 41, -112}},
    {"1x3",
     "GGACAGCCAAAAGGCAAAGUCC",
     "((((.(((....)))...))))",
     {0, -33, -24, -22, 23, -34, -33, 41, -82}},
    {"2x3",
     "GGACAAGCCAAAAGGCAAAGUCC",
     "((((..(((....)))...))))",
     {0, -33, -24, -22, 26, -34, -33, 41, -79}},
    {"3x3",
     "GGACAAAGCCAAAAGGCAAAGUCC",
     "((((...(((....)))...))))",
     {0, -33, -24, -22, 20, -34, -33, 41, -85}},
    {"2x5",
     "GGACAAGCCAAAAGGCAAAAAGUCC",
     "((((..(((....))).....))))",
     {0, -33, -24, -22, 39, -34, -33, 41, -66}},
    {"4x4 with a G-G mismatch",
     "GGACGAAAGCCAAAAGGCAAAGGUCC",
     "((((....(((....)))....))))",
     {0, -33, -24, -22, 13, -34, -33, 41, -92}},
    {"16x16",
     "GGACAAAAAAAAAAAAAAAAGCCAAAAGGCAAAAAAAAAAAAAAAAGUCC",
     "((((................(((....)))................))))",
     {0, -33, -24, -22, 38, -34, -33, 41, -67}},
    {"1x32",
     "GGACAGCCAAAAGGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGUCC",
     "((((.(((....)))................................))))",
     {0, -33, -24, -22, 68, -34, -33, 41, -37}},
    {"closing pair at both sequence ends",
     "CAGCCAAAAGGCAG",
     "(.(((....))).)",
     {0, 9, -34, -33, 41, -17}},
    {"2x3 with N",
     "GGACANGCCAAAAGGCAAAGUCC",
     "((((..(((....)))...))))",
     {0, -33, -24, -22, 26, -34, -33, 41, -79}},
};

const std::vector<Efn2Case> kBulges = {
    {"one nucleotide",
     "GGACAGCCAAAAGGCGUCC",
     "((((.(((....)))))))",
     {0, -33, -24, -22, 14, -34, -33, 41, -91}},
    {"one nucleotide on the 3' side",
     "GGACGCCAAAAGGCAGUCC",
     "(((((((....))).))))",
     {0, -33, -24, -22, 14, -34, -33, 41, -91}},
    {"G before a G: 2 positions",
     "GGACGGCCAAAAGGCGUCC",
     "((((.(((....)))))))",
     {0, -33, -24, -22, 10, -34, -33, 41, -95}},
    {"C after a C: 2 positions, C bonus",
     "GGACCGCCAAAAGGCGUCC",
     "((((.(((....)))))))",
     {0, -33, -24, -22, 1, -34, -33, 41, -104}},
    {"C between Gs: no bonus",
     "GGAGCGCCAAAAGGCCUCC",
     "((((.(((....)))))))",
     {0, -33, -24, -21, 5, -34, -33, 41, -99}},
    {"Gs up to the 5' end: 3 positions",
     "GGGCCAAAAGGCC",
     "(.(((....))))",
     {0, -2, -34, -33, 41, -28}},
    {"Cs up to the 3' end: 3 positions",
     "GGCCAAAAGGCCC",
     "((((....))).)",
     {0, -11, -34, -33, 41, -37}},
    {"N", "GGACNGCCAAAAGGCGUCC", "((((.(((....)))))))", {0, -33, -24, -22, 14, -34, -33, 41, -91}},
    {"two nucleotides",
     "GGACAAGCCAAAAGGCGUCC",
     "((((..(((....)))))))",
     {0, -33, -24, -22, 28, -34, -33, 41, -77}},
    {"three nucleotides",
     "GGACAAAGCCAAAAGGCGUCC",
     "((((...(((....)))))))",
     {0, -33, -24, -22, 32, -34, -33, 41, -73}},
    {"35 nucleotides",
     "GGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCCAAAAGGCGUCC",
     "((((...................................(((....)))))))",
     {0, -33, -24, -22, 63, -34, -33, 41, -42}},
};

const std::vector<Efn2Case> kExteriorLoops = {
    {"no helix", "AAAAAA", "......", {0, 0}},
    {"G-C helix", "GGACAAAAGUCC", "((((....))))", {0, -33, -24, -22, 41, -38}},
    {"A-U helix", "AGACAAAAGUCU", "((((....))))", {5, -21, -24, -22, 41, -21}},
    {"G-U helix", "GGACAAAAGUCU", "((((....))))", {5, -21, -24, -22, 41, -21}},
    {"3' dangle", "GGACAAAAGUCCG", "((((....)))).", {-17, -33, -24, -22, 41, -55}},
    {"5' dangle", "AGGACAAAAGUCC", ".((((....))))", {-2, -33, -24, -22, 41, -40}},
    {"mismatch not better than the dangle",
     "AGGACAAAAGUCCG",
     ".((((....)))).",
     {-17, -33, -24, -22, 41, -55}},
    {"terminal mismatch", "GCGACAAAAGUCGA", ".((((....)))).", {-13, -24, -24, -22, 41, -42}},
    {"terminal mismatch among more nucleotides",
     "AAGCGACAAAAGUCGAAA",
     "...((((....))))...",
     {-13, -24, -24, -22, 41, -42}},
    {"flush coaxial stack",
     "GGACAAAAGUCCGGACAAAAGUCC",
     "((((....))))((((....))))",
     {-24, -33, -24, -22, 41, -33, -24, -22, 41, -100}},
    {"flush coaxial stack of A-U ends",
     "AGACAAAAGUCUAGACAAAAGUCU",
     "((((....))))((((....))))",
     {-3, -21, -24, -22, 41, -21, -24, -22, 41, -55}},
    {"flush coaxial stack with flanks",
     "AGGACAAAAGUCCGGACAAAAGUCCA",
     ".((((....))))((((....)))).",
     {-24, -33, -24, -22, 41, -33, -24, -22, 41, -100}},
    {"coaxial stack across one nucleotide, mismatch after",
     "GGACAAAAGUCCAGGACAAAAGUCCA",
     "((((....)))).((((....)))).",
     {-37, -33, -24, -22, 41, -33, -24, -22, 41, -113}},
    {"coaxial stack across one nucleotide, mismatch before",
     "AGGACAAAAGUCCAGGACAAAAGUCC",
     ".((((....)))).((((....))))",
     {-37, -33, -24, -22, 41, -33, -24, -22, 41, -113}},
    {"one nucleotide between helices, no other",
     "GGACAAAAGUCCAGGACAAAAGUCC",
     "((((....)))).((((....))))",
     {-17, -33, -24, -22, 41, -33, -24, -22, 41, -93}},
    {"two nucleotides between helices",
     "GGACAAAAGUCCAAGGACAAAAGUCC",
     "((((....))))..((((....))))",
     {-19, -33, -24, -22, 41, -33, -24, -22, 41, -95}},
    {"three nucleotides between helices",
     "GGACAAAAGUCCAGAGGACAAAAGUCC",
     "((((....))))...((((....))))",
     {-19, -33, -24, -22, 41, -33, -24, -22, 41, -95}},
    {"three helices in a row",
     "GGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCC",
     "((((....))))((((....))))((((....))))",
     {-24, -33, -24, -22, 41, -33, -24, -22, 41, -33, -24, -22, 41, -138}},
    {"three helices, both kinds of coaxial stack",
     "AGGACAAAAGUCCAGGACAAAAGUCCGGACAAAAGUCCA",
     ".((((....)))).((((....))))((((....)))).",
     {-54, -33, -24, -22, 41, -33, -24, -22, 41, -33, -24, -22, 41, -168}},
    {"four different helices",
     "CAGACAAAAGUCUGGCGGAAACGCCGACAAAAGUCGAAGGACAAAAGUCUU",
     ".((((....)))).(((....)))((((....))))..((((....)))).",
     {-41, -21, -24, -22, 41, -34, -24, 31, -24, -24, -22, 41, -21, -24, -22, 41, -149}},
    {"helix of one pair", "GCAAAAGCA", ".(....)..", {-11, 41, 30}},
    {"N next to the helix", "NGGACAAAAGUCCN", ".((((....)))).", {0, -33, -24, -22, 41, -38}},
    {"N inside a helix of one pair", "GCNAAAGCA", ".(....)..", {kInf, 56, kInf}},
    {"pair that cannot form", "GGACAAAAGUCAA", "((((....)))).", {kInf, kInf, -24, -22, 41, kInf}},
    {"two hairpins",
     "AGGACAAAAGUCCAGCGGAAACGCG",
     ".((((....)))).(((....))).",
     {-37, -33, -24, -22, 41, -34, -24, 31, -102}},
    {"branch with a bulge and an internal loop",
     "CGCGAAGGACGUUCGCGAGCCUCGCAGGACAAAAGUCCAA",
     ".((((.((.((....))..)))))).((((....))))..",
     {-37, -34, -24, -24, 13, -33, 23, -24, 45, -33, -24, -22, 41, -133}},
    {"three branches",
     "GCGGAAACGCGCGAAGGACGUUCGCGAGCCUCGCUAGACAAAAGUCUC",
     "(((....)))((((.((.((....))..)))))).((((....)))).",
     {-24, -34, -24, 31, -34, -24, -24, 13, -33, 23, -24, 45, -21, -24, -22, 41, -135}},
};

const std::vector<Efn2Case> kMultibranchLoops = {
    {"three helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))))",
     {0, -34, 82, -33, -24, -22, 41, -33, -24, -22, 41, -28}},
    {"three helices, one nucleotide",
     "GCGGACAAAAGUCCAGGACAAAAGUCCGC",
     "((((((....)))).((((....))))))",
     {0, -34, 71, -33, -24, -22, 41, -33, -24, -22, 41, -39}},
    {"three helices, two nucleotides",
     "GCGGACAAAAGUCCAAGGACAAAAGUCCGC",
     "((((((....))))..((((....))))))",
     {0, -34, 46, -33, -24, -22, 41, -33, -24, -22, 41, -64}},
    {"three helices, a nucleotide in every gap",
     "GCAGGACAAAAGUCCAGGACAAAAGUCCAGC",
     "((.((((....)))).((((....)))).))",
     {0, -34, 21, -33, -24, -22, 41, -33, -24, -22, 41, -89}},
    {"four helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))((((....))))))",
     {0, -34, 21, -33, -24, -22, 41, -33, -24, -22, 41, -33, -24, -22, 41, -127}},
    {"four helices, one nucleotide",
     "GCGGACAAAAGUCCGGACAAAAGUCCAGGACAAAAGUCCGC",
     "((((((....))))((((....)))).((((....))))))",
     {0, -34, 26, -33, -24, -22, 41, -33, -24, -22, 41, -33, -24, -22, 41, -122}},
    {"five helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))((((....))))((((....))))))",
     {0,  -34, 46,  -33, -24, -22, 41,  -33, -24, -22,
      41, -33, -24, -22, 41,  -33, -24, -22, 41,  -140}},
    {"five helices, one nucleotide",
     "GCAGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((.((((....))))((((....))))((((....))))((((....))))))",
     {0,  -34, 33,  -33, -24, -22, 41,  -33, -24, -22,
      41, -33, -24, -22, 41,  -33, -24, -22, 41,  -153}},
    {"all the nucleotides in one gap",
     "GCAAAAAGGACAAAAGUCCGGACAAAAGUCCGC",
     "((.....((((....))))((((....))))))",
     {0, -34, 52, -33, -24, -22, 41, -33, -24, -22, 41, -58}},
    {"gaps of 1, 0 and 1",
     "GCAGGACAAAAGUCCGGACAAAAGUCCAGC",
     "((.((((....))))((((....)))).))",
     {0, -34, 40, -33, -24, -22, 41, -33, -24, -22, 41, -70}},
    {"gaps of 2, 1 and 0",
     "GCAGGGACAAAAGUCCAGGACAAAAGUCCGC",
     "((..((((....)))).((((....))))))",
     {0, -34, 34, -33, -24, -22, 41, -33, -24, -22, 41, -76}},
    {"eight nucleotides",
     "GCAAAGGACAAAAGUCCAAGGACAAAAGUCCAAAGC",
     "((...((((....))))..((((....))))...))",
     {0, -34, 30, -33, -24, -22, 41, -33, -24, -22, 41, -80}},
    {"nine nucleotides",
     "GCAAAGGACAAAAGUCCAAAGGACAAAAGUCCAAAGC",
     "((...((((....))))...((((....))))...))",
     {0, -34, 25, -33, -24, -22, 41, -33, -24, -22, 41, -85}},
    {"twenty nucleotides",
     "GCAAAAAAAGGACAAAAGUCCAAAAAAGGACAAAAGUCCAAAAAAAGC",
     "((.......((((....))))......((((....)))).......))",
     {0, -34, 40, -33, -24, -22, 41, -33, -24, -22, 41, -70}},
    {"A-U closing pair",
     "GAGGACAAAAGUCCAGGACAAAAGUCCUC",
     "((((((....)))).((((....))))))",
     {0, -24, 79, -33, -24, -22, 41, -33, -24, -22, 41, -21}},
    {"U-G closing pair",
     "GUAGGACAAAAGUCCAGGACAAAAGUCCAGC",
     "((.((((....)))).((((....)))).))",
     {0, -25, 26, -33, -24, -22, 41, -33, -24, -22, 41, -75}},
    {"branches that end in A-U and G-U",
     "GCAGACAAAAGUCUAGGACAAAAGUCUAGC",
     "((((((....)))).((((....)))).))",
     {0, -34, 59, -21, -24, -22, 41, -21, -24, -22, 41, -27}},
    {"G and A around the branches",
     "GCGCGACAAAAGUCGAGCGACAAAAGUCGAGC",
     "((.((((....))))..((((....)))).))",
     {0, -34, 33, -24, -24, -22, 41, -24, -24, -22, 41, -59}},
    {"a branch of one pair",
     "GCACAAAAGAGGACAAAAGUCCGC",
     "((.(....).((((....))))))",
     {0, -34, 46, 41, -33, -24, -22, 41, 15}},
    {"five different helices",
     "GUCAGACAAAAGUCUGGCGGAAACGCCGACAAAAGUCGAAGGACAAAAGUCUUAC",
     "((.((((....)))).(((....)))((((....))))..((((....)))).))",
     {0, -22, 9, -21, -24, -22, 41, -34, -24, 31, -24, -24, -22, 41, -21, -24, -22, 41, -121}},
    {"a branch with loops of its own",
     "GCAGCGAAGGACGUUCGCGAGCCUCGCGCGGAAACGCGGC",
     "((.((((.((.((....))..))))))(((....))).))",
     {0, -34, 40, -34, -24, -24, 13, -33, 23, -24, 45, -34, -24, 31, -79}},
    {"N in every gap",
     "GCNGGACAAAAGUCCNGGACAAAAGUCCNGC",
     "((.((((....)))).((((....)))).))",
     {0, -34, 75, -33, -24, -22, 41, -33, -24, -22, 41, -35}},
    {"pair that cannot form",
     "GCGGACAAAAGUCAAGGACAAAAGUCCGC",
     "((((((....)))).((((....))))))",
     {0, -34, kInf, kInf, -24, -22, 41, -33, -24, -22, 41, kInf}},
    {"N inside a branch of one pair",
     "GCACNAAAGAGGACAAAAGUCCGC",
     "((.(....).((((....))))))",
     {0, -34, kInf, 56, -33, -24, -22, 41, kInf}},
    {"N beyond a closing helix of one pair",
     "NCGGACAAAAGUCCAGGACAAAAGUCCGA",
     ".(((((....)))).((((....))))).",
     {-11, kInf, -33, -24, -22, 41, -33, -24, -22, 41, kInf}},
    {"loop inside a loop",
     "GCAAGGACAAAAGUCCGCGGAAACGCAUAGGACAAAAGUCCGC",
     "((.(((((....))))(((....))).).((((....))))))",
     {0, -34, 52, 79, -33, -24, -22, 41, -34, -24, 31, -33, -24, -22, 41, -6}},
    {"one multibranch loop",
     "AGCGGACAAAAGUCCAGGACAAAAGUCCGCG",
     ".((((((....)))).((((....)))))).",
     {-17, -34, 71, -33, -24, -22, 41, -33, -24, -22, 41, -56}},
    {"three-way junction, nothing unpaired",
     "GACGCGGAAACGCCGACAAAAGUCGGUC",
     "((((((....)))((((....)))))))",
     {0, -24, -22, 73, -34, -24, 31, -24, -24, -22, 41, -29}},
    {"a loop inside a loop",
     "GGCAAGGACAAAAGUCCGCGGAAACGCAUAGGACAAAAGUCCGCAGCGGAAACGC",
     ".((.(((((....))))(((....))).).((((....)))))).(((....)))",
     {-36, -34, 52, 79, -33, -24, -22, 41, -34, -24, 31, -33, -24, -22, 41, -34, -24, 31, -69}},
    {"two multibranch loops and a branch with loops",
     "GCGGACAAAAGUCCAAGCGGAAACGCAGCUCUGCGAAGGACGUUCGCGAGCCUCGCAAGACAAAAGUCUCAAAAGAGC",
     "((((((....))))..(((....))).)).((((((.((.((....))..)))))).((((....))))(....))).",
     {-27, -34, 33, -33, -24, -22, 41, -34, -24, 31,  -21, 39, -34,
      -24, -24, 13, -33, 23,  -24, 45, -21, -24, -22, 41,  41, -118}},
};

const std::vector<Efn2Case> kMultibranchLoopsSimple = {
    {"three helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))))",
     {0, -34, 51, -33, -24, -22, 41, -33, -24, -22, 41, -59}},
    {"three helices, one nucleotide",
     "GCGGACAAAAGUCCAGGACAAAAGUCCGC",
     "((((((....)))).((((....))))))",
     {0, -34, 34, -33, -24, -22, 41, -33, -24, -22, 41, -76}},
    {"three helices, two nucleotides",
     "GCGGACAAAAGUCCAAGGACAAAAGUCCGC",
     "((((((....))))..((((....))))))",
     {0, -34, 34, -33, -24, -22, 41, -33, -24, -22, 41, -76}},
    {"three helices, a nucleotide in every gap",
     "GCAGGACAAAAGUCCAGGACAAAAGUCCAGC",
     "((.((((....)))).((((....)))).))",
     {0, -34, 21, -33, -24, -22, 41, -33, -24, -22, 41, -89}},
    {"four helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))((((....))))))",
     {0, -34, 21, -33, -24, -22, 41, -33, -24, -22, 41, -33, -24, -22, 41, -127}},
    {"four helices, one nucleotide",
     "GCGGACAAAAGUCCGGACAAAAGUCCAGGACAAAAGUCCGC",
     "((((((....))))((((....)))).((((....))))))",
     {0, -34, 21, -33, -24, -22, 41, -33, -24, -22, 41, -33, -24, -22, 41, -127}},
    {"five helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))((((....))))((((....))))))",
     {0,  -34, 15,  -33, -24, -22, 41,  -33, -24, -22,
      41, -33, -24, -22, 41,  -33, -24, -22, 41,  -171}},
    {"five helices, one nucleotide",
     "GCAGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((.((((....))))((((....))))((((....))))((((....))))))",
     {0,  -34, -2,  -33, -24, -22, 41,  -33, -24, -22,
      41, -33, -24, -22, 41,  -33, -24, -22, 41,  -188}},
    {"all the nucleotides in one gap",
     "GCAAAAAGGACAAAAGUCCGGACAAAAGUCCGC",
     "((.....((((....))))((((....))))))",
     {0, -34, 34, -33, -24, -22, 41, -33, -24, -22, 41, -76}},
    {"gaps of 1, 0 and 1",
     "GCAGGACAAAAGUCCGGACAAAAGUCCAGC",
     "((.((((....))))((((....)))).))",
     {0, -34, 34, -33, -24, -22, 41, -33, -24, -22, 41, -76}},
    {"gaps of 2, 1 and 0",
     "GCAGGGACAAAAGUCCAGGACAAAAGUCCGC",
     "((..((((....)))).((((....))))))",
     {0, -34, 22, -33, -24, -22, 41, -33, -24, -22, 41, -88}},
    {"eight nucleotides",
     "GCAAAGGACAAAAGUCCAAGGACAAAAGUCCAAAGC",
     "((...((((....))))..((((....))))...))",
     {0, -34, 24, -33, -24, -22, 41, -33, -24, -22, 41, -86}},
    {"nine nucleotides",
     "GCAAAGGACAAAAGUCCAAAGGACAAAAGUCCAAAGC",
     "((...((((....))))...((((....))))...))",
     {0, -34, 24, -33, -24, -22, 41, -33, -24, -22, 41, -86}},
    {"twenty nucleotides",
     "GCAAAAAAAGGACAAAAGUCCAAAAAAGGACAAAAGUCCAAAAAAAGC",
     "((.......((((....))))......((((....)))).......))",
     {0, -34, 24, -33, -24, -22, 41, -33, -24, -22, 41, -86}},
    {"A-U closing pair",
     "GAGGACAAAAGUCCAGGACAAAAGUCCUC",
     "((((((....)))).((((....))))))",
     {0, -24, 42, -33, -24, -22, 41, -33, -24, -22, 41, -58}},
    {"U-G closing pair",
     "GUAGGACAAAAGUCCAGGACAAAAGUCCAGC",
     "((.((((....)))).((((....)))).))",
     {0, -25, 26, -33, -24, -22, 41, -33, -24, -22, 41, -75}},
    {"branches that end in A-U and G-U",
     "GCAGACAAAAGUCUAGGACAAAAGUCUAGC",
     "((((((....)))).((((....)))).))",
     {0, -34, 53, -21, -24, -22, 41, -21, -24, -22, 41, -33}},
    {"G and A around the branches",
     "GCGCGACAAAAGUCGAGCGACAAAAGUCGAGC",
     "((.((((....))))..((((....)))).))",
     {0, -34, 27, -24, -24, -22, 41, -24, -24, -22, 41, -65}},
    {"a branch of one pair",
     "GCACAAAAGAGGACAAAAGUCCGC",
     "((.(....).((((....))))))",
     {0, -34, 40, 41, -33, -24, -22, 41, 9}},
    {"five different helices",
     "GUCAGACAAAAGUCUGGCGGAAACGCCGACAAAAGUCGAAGGACAAAAGUCUUAC",
     "((.((((....)))).(((....)))((((....))))..((((....)))).))",
     {0, -22, 2, -21, -24, -22, 41, -34, -24, 31, -24, -24, -22, 41, -21, -24, -22, 41, -128}},
    {"a branch with loops of its own",
     "GCAGCGAAGGACGUUCGCGAGCCUCGCGCGGAAACGCGGC",
     "((.((((.((.((....))..))))))(((....))).))",
     {0, -34, 34, -34, -24, -24, 13, -33, 23, -24, 45, -34, -24, 31, -85}},
    {"N in every gap",
     "GCNGGACAAAAGUCCNGGACAAAAGUCCNGC",
     "((.((((....)))).((((....)))).))",
     {0, -34, 75, -33, -24, -22, 41, -33, -24, -22, 41, -35}},
    {"pair that cannot form",
     "GCGGACAAAAGUCAAGGACAAAAGUCCGC",
     "((((((....)))).((((....))))))",
     {0, -34, kInf, kInf, -24, -22, 41, -33, -24, -22, 41, kInf}},
    {"N inside a branch of one pair",
     "GCACNAAAGAGGACAAAAGUCCGC",
     "((.(....).((((....))))))",
     {0, -34, kInf, 56, -33, -24, -22, 41, kInf}},
    {"N beyond a closing helix of one pair",
     "NCGGACAAAAGUCCAGGACAAAAGUCCGA",
     ".(((((....)))).((((....))))).",
     {-11, kInf, -33, -24, -22, 41, -33, -24, -22, 41, kInf}},
    {"loop inside a loop",
     "GCAAGGACAAAAGUCCGCGGAAACGCAUAGGACAAAAGUCCGC",
     "((.(((((....))))(((....))).).((((....))))))",
     {0, -34, 46, 42, -33, -24, -22, 41, -34, -24, 31, -33, -24, -22, 41, -49}},
    {"one multibranch loop",
     "AGCGGACAAAAGUCCAGGACAAAAGUCCGCG",
     ".((((((....)))).((((....)))))).",
     {-17, -34, 34, -33, -24, -22, 41, -33, -24, -22, 41, -93}},
    {"three-way junction, nothing unpaired",
     "GACGCGGAAACGCCGACAAAAGUCGGUC",
     "((((((....)))((((....)))))))",
     {0, -24, -22, 42, -34, -24, 31, -24, -24, -22, 41, -60}},
    {"a loop inside a loop",
     "GGCAAGGACAAAAGUCCGCGGAAACGCAUAGGACAAAAGUCCGCAGCGGAAACGC",
     ".((.(((((....))))(((....))).).((((....)))))).(((....)))",
     {-36, -34, 46, 42, -33, -24, -22, 41, -34, -24, 31, -33, -24, -22, 41, -34, -24, 31, -112}},
    {"two multibranch loops and a branch with loops",
     "GCGGACAAAAGUCCAAGCGGAAACGCAGCUCUGCGAAGGACGUUCGCGAGCCUCGCAAGACAAAAGUCUCAAAAGAGC",
     "((((((....))))..(((....))).)).((((((.((.((....))..)))))).((((....))))(....))).",
     {-27, -34, 21, -33, -24, -22, 41, -34, -24, 31,  -21, 34, -34,
      -24, -24, 13, -33, 23,  -24, 45, -21, -24, -22, 41,  41, -135}},
};

const std::vector<Efn2Case> kForbidden = {
    {"C-A pair closing a hairpin", "GGACAAAAAUCC", "((((....))))", {0, -33, -24, kInf, kInf, kInf}},
    {"hairpin of two nucleotides", "GGACAAGUCC", "((((..))))", {0, -33, -24, -22, kInf, kInf}},
    {"N in a 1x1 loop",
     "GGACNGCCAAAAGGCAGUCC",
     "((((.(((....))).))))",
     {0, -33, -24, -22, kInf, -34, -33, 41, kInf}},
    {"N next to a hairpin's closing helix of one pair",
     "ANCAAAAGNA",
     "..(....)..",
     {0, kInf, kInf}},
};

const std::vector<Efn2Case> kM6A = {
    {"m6A after a G-C helix", "GGACAAAAGUCC6", "((((....)))).", {-13, -33, -24, -22, 41, -51}},
    {"m6A on both sides", "6AGACAAAAGUCU6", ".((((....)))).", {-8, -21, -24, -22, 41, -34}},
    {"m6A-U helix end", "A6GACAAAAGUCUA", ".((((....)))).", {-13, -16, -24, -22, 41, -34}},
    {"flush coaxial stack on an m6A-U pair",
     "GGACAAAAGUCC6GACAAAAGUCU",
     "((((....))))((((....))))",
     {-12, -33, -24, -22, 41, -16, -24, -22, 41, -71}},
    {"coaxial stack across m6A",
     "GGACAAAAGUCC6GGACAAAAGUCCA",
     "((((....)))).((((....)))).",
     {-40, -33, -24, -22, 41, -33, -24, -22, 41, -116}},
    {"total with m6A",
     "6GGACAAAAGUCC6G6GGA6ACUCA",
     ".((((....)))).(((....))).",
     {-40, -33, -24, -22, 41, -18, -16, 31, -81}},
    {"m6A in every gap",
     "GC6GGACAAAAGUCC6GGACAAAAGUCC6GC",
     "((.((((....)))).((((....)))).))",
     {0, -34, 22, -33, -24, -22, 41, -33, -24, -22, 41, -88}},
    {"m6A next to the closing pair",
     "GC6GGACAAAAGUCCGGACAAAAGUCC6GC",
     "((.((((....))))((((....)))).))",
     {0, -34, 39, -33, -24, -22, 41, -33, -24, -22, 41, -71}},
    {"a branch that ends in m6A-U",
     "GCA6GACAAAAGUCUAGGACAAAAGUCCGC",
     "((.((((....)))).((((....))))))",
     {0, -34, 44, -16, -24, -22, 41, -33, -24, -22, 41, -49}},
    {"m6A-U closing pair",
     "G6GGACAAAAGUCCAGGACAAAAGUCCUC",
     "((((((....)))).((((....))))))",
     {0, -18, 79, -33, -24, -22, 41, -33, -24, -22, 41, -15}},
    {"m6A between two helices, nothing else",
     "GCGGACAAAAGUCC6GGACAAAAGUCCGC",
     "((((((....)))).((((....))))))",
     {0, -34, 75, -33, -24, -22, 41, -33, -24, -22, 41, -35}},
    {"total with m6A",
     "6GCGGACAAAAGUCC6G6GGA6ACUCAGCA",
     ".((((((....)))).(((....))).)).",
     {-18, -34, 39, -33, -24, -22, 41, -18, -16, 31, -54}},
};

const std::vector<Efn2Case> kPseudouridine = {
    {"pseudouridine after a G-C helix",
     "GGACAAAAGUCCP",
     "((((....)))).",
     {-15, -33, -24, -22, 41, -53}},
    {"flush coaxial stack on an A-pseudouridine pair",
     "GGACAAAAGUCCAGACAAAAGUCP",
     "((((....))))((((....))))",
     {-20, -33, -24, -22, 41, -27, -24, -22, 41, -90}},
    {"coaxial stack across pseudouridine",
     "GGACAAAAGUCCPGGACAAAAGUCCA",
     "((((....)))).((((....)))).",
     {-40, -33, -24, -22, 41, -33, -24, -22, 41, -116}},
    {"Ψ in every gap",
     "GCPGGACAAAAGUCCPGGACAAAAGUCCPGC",
     "((.((((....)))).((((....)))).))",
     {0, -34, 17, -33, -24, -22, 41, -33, -24, -22, 41, -93}},
    {"a branch that ends in Ψ-A",
     "GCAPGACAAAAGUCAAGGACAAAAGUCCGC",
     "((.((((....)))).((((....))))))",
     {0, -34, 50, -24, -24, -22, 41, -33, -24, -22, 41, -51}},
    {"A-Ψ closing pair",
     "GAGGACAAAAGUCCAGGACAAAAGUCCPC",
     "((((((....)))).((((....))))))",
     {0, -24, 72, -33, -24, -22, 41, -33, -24, -22, 41, -28}},
    {"Ψ between two helices, nothing else",
     "GCGGACAAAAGUCCPGGACAAAAGUCCGC",
     "((((((....)))).((((....))))))",
     {0, -34, 73, -33, -24, -22, 41, -33, -24, -22, 41, -37}},
    {"total with Ψ",
     "PGCGGACAAAAGUCCPGAGGAPACPCAGCA",
     ".((((((....)))).(((....))).)).",
     {-17, -34, 40, -33, -24, -22, 41, -24, -27, 31, -69}},
};

const std::vector<Efn2Case> kM1Pseudouridine = {
    {"m1-pseudouridine after a G-C helix",
     "GGACAAAAGUCC1",
     "((((....)))).",
     {-12, -33, -24, -22, 41, -50}},
    {"flush coaxial stack on an A-m1-pseudouridine pair",
     "GGACAAAAGUCCAGACAAAAGUC1",
     "((((....))))((((....))))",
     {-26, -33, -24, -22, 41, -22, -24, -22, 41, -91}},
    {"coaxial stack across m1-pseudouridine",
     "GGACAAAAGUCC1GGACAAAAGUCCA",
     "((((....)))).((((....)))).",
     {-40, -33, -24, -22, 41, -33, -24, -22, 41, -116}},
    {"m1Ψ in every gap",
     "GC1GGACAAAAGUCC1GGACAAAAGUCC1GC",
     "((.((((....)))).((((....)))).))",
     {0, -34, 26, -33, -24, -22, 41, -33, -24, -22, 41, -84}},
    {"a branch that ends in m1Ψ-A",
     "GCA1GACAAAAGUCAAGGACAAAAGUCCGC",
     "((.((((....)))).((((....))))))",
     {0, -34, 47, -26, -24, -22, 41, -33, -24, -22, 41, -56}},
    {"A-m1Ψ closing pair",
     "GAGGACAAAAGUCCAGGACAAAAGUCC1C",
     "((((((....)))).((((....))))))",
     {0, -25, 73, -33, -24, -22, 41, -33, -24, -22, 41, -28}},
    {"m1Ψ between two helices, nothing else",
     "GCGGACAAAAGUCC1GGACAAAAGUCCGC",
     "((((((....)))).((((....))))))",
     {0, -34, 76, -33, -24, -22, 41, -33, -24, -22, 41, -34}},
    {"total with m1Ψ",
     "1GCGGACAAAAGUCC1GAGGA1AC1CAGCA",
     ".((((((....)))).(((....))).)).",
     {-17, -34, 40, -33, -24, -22, 41, -25, -22, 31, -65}},
};

const std::vector<Efn2Case> kDna = {
    {"one multibranch loop",
     "AGCGGACAAAAGUCCAGGACAAAAGUCCGCG",
     ".((((((....)))).((((....)))))).",
     {-11, -22, 5, -18, -13, -14, 32, -18, -13, -14, 32, -54}},
    {"three helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))))",
     {0, -22, 14, -18, -13, -14, 32, -18, -13, -14, 32, -34}},
    {"all the nucleotides in one gap",
     "GCAAAAAGGACAAAAGUCCGGACAAAAGUCCGC",
     "((.....((((....))))((((....))))))",
     {0, -22, 13, -18, -13, -14, 32, -18, -13, -14, 32, -35}},
    {"twenty nucleotides",
     "GCAAAAAAAGGACAAAAGUCCAAAAAAGGACAAAAGUCCAAAAAAAGC",
     "((.......((((....))))......((((....)))).......))",
     {0, -22, 29, -18, -13, -14, 32, -18, -13, -14, 32, -19}},
    {"two multibranch loops, a bulge and an internal loop",
     "GCGGACAAAAGUCCAAGCGGAAACGCAGCUCUGCGAAGGACGUUCGCGAGCCUCGCAAGACAAAAGUCUCAAAAGAGC",
     "((((((....))))..(((....))).)).((((((.((.((....))..)))))).((((....))))(....))).",
     {-28, -22, 0,  -18, -13, -14, 32, -22, -22, 24,  -13, 12, -22,
      -22, -13, 12, -18, 10,  -22, 25, -13, -13, -14, 32,  32, -110}},
    {"2x3 internal loop",
     "GGACAAGCCAAAAGGCAAAGUCC",
     "((((..(((....)))...))))",
     {0, -18, -13, -14, 19, -22, -18, 32, -34}},
    {"C bulge after a C",
     "GGACCGCCAAAAGGCGUCC",
     "((((.(((....)))))))",
     {0, -18, -13, -14, 3, -22, -18, 32, -50}},
    {"35 nucleotides in a hairpin",
     "GGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGUCC",
     "((((...................................))))",
     {0, -18, -13, -14, 57, 12}},
    {"helices next to each other",
     "AGGACAAAAGUCCAGGACAAAAGUCCGGACAAAAGUCCA",
     ".((((....)))).((((....))))((((....)))).",
     {-42, -18, -13, -14, 32, -18, -13, -14, 32, -18, -13, -14, 32, -81}},
};

const std::vector<Efn2Case> kAt10C = {
    {"one multibranch loop",
     "AGCGGACAAAAGUCCAGGACAAAAGUCCGCG",
     ".((((((....)))).((((....)))))).",
     {-23, -44, 101, -42, -33, -30, 33, -42, -33, -30, 33, -110}},
    {"three helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))))",
     {0, -44, 111, -42, -33, -30, 33, -42, -33, -30, 33, -77}},
    {"all the nucleotides in one gap",
     "GCAAAAAGGACAAAAGUCCGGACAAAAGUCCGC",
     "((.....((((....))))((((....))))))",
     {0, -44, 74, -42, -33, -30, 33, -42, -33, -30, 33, -114}},
    {"twenty nucleotides",
     "GCAAAAAAAGGACAAAAGUCCAAAAAAGGACAAAAGUCCAAAAAAAGC",
     "((.......((((....))))......((((....)))).......))",
     {0, -44, 44, -42, -33, -30, 33, -42, -33, -30, 33, -144}},
    {"two multibranch loops, a bulge and an internal loop",
     "GCGGACAAAAGUCCAAGCGGAAACGCAGCUCUGCGAAGGACGUUCGCGAGCCUCGCAAGACAAAAGUCUCAAAAGAGC",
     "((((((....))))..(((....))).)).((((((.((.((....))..)))))).((((....))))(....))).",
     {-35, -44, 43, -42, -33, -30, 33, -44, -31, 21,  -28, 48, -44,
      -31, -33, 12, -42, 24,  -31, 39, -28, -33, -30, 33,  33, -273}},
    {"2x3 internal loop",
     "GGACAAGCCAAAAGGCAAAGUCC",
     "((((..(((....)))...))))",
     {0, -42, -33, -30, 20, -44, -42, 33, -138}},
    {"C bulge after a C",
     "GGACCGCCAAAAGGCGUCC",
     "((((.(((....)))))))",
     {0, -42, -33, -30, 1, -44, -42, 33, -157}},
    {"35 nucleotides in a hairpin",
     "GGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGUCC",
     "((((...................................))))",
     {0, -42, -33, -30, 55, -50}},
    {"helices next to each other",
     "AGGACAAAAGUCCAGGACAAAAGUCCGGACAAAAGUCCA",
     ".((((....)))).((((....))))((((....)))).",
     {-72, -42, -33, -30, 33, -42, -33, -30, 33, -42, -33, -30, 33, -288}},
};

const std::vector<Efn2Case> kAt60C = {
    {"one multibranch loop",
     "AGCGGACAAAAGUCCAGGACAAAAGUCCGCG",
     ".((((((....)))).((((....)))))).",
     {-12, -25, 46, -26, -17, -15, 48, -26, -17, -15, 48, -11}},
    {"three helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))))",
     {0, -25, 58, -26, -17, -15, 48, -26, -17, -15, 48, 13}},
    {"all the nucleotides in one gap",
     "GCAAAAAGGACAAAAGUCCGGACAAAAGUCCGC",
     "((.....((((....))))((((....))))))",
     {0, -25, 33, -26, -17, -15, 48, -26, -17, -15, 48, -12}},
    {"twenty nucleotides",
     "GCAAAAAAAGGACAAAAGUCCAAAAAAGGACAAAAGUCCAAAAAAAGC",
     "((.......((((....))))......((((....)))).......))",
     {0, -25, 37, -26, -17, -15, 48, -26, -17, -15, 48, -8}},
    {"two multibranch loops, a bulge and an internal loop",
     "GCGGACAAAAGUCCAAGCGGAAACGCAGCUCUGCGAAGGACGUUCGCGAGCCUCGCAAGACAAAAGUCUCAAAAGAGC",
     "((((((....))))..(((....))).)).((((((.((.((....))..)))))).((((....))))(....))).",
     {-19, -25, 25, -26, -17, -15, 48, -25, -18, 40,  -15, 32, -25,
      -18, -17, 13, -26, 22,  -18, 51, -15, -17, -15, 48,  48, 16}},
    {"2x3 internal loop",
     "GGACAAGCCAAAAGGCAAAGUCC",
     "((((..(((....)))...))))",
     {0, -26, -17, -15, 31, -25, -26, 48, -30}},
    {"C bulge after a C",
     "GGACCGCCAAAAGGCGUCC",
     "((((.(((....)))))))",
     {0, -26, -17, -15, 0, -25, -26, 48, -61}},
    {"35 nucleotides in a hairpin",
     "GGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGUCC",
     "((((...................................))))",
     {0, -26, -17, -15, 72, 14}},
    {"helices next to each other",
     "AGGACAAAAGUCCAGGACAAAAGUCCGGACAAAAGUCCA",
     ".((((....)))).((((....))))((((....)))).",
     {-38, -26, -17, -15, 48, -26, -17, -15, 48, -26, -17, -15, 48, -68}},
};

const std::vector<Efn2Case> kAt90C = {
    {"one multibranch loop",
     "AGCGGACAAAAGUCCAGGACAAAAGUCCGCG",
     ".((((((....)))).((((....)))))).",
     {-5, -14, 17, -16, -7, -6, 55, -16, -7, -6, 55, 50}},
    {"three helices, nothing unpaired",
     "GCGGACAAAAGUCCGGACAAAAGUCCGC",
     "((((((....))))((((....))))))",
     {0, -14, 29, -16, -7, -6, 55, -16, -7, -6, 55, 67}},
    {"all the nucleotides in one gap",
     "GCAAAAAGGACAAAAGUCCGGACAAAAGUCCGC",
     "((.....((((....))))((((....))))))",
     {0, -14, 11, -16, -7, -6, 55, -16, -7, -6, 55, 49}},
    {"twenty nucleotides",
     "GCAAAAAAAGGACAAAAGUCCAAAAAAGGACAAAAGUCCAAAAAAAGC",
     "((.......((((....))))......((((....)))).......))",
     {0, -14, 37, -16, -7, -6, 55, -16, -7, -6, 55, 75}},
    {"two multibranch loops, a bulge and an internal loop",
     "GCGGACAAAAGUCCAAGCGGAAACGCAGCUCUGCGAAGGACGUUCGCGAGCCUCGCAAGACAAAAGUCUCAAAAGAGC",
     "((((((....))))..(((....))).)).((((((.((.((....))..)))))).((((....))))(....))).",
     {-10, -14, 16, -16, -7, -6,  55, -14, -10, 50, -7, 27, -14,
      -10, -7,  14, -16, 21, -10, 56, -7,  -7,  -6, 55, 55, 188}},
    {"2x3 internal loop",
     "GGACAAGCCAAAAGGCAAAGUCC",
     "((((..(((....)))...))))",
     {0, -16, -7, -6, 37, -14, -16, 55, 33}},
    {"C bulge after a C",
     "GGACCGCCAAAAGGCGUCC",
     "((((.(((....)))))))",
     {0, -16, -7, -6, 0, -14, -16, 55, -4}},
    {"35 nucleotides in a hairpin",
     "GGACAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGUCC",
     "((((...................................))))",
     {0, -16, -7, -6, 82, 53}},
    {"helices next to each other",
     "AGGACAAAAGUCCAGGACAAAAGUCCGGACAAAAGUCCA",
     ".((((....)))).((((....))))((((....)))).",
     {-18, -16, -7, -6, 55, -16, -7, -6, 55, -16, -7, -6, 55, 60}},
};

// RNAstructure's folding recursion without coaxial stacking (see FoldCase), for the structures
// of the tables above.
const std::vector<FoldCase> kExteriorLoopsWithoutCoaxialStacking = {
    {"no helix", 0},
    {"G-C helix", -38},
    {"A-U helix", -21},
    {"G-U helix", -21},
    {"3' dangle", -55},
    {"5' dangle", -40},
    {"mismatch not better than the dangle", -55},
    {"terminal mismatch", -42},
    {"terminal mismatch among more nucleotides", -42},
    {"flush coaxial stack", -76},
    {"flush coaxial stack of A-U ends", -42},
    {"flush coaxial stack with flanks", -95},
    {"coaxial stack across one nucleotide, mismatch after", -110},
    {"coaxial stack across one nucleotide, mismatch before", -93},
    {"one nucleotide between helices, no other", -93},
    {"two nucleotides between helices", -95},
    {"three nucleotides between helices", -95},
    {"three helices in a row", -114},
    {"three helices, both kinds of coaxial stack", -148},
    {"four different helices", -124},
    {"helix of one pair", 30},
    {"N next to the helix", -38},
    {"two hairpins", -99},
    {"branch with a bulge and an internal loop", -130},
    {"three branches", -119},
};

// It cannot build the structures with a pair that does not form or a forbidden table entry.
const std::vector<FoldCase> kMultibranchLoopsWithoutCoaxialStacking = {
    {"three helices, nothing unpaired", -35},
    {"three helices, one nucleotide", -52},
    {"three helices, two nucleotides", -54},
    {"three helices, a nucleotide in every gap", -86},
    {"four helices, nothing unpaired", -79},
    {"four helices, one nucleotide", -96},
    {"five helices, nothing unpaired", -123},
    {"five helices, one nucleotide", -140},
    {"all the nucleotides in one gap", -54},
    {"gaps of 1, 0 and 1", -69},
    {"gaps of 2, 1 and 0", -69},
    {"eight nucleotides", -86},
    {"nine nucleotides", -86},
    {"twenty nucleotides", -86},
    {"A-U closing pair", -37},
    {"U-G closing pair", -62},
    {"branches that end in A-U and G-U", -15},
    {"G and A around the branches", -58},
    {"a branch of one pair", 16},
    {"five different helices", -78},
    {"a branch with loops of its own", -78},
    {"N in every gap", -35},
    {"loop inside a loop", -18},
    {"one multibranch loop", -69},
    {"three-way junction, nothing unpaired", -27},
    {"a loop inside a loop", -62},
    {"two multibranch loops and a branch with loops", -76},
};

const std::vector<FoldCase> kM6AWithoutCoaxialStacking = {
    {"m6A after a G-C helix", -51},
    {"m6A on both sides", -34},
    {"m6A-U helix end", -34},
    {"flush coaxial stack on an m6A-U pair", -59},
    {"coaxial stack across m6A", -106},
    {"total with m6A", -76},
    {"m6A in every gap", -74},
    {"m6A next to the closing pair", -61},
    {"a branch that ends in m6A-U", -44},
    {"m6A-U closing pair", -36},
    {"m6A between two helices, nothing else", -48},
    {"total with m6A", -48},
};

const std::vector<FoldCase> kDnaWithoutCoaxialStacking = {
    {"one multibranch loop", -32},
    {"three helices, nothing unpaired", -12},
    {"all the nucleotides in one gap", -21},
    {"twenty nucleotides", -5},
    {"two multibranch loops, a bulge and an internal loop", -60},
    {"2x3 internal loop", -34},
    {"C bulge after a C", -50},
    {"35 nucleotides in a hairpin", 12},
    {"helices next to each other", -66},
};

const std::vector<FoldCase> kAt90CWithoutCoaxialStacking = {
    {"one multibranch loop", 77},
    {"three helices, nothing unpaired", 87},
    {"all the nucleotides in one gap", 82},
    {"twenty nucleotides", 72},
    {"two multibranch loops, a bulge and an internal loop", 233},
    {"2x3 internal loop", 33},
    {"C bulge after a C", -4},
    {"35 nucleotides in a hairpin", 53},
    {"helices next to each other", 68},
};

}  // namespace

TEST(Eval, HairpinsMatchEfn2) { expect_efn2(kHairpins, load("turner2004")); }

TEST(Eval, InternalLoopsMatchEfn2) { expect_efn2(kInternalLoops, load("turner2004")); }

TEST(Eval, BulgesMatchEfn2) { expect_efn2(kBulges, load("turner2004")); }

TEST(Eval, ExteriorLoopsMatchEfn2) { expect_efn2(kExteriorLoops, load("turner2004")); }

TEST(Eval, MultibranchLoopsMatchEfn2) { expect_efn2(kMultibranchLoops, load("turner2004")); }

// efn2 --simple: the multibranch rule of RNAstructure's folding algorithms, without the terms
// that depend on the loop as a whole. The same structures as above.
TEST(Eval, SimpleMultibranchRuleMatchesEfn2) {
    expect_efn2(kMultibranchLoopsSimple, load("turner2004"), rs::MultibranchRule::Simple);
}

// A pair RNAstructure cannot form, a hairpin of fewer than 3 nucleotides, N where a table has
// no entry for it: efn2 prints sums around 1400 kcal/mol.
TEST(Eval, ForbiddenLoopsAreInfinite) { expect_efn2(kForbidden, load("turner2004")); }

// efn2 --alphabet m6A, PSU and 1
TEST(Eval, ModifiedBasesMatchEfn2) {
    expect_efn2(kM6A, load("m6A"));
    expect_efn2(kPseudouridine, load("Ψ"));
    expect_efn2(kM1Pseudouridine, load("m1Ψ"));

    // All three merged into turner2004: these structures get the same energies, and those
    // without modified bases the energies of turner2004 alone.
    const rs::Datatable merged = load("turner2004", {"m6A", "Ψ", "m1Ψ"});
    expect_efn2(kM6A, merged);
    expect_efn2(kPseudouridine, merged);
    expect_efn2(kM1Pseudouridine, merged);
    expect_efn2(kMultibranchLoops, merged);
    expect_efn2(kInternalLoops, merged);
}

// efn2 --alphabet dna: the one shipped set with a penalty per unpaired nucleotide of a
// multibranch loop
TEST(Eval, DnaMatchesEfn2) { expect_efn2(kDna, load("dna")); }

// efn2 --temperature 283.15, 333.15 and 363.15. At 90 degrees C the asymmetry and strain
// parameters of multibranch loops are negative.
TEST(Eval, OtherTemperaturesMatchEfn2) {
    expect_efn2(kAt10C, load("turner2004", {}, 10.0));
    expect_efn2(kAt60C, load("turner2004", {}, 60.0));
    expect_efn2(kAt90C, load("turner2004", {}, 90.0));
}

// The functions work from encoded nucleotides alone: there is no sequence here, only codes.
TEST(Eval, EncodedValuesWithoutASequence) {
    const rs::Datatable d = load("turner2004");
    const rs::Code A = d.code("A");
    const rs::Code C = d.code("C");
    const rs::Code G = d.code("G");
    const rs::Code U = d.code("U");
    const rs::Code none = rs::kNoNucleotide;

    // Two adjacent pairs are a stack: the table entry
    EXPECT_EQ(rs::internal_loop_energy(0, 0, G, C, A, U, none, none, none, none, none, none, none,
                                       none, 1, d),
              d.stack(G, C, A, U));
    EXPECT_EQ(d.stack(G, C, A, U), -24);

    // The hairpin C(AAAA)G after an A-U pair, as in GGACAAAAGUCC: efn2 -w has Hairpin = 4.1
    const rs::Code loop[] = {C, A, A, A, A, G};
    EXPECT_EQ(rs::hairpin_energy(4, loop, A, G, U, d), 41);
    // The same loop at the two ends of a strand, as in CAAAAG: also 4.1
    EXPECT_EQ(rs::hairpin_energy(4, loop, none, none, none, d), 41);

    // A C-G helix that goes on with G-C, then an A: efn2 -w has Exterior loop = -1.1 for
    // CGACAAAAGUCGA, and -1.3 with a G before the helix as well (GCGACAAAAGUCGA).
    using Element = rs::LoopElement;
    EXPECT_EQ(rs::exterior_energy({Element::helix(C, G, G, C), Element::unpaired(A)}, d), -11);
    EXPECT_EQ(rs::exterior_energy(
                  {Element::unpaired(G), Element::helix(C, G, G, C), Element::unpaired(A)}, d),
              -13);
    // An A-U helix that goes on with G-C, alone: the end of the helix, 0.5 (AGACAAAAGUCU)
    EXPECT_EQ(rs::exterior_energy({Element::helix(A, U, G, C)}, d), 5);

    // The two sides of an internal loop can come from different strands. efn2 -w has 0.9 for
    // the 1x1 loop of GGACAGCCAAAAGGCAGUCC, ((((.(((....))).)))): here its 5' side ACAGC and
    // its 3' side GCAGU are two arrays.
    const rs::Code five[] = {A, C, A, G, C};   // beyond_i, i, si1 = sp1, p, beyond_p
    const rs::Code three[] = {G, C, A, G, U};  // beyond_q, q, sq1 = sj1, j, beyond_j
    EXPECT_EQ(
        rs::internal_loop_energy(1, 1, five[1], three[3], five[3], three[1], five[2], three[2],
                                 five[2], three[2], five[0], three[4], five[4], three[0], 1, d),
        9);
}

// The stem-level functions: one helix with what stacks on its end, and two helices that stack
// on each other. The values are efn2 -w's "Exterior loop" and "Multiloop" for structures with
// exactly that in them.
TEST(Eval, StemsAndCoaxialStacks) {
    const rs::Datatable d = load("turner2004");
    const rs::Code A = d.code("A");
    const rs::Code C = d.code("C");
    const rs::Code G = d.code("G");
    const rs::Code U = d.code("U");
    const rs::Code none = rs::kNoNucleotide;
    using Stacking = rs::Stacking;

    // The helix GGAC/GUCC in an exterior loop, between an A and a G
    const auto stem = [&](Stacking stacking) {
        return rs::exterior_stem_energy(G, C, G, C, A, G, stacking, d);
    };
    EXPECT_EQ(stem(Stacking::None), 0);       // GGACAAAAGUCC
    EXPECT_EQ(stem(Stacking::Dangle5), -2);   // AGGACAAAAGUCC
    EXPECT_EQ(stem(Stacking::Dangle3), -17);  // GGACAAAAGUCCG
    // efn2 takes the lowest: for AGGACAAAAGUCCG the mismatch is not better than the 3' dangle.
    EXPECT_GT(stem(Stacking::Mismatch), stem(Stacking::Dangle3));
    using Element = rs::LoopElement;
    EXPECT_EQ(rs::exterior_energy(
                  {Element::unpaired(A), Element::helix(G, C, G, C), Element::unpaired(G)}, d),
              -17);
    // For GCGACAAAAGUCGA it is: -1.3
    EXPECT_EQ(rs::exterior_stem_energy(C, G, G, C, G, A, Stacking::Mismatch, d), -13);
    // The end of an A-U helix costs 0.5: AGACAAAAGUCU
    EXPECT_EQ(rs::exterior_stem_energy(A, U, G, C, none, none, Stacking::None, d), 5);

    // A nucleotide that is not there cannot stack.
    EXPECT_EQ(rs::exterior_stem_energy(G, C, G, C, none, G, Stacking::Dangle5, d), kInf);
    EXPECT_EQ(rs::exterior_stem_energy(G, C, G, C, A, none, Stacking::Mismatch, d), kInf);
    // Nor can anything on a pair RNAstructure does not form.
    EXPECT_EQ(rs::exterior_stem_energy(A, C, G, C, A, G, Stacking::Dangle3, d), kInf);

    // Two helices next to each other: GGACAAAAGUCCGGACAAAAGUCC has -2.4 and
    // AGACAAAAGUCUAGACAAAAGUCU -0.3, which is 0.5 for each of the helix ends and -1.3.
    EXPECT_EQ(rs::coaxial_stack_energy(G, C, G, C, d), -24);
    EXPECT_EQ(rs::coaxial_stack_energy(A, U, A, U, d), -13);
    // One nucleotide between them, and one after or before: GGACAAAAGUCCAGGACAAAAGUCCA and
    // AGGACAAAAGUCCAGGACAAAAGUCC both have -3.7.
    EXPECT_EQ(rs::coaxial_stack_mismatch3_energy(G, C, G, C, A, A, d), -37);
    EXPECT_EQ(rs::coaxial_stack_mismatch5_energy(G, C, G, C, A, A, d), -37);

    // A multibranch loop by the simple rule is a sum: GCGGACAAAAGUCCGGACAAAAGUCCGC,
    // ((((((....))))((((....)))))), has 5.1 with efn2 -s. The pair C-G that closes it is a helix
    // the loop reaches at G and leaves at C, whose next pair, G-C, is outside it.
    const int closing = rs::multibranch_stem_energy(G, C, C, G, C, G, Stacking::None, d);
    const int branch = rs::multibranch_stem_energy(G, C, G, C, C, G, Stacking::None, d);
    EXPECT_EQ(closing, d.misc().efn2c);  // these helix ends are 0
    EXPECT_EQ(d.misc().efn2a + closing + 2 * branch + rs::coaxial_stack_energy(G, C, G, C, d), 51);
    // efn2's default energy for the loop is 8.2: three helices with nothing unpaired are
    // strained.
    EXPECT_EQ(51 + d.misc().strain, 82);
}

// The loops of the tables above, put together from the stem-level functions.
TEST(Eval, StemEnergiesAddUpToTheLoops) {
    const rs::Datatable rna = load("turner2004");
    expect_stem_sums(kExteriorLoops, rna);
    expect_stem_sums(kMultibranchLoops, rna);
    expect_stem_sums(kForbidden, rna);
    expect_stem_sums(kDna, load("dna"));
    expect_stem_sums(kAt90C, load("turner2004", {}, 90.0));

    const rs::Datatable merged = load("turner2004", {"m6A", "Ψ", "m1Ψ"});
    expect_stem_sums(kM6A, merged);
    expect_stem_sums(kPseudouridine, merged);
    expect_stem_sums(kM1Pseudouridine, merged);
}

// CoaxialStacking::Off: the helices of a loop do not stack on each other, only unpaired
// nucleotides stack on them. The loop values are what is left of the structures' totals in
// kExteriorLoopsWithoutCoaxialStacking and kMultibranchLoopsWithoutCoaxialStacking.
TEST(Eval, WithoutCoaxialStacking) {
    const rs::Datatable d = load("turner2004");
    const rs::MiscLoop& m = d.misc();
    const rs::Code A = d.code("A");
    const rs::Code C = d.code("C");
    const rs::Code G = d.code("G");
    using Element = rs::LoopElement;
    const Element a = Element::unpaired(A);
    const Element helix = Element::helix(G, C, G, C);  // as GGACAAAAGUCC
    const rs::CoaxialStacking off = rs::CoaxialStacking::Off;

    // GGACAAAAGUCCGGACAAAAGUCC: two helices and nothing unpaired. All that can stack is one
    // helix on the other.
    EXPECT_EQ(rs::exterior_energy({helix, helix}, d), -24);
    EXPECT_EQ(rs::exterior_energy({helix, helix}, d, off), 0);

    // GGACAAAAGUCCAGGACAAAAGUCCA: with coaxial stacking the helices stack across the A between
    // them, -3.7. Without, each takes the A after it as a 3' dangle, twice -1.7.
    const int dangle3 = rs::exterior_stem_energy(G, C, G, C, A, A, rs::Stacking::Dangle3, d);
    EXPECT_EQ(dangle3, -17);
    EXPECT_EQ(rs::exterior_energy({helix, a, helix, a}, d), -37);
    EXPECT_EQ(rs::exterior_energy({helix, a, helix, a}, d, off), 2 * dangle3);
    // AGGACAAAAGUCCAGGACAAAAGUCC: the second helix has no nucleotide after it, and the A
    // between them stacks once, on the first helix.
    EXPECT_EQ(rs::exterior_energy({a, helix, a, helix}, d), -37);
    EXPECT_EQ(rs::exterior_energy({a, helix, a, helix}, d, off), dangle3);

    // A loop with one helix has no coaxial stack to leave out.
    EXPECT_EQ(rs::exterior_energy({a, helix, a}, d, off), rs::exterior_energy({a, helix, a}, d));

    // The multibranch loop of GCGGACAAAAGUCCGGACAAAAGUCCGC, ((((((....))))((((....)))))): its
    // three helices have no unpaired nucleotide to stack on them, so what is left is the
    // offset and the penalty per helix (these helix ends are 0). With coaxial stacking it is
    // 5.1 by the simple rule and 8.2 by efn2's.
    const auto loop = [&](rs::MultibranchRule rule, rs::CoaxialStacking coaxial) {
        return rs::multibranch_energy(C, G, G, C, {helix, helix}, d, rule, coaxial);
    };
    EXPECT_EQ(loop(rs::MultibranchRule::Simple, rs::CoaxialStacking::On), 51);
    EXPECT_EQ(loop(rs::MultibranchRule::Simple, off), m.efn2a + 3 * m.efn2c);
    EXPECT_EQ(loop(rs::MultibranchRule::Simple, off), 75);
    // The rule and the stacking are independent: efn2's rule adds the strain of three helices
    // with nothing unpaired either way.
    EXPECT_EQ(loop(rs::MultibranchRule::Efn2, rs::CoaxialStacking::On), 82);
    EXPECT_EQ(loop(rs::MultibranchRule::Efn2, off), 75 + m.strain);
}

// RNAstructure's own folding recursion with its coaxial stacking disabled gives these
// structures the same energies.
TEST(Eval, WithoutCoaxialStackingMatchesRNAstructureFolding) {
    const rs::Datatable rna = load("turner2004");
    expect_fold_without_coaxial_stacking(kExteriorLoops, kExteriorLoopsWithoutCoaxialStacking, rna);
    expect_fold_without_coaxial_stacking(kMultibranchLoops, kMultibranchLoopsWithoutCoaxialStacking,
                                         rna);
    expect_fold_without_coaxial_stacking(kM6A, kM6AWithoutCoaxialStacking, load("m6A"));
    expect_fold_without_coaxial_stacking(kDna, kDnaWithoutCoaxialStacking, load("dna"));
    expect_fold_without_coaxial_stacking(kAt90C, kAt90CWithoutCoaxialStacking,
                                         load("turner2004", {}, 90.0));
}

// A branch that is not a helix (a pseudoknot, for a program that has them): no helix end,
// nothing stacks on it or across it, and a multibranch loop counts it as a helix. efn2 has no
// such rule, so the expected values are put together from other results.
TEST(Eval, BranchesThatAreNotHelices) {
    const rs::Datatable d = load("turner2004");
    const rs::MiscLoop& m = d.misc();
    using Element = rs::LoopElement;
    const rs::Code A = d.code("A");
    const rs::Code C = d.code("C");
    const rs::Code G = d.code("G");
    const Element a = Element::unpaired(A);
    const Element g = Element::unpaired(G);
    const Element helix = Element::helix(G, C, G, C);  // as GGACAAAAGUCC
    const Element other = Element::other(C, C);

    // Exterior loop: the helices on its two sides are scored apart.
    EXPECT_EQ(rs::exterior_energy({a, helix, a, other, a, helix, g}, d),
              rs::exterior_energy({a, helix, a}, d) + rs::exterior_energy({a, helix, g}, d));
    EXPECT_EQ(rs::exterior_energy({helix, other, helix}, d), 0);
    EXPECT_LT(rs::exterior_energy({helix, helix}, d), 0);

    // Multibranch loop closed by a C-G pair whose helix goes on with G-C
    const int closing_end = d.helix_ends(G, C, C, G);
    // Around one such branch: two helices
    EXPECT_EQ(rs::multibranch_energy(C, G, G, C, {other}, d), m.efn2a + 2 * m.efn2c + closing_end);
    // Around two: three helices with nothing unpaired, which are strained
    EXPECT_EQ(rs::multibranch_energy(C, G, G, C, {other, other}, d),
              m.efn2a + 3 * m.efn2c + m.strain + closing_end);
    // The nucleotides next to it can only stack on the closing pair: a 3' dangle, a 5' dangle
    // or a mismatch (Datatable::dangle and tstkm, with the pair as the loop sees it).
    const int stacking = std::min(
        {0, int{d.dangle(C, G, A, 1)}, int{d.dangle(C, G, A, 2)}, int{d.tstkm(C, G, A, A)}});
    EXPECT_LT(stacking, 0);
    EXPECT_EQ(rs::multibranch_energy(C, G, G, C, {a, other, a}, d),
              m.efn2a + 2 * m.efn2c + closing_end + stacking);
}

// Nothing in the structure: the exterior loop is all there is.
TEST(Eval, UnpairedSequence) {
    const rs::Datatable d = load("turner2004");
    EXPECT_EQ(rs::exterior_energy(elements(d.encode("ACGUACGU"), pair_table("........"), 0, 8), d),
              0);
    EXPECT_EQ(rs::exterior_energy({}, d), 0);
}
