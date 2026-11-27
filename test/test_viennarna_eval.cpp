// ViennaRNA loop energies: vrna_E_internal() applies the salt correction and noGUclosure on top
// of the loop energy, the stem energies pick a mismatch or dangle by the neighbours given, and
// hairpins use the special-loop tables.
#include <ViennaRNA/eval/exterior.hpp>
#include <ViennaRNA/eval/hairpin.hpp>
#include <ViennaRNA/eval/internal.hpp>
#include <ViennaRNA/eval/multibranch.hpp>
#include <ViennaRNA/model.hpp>
#include <ViennaRNA/params/basic.hpp>
#include <ViennaRNA/params/salt.hpp>
#include <gtest/gtest.h>

#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>

namespace v = thermorna::viennarna;

namespace {

struct FreeParams {
    void operator()(v::vrna_param_t* P) const { std::free(P); }
};
using Params = std::unique_ptr<v::vrna_param_t, FreeParams>;

Params make_params(double salt, int no_gu_closure) {
    v::vrna_md_t md;
    v::vrna_md_set_default(&md);
    md.salt = salt;
    md.noGUclosure = no_gu_closure;
    return Params(v::vrna_params(&md));
}

// Pair types: 1 CG, 2 GC, 3 GU, 4 UG, 5 AU, 6 UA. Bases: 1 A, 2 C, 3 G, 4 U.
constexpr unsigned int CG = 1, GC = 2, GU = 3, UG = 4, AU = 5;

bool closed_by_gu(unsigned int t1, unsigned int t2) {
    return t1 == GU || t1 == UG || t2 == GU || t2 == UG;
}

}  // namespace

// The model options change nothing but their own corrections.
TEST(ViennaRNAEval, InternalModelOptionsOnlyAddTheirCorrections) {
    Params P = make_params(VRNA_MODEL_DEFAULT_SALT, 0);
    Params G = make_params(VRNA_MODEL_DEFAULT_SALT, 1);
    Params S = make_params(0.05, 0);
    for (unsigned int n1 = 0; n1 <= 33; ++n1)
        for (unsigned int n2 = 0; n2 <= 33; ++n2) {
            const unsigned int backbones = n1 + n2 + 2;
            const int salt =
                (n1 == 0 && n2 == 0) ? S->SaltStack
                : (backbones <= MAXLOOP + 1)
                    ? S->SaltLoop[backbones]
                    : v::vrna_salt_loop_int(backbones, S->model_details.salt, S->temperature + K0,
                                            S->model_details.backbone_length);
            for (unsigned int t1 = 1; t1 <= NBPAIRS; ++t1)
                for (unsigned int t2 = 1; t2 <= NBPAIRS; ++t2)
                    for (int a = 0; a <= 4; a += 2)
                        for (int b = 0; b <= 4; b += 2)
                            for (int c = 1; c <= 4; c += 3)
                                for (int d = 1; d <= 4; d += 3) {
                                    const int e =
                                        v::vrna_E_internal(n1, n2, t1, t2, a, b, c, d, P.get());
                                    const bool forbidden = (n1 + n2 > 0) && closed_by_gu(t1, t2);
                                    ASSERT_EQ(
                                        v::vrna_E_internal(n1, n2, t1, t2, a, b, c, d, G.get()),
                                        forbidden ? INF : e)
                                        << n1 << "x" << n2 << " types " << t1 << "," << t2;
                                    ASSERT_EQ(
                                        v::vrna_E_internal(n1, n2, t1, t2, a, b, c, d, S.get()),
                                        e + salt)
                                        << n1 << "x" << n2 << " types " << t1 << "," << t2;
                                }
        }
}

TEST(ViennaRNAEval, InternalForbidsGUClosingPairsWithNoGUClosure) {
    Params P = make_params(VRNA_MODEL_DEFAULT_SALT, 0);
    Params G = make_params(VRNA_MODEL_DEFAULT_SALT, 1);
    EXPECT_EQ(v::vrna_E_internal(1, 2, GU, CG, 1, 2, 3, 4, G.get()), INF);
    EXPECT_EQ(v::vrna_E_internal(2, 2, CG, UG, 1, 2, 3, 4, G.get()), INF);
    EXPECT_EQ(v::vrna_E_internal(0, 3, GU, GC, 1, 2, 3, 4, G.get()), INF);  // bulge
    EXPECT_EQ(v::vrna_E_internal(1, 2, CG, AU, 1, 2, 3, 4, G.get()),
              v::vrna_E_internal(1, 2, CG, AU, 1, 2, 3, 4, P.get()));
    // a stack is not a loop
    EXPECT_EQ(v::vrna_E_internal(0, 0, GU, CG, 1, 2, 3, 4, G.get()),
              v::vrna_E_internal(0, 0, GU, CG, 1, 2, 3, 4, P.get()));
}

TEST(ViennaRNAEval, InternalAddsSaltCorrections) {
    Params P = make_params(VRNA_MODEL_DEFAULT_SALT, 0);
    Params S = make_params(0.05, 0);
    ASSERT_NE(S->SaltLoop[7], 0);
    ASSERT_NE(S->SaltStack, 0);
    EXPECT_EQ(v::vrna_E_internal(2, 3, CG, GC, 1, 2, 3, 4, S.get()),
              v::vrna_E_internal(2, 3, CG, GC, 1, 2, 3, 4, P.get()) + S->SaltLoop[7]);
    EXPECT_EQ(v::vrna_E_internal(0, 0, CG, GC, 1, 2, 3, 4, S.get()),
              v::vrna_E_internal(0, 0, CG, GC, 1, 2, 3, 4, P.get()) + S->SaltStack);
    // loops longer than the table: the correction is computed
    EXPECT_NE(v::vrna_E_internal(20, 20, CG, GC, 1, 2, 3, 4, S.get()),
              v::vrna_E_internal(20, 20, CG, GC, 1, 2, 3, 4, P.get()));
}

TEST(ViennaRNAEval, StemEnergies) {
    Params P = make_params(VRNA_MODEL_DEFAULT_SALT, 0);
    EXPECT_EQ(v::vrna_E_multibranch_stem(CG, -1, -1, P.get()), P->MLintern[CG]);
    EXPECT_EQ(v::vrna_E_multibranch_stem(AU, -1, -1, P.get()), P->MLintern[AU] + P->TerminalAU);
    EXPECT_EQ(v::vrna_E_multibranch_stem(CG, 2, -1, P.get()), P->MLintern[CG] + P->dangle5[CG][2]);
    EXPECT_EQ(v::vrna_E_multibranch_stem(CG, -1, 3, P.get()), P->MLintern[CG] + P->dangle3[CG][3]);
    EXPECT_EQ(v::vrna_E_multibranch_stem(GU, 2, 3, P.get()),
              P->MLintern[GU] + P->mismatchM[GU][2][3] + P->TerminalAU);
    EXPECT_EQ(v::vrna_E_exterior_stem(CG, -1, -1, P.get()), 0);
    EXPECT_EQ(v::vrna_E_exterior_stem(UG, 2, 3, P.get()), P->mismatchExt[UG][2][3] + P->TerminalAU);
    EXPECT_EQ(v::vrna_E_exterior_stem(CG, 2, -1, P.get()), P->dangle5[CG][2]);
}

TEST(ViennaRNAEval, HairpinEnergies) {
    Params P = make_params(VRNA_MODEL_DEFAULT_SALT, 0);
    // no special loops of length 5: the loop term plus the terminal mismatch
    EXPECT_EQ(v::vrna_E_hairpin(5, CG, 1, 2, "CAAAAAG", P.get()),
              P->hairpin[5] + P->mismatchH[CG][1][2]);
    // a tabulated tetraloop has its own energy
    ASSERT_GE(std::strlen(P->Tetraloops), 6u);
    const std::string tetraloop(P->Tetraloops, 6);
    EXPECT_EQ(v::vrna_E_hairpin(4, CG, 1, 2, tetraloop.c_str(), P.get()), P->Tetraloop_E[0]);

    Params S = make_params(0.05, 0);
    EXPECT_EQ(v::vrna_E_hairpin(5, CG, 1, 2, "CAAAAAG", S.get()),
              S->hairpin[5] + S->mismatchH[CG][1][2] + S->SaltLoop[6]);
}

TEST(ViennaRNAEval, NullParametersGiveINF) {
    EXPECT_EQ(v::vrna_E_internal(2, 3, CG, GC, 1, 2, 3, 4, nullptr), INF);
    EXPECT_EQ(v::vrna_E_multibranch_stem(CG, 2, 3, nullptr), INF);
    EXPECT_EQ(v::vrna_E_hairpin(5, CG, 1, 2, "CAAAAAG", nullptr), INF);
}
