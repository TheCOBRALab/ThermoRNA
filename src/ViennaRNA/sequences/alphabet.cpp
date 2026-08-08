/*
 *  alphabet.c
 *
 *  Code for handling nucleotide and base pair alphabet
 *
 *  Part of the ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/sequences/alphabet.hpp"
#include "ViennaRNA/utils/basic.hpp"
#include "ViennaRNA/utils/utils.hpp"

#include <ctype.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cctype>
#include <cstring>
#include <limits>
#include <stdexcept>

/*
 * For now, we neglect all non-standard nucleotides in an input sequence, i.e. only
 * ACGTUN is allowed.
 *
 * However, the standard nucleotide ambiguity code table would allow for many more:
 *
 * A = Adenylic acid
 * C = Cytidylic acid
 * G = Guanylic acid
 * T = Thymidylic acid
 * U = Uridylic acid
 * I = Inosylic acid
 * R = A or G = puRine
 * Y = C or T = pYrimidine
 * K = G or T = Keto
 * M = A or C = aMino
 * S = G or C = Strong base pair
 * W = A or T = Weak base pair
 * B = not A (G or C or T)
 * D = not C (A or G or T)
 * H = not G (A or C or T)
 * V = not T/U (A or C or G)
 * N = aNy base  (by convention, X is used for unknown amino acids, N for unknown nucleotides)
 *
 * For the future, we aim to accept all of the above codes.
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

namespace thermorna::viennarna {
PRIVATE const char Law_and_Order[] = "_ACGUTXKI";

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC int vrna_nucleotide_IUPAC_identity(char nt, char mask) {
    char n1, n2;
    const char* p = nullptr;

    n1 = static_cast<char>(toupper(static_cast<unsigned char>(nt)));
    n2 = static_cast<char>(toupper(static_cast<unsigned char>(mask)));

    switch (n1) {
        case 'A': p = strchr("ARMWDHVN", n2); break;
        case 'C': p = strchr("CYMSBHVN", n2); break;
        case 'G': p = strchr("GRKSBDVN", n2); break;
        case 'T': p = strchr("TYKWBDHN", n2); break;
        case 'U': p = strchr("UYKWBDHN", n2); break;
        case 'I': p = strchr("IN", n2); break;
        case 'R': p = strchr("AGR", n2); break;
        case 'Y': p = strchr("CTUY", n2); break;
        case 'K': p = strchr("GTUK", n2); break;
        case 'M': p = strchr("ACM", n2); break;
        case 'S': p = strchr("GCS", n2); break;
        case 'W': p = strchr("ATUW", n2); break;
        case 'B': p = strchr("GCTBU", n2); break;
        case 'D': p = strchr("AGTUD", n2); break;
        case 'H': p = strchr("ACTUH", n2); break;
        case 'V': p = strchr("ACGV", n2); break;
        case 'N': p = strchr("ACGTUN", n2); break;
    }

    return p ? 1 : 0;
}

PUBLIC short* vrna_seq_encode(const char* sequence, const vrna_md_t* md) {
    unsigned int i, l;
    short* S = NULL;

    if (sequence && md) {
        S = vrna_seq_encode_simple(sequence, md);

        l = (unsigned int) strlen(sequence);

        for (i = 1; i <= l; i++) S[i] = md->alias[S[i]];

        S[l + 1] = S[1];
        S[0] = S[l];
    }

    return S;
}

PUBLIC short* vrna_seq_encode_simple(const char* sequence, const vrna_md_t* md) {
    unsigned int i, l;
    short* S = NULL;

    if (sequence && md) {
        l = (unsigned int) strlen(sequence);
        S = (short*) vrna_alloc(sizeof(short) * (l + 2));

        for (i = 1; i <= l; i++) /* make numerical encoding of sequence */
            S[i] = (short) vrna_nucleotide_encode(sequence[i - 1], md);

        S[l + 1] = S[1];
        S[0] = (short) l;
    }

    return S;
}

PUBLIC int vrna_nucleotide_encode(char c, const vrna_md_t* md) {
    /* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */
    int code = -1;

    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));

    if (md) {
        if (md->energy_set > 0) {
            code = (int) (c - 'A') + 1;
        } else {
            const char* pos;
            pos = strchr(Law_and_Order, c);
            if (pos == NULL)
                code = 0;
            else
                code = (int) (pos - Law_and_Order);

            if (code > 5) code = 0;

            if (code > 4) code--; /* make T and U equivalent */
        }
    }

    return code;
}

PUBLIC char vrna_nucleotide_decode(int enc, const vrna_md_t* md) {
    if (md) {
        if (md->energy_set > 0)
            return static_cast<char>(enc + 'A' - 1);
        else
            return (char) Law_and_Order[enc];
    } else {
        return (char) 0;
    }
}

PUBLIC void vrna_aln_encode(const char* sequence, short** S_p, short** s5_p, short** s3_p,
                            char** ss_p, unsigned int** as_p, const vrna_md_t* md) {
    std::size_t i, l, p;

    l = strlen(sequence);
    if (l > std::numeric_limits<unsigned int>::max()) {
        throw std::length_error("Sequence too long");
    }
        
    (*s5_p) = (short*) vrna_alloc((l + 2) * sizeof(short));
    (*s3_p) = (short*) vrna_alloc((l + 2) * sizeof(short));
    (*as_p) = (unsigned int*) vrna_alloc((l + 2) * sizeof(unsigned int));
    (*ss_p) = (char*) vrna_alloc((l + 2) * sizeof(char));

    /* make numerical encoding of sequence */
    (*S_p) = vrna_seq_encode_simple(sequence, md);

    (*s5_p)[0] = (*s5_p)[1] = 0;

    if (md->oldAliEn) {
        /* use alignment sequences in all energy evaluations */
        (*ss_p)[0] = sequence[0];
        for (i = 1; i < l; i++) {
            (*s5_p)[i] = (*S_p)[i - 1];
            (*s3_p)[i] = (*S_p)[i + 1];
            (*ss_p)[i] = sequence[i];
            (*as_p)[i] = static_cast<unsigned int>(i);
        }
        (*ss_p)[l] = sequence[l];
        (*as_p)[l] = static_cast<unsigned int>(l);
        (*s5_p)[l] = (*S_p)[l - 1];
        (*s3_p)[l] = 0;
        (*S_p)[l + 1] = (*S_p)[1];
        (*s5_p)[1] = 0;
        if (md->circ) {
            (*s5_p)[1] = (*S_p)[l];
            (*s3_p)[l] = (*S_p)[1];
            (*ss_p)[l + 1] = static_cast<char>((*S_p)[1]);
        }
    } else {
        if (md->circ) {
            for (i = l; i > 0; i--) {
                char c5;
                c5 = sequence[i - 1];
                if ((c5 == '-') || (c5 == '_') || (c5 == '~') || (c5 == '.')) continue;

                (*s5_p)[1] = (*S_p)[i];
                break;
            }
            for (i = 1; i <= l; i++) {
                char c3;
                c3 = sequence[i - 1];
                if ((c3 == '-') || (c3 == '_') || (c3 == '~') || (c3 == '.')) continue;

                (*s3_p)[l] = (*S_p)[i];
                break;
            }
        } else {
            (*s5_p)[1] = (*s3_p)[l] = 0;
        }

        for (i = 1, p = 0; i <= l; i++) {
            char c5;
            c5 = sequence[i - 1];
            if ((c5 == '-') || (c5 == '_') || (c5 == '~') || (c5 == '.')) {
                (*s5_p)[i + 1] = (*s5_p)[i];
            } else {
                /* no gap */
                (*ss_p)[p++] = sequence[i - 1]; /*start at 0!!*/
                (*s5_p)[i + 1] = (*S_p)[i];
            }

            (*as_p)[i] = static_cast<unsigned int>(p);
        }
        for (i = l; i >= 1; i--) {
            char c3;
            c3 = sequence[i - 1];
            if ((c3 == '-') || (c3 == '_') || (c3 == '~') || (c3 == '.'))
                (*s3_p)[i - 1] = (*s3_p)[i];
            else
                (*s3_p)[i - 1] = (*S_p)[i];
        }
    }
}

PUBLIC unsigned int vrna_get_ptype_md(int i, int j, const vrna_md_t* md) {
    unsigned int tt = (unsigned int) md->pair[i][j];

    return (tt == 0) ? 7 : tt;
}

PUBLIC unsigned int vrna_get_ptype(int ij, const char* ptype) {
    unsigned int tt = (unsigned int) ptype[ij];

    return (tt == 0) ? 7 : tt;
}

PUBLIC unsigned int vrna_get_ptype_window(int i, int j, const char** ptype) {
    unsigned int tt = (unsigned int) ptype[i][j - i];

    return (tt == 0) ? 7 : tt;
}

}  // namespace thermorna::viennarna