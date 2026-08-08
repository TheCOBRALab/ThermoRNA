#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/eval/hairpin.hpp"
#include "ViennaRNA/params/basic.hpp"
#include "ViennaRNA/params/constants.hpp"
#include "ViennaRNA/params/salt.hpp"

#include <cstring>

#ifdef __GNUC__
#define INLINE inline
#else
#define INLINE
#endif

namespace thermorna::viennarna {
PUBLIC int vrna_E_hairpin(unsigned int size, unsigned int type, int si1, int sj1,
                          const char* sequence, vrna_param_t* P) {
    int energy, salt_correction;

    energy = INF;

    if (P) {
        salt_correction = 0;

        if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
            if (size <= MAXLOOP)
                salt_correction = P->SaltLoop[size + 1];
            else
                salt_correction =
                    vrna_salt_loop_int(size + 1, P->model_details.salt, P->temperature + K0,
                                       P->model_details.backbone_length);
        }

        if (size <= 30)
            energy = P->hairpin[size];
        else
            energy = P->hairpin[30] + (int) (P->lxc * log((size) / 30.));

        energy += salt_correction;

        if (size < 3) return energy; /* should only be the case when folding alignments */

        if ((sequence) && (P->model_details.special_hp)) {
            if (size == 4) {
                /* check for tetraloop bonus */
                char tl[7] = {0}, *ts;
                memcpy(tl, sequence, sizeof(char) * 6);
                tl[6] = '\0';
                if ((ts = strstr(P->Tetraloops, tl)))
                    return P->Tetraloop_E[(ts - P->Tetraloops) / 7] + salt_correction;
            } else if (size == 6) {
                char tl[9] = {0}, *ts;
                memcpy(tl, sequence, sizeof(char) * 8);
                tl[8] = '\0';
                if ((ts = strstr(P->Hexaloops, tl)))
                    return P->Hexaloop_E[(ts - P->Hexaloops) / 9] + salt_correction;
            } else if (size == 3) {
                char tl[6] = {0}, *ts;
                memcpy(tl, sequence, sizeof(char) * 5);
                tl[5] = '\0';
                if ((ts = strstr(P->Triloops, tl)))
                    return P->Triloop_E[(ts - P->Triloops) / 6] + salt_correction;

                return energy + (type > 2 ? P->TerminalAU : 0);
            }
        }

        if ((si1 >= 0) && (sj1 >= 0)) energy += P->mismatchH[type][si1][sj1];
    }

    return energy;
}

}  // namespace thermorna::viennarna