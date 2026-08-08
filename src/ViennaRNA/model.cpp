#include "ViennaRNA/model.hpp"
#include "ViennaRNA/params/constants.hpp"
#include "ViennaRNA/utils/utils.hpp"


#include <cstring>

namespace thermorna::viennarna {
/*  below are the evil global variables that will vanish
 *  as soon as we drop backward compatibility in ViennaRNA
 *  Package v3
 */

double temperature = VRNA_MODEL_DEFAULT_TEMPERATURE;
double pf_scale = VRNA_MODEL_DEFAULT_PF_SCALE;
int dangles = VRNA_MODEL_DEFAULT_DANGLES;
int tetra_loop = VRNA_MODEL_DEFAULT_SPECIAL_HP;
int noGU = VRNA_MODEL_DEFAULT_NO_GU;
int no_closingGU = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE;
int energy_set = VRNA_MODEL_DEFAULT_ENERGY_SET;
char *nonstandards = nullptr;
int max_bp_span = VRNA_MODEL_DEFAULT_MAX_BP_SPAN;
int logML = VRNA_MODEL_DEFAULT_LOG_ML;

/* below are some more deprecated global symbols we need to get rid off */

int james_rule = 1; /* interior loops of size 2 get energy 0.8Kcal and
                     * no mismatches (no longer used) */

#define   BP_REV_DEFAULT        { 0, 2, 1, 4, 3, 6, 5, 7 }

#define   BP_ALIAS_DEFAULT      { 0, 1, 2, 3, 4, 3, 2, 0 }

#define   BP_ENCODING_DEFAULT \
  /*  _  A  C  G  U  X  K  I */ \
  { { 0, 0, 0, 0, 0, 0, 0, 0 }, \
    { 0, 0, 0, 0, 5, 0, 0, 5 }, \
    { 0, 0, 0, 1, 0, 0, 0, 0 }, \
    { 0, 0, 2, 0, 3, 0, 0, 0 }, \
    { 0, 6, 0, 4, 0, 0, 0, 6 }, \
    { 0, 0, 0, 0, 0, 0, 2, 0 }, \
    { 0, 0, 0, 0, 0, 1, 0, 0 }, \
    { 0, 6, 0, 0, 5, 0, 0, 0 } }

#define   DM_DEFAULT \
  { { 0, 0, 0, 0, 0, 0, 0 }, /* hamming distance between pairs */ \
    { 0, 0, 2, 2, 1, 2, 2 } /* CG */, \
    { 0, 2, 0, 1, 2, 2, 2 } /* GC */, \
    { 0, 2, 1, 0, 2, 1, 2 } /* GU */, \
    { 0, 1, 2, 2, 0, 2, 1 } /* UG */, \
    { 0, 2, 2, 1, 2, 0, 2 } /* AU */, \
    { 0, 2, 2, 2, 1, 2, 0 } /* UA */ }

static vrna_md_t defaults = {
    VRNA_MODEL_DEFAULT_TEMPERATURE,       // temperature
    1.,                                   // betaScale
    VRNA_MODEL_DEFAULT_PF_SMOOTH,         // pf_smooth
    VRNA_MODEL_DEFAULT_DANGLES,           // dangles
    VRNA_MODEL_DEFAULT_SPECIAL_HP,        // special_hp
    VRNA_MODEL_DEFAULT_NO_LP,             // noLP
    VRNA_MODEL_DEFAULT_NO_GU,             // noGU
    VRNA_MODEL_DEFAULT_NO_GU_CLOSURE,     // noGUclosure
    VRNA_MODEL_DEFAULT_LOG_ML,             // logML
    VRNA_MODEL_DEFAULT_CIRC,               // circ
    VRNA_MODEL_DEFAULT_CIRC_PENALTY,       // circ_penalty
    VRNA_MODEL_DEFAULT_GQUAD,              // gquad
    VRNA_MODEL_DEFAULT_UNIQ_ML,            // uniq_ML
    VRNA_MODEL_DEFAULT_ENERGY_SET,         // energy_set
    VRNA_MODEL_DEFAULT_BACKTRACK,          // backtrack
    VRNA_MODEL_DEFAULT_BACKTRACK_TYPE,     // backtrack_type
    VRNA_MODEL_DEFAULT_COMPUTE_BPP,        // compute_bpp
    { 0 },                                 // nonstandards[64]
    VRNA_MODEL_DEFAULT_MAX_BP_SPAN,        // max_bp_span
    TURN,                                  // min_loop_size
    VRNA_MODEL_DEFAULT_WINDOW_SIZE,        // window_size
    VRNA_MODEL_DEFAULT_ALI_OLD_EN,         // oldAliEn
    VRNA_MODEL_DEFAULT_ALI_RIBO,           // ribo
    VRNA_MODEL_DEFAULT_ALI_CV_FACT,        // cv_fact
    VRNA_MODEL_DEFAULT_ALI_NC_FACT,        // nc_fact
    1.07,                                  // sfact
    BP_REV_DEFAULT,                        // rtype[8]
    BP_ALIAS_DEFAULT,                      // alias[MAXALPHA + 1]
    BP_ENCODING_DEFAULT,                   // pair[MAXALPHA + 1][MAXALPHA + 1]
    DM_DEFAULT,                            // pair_dist[7][7]
    VRNA_MODEL_DEFAULT_SALT,               // salt
    VRNA_MODEL_DEFAULT_SALT_MLLOWER,       // saltMLLower
    VRNA_MODEL_DEFAULT_SALT_MLUPPER,       // saltMLUpper
    VRNA_MODEL_DEFAULT_SALT_DPXINIT,       // saltDPXInit
    VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT,  // saltDPXInitFact
    VRNA_MODEL_DEFAULT_HELICAL_RISE,       // helical_rise
    VRNA_MODEL_DEFAULT_BACKBONE_LENGTH,    // backbone_length
    VRNA_MODEL_DEFAULT_CIRC_ALPHA0         // circ_alpha0
};

void vrna_md_set_default(vrna_md_t *md) {
    if (md) /* copy defaults */
        vrna_md_copy(md, &defaults);
}

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
vrna_md_t *vrna_md_copy(vrna_md_t *md_to, const vrna_md_t *md_from) {
    int i;
    vrna_md_t *md;

    md = nullptr;

    /* only process if md_from is non-NULL */
    if (md_from) {
        if (!md_to) /* create container to be filled */
            md = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
        else
            /* or directly write to target */
            md = md_to;

        /* check if not the same object */
        if (md_to != md_from) {
            /* copy simple members */
            memcpy(md, md_from, sizeof(vrna_md_t));
            /* copy arrays */
            memcpy(md->rtype, &(md_from->rtype[0]), 8 * sizeof(int));
            memcpy(md->alias, &(md_from->alias[0]), (MAXALPHA + 1) * sizeof(short));
            memcpy(md->nonstandards, &(md_from->nonstandards[0]), 64 * sizeof(char));
            /* copy matrices */
            for (i = 0; i <= MAXALPHA; i++)
                memcpy(md->pair[i], (md_from->pair[i]), (MAXALPHA + 1) * sizeof(int));
        }
    }

    return md;
}
} // namespace thermorna::viennarna