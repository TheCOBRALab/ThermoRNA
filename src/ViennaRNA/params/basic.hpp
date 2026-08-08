#pragma once

#ifdef VRNA_WARN_DEPRECATED
#if defined(__clang__)
#define DEPRECATED(func, msg) func __attribute__((deprecated("", msg)))
#elif defined(__GNUC__)
#define DEPRECATED(func, msg) func __attribute__((deprecated(msg)))
#else
#define DEPRECATED(func, msg) func
#endif
#else
#define DEPRECATED(func, msg) func
#endif

#include "ViennaRNA/model.hpp"
#include "ViennaRNA/params/constants.hpp"

namespace thermorna::viennarna {
/** @brief Typename for the free energy parameter data structure #vrna_params */
typedef struct vrna_param_s vrna_param_t;
/** @brief Typename for the Boltzmann factor data structure #vrna_exp_params */
typedef struct vrna_exp_param_s vrna_exp_param_t;

#define VRNA_GQUAD_MAX_STACK_SIZE 7
#define VRNA_GQUAD_MIN_STACK_SIZE 2
#define VRNA_GQUAD_MAX_LINKER_LENGTH 15
#define VRNA_GQUAD_MIN_LINKER_LENGTH 1
#define VRNA_GQUAD_MIN_BOX_SIZE \
    ((4 * VRNA_GQUAD_MIN_STACK_SIZE) + (3 * VRNA_GQUAD_MIN_LINKER_LENGTH))
#define VRNA_GQUAD_MAX_BOX_SIZE \
    ((4 * VRNA_GQUAD_MAX_STACK_SIZE) + (3 * VRNA_GQUAD_MAX_LINKER_LENGTH))

/**
 *  @brief The datastructure that contains temperature scaled energy parameters.
 */
struct vrna_param_s {
    int id;
    int stack[NBPAIRS + 1][NBPAIRS + 1];
    int hairpin[31];
    int bulge[MAXLOOP + 1];
    int internal_loop[MAXLOOP + 1];
    int mismatchExt[NBPAIRS + 1][5][5];
    int mismatchI[NBPAIRS + 1][5][5];
    int mismatch1nI[NBPAIRS + 1][5][5];
    int mismatch23I[NBPAIRS + 1][5][5];
    int mismatchH[NBPAIRS + 1][5][5];
    int mismatchM[NBPAIRS + 1][5][5];
    int dangle5[NBPAIRS + 1][5];
    int dangle3[NBPAIRS + 1][5];
    int int11[NBPAIRS + 1][NBPAIRS + 1][5][5];
    int int21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
    int int22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
    int ninio[5];
    double lxc;
    int MLbase;
    int MLintern[NBPAIRS + 1];
    int MLclosing;
    int TerminalAU;
    int DuplexInit;
    int Tetraloop_E[200];
    char Tetraloops[1401];
    int Triloop_E[40];
    char Triloops[241];
    int Hexaloop_E[40];
    char Hexaloops[1801];
    int TripleC;
    int MultipleCA;
    int MultipleCB;
    int gquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
    int gquadLayerMismatch;
    int gquadLayerMismatchMax;

    double temperature; /**<  @brief  Temperature used for loop contribution scaling */

    vrna_md_t model_details; /**<  @brief  Model details to be used in the recursions */
    char param_file[256];    /**<  @brief  The filename the parameters were derived from, or empty
                                string if they represent the default */
    int SaltStack;
    int SaltLoop[MAXLOOP + 2];
    double SaltLoopDbl[MAXLOOP + 2];
    int SaltMLbase;
    int SaltMLintern;
    int SaltMLclosing;
    int SaltDPXInit;
};

/**
 *  @brief  The data structure that contains temperature scaled Boltzmann weights of the energy
 * parameters.
 */
struct vrna_exp_param_s {
    int id; /**<  @brief  An identifier for the data structure
             *    @deprecated This attribute will be removed in version 3
             */
    double expstack[NBPAIRS + 1][NBPAIRS + 1];
    double exphairpin[31];
    double expbulge[MAXLOOP + 1];
    double expinternal[MAXLOOP + 1];
    double expmismatchExt[NBPAIRS + 1][5][5];
    double expmismatchI[NBPAIRS + 1][5][5];
    double expmismatch23I[NBPAIRS + 1][5][5];
    double expmismatch1nI[NBPAIRS + 1][5][5];
    double expmismatchH[NBPAIRS + 1][5][5];
    double expmismatchM[NBPAIRS + 1][5][5];
    double expdangle5[NBPAIRS + 1][5];
    double expdangle3[NBPAIRS + 1][5];
    double expint11[NBPAIRS + 1][NBPAIRS + 1][5][5];
    double expint21[NBPAIRS + 1][NBPAIRS + 1][5][5][5];
    double expint22[NBPAIRS + 1][NBPAIRS + 1][5][5][5][5];
    double expninio[5][MAXLOOP + 1];
    double lxc;
    double expMLbase;
    double expMLintern[NBPAIRS + 1];
    double expMLclosing;
    double expTermAU;
    double expDuplexInit;
    double exptetra[40];
    double exptri[40];
    double exphex[40];
    char Tetraloops[1401];
    double expTriloop[40];
    char Triloops[241];
    char Hexaloops[1801];
    double expTripleC;
    double expMultipleCA;
    double expMultipleCB;
    double expgquad[VRNA_GQUAD_MAX_STACK_SIZE + 1][3 * VRNA_GQUAD_MAX_LINKER_LENGTH + 1];
    double expgquadLayerMismatch;
    int gquadLayerMismatchMax;

    double kT;
    double pf_scale; /**<  @brief    Scaling factor to avoid over-/underflows */

    double temperature; /**<  @brief    Temperature used for loop contribution scaling */
    double alpha;       /**<  @brief    Scaling factor for the thermodynamic temperature
                         *    @details  This allows for temperature scaling in Boltzmann
                         *              factors independently from the energy contributions.
                         *              The resulting Boltzmann factors are then computed by
                         *              @f$ e^{-E/(\alpha \cdot K \cdot T)} @f$
                         */

    vrna_md_t model_details; /**<  @brief  Model details to be used in the recursions */
    char param_file[256];    /**<  @brief  The filename the parameters were derived from, or empty
                                string if they represent the default */

    double expSaltStack;
    double expSaltLoop[MAXLOOP + 2];
    double SaltLoopDbl[MAXLOOP + 2];
    int SaltMLbase;
    int SaltMLintern;
    int SaltMLclosing;
    int SaltDPXInit;
};

/**
 *  @brief  Get a data structure containing prescaled free energy parameters
 *
 *  If a NULL pointer is passed for the model details parameter, the default
 *  model parameters are stored within the requested #vrna_param_t structure.
 *
 *  @see #vrna_md_t, vrna_md_set_default(), vrna_exp_params()
 *
 *  @param  md  A pointer to the model details to store inside the structure (Maybe NULL)
 *  @return     A pointer to the memory location where the requested parameters are stored
 */
vrna_param_t* vrna_params(vrna_md_t* md);

/**
 *  @brief Get a copy of the provided free energy parameters
 *
 *  If NULL is passed as parameter, a default set of energy parameters is created
 *  and returned.
 *
 *  @see vrna_params(), #vrna_param_t
 *
 *  @param  par   The free energy parameters that are to be copied (Maybe NULL)
 *  @return       A copy or a default set of the (provided) parameters
 */
vrna_param_t* vrna_params_copy(vrna_param_t* par);

/**
 *  @brief  Get a data structure containing prescaled free energy parameters
 *          already transformed to Boltzmann factors
 *
 *  This function returns a data structure that contains all necessary precomputed
 *  energy contributions for each type of loop.
 *
 *  In contrast to vrna_params(), the free energies within this data structure
 *  are stored as their Boltzmann factors, i.e.
 *
 *  @f$ exp(-E / kT) @f$
 *
 *  where @f$ E @f$ is the free energy.
 *
 *  If a NULL pointer is passed for the model details parameter, the default
 *  model parameters are stored within the requested #vrna_exp_param_t structure.
 *
 *  @see #vrna_md_t, vrna_md_set_default(), vrna_params(), vrna_rescale_pf_params()
 *
 *  @param  md  A pointer to the model details to store inside the structure (Maybe NULL)
 *  @return     A pointer to the memory location where the requested parameters are stored
 */
vrna_exp_param_t* vrna_exp_params(vrna_md_t* md);

/**
 *  @brief  Get a data structure containing prescaled free energy parameters
 *          already transformed to Boltzmann factors (alifold version)
 *
 *  If a NULL pointer is passed for the model details parameter, the default
 *  model parameters are stored within the requested #vrna_exp_param_t structure.
 *
 *  @see #vrna_md_t, vrna_md_set_default(), vrna_exp_params(), vrna_params()
 *
 *  @param  n_seq   The number of sequences in the alignment
 *  @param  md      A pointer to the model details to store inside the structure (Maybe NULL)
 *  @return         A pointer to the memory location where the requested parameters are stored
 */
vrna_exp_param_t* vrna_exp_params_comparative(unsigned int n_seq, vrna_md_t* md);

/**
 *  @brief Get a copy of the provided free energy parameters (provided as Boltzmann factors)
 *
 *  If NULL is passed as parameter, a default set of energy parameters is created
 *  and returned.
 *
 *  @see vrna_exp_params(), #vrna_exp_param_t
 *
 *  @param  par   The free energy parameters that are to be copied (Maybe NULL)
 *  @return       A copy or a default set of the (provided) parameters
 */
vrna_exp_param_t* vrna_exp_params_copy(vrna_exp_param_t* par);

}  // namespace thermorna::viennarna
