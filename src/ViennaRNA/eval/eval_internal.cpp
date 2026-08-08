#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/loading_params/params/basic.hpp"
#include "ViennaRNA/eval/internal.hpp"
#include "ViennaRNA/loading_params/params/default.hpp"
#include "ViennaRNA/loading_params/params/salt.hpp"

namespace thermorna::viennarna {
PUBLIC int vrna_E_internal(unsigned int  n1,
                unsigned int  n2,
                unsigned int  type,
                unsigned int  type_2,
                int           si1,
                int           sj1,
                int           sp1,
                int           sq1,
                vrna_param_t  *P)
{
  /* compute energy of degree 2 loop (stack bulge or internal) */
  unsigned int  nl, ns, u, backbones, no_close;
  int           energy, salt_stack_correction, salt_loop_correction;

  energy = INF;

  if (P) {
    no_close              = 0;
    salt_stack_correction = P->SaltStack;
    salt_loop_correction  = 0;

    if ((P->model_details.noGUclosure) &&
        ((type_2 == 3) || (type_2 == 4) || (type == 3) || (type == 4)))
      no_close = 1;

    if (n1 > n2) {
      nl  = n1;
      ns  = n2;
    } else {
      nl  = n2;
      ns  = n1;
    }

    if (nl == 0)  /* stack */
      return P->stack[type][type_2] +
             salt_stack_correction;

    if (no_close)
      return INF;

    if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
      /* salt correction for loop */
      backbones = nl + ns + 2;
      if (backbones <= MAXLOOP + 1)
        salt_loop_correction = P->SaltLoop[backbones];
      else
        salt_loop_correction = vrna_salt_loop_int(backbones,
                                                  P->model_details.salt,
                                                  P->temperature + K0,
                                                  P->model_details.backbone_length);
    }

    energy = 0;

    switch (ns) {
      case 0:
        /* bulge */
        energy = (nl <= MAXLOOP) ?
                 P->bulge[nl] :
                 (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
        if (nl == 1) {
          /* add stacking energy for 1-bulges */
          energy += P->stack[type][type_2];
        } else {
          if (type > 2)
            energy += P->TerminalAU;

          if (type_2 > 2)
            energy += P->TerminalAU;
        }

        break;

      case 1:
        if (nl == 1) {
          /* 1x1 loop */
          energy = P->int11[type][type_2][si1][sj1];
        } else if (nl == 2) {
          /* 2x1 loop */
          if (n1 == 1)
            energy = P->int21[type][type_2][si1][sq1][sj1];
          else
            energy = P->int21[type_2][type][sq1][si1][sp1];
        } else {
          /* 1xn loop */
          energy = (nl + 1 <= MAXLOOP) ?
                   (P->internal_loop[nl + 1]) :
                   (P->internal_loop[30] + (int)(P->lxc * log((nl + 1) / 30.)));
          energy  += MIN2(MAX_NINIO, (int)(nl - ns) * P->ninio[2]);
          energy  += P->mismatch1nI[type][si1][sj1] +
                     P->mismatch1nI[type_2][sq1][sp1];
        }

        break;

      case 2:
        if (nl == 2) {
          /* 2x2 loop */
          energy = P->int22[type][type_2][si1][sp1][sq1][sj1];
          break;
        } else if (nl == 3) {
          /* 2x3 loop */
          energy = P->internal_loop[5] +
                   P->ninio[2];
          energy += P->mismatch23I[type][si1][sj1] +
                    P->mismatch23I[type_2][sq1][sp1];
          break;
        }

      /* fall through */

      default:
        /* generic internal loop (no else here!)*/
        u       = nl + ns;
        energy  += (u <= MAXLOOP) ?
                   (P->internal_loop[u]) :
                   (P->internal_loop[30] + (int)(P->lxc * log((u) / 30.)));

        energy += std::min(MAX_NINIO, (int)(nl - ns) * P->ninio[2]);

        energy += P->mismatchI[type][si1][sj1] +
                  P->mismatchI[type_2][sq1][sp1];

        break;
    }

    energy += salt_loop_correction;
  }

  return energy;
}
}