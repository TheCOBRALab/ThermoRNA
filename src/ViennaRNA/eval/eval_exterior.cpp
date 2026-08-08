#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "exterior.hpp"
#include "ViennaRNA/loading_params/params/default.hpp"
#include "ViennaRNA/loading_params/params/basic.hpp"


namespace thermorna::viennarna {
PUBLIC int
vrna_E_exterior_stem(unsigned int type,
                     int          n5d,
                     int          n3d,
                     vrna_param_t *p)
{
  int energy = 0;

  if (n5d >= 0 && n3d >= 0)
    energy += p->mismatchExt[type][n5d][n3d];
  else if (n5d >= 0)
    energy += p->dangle5[type][n5d];
  else if (n3d >= 0)
    energy += p->dangle3[type][n3d];

  if (type > 2)
    energy += p->TerminalAU;

  return energy;
}
}