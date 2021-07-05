#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>  // optional
#include "cartogram.h"
#include "checkring.h"
#include "dcn.h"
#include "invcarto.h"
#include "geomcarto.h"
#include "gridanalysis.h"
#include "makegraticule.h"
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
  CALLDEF(makefinalgraticule, 5),
  CALLDEF(makeoriggraticule, 3),
  CALLDEF(invcarto, 7),
  CALLDEF(gridanalysis, 4),
  CALLDEF(geomcarto, 8),
  CALLDEF(dcn, 12),
  CALLDEF(checkring, 3),
  CALLDEF(cartogramR, 12),
   {NULL, NULL, 0}
};

void
attribute_visible  // optional
R_init_cartogramR(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

