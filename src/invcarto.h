#ifndef __INVCARTO_H_
#define __INVCARTO_H_
#include "typepoint.h"
/********************************** Macros. **********************************/
#define FREEIC1 free(xdisp); \
  free(ydisp); \
  free(invproj); \
  for (i=0; i<=lx; i++) free(projgrid[i]); \
  free(projgrid); \
  for (i=0; i<4*lx*ly; i++) free(tri[i]); \
  free(tri); \
  for (i=0; i<lx; i++) free(xyhalfshift2tri[i]); \
  free(xyhalfshift2tri)
#define FREEIC2 free(proj); \
  free(invproj2)
/**************************** Function prototypes. ***************************/
SEXP invcarto (SEXP rgridx, SEXP rgridy, SEXP rpadding,
		 SEXP rLL, SEXP rbbox, SEXP roptions, SEXP rdensityct);
#endif // __INVCARTO_H_
