#ifndef __GEOMCARTO_H__
#define __GEOMCARTO_H__
/********************************** Macros. **********************************/
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)>(b)) ? (b) : (a))
/**************************** Function prototypes. ***************************/
SEXP geomcarto (SEXP rygeom, SEXP rmultipoly, SEXP rgridx, SEXP rgridy, SEXP rpadding,
		 SEXP rLL, SEXP rbbox, SEXP roptions) ;
#endif
