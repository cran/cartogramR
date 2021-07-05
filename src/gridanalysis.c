/** Calculate the grid used in cartogramR for a given L and padding
 * return a list of sfg POINT (the grid)
 *
 *******************************************************************/

/******************************** Inclusions. **********************/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "caracmap.h"
#include "gridanalysis.h"


/**************************** Function prototypes. ***************************/
/********************************** Functions *******************************/

/** \fn gridanalysis
 *  \brief Calculate the grid used in cartogramR for a given L and padding
 * return a list of sfg POINT (the grid)
 *
 * \param  rpadding : SEXP, the padding used in cartogram
 *         Determines space between map and boundary (default to 1.5)
 * \param  rLL : SEXP, the value of L in cartogram  (default is 512),
 *         must be a power of two (for fftw)
 * \param  rbbox: SEXP, the bounding box in cartogram
 * \param  roptions Integer vector of options
 * \return rans : SEXP, The R list of sfg POINT (the grid)
 *******************************************************************/

SEXP gridanalysis (SEXP rpadding, SEXP rLL, SEXP rbbox, SEXP roptions)
{
  /*****************************************************************************/
  /* input and output from/to R */
  /*****************************************************************************/
  /* list output is an R object */
  SEXP  rans, rcoord, rclass, rclassans, rcrs, rclasscrs, rnames,
    rnamesbbox, rclassbbox, rbbox2;
  /*****************************************************************************/
  /* processing input  from R */
  /*****************************************************************************/
  /* double */
  rbbox = PROTECT(rbbox);
  rpadding = PROTECT(rpadding);
  double *bbox, padding;
  bbox = REAL(rbbox);
  padding = REAL(rpadding)[0];
  /* integer : option(s)   */
  rLL = PROTECT(rLL);
  roptions = PROTECT(roptions);
  int LL, *options;
  LL = INTEGER(rLL)[0];
  options = INTEGER(roptions);
  /************************************************************************/
  /* local variables */
  /************************************************************************/
  int i, j, lx, ly, iter;
  /* number of rows in y_geom */
  int n_comp, *caracmapi;
  double map_minx, map_maxx, map_miny, map_maxy, *caracmapd,
    latt_const, new_minx, new_miny;
  /* bounding box (SF order)*/
  map_minx= bbox[0];
  map_miny = bbox[1];
  map_maxx = bbox[2];
  map_maxy = bbox[3];
  /************************************************************************/
  /* Map  */
  /************************************************************************/
  caracmapd = (double *) R_alloc(3, sizeof(double));
  caracmapi = (int *) R_alloc(2, sizeof(int));
  caract_map(caracmapd, caracmapi, padding, LL, map_maxx, map_maxy, map_minx,map_miny);
  lx = caracmapi[0];
  ly = caracmapi[1];
  latt_const = caracmapd[0];
  new_minx = caracmapd[1];
  new_miny = caracmapd[2];
  /************************************************************************/
  /* result  R list of n_comp components */
  /************************************************************************/
  n_comp = lx * ly;
  rans  = PROTECT(allocVector(VECSXP, n_comp));
   /************************************************************************/
   /* class and other attributes of rans */
   /************************************************************************/
   rclassans = PROTECT(allocVector(STRSXP, 2));
   SET_STRING_ELT(rclassans, 0, mkChar("sfc_POINT"));
   SET_STRING_ELT(rclassans, 1, mkChar("sfc"));
   classgets(rans, rclassans);
   setAttrib(rans, install("precision"), ScalarReal(0));
   setAttrib(rans, install("n_empty"), ScalarInteger(0));
  /************************************************************************/
  /* attribute for ans : crs */
  /* list of two components : 2 vectors of character, NA and NA */
  /************************************************************************/
   rcrs  = PROTECT(allocVector(VECSXP, 2));
   SEXP rinput = PROTECT(allocVector(STRSXP, 1));
   SET_STRING_ELT(rinput, 0, NA_STRING);
   SEXP rwkt = PROTECT(allocVector(STRSXP, 1));
   SET_STRING_ELT(rwkt, 0, NA_STRING);
   /* set the vectors in rcrs */
   SET_VECTOR_ELT(rcrs, 0, rinput);
   SET_VECTOR_ELT(rcrs, 1, rwkt);
   /* names of components of the list rcrs */
   rnames = PROTECT(allocVector(STRSXP, 2));
   SET_STRING_ELT(rnames, 0, mkChar("input"));
   SET_STRING_ELT(rnames, 1, mkChar("wkt"));
   /* assign names to list  rcrs */
   setAttrib(rcrs, R_NamesSymbol, rnames);
   /* class crs */
   rclasscrs = PROTECT(allocVector(STRSXP, 1));
   SET_STRING_ELT(rclasscrs, 0, mkChar("crs"));
   classgets(rcrs, rclasscrs);
   /* attribute crs */
   setAttrib(rans, install("crs"), rcrs);
  /************************************************************************/
  /* class of  each component of rans  */
  /************************************************************************/
  rclass = PROTECT(allocVector(STRSXP, 3));
  SET_STRING_ELT(rclass, 0, mkChar("XY"));
  SET_STRING_ELT(rclass, 1, mkChar("POINT"));
  SET_STRING_ELT(rclass, 2, mkChar("sfg"));
  /************************************************************************/
  /* Grid  */
  /************************************************************************/
  double coordx, coordy, minx, miny, maxx, maxy, *bbox2;
   iter=0;
   for (i=0; i<lx; i++) {
      for (j=0; j<ly; j++) {
	/* grid */
	coordx = i + 0.5;
	coordy = j + 0.5;
	/* rescale and set */
	coordx = coordx*latt_const + new_minx;
	coordy = coordy*latt_const + new_miny;
	if (iter==0) {
	  minx=coordx;
	  miny=coordy;
	  maxx=coordx;
	  maxy=coordy;
	} else {
	  minx=fmin2(minx, coordx);
	  miny=fmin2(miny, coordy);
	  maxx=fmax2(maxx, coordx);
	  maxy=fmax2(maxy, coordy);
	}
	/* save  */
	rcoord  = PROTECT(allocVector(REALSXP, 2));
	REAL(rcoord)[0] = coordx;
	REAL(rcoord)[1] = coordy;
	classgets(rcoord, rclass);
	SET_VECTOR_ELT(rans, iter, rcoord);
	UNPROTECT(1);
	iter++;
      }
    }
   /* bbox : last attribute for ans */
   rbbox2 = PROTECT(allocVector(REALSXP, 4));
   /* assign values to vector rbbox2 */
   bbox2 = REAL(rbbox2);
   bbox2[0]=minx;
   bbox2[1]=miny;
   bbox2[2]=maxx;
   bbox2[3]=maxy;
   /* names of components of the vector rbbox2 */
   rnamesbbox = PROTECT(allocVector(STRSXP, 4));
   SET_STRING_ELT(rnamesbbox, 0, mkChar("xmin"));
   SET_STRING_ELT(rnamesbbox, 1, mkChar("ymin"));
   SET_STRING_ELT(rnamesbbox, 2, mkChar("xmax"));
   SET_STRING_ELT(rnamesbbox, 3, mkChar("ymax"));
   /* assign names to vector rbbox2 */
   setAttrib(rbbox2, R_NamesSymbol, rnamesbbox);
   /* class bbox */
   rclassbbox = PROTECT(allocVector(STRSXP, 1));
   SET_STRING_ELT(rclassbbox, 0, mkChar("bbox"));
   classgets(rbbox2, rclassbbox);
   /* set bbox attribute to rans */
   setAttrib(rans, install("bbox"), rbbox2);
   /* unprotect and return */
   UNPROTECT(6);
   UNPROTECT(9); /* class and attributes */
   return rans;
   /* class for ans */
}
