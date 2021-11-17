/** Calculate the grid used in cartogramR for a given L and padding
 * return a list of sfg POINT (the grid)
 *
 *******************************************************************/

/******************************** Inclusions. **********************/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "caracmap.h"
#include "makegraticule.h"
#include "interpol2.h"


/**************************** Function prototypes. ***************************/
/********************************** Functions *******************************/

/** \fn makeoriggraticule
 *  \brief Calculate the graticule used in cartogramR for a given L and padding
 * return a list of sfg POLYGON (the graticule)
 *
 * \param  rpadding : SEXP, the padding used in cartogram
 *         Determines space between map and boundary (default to 1.5)
 * \param  rLL : SEXP, the value of L in cartogram  (default is 512),
 *         must be a power of two (for fftw)
 * \param  rbbox: SEXP, the bounding box in cartogram
 * \return rans : SEXP, The R list of sfg POLYGON
 *******************************************************************/

SEXP makeoriggraticule (SEXP rpadding, SEXP rLL, SEXP rbbox)
{
  /*****************************************************************************/
  /* input and output from/to R */
  /*****************************************************************************/
  /* list output is an R object */
  SEXP  rans, rcoord, rclass, rclassans, rcrs, rclasscrs, rnames,
    rnamesbbox, rclassbbox, rbbox2, rcolans, rpolygext, rscalarreal,
    rscalarinteger;
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
  int LL;
  LL = INTEGER(rLL)[0];
  /************************************************************************/
  /* local variables */
  /************************************************************************/
  int i, j, lx, ly;
  /* number of rows in y_geom */
  int *caracmapi;
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
  /* result  R list of lx components */
  /************************************************************************/
  rans  = PROTECT(allocVector(VECSXP, lx));
   /************************************************************************/
   /* class and other attributes of rans */
   /************************************************************************/
   rclassans = PROTECT(allocVector(STRSXP, 2));
   SET_STRING_ELT(rclassans, 0, mkChar("sfc_MULTIPOLYGON"));
   SET_STRING_ELT(rclassans, 1, mkChar("sfc"));
   classgets(rans, rclassans);
   rscalarreal = PROTECT(ScalarReal(0));
   rscalarinteger = PROTECT(ScalarInteger(0));
   setAttrib(rans, install("precision"), rscalarreal);
   setAttrib(rans, install("n_empty"), rscalarinteger);
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
  SET_STRING_ELT(rclass, 1, mkChar("MULTIPOLYGON"));
  SET_STRING_ELT(rclass, 2, mkChar("sfg"));
  /************************************************************************/
  /* Graticule  */
  /*  multipolygon/list with (lx) components */
  /*  each component i is a list of 1 component (EXT, no hole) */
  /*  this component EXT is a matrix 5x2 of vertices*/
  /************************************************************************/
  double coordx, coordy, coordxp1, coordyp1, minx=0.0, miny=0.0,
    maxx=0.0, maxy=0.0, *bbox2, *coord;
   for (i=0; i<(lx); i++) {
    rcolans  = PROTECT(allocVector(VECSXP, ly));
      for (j=0; j<(ly); j++) {
         rpolygext  = PROTECT(allocVector(VECSXP,1));
	/* grid */
	coordx = (double) i ;
	coordy = (double) j ;
	coordxp1 = (double) (i + 1);
	coordyp1 = (double) (j + 1);
	/* rescale and set */
	coordx = coordx*latt_const + new_minx;
	coordy = coordy*latt_const + new_miny;
	coordxp1 = coordxp1*latt_const + new_minx;
	coordyp1 = coordyp1*latt_const + new_miny;
	if (i==0 && j==0) {
	  minx=coordx;
	  miny=coordy;
	  maxx=coordxp1;
	  maxy=coordyp1;
	} else {
	  minx=fmin2(minx, coordx);
	  miny=fmin2(miny, coordy);
	  maxx=fmax2(maxx, coordxp1);
	  maxy=fmax2(maxy, coordyp1);
	}
	/* save  */
	rcoord  = PROTECT(allocMatrix(REALSXP, 5, 2));
    coord = REAL(rcoord);
	coord[0] = coordx;
	coord[1] = coordxp1;
	coord[2] = coordxp1;
	coord[3] = coordx;
	coord[4] = coordx;
	coord[5] = coordy;
	coord[6] = coordy;
	coord[7] = coordyp1;
	coord[8] = coordyp1;
	coord[9] = coordy;
	SET_VECTOR_ELT(rpolygext, 0, rcoord);
	SET_VECTOR_ELT(rcolans, j, rpolygext);
	UNPROTECT(2); /*rcoord + rpolygext*/
      }
	classgets(rcolans, rclass);
    SET_VECTOR_ELT(rans, i, rcolans);
	UNPROTECT(1); /*rcolans*/
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
   UNPROTECT(2);/*rscalarS*/
   UNPROTECT(6);
   UNPROTECT(8); /* class and attributes */
   return rans;
   /* class for ans */
}


/** \fn makefinalgraticule
 *  \brief Calculate the final graticule obtained in cartogramR
 * return a list of sfg POLYGON (the graticule)
 *
 * \param  rpadding : SEXP, the padding used in cartogram
 *         Determines space between map and boundary (default to 1.5)
 * \param  rLL : SEXP, the value of L in cartogram  (default is 512),
 *         must be a power of two (for fftw)
 * \param  rbbox: SEXP, the bounding box in cartogram
 * \param  rgridx : SEXP, The double SEXP matrix which contains the
 *           transformed grid (final grid after deformation ->
 *           will give  a discrete representation of the deformation on x axis)
 * \param  rgridy : SEXP, The double SEXP matrix which contains the
 *           transformed grid (final grid after deformation ->
 *           will give  a discrete representation of the deformation on y axis)
  *
 * \return rans : SEXP, The R list of sfg POLYGON
 *******************************************************************/

SEXP makefinalgraticule (SEXP rpadding, SEXP rLL, SEXP rbbox,
                         SEXP rgridx, SEXP rgridy)
{
  /*****************************************************************************/
  /* input and output from/to R */
  /*****************************************************************************/
  /* list output is an R object */
  SEXP  rans, rcoord, rclass, rclassans, rcrs, rclasscrs, rnames,
    rnamesbbox, rclassbbox, rbbox2, rcolans, rpolygext,
    rscalarreal, rscalarinteger;
  /*****************************************************************************/
  /* processing input  from R */
  /*****************************************************************************/
  /* double */
  rgridx = PROTECT(coerceVector(rgridx, REALSXP));
  rgridy = PROTECT(coerceVector(rgridy, REALSXP));
  rbbox = PROTECT(rbbox);
  rpadding = PROTECT(rpadding);
  double  *gridx, *gridy, *bbox, padding;
  gridx = REAL(rgridx);
  gridy = REAL(rgridy);
  bbox = REAL(rbbox);
  padding = REAL(rpadding)[0];
  /* integer : option(s)   */
  rLL = PROTECT(rLL);
  int LL;
  LL = INTEGER(rLL)[0];
  /************************************************************************/
  /* local variables */
  /************************************************************************/
  int i, j, lx, ly;
  /* number of rows in y_geom */
  int *caracmapi;
  double map_minx, map_maxx, map_miny, map_maxy, *caracmapd,
    latt_const, new_minx, new_miny, *xdisp, *ydisp;
  /* bounding box (SF order)*/
  map_minx= bbox[0];
  map_miny = bbox[1];
  map_maxx = bbox[2];
  map_maxy = bbox[3];
  /* Map  */
  caracmapd = (double *) R_alloc(3, sizeof(double));
  caracmapi = (int *) R_alloc(2, sizeof(int));
  caract_map(caracmapd, caracmapi, padding, LL, map_maxx, map_maxy, map_minx,map_miny);
  lx = caracmapi[0];
  ly = caracmapi[1];
  latt_const = caracmapd[0];
  new_minx = caracmapd[1];
  new_miny = caracmapd[2];
  /* Grid -> displacement vector*/
  xdisp = (double *) R_alloc(lx*ly, sizeof(double));
  ydisp = (double *) R_alloc(lx*ly, sizeof(double));
  /* error */
  int errorloc=0, valoption=0, *options;
  options = &valoption;
  /************************************************************************/
  /* result  R list of lx components */
  /************************************************************************/
  rans  = PROTECT(allocVector(VECSXP, lx));
   /************************************************************************/
   /* class and other attributes of rans */
   /************************************************************************/
   rclassans = PROTECT(allocVector(STRSXP, 2));
   SET_STRING_ELT(rclassans, 0, mkChar("sfc_MULTIPOLYGON"));
   SET_STRING_ELT(rclassans, 1, mkChar("sfc"));
   classgets(rans, rclassans);
   rscalarreal = PROTECT(ScalarReal(0));
   rscalarinteger = PROTECT(ScalarInteger(0));
   setAttrib(rans, install("precision"), rscalarreal);
   setAttrib(rans, install("n_empty"), rscalarinteger);
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
  SET_STRING_ELT(rclass, 1, mkChar("MULTIPOLYGON"));
  SET_STRING_ELT(rclass, 2, mkChar("sfg"));
   /************************************************************************/
  /* Grid -> displacement vector*/
  /************************************************************************/
  for (i=0; i<lx; i++) {
    for (j=0; j<ly; j++) {
      xdisp[i*ly + j] = gridx[i*ly + j] - i - 0.5;
      ydisp[i*ly + j] = gridy[i*ly + j] - j - 0.5;
    }
  }
  /************************************************************************/
  /* Graticule  */
  /*  multipolygon/list with (lx) components */
  /*  each component i is a list of 1 component (EXT, no hole) */
  /*  this component EXT is a matrix 5x2 of vertices*/
  /************************************************************************/
  double coordx00, coordx10, coordx11, coordx01,
    coordy00, coordy10, coordy11, coordy01,
    minxp, minyp, maxxp, maxyp, ctx, ctxp1, cty, ctyp1,
    minx=0.0, miny=0.0, maxx=0.0, maxy=0.0, *bbox2, *coord;
   for (i=0; i<(lx); i++) {
    rcolans  = PROTECT(allocVector(VECSXP, ly));
      for (j=0; j<(ly); j++) {
         rpolygext  = PROTECT(allocVector(VECSXP,1));
         /* grid */
         ctx = (double) i;
         cty = (double) j;
         ctxp1 = (double) (i + 1);
         ctyp1 = (double) (j + 1);
    /* interpolate */
          coordx00=interpol2(ctx, cty, xdisp, 'x',
                             options, &errorloc, lx, ly)  + ctx;
          coordx10=interpol2(ctxp1, cty, xdisp, 'x',
                             options, &errorloc, lx, ly) + ctxp1;
          coordx11=interpol2(ctxp1, ctyp1, xdisp, 'x',
                             options, &errorloc, lx, ly) + ctxp1;
          coordx01=interpol2(ctx, ctyp1, xdisp, 'x',
                             options, &errorloc, lx, ly) + ctx;
          coordy00=interpol2(ctx, cty, ydisp, 'y',
                             options, &errorloc, lx, ly) + cty;
          coordy10=interpol2(ctxp1, cty, ydisp, 'y',
                             options, &errorloc, lx, ly) + cty;
          coordy11=interpol2(ctxp1, ctyp1, ydisp, 'y',
                             options, &errorloc, lx, ly) + ctyp1;
          coordy01=interpol2(ctx, ctyp1, ydisp, 'y',
                             options, &errorloc, lx, ly) + ctyp1;
 	/* rescale and set */
          coordx00 = coordx00*latt_const + new_minx;
          coordx10 = coordx10*latt_const + new_minx;
          coordx11 = coordx11*latt_const + new_minx;
          coordx01 = coordx01*latt_const + new_minx;
          coordy00 = coordy00*latt_const + new_miny;
          coordy10 = coordy10*latt_const + new_miny;
          coordy11 = coordy11*latt_const + new_miny;
          coordy01 = coordy01*latt_const + new_miny;
         if (i==0 && j==0) {
            minx=coordx00;
            miny=coordy00;
            maxx=fmax2(coordx10,coordx11);
            maxy=fmax2(coordy01,coordy11);
          } else {
            minxp =fmin2(coordx00,coordx01);
            minyp =fmin2(coordy00,coordy10);
            maxxp = fmax2(coordx10,coordx11);
            maxyp = fmax2(coordy01,coordy11);
            minx=fmin2(minx, minxp);
            miny=fmin2(miny, minyp);
            maxx=fmax2(maxx, maxxp);
            maxy=fmax2(maxy, maxyp);
          }
          /* save  */
          rcoord  = PROTECT(allocMatrix(REALSXP, 5, 2));
          coord = REAL(rcoord);
          coord[0] = coordx00;
          coord[1] = coordx10;
          coord[2] = coordx11;
          coord[3] = coordx01;
          coord[4] = coordx00;
          coord[5] = coordy00;
          coord[6] = coordy10;
          coord[7] = coordy11;
          coord[8] = coordy01;
          coord[9] = coordy00;
          SET_VECTOR_ELT(rpolygext, 0, rcoord);
          SET_VECTOR_ELT(rcolans, j, rpolygext);
          UNPROTECT(2); /*rcoord + rpolygext*/
      }
      classgets(rcolans, rclass);
      SET_VECTOR_ELT(rans, i, rcolans);
      UNPROTECT(1); /*rcolans*/
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
   UNPROTECT(2);/*rscalarS*/
   UNPROTECT(8);
   UNPROTECT(8); /* class and attributes */
   return rans;
   /* class for ans */
}


