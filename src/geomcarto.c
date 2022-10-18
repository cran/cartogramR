/** Apply the deformation used to build a cartogram to a set of
 * simple geometry
 *
 *******************************************************************/

/******************************** Inclusions. **********************/

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> /* for fmin2 */
#include "geomcarto.h"
#include "caracmap.h"
#include "interpol2.h"


/********************************** Functions *******************************/

/** \fn geomcarto
 *  \brief Apply the deformation used to build a cartogram to a set of
 * simple geometry
 *
 * \param  rygeom : SEXP, The R list of simple geometry to which the deformation
 *          used to build a cartogram is applied
 * \param  rmultipoly : SEXP, The integer SEXP vector that define the type
 *          of geometry
 * \param  rgridx : SEXP, The double SEXP matrix which contains the
 *           transformed grid (final grid after deformation ->
 *           will give  a discrete representation of the deformation on x axis)
 * \param  rgridy : SEXP, The double SEXP matrix which contains the
 *           transformed grid (final grid after deformation ->
 *           will give  a discrete representation of the deformation on y axis)
 * \param  rpadding : SEXP, the padding used in cartogram
 *         Determines space between map and boundary (default to 1.5)
 * \param  rLL : SEXP, the value of L in cartogram  (default is 512),
 *         must be a power of two (for fftw)
 * \param  rbbox: SEXP, the bounding box in cartogram
 * \param  roptions Integer vector of options
 * \return rygeom : SEXP, The R list of simple geometry
 *******************************************************************/

SEXP geomcarto (SEXP rygeomd, SEXP rmultipoly, SEXP rgridx, SEXP rgridy,
                SEXP rpadding, SEXP rLL, SEXP rbboxd, SEXP roptions)
{
  /*****************************************************************************/
  /* input and output from/to R */
  /*****************************************************************************/
  /* list output is an R object */
  SEXP rygeom = PROTECT(duplicate(rygeomd));
  SEXP rbbox=PROTECT(duplicate(rbboxd)), rnamesbbox, rclassbbox;
  /*****************************************************************************/
  /* processing input  from R */
  /*****************************************************************************/
  /* double */
  rgridx = PROTECT(coerceVector(rgridx, REALSXP));
  rgridy = PROTECT(coerceVector(rgridy, REALSXP));
  rpadding = PROTECT(coerceVector(rpadding, REALSXP));
  double *gridx, *gridy, *bbox, padding;
  gridx = REAL(rgridx);
  gridy = REAL(rgridy);
  bbox = REAL(rbbox);
  padding = REAL(rpadding)[0];
  /* integer : option(s)   */
  rLL = PROTECT(coerceVector(rLL, INTSXP));
  roptions = PROTECT(coerceVector(roptions, INTSXP));
  int LL, *options, *multipoly;
  LL = INTEGER(rLL)[0];
  options = INTEGER(roptions);
  /* rmultipoly is a vector of 0/1 which indicates if row i is
   * - 0=point: list of length 2 with coord: (x , y)
   * - 1=multipoint: liste of length 2n with coords :
   *  (x1, x2, ..., xn, y1, y2, ..., yn)
   * - 2=line: liste of length 2n with coords : ls=(x1, x2, ..., xn, y1, y2, ..., yn)
   * - 3=multiline: list of lines (ls1, ls2, ...)
   * - 4=polygon: list  pol=(EXT, HOLE1, HOLE2, HOLE3,...) where each component
   *            is a matrix of coords
   * - 5=multipolygon: list of polygons: (pol1, pol22, ...) */
  rmultipoly = PROTECT(coerceVector(rmultipoly, INTSXP));
  multipoly  = INTEGER(rmultipoly);
  /************************************************************************/
  /* local variables */
  /************************************************************************/
  int i, j, jj, k, lx, ly, iter;
  /* number of rows in y_geom */
  int n_rows, *caracmapi;
  double map_minx, map_maxx, map_miny, map_maxy, *caracmapd, *xdisp, *ydisp;
  n_rows = length(rygeom);
  /* bounding box */
  map_minx= bbox[0];
  map_maxx = bbox[1];
  map_miny = bbox[2];
  map_maxy = bbox[3];
  /* error */
  int errorloc=0;

   /************************************************************************/
  /* Map  */
  /************************************************************************/
  caracmapd = (double *) R_alloc(3, sizeof(double));
  caracmapi = (int *) R_alloc(2, sizeof(int));
  caract_map(caracmapd, caracmapi, padding, LL, map_maxx, map_maxy, map_minx,map_miny);
  lx = caracmapi[0];
  ly = caracmapi[1];

  xdisp = (double *) R_alloc(lx*ly, sizeof(double));
  ydisp = (double *) R_alloc(lx*ly, sizeof(double));
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
  /* Result: rygeom  */
  /************************************************************************/
  /* Read polygon from R list to n_polycorn (and assign memory)*/
  /************************************************************************/
  int nbinpoly, nbpts,  nblistmulti;

  SEXP rlistcoord2,  rcoordvert2, rlistmulti2;
  double *coordvert2, coordx, coordy, coordxx, coordyy, minx=0.0, miny=0.0, maxx=0.0, maxy=0.0;
  iter=0;
  for (i=0; i<n_rows; i++) {
    rlistcoord2 = PROTECT(VECTOR_ELT(rygeom, i));
    nbinpoly = length(rlistcoord2);
    if (multipoly[i] < 3) {
      /* R vector of "points" (x1, x2, ..., xnbpts, y1, y2, ..., ynbpts) */
      coordvert2 = REAL(rlistcoord2);
      nbpts = (int) length(rlistcoord2)/2;
      for (j=0; j<nbpts; j++) {
	coordx =  coordvert2[j];
	coordy = coordvert2[j + nbpts];
	/* chgt coord -> */
	coordxx = (coordx - caracmapd[1])/caracmapd[0];
	coordyy = (coordy - caracmapd[2])/caracmapd[0];
	/* calculus */
	coordx = interpol2(coordxx, coordyy, xdisp, 'x', options,
			   &errorloc, lx, ly) + coordxx;
	if (errorloc>0) break;
	coordy = interpol2(coordxx, coordyy, ydisp, 'y', options,
			   &errorloc, lx, ly) + coordyy;
	if (errorloc>0) break;
	/* chgt coord <- */
	coordx = coordx * caracmapd[0] +caracmapd[1];
	coordy = coordy * caracmapd[0] +caracmapd[2];
	/* final */
	coordvert2[j] = coordx;
	coordvert2[j + nbpts] = coordy;
	if (iter==0) {
	  minx=coordx;
	  miny=coordy;
	  maxx=coordx;
	  maxy=coordy;
	  iter=1;
	} else {
	  minx=fmin2(minx, coordx);
	  miny=fmin2(miny, coordy);
	  maxx=fmax2(maxx, coordx);
	  maxy=fmax2(maxy, coordy);
	}
      }
    } else {
     if (multipoly[i]==3 || multipoly[i]==4 ) {
	/********************************************/
	/* simple sf polygons  + multilines*/
	/********************************************/
	for (j=0; j<nbinpoly; j++) {
	  /* R list of matrices/vector of coordinates */
	  rcoordvert2 = PROTECT(VECTOR_ELT(rlistcoord2, j));
	  coordvert2 = REAL(rcoordvert2);
	  nbpts = (int) length(rcoordvert2)/2;
	  /********************************************/
	  for (k=0; k<nbpts; k++) {
	    /* loop on coordinates */
	    coordx = coordvert2[k];
	    coordy = coordvert2[k + nbpts];
	    /* chgt coord -> */
	    coordxx = (coordx - caracmapd[1])/caracmapd[0];
	    coordyy = (coordy - caracmapd[2])/caracmapd[0];
	    /* calculus */
	    coordx = interpol2(coordxx, coordyy, xdisp, 'x',
			       options, &errorloc, lx, ly) + coordxx;
	    if (errorloc>0) break;
	    coordy = interpol2(coordxx, coordyy, ydisp, 'y',
			       options, &errorloc, lx, ly) + coordyy;
	    if (errorloc>0) break;
	    /* chgt coord <- */
	    coordy = coordy * caracmapd[0] +caracmapd[2];
	    coordx = coordx * caracmapd[0] +caracmapd[1];
	    /* final */
	    coordvert2[k]=coordx;
	    coordvert2[k + nbpts]=coordy;
	    	      if (iter==0) {
		minx=coordx;
		miny=coordy;
		maxx=coordx;
		maxy=coordy;
		iter=1;
	      } else {
		minx=fmin2(minx, coordx);
		miny=fmin2(miny, coordy);
		maxx=fmax2(maxx, coordx);
		maxy=fmax2(maxy, coordy);
	      }

	  }
	  UNPROTECT(1); /* rcoordvert2 */
	  if (errorloc>0)  break;
	}
      } else if  (multipoly[i]==5){
	/********************************************/
	/* multi sf polygons */
	/********************************************/
	for (j=0; j<nbinpoly; j++) {
	  /* R list of polygons */
	  rlistmulti2 = PROTECT(coerceVector(VECTOR_ELT(rlistcoord2, j), VECSXP));
	  nblistmulti = length(rlistmulti2);
	  for (jj=0; jj<nblistmulti; jj++) {
	    /********************************************/
	    /* get values it is an sf POLYGON*/
	    /********************************************/
	    /* R list of matrices/vector of coordinates */
	    rcoordvert2 = PROTECT(coerceVector(VECTOR_ELT(rlistmulti2, jj), REALSXP));
	    coordvert2 = REAL(rcoordvert2);
	    nbpts = (int) length(rcoordvert2)/2;
	    /********************************************/
	    /* values */
	    /********************************************/
	    for (k=0; k<nbpts; k++) {
	      /* loop on coordinates */
	      coordx = coordvert2[k];
	      coordy = coordvert2[k + nbpts];
	      /* chgt coord -> */
	      coordxx = (coordx - caracmapd[1])/caracmapd[0];
	      coordyy = (coordy - caracmapd[2])/caracmapd[0];
	      /* calculus */
	      coordx = interpol2(coordxx, coordyy, xdisp, 'x', options,
				 &errorloc, lx, ly) + coordxx;
	      if (errorloc>0) break;
	      coordy = interpol2(coordxx, coordyy, ydisp, 'y', options,
				 &errorloc, lx, ly) + coordyy;
	      if (errorloc>0) break;
	      /* chgt coord <- */
	      coordy = coordy * caracmapd[0] +caracmapd[2];
	      coordx = coordx * caracmapd[0] +caracmapd[1];
	      /* final */
	      coordvert2[k]=coordx;
	      coordvert2[k + nbpts]=coordy;
	      if (iter==0) {
		minx=coordx;
		miny=coordy;
		maxx=coordx;
		maxy=coordy;
		iter=1;
	      } else {
		minx=fmin2(minx, coordx);
		miny=fmin2(miny, coordy);
		maxx=fmax2(maxx, coordx);
		maxy=fmax2(maxy, coordy);
	      }
	    }
	    UNPROTECT(1); /* rcoordvert2 */
	    if (errorloc>0)  break;
	  }
	  UNPROTECT(1); /* rlistmulti2 */
	  if (errorloc>0)  break;
	}
      }
    }
    /*SET_VECTOR_ELT(rygeom, i, rlistcoord2);*/
    UNPROTECT(1); /* rlistcoord2 */
    if (errorloc>0)  break;
  }
  bbox[0]=minx;
  bbox[1]=miny;
  bbox[2]=maxx;
  bbox[3]=maxy;
  /* nases of components of the vector rbbox */
  rnamesbbox = PROTECT(allocVector(STRSXP, 4));
  SET_STRING_ELT(rnamesbbox, 0, mkChar("xmin"));
  SET_STRING_ELT(rnamesbbox, 1, mkChar("ymin"));
  SET_STRING_ELT(rnamesbbox, 2, mkChar("xmax"));
  SET_STRING_ELT(rnamesbbox, 3, mkChar("ymax"));
  /* assign names to vector rbbox */
  setAttrib(rbbox, R_NamesSymbol, rnamesbbox);
  /* class bbox */
  rclassbbox = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(rclassbbox, 0, mkChar("bbox"));
  classgets(rbbox, rclassbbox);
  /* set bbox attribute to rygeom */
  setAttrib(rygeom, install("bbox"), rbbox);
  UNPROTECT(10); /* rygeom + rbbox + rgridx + rgridy + rpadding + rLL
                  + roptions + rmultipoly + rnamesbbox +
                  rclassbbox*/
  if (errorloc>0)  error("interpolation error");
  return rygeom;
}

