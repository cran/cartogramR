/** \brief Program for the construction of flow based cartograms.
 *
 * Program for the construction of cartograms. The program needs polygon
 * coordinates and a value for the (relative) target area for each region as
 * input, calculates the new coordinates.

 * If you use output created by this program please acknowledge the use of
 * this code and its first publication in:
 * Michael T. Gastner, Vivien Seguy and Pratyush More, "A fast flow-based
 * algorithm for creating density-equalizing map  projections".
 ******************************---------------------------- */



/* -------****************** Inclusions. *****------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Rmath.h> /* for fmin2 */
#include "cartogram.h"

/* -------------------------------- Main. ------------------------------- */
/** \fn cartogramR
 * \param  rcentroidx Real vector of x-coordinates of centroids
 * \param  rcentroidy Real vector of y-coordinates of centroids
 * \param  rygeomd  The R list of polygons (each component i is either a
 *                     multipolygon or a polygon) of length n_rows
 * \param rvarregion
 * \param rdimpoly
 * \param rbboxd Real vector of bounding box
 * \param rparamsdouble Real vector of real parameters
 * \param rparamsint Integer vector of integer parameters
 * \param roptions Integer vector of options
 * \return rans a R list
 * [[1]] rygeom the sf object (aka the cartogram)
 * [[2]] original area
 * [[3]] final area
 * [[4]] final x-coordinates of centroids
 * [[5]] final y-coordinates of centroids
 * [[6]] grid x matrix (optionnal)
 * [[7]] grid y matrix (optionnal)
 ******************************---------------------------- */

/*------------------------ Global variables. ----------------------*/
int L;
double MAX_PERMITTED_AREA_ERROR, MIN_POP_FAC, PADDING, BLUR_WIDTH;
/* Variables for map. */
double *area_err, *cart_area, map_maxx, map_maxy, map_minx, map_miny,
  *target_area, *bbox, *coordvertices;
int max_id, min_id, n_poly, *n_polyinreg, n_reg, **polyinreg;
POINT **cartcorn, **origcorn, **polycorn, *proj, *proj2, *proj3, *projinit, *projtmp;
size_t projsize;
/* Variables for digitizing the density. */
double *rho_ft, *rho_init;
fftw_plan plan_fwd;
int lx, ly;
/*----------------- end of Global variables. ----------------------*/


SEXP cartogramR (SEXP rcentroidx, SEXP rcentroidy, SEXP rygeomd,
    SEXP rvarregion, SEXP rnb_polyinreg, SEXP rn_polycorn,
		 SEXP rdimpoly, SEXP rbboxd, SEXP rparamsdouble, SEXP rparamsint,
		 SEXP roptions, SEXP rmultipoly)
		/* (int argc, char* argv[]) */
{
  /*---------------------------------------*/
  /* input and output from/to R */
  /*---------------------------------------*/
  /*   - the initial centroid coordinates (of the regions): input (rcentroidx) */
  /*   - the transformed centroid coordinates: ouput (rcentroidx2) */
  SEXP rcentroidx2,  rcentroidy2;
  rcentroidx2 = PROTECT(duplicate(rcentroidx));
  rcentroidy2 = PROTECT(duplicate(rcentroidy));
  double *centroidx, *centroidy;
  centroidx = REAL(rcentroidx2);
  centroidy = REAL(rcentroidy2);
  /* list output */
  /* rygeom/rygeomd is an R list that contains  */
  /*   - the initial polygons (the regions): input (rygeomd) */
  /*   - the result/the cartogram: output (rygeom) */
  SEXP rygeom = PROTECT(duplicate(rygeomd));
  SEXP rbbox=PROTECT(duplicate(rbboxd)), rnamesbbox, rclassbbox;
  /* double input/output */
  /* rcentroidx and rcentroidy contains */
  /*   - initial regions centroids coordinates: input */
  /*--------------------------------------------*/
  /* input  from R */
  /*--------------------------------------------*/
  /* double */
  rvarregion = PROTECT(rvarregion);
  rparamsdouble = PROTECT(rparamsdouble);
  double *varregion, *bbox, *paramsdouble, *final_area, *original_area;
  varregion = REAL(rvarregion);
  bbox = REAL(rbbox);
  paramsdouble = REAL(rparamsdouble);
  /* integer   */
  rnb_polyinreg = PROTECT(rnb_polyinreg);
  rn_polycorn = PROTECT(rn_polycorn);
  rdimpoly = PROTECT(rdimpoly);
  rparamsint = PROTECT(rparamsint);
  roptions = PROTECT(roptions);
  rmultipoly = PROTECT(rmultipoly);
  int *nb_polyinreg, *n_polycorn, *dimpoly, *paramsint, *options, *multipoly;
  nb_polyinreg = INTEGER(rnb_polyinreg);
  n_polycorn = INTEGER(rn_polycorn);
  dimpoly = INTEGER(rdimpoly);
  paramsint = INTEGER(rparamsint);
  options = INTEGER(roptions);
  multipoly = INTEGER(rmultipoly);
  /* error and max iterations*/
  int errorloc=0, maxit;
  /* number of rows in y_geom and total number of polgons*/
  int n_rows, nbpoly;
  n_rows = length(rygeom);
  nbpoly = length(rn_polycorn);

  /* set parameters to their given value*/
  MAX_PERMITTED_AREA_ERROR=paramsdouble[0];
  MIN_POP_FAC=paramsdouble[1];
  PADDING=paramsdouble[2];
  BLUR_WIDTH=paramsdouble[3];
  L=paramsint[0];
  maxit=paramsint[1];

  /* number of polygons */
  n_poly = dimpoly[0];
  /* number of regions and id of region */
  n_reg = dimpoly[1];
  min_id = dimpoly[2];
  max_id = dimpoly[3];
  /* bounding box */
  map_minx= bbox[0];
  map_maxx = bbox[1];
  map_miny = bbox[2];
  map_maxy = bbox[3];


  /*-----------------------------------------------------------*/
  /* local variables */
  /*-----------------------------------------------------------*/
  Rboolean diff, gridexport, absrel;
  double  curcrit, cart_tot_area, correction_factor=1.0, init_tot_area, scale_map;
  int i, integration, j;
  /* process options */
  if (options[1]>0)   diff = TRUE; else diff = FALSE;
  if (options[2]>0)   gridexport = TRUE; else gridexport = FALSE;
  if (options[3]>0)   absrel= TRUE; else absrel = FALSE;

  /*------------------------------------------------------------*/
  /* Result: rans  */
  /* [[1]] rygeom the sf object */
  /* [[2]] original area */
  /* [[3]] final area */
  /* [[4]] final x-coordinates of centroids */
  /* [[5]] final y-coordinates of centroids */
  /* [[6]] grid x matrix*/
  /* [[7]] grid y matrix */
  /*------------------------------------------------------------*/
  SEXP rans;
  if (gridexport) {
    rans = PROTECT(allocVector(VECSXP, 7));
  } else {
    rans = PROTECT(allocVector(VECSXP, 5));
  }
  SEXP rfinal_area, roriginal_area;
  roriginal_area = PROTECT(allocVector(REALSXP, n_reg));
  rfinal_area = PROTECT(allocVector(REALSXP, n_reg));
  final_area = REAL(rfinal_area);
  original_area = REAL(roriginal_area);

  /*-------------------------------------------------------------*/
  /* Read polygon from R list to n_polycorn (and assign memory)*/
  /*-------------------------------------------------------------*/
  int jj,k, nbinpoly, nbvert, nblistmulti, ctrpoly=0;

  /* allocation step 1*/
  polycorn = (POINT**) malloc(nbpoly * sizeof(POINT*));
  /* allocation polygon 1 */
  polycorn[ctrpoly] = (POINT*) malloc(n_polycorn[ctrpoly] * sizeof(POINT));

  SEXP rlistcoord,  rcoordvert, rlistmulti;
  double *coordvert;
  for (i=0; i<n_rows; i++) {
    rlistcoord = PROTECT(VECTOR_ELT(rygeomd, i));
    nbinpoly = length(rlistcoord);
    for (j=0; j<nbinpoly; j++) {
      if (multipoly[i]==0) {
        /* simple sf polygons */
        /* get values directly it is an sf POLYGON*/
        rcoordvert = PROTECT(VECTOR_ELT(rlistcoord, j));
        coordvert = REAL(rcoordvert);
        nbvert = n_polycorn[ctrpoly];
        /* setting values */
        for (k=0; k<nbvert; k++) {
          polycorn[ctrpoly][k].x=coordvert[k];
          polycorn[ctrpoly][k].y=coordvert[nbvert+k];
        }
        UNPROTECT(1); /* rcoordvert */
        ctrpoly++;
        /* allocate next polycorn[ctrpoly] except for the last */
        if (ctrpoly<nbpoly) {
          polycorn[ctrpoly] = (POINT*) malloc(n_polycorn[ctrpoly] * sizeof(POINT));
        }
      } else {
        /* multi sf polygons */
        rlistmulti = PROTECT(VECTOR_ELT(rlistcoord, j));
        nblistmulti = length(rlistmulti);
        for (jj=0; jj<nblistmulti; jj++) {
          /* get values it is an sf POLYGON*/
          rcoordvert = PROTECT(VECTOR_ELT(rlistmulti, jj));
          coordvert = REAL(rcoordvert);
          nbvert = n_polycorn[ctrpoly];
          /* setting values */
          for (k=0; k<nbvert; k++) {
            polycorn[ctrpoly][k].x=coordvert[k];
            polycorn[ctrpoly][k].y=coordvert[nbvert+k];
          }
          UNPROTECT(1); /* rcoordvert */
          ctrpoly++;
          /* allocate next polycorn[ctrpoly] except for the last */
          if (ctrpoly<nbpoly) {
            polycorn[ctrpoly] = (POINT*) malloc(n_polycorn[ctrpoly] * sizeof(POINT));
          }
        }
        UNPROTECT(1); /* rlistmulti */
      }
    }
    UNPROTECT(1); /* rlistcoord */
  }


  /*-------------------------------------------------------------------*/
  /* Fill the lx-times-ly grid with   */
  /* density and print a map. rho_ft[] will be filled with the Fourier */
  /* transform of the initial density.                                 */
  /*-------------------------------------------------------------------*/
  fill_with_density1(centroidx, centroidy, n_polycorn, varregion,
                     nb_polyinreg, options, original_area);
  SET_VECTOR_ELT(rans, 1, roriginal_area);
  /*-------------------------------------------------------------------*/
  /* Allocate memory for the projected positions. */
  /*-------------------------------------------------------------------*/
  projsize = lx * ly * sizeof(POINT);
  proj = (POINT*) malloc(projsize);
  proj3 = (POINT*) malloc(projsize);
  cartcorn = (POINT**) malloc(n_poly * sizeof(POINT*));
  for (i=0; i<n_poly; i++)
    cartcorn[i] = (POINT*) malloc(n_polycorn[i] * sizeof(POINT));
  area_err = (double*) malloc(n_reg * sizeof(double));
  cart_area = (double*) malloc(n_reg * sizeof(double));
  /*-------------------------------------------------------------------*/
  /* Are we already at final point ?  */
  /*-------------------------------------------------------------------*/
  if (absrel) {
    curcrit=max_area_err(area_err, cart_area, n_polycorn, polycorn, &init_tot_area);
  } else {
    /* the first time we need to rescale the abserror criterion on cartogram scale */
    scale_map = scale_map_factor();
    MAX_PERMITTED_AREA_ERROR /= scale_map;
    curcrit = max_absarea_err(area_err, cart_area, n_polycorn, polycorn, &init_tot_area);
  }
  if (curcrit <= MAX_PERMITTED_AREA_ERROR) {
    /* No need to compute anything. */
    /* -------------------------------------------------------- */
    /* Export the grid if needed */
    /* -------------------------------------------------------- */
    if (gridexport) {
      /* export grid */
      SEXP ransx = PROTECT(allocMatrix(REALSXP, ly, lx));
      SEXP ransy = PROTECT(allocMatrix(REALSXP, ly, lx));
      double *ansx = REAL(ransx);
      double *ansy = REAL(ransy);
      for (i=0; i<lx; i++) {
        for (j=0; j<ly; j++) {
          /* -------**** Rescale grid coordinate *******------- */
          ansx[i*ly + j] = correction_factor * (proj3[i*ly + j].x - 0.5*lx) + 0.5*lx;
          ansy[i*ly + j] = correction_factor * (proj3[i*ly + j].y - 0.5*ly) + 0.5*ly;
          /* please stay on grid (TODO maybe a warning if outside of the grid)*/
          if (ansx[i*ly + j] <0)  ansx[i*ly + j]=0;
          if (ansy[i*ly + j] <0)  ansy[i*ly + j]=0;
          if (ansx[i*ly + j] >lx) ansx[i*ly + j]=lx;
          if (ansy[i*ly + j] >ly) ansy[i*ly + j]=ly;
        }
      }
      SET_VECTOR_ELT(rans, 5, ransx);
      SET_VECTOR_ELT(rans, 6, ransy);
    }
    /*--------------------------------------------------------------
     * export polycorn (ugly global variable) to rygeom an R sf list
     * four  steps
     ----------------------------------------------------------------*/
    /* first return to the original scale */
    inv_rescale_map (centroidx, centroidy, n_polycorn, options);
    /* second: export centroid */
    SET_VECTOR_ELT(rans, 3, rcentroidx2);
    SET_VECTOR_ELT(rans, 4, rcentroidy2);

    /* third: export the cartogram in the R object sf geometry */
    ctrpoly=0;
    double minx, maxx, miny, maxy;
    for (i=0; i<n_rows; i++) {
      rlistcoord = PROTECT(VECTOR_ELT(rygeom, i));
      nbinpoly = length(rlistcoord);
      for (j=0; j<nbinpoly; j++) {
        if (multipoly[i]==0) {
          /* simple sf polygons */
          /* R setup */
          rcoordvert = PROTECT(VECTOR_ELT(rlistcoord, j));
          coordvert = REAL(rcoordvert);
          nbvert = n_polycorn[ctrpoly];
          /* setting values */
          for (k=0; k<nbvert; k++) {
            coordvert[k]=polycorn[ctrpoly][k].x;
            coordvert[nbvert+k]=polycorn[ctrpoly][k].y;
            if ((k==0)&&(i==0)) {
              minx=polycorn[ctrpoly][k].x;
              miny=polycorn[ctrpoly][k].y;
              maxx=polycorn[ctrpoly][k].x;
              maxy=polycorn[ctrpoly][k].y;
            } else {
              minx=fmin2(minx, polycorn[ctrpoly][k].x);
              miny=fmin2(miny, polycorn[ctrpoly][k].y);
              maxx=fmax2(maxx, polycorn[ctrpoly][k].x);
              maxy=fmax2(maxy, polycorn[ctrpoly][k].y);
            }
          }
          UNPROTECT(1); /* rcoordvert */
          ctrpoly++;
        } else {
          /* multi sf polygons */
          rlistmulti = PROTECT(VECTOR_ELT(rlistcoord, j));
          nblistmulti = length(rlistmulti);
          for (jj=0; jj<nblistmulti; jj++) {
            /* R setup  */
            rcoordvert = PROTECT(VECTOR_ELT(rlistmulti, jj));
            coordvert = REAL(rcoordvert);
            nbvert = n_polycorn[ctrpoly];
            /* setting values */
            for (k=0; k<nbvert; k++) {
              coordvert[k]=polycorn[ctrpoly][k].x;
              coordvert[nbvert+k]=polycorn[ctrpoly][k].y;
              if ((k==0)&&(i==0)) {
                minx=polycorn[ctrpoly][k].x;
                miny=polycorn[ctrpoly][k].y;
                maxx=polycorn[ctrpoly][k].x;
                maxy=polycorn[ctrpoly][k].y;
              } else {
                minx=fmin2(minx, polycorn[ctrpoly][k].x);
                miny=fmin2(miny, polycorn[ctrpoly][k].y);
                maxx=fmax2(maxx, polycorn[ctrpoly][k].x);
                maxy=fmax2(maxy, polycorn[ctrpoly][k].y);
              }
            }
            UNPROTECT(1); /* rcoordvert */
            ctrpoly++;
          }
          UNPROTECT(1); /* rlistmulti */
        }
      }
      UNPROTECT(1); /* rlistcoord */
    }
    bbox[0]=minx;
    bbox[1]=miny;
    bbox[2]=maxx;
    bbox[3]=maxy;
    /* names of components of the vector rbbox */
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
    SET_VECTOR_ELT(rans, 0, rygeom);
    /* four:  Areas (on initial scale) */
    for (i=0; i<n_reg; i++) {
      final_area[i]=0;
      for (j=0; j<n_polyinreg[i]; j++)
        final_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
                                      polycorn[polyinreg[i][j]]);
    }
    SET_VECTOR_ELT(rans, 2, rfinal_area);
    /*------------------------- end export------------------------*/
    /*------------------------ Free memory. ----------------------*/
    UNPROTECT(12);
    UNPROTECT(3); /* rygeom + rcentroidx2 + rcentroidy2 */
    FREEG1 ;
    Rprintf("max. rel. area error: %f\n", curcrit);
    warning("at the init of the program, we are already at the end, exiting...");
  } else {
    /* -------------------------------------------------------- */
    /*  the program starts */
    /* -------------------------------------------------------- */
    /* proj[i*ly+j] will store the current position of the point that started  */
    /* at (i+0.5, j+0.5).                                                      */

    projinit = (POINT*) malloc(projsize);
    for (i = 0; i < lx; i++) {
      for (j = 0; j < ly; j++) {
        projinit[i * ly + j].x = i + 0.5;
        projinit[i * ly + j].y = j + 0.5;
      }
    }
    memcpy(proj, projinit, projsize);
    memcpy(proj3, projinit, projsize);
    for (i=0; i<lx; i++)
      for (j=0; j<ly; j++) {
        proj[i*ly + j].x = i + 0.5;
        proj[i*ly + j].y = j + 0.5;
        proj3[i*ly + j].x = i + 0.5;
        proj3[i*ly + j].y = j + 0.5;
      }

    /* -------------------------------------------------------- */
    /* ------- First integration of the equations of motion.*/
    /* -------------------------------------------------------- */
    if (options[0]>0) Rprintf("Starting iteration 1 ");
    if (options[0]>1) Rprintf("\n");
    if (!diff)
      ffb_integrate(options, errorloc);
    else
      diff_integrate(options, errorloc);
    if (errorloc>0) {
      /* ------- free memory on error ------- */
      FREEG1 ;
      UNPROTECT(12);
      UNPROTECT(3); /* rygeom + rcentroidx2 + rcentroidy2*/
      error("error in ffb_integrate/diff_integrate");
      return rans;
    }
    project(centroidx, centroidy, FALSE, options, errorloc, n_polycorn, gridexport);  /* FALSE because we do not need to project the graticule. */
    if (errorloc>0) {
      /* ------- free memory on error ------- */
      FREEG1 ;
      UNPROTECT(12);
      UNPROTECT(3); /* rygeom + rcentroidx2 + rcentroidy2*/
      error("error in project");
      return rans;
    }

    if (absrel) {
      curcrit = max_area_err(area_err, cart_area, n_polycorn,
                             cartcorn, &cart_tot_area);
      if (options[0]>0) Rprintf("max. rel. area error: %f\n", curcrit);
    } else {
      curcrit = max_absarea_err(area_err, cart_area, n_polycorn,
                                cartcorn, &cart_tot_area);
      if (options[0]>0) Rprintf("max. abs. area error: %f\n", curcrit);
    }

    /* -------------------------------------------------------- */
    /* Additional integrations to come closer to target areas.*/
    /* -------------------------------------------------------- */
    proj2 = (POINT*) malloc(projsize);
    integration = 1;
    while (curcrit > MAX_PERMITTED_AREA_ERROR && integration<=maxit) {

      fill_with_density2(n_polycorn);

      /* Copy the current graticule before resetting. We will construct the    */
      /* final graticule by interpolating proj2 on the basis of proj.          */

      projtmp = proj2;
      proj2 = proj;
      proj = projtmp;
      memcpy(proj, projinit, projsize);
 /* for (i=0; i<lx; i++) */
 /*        for (j=0; j<ly; j++) { */
 /*          proj2[i*ly + j].x = proj[i*ly + j].x; */
 /*          proj2[i*ly + j].y = proj[i*ly + j].y; */
 /*        } */
 /*      for (i=0; i<lx; i++) */
 /*        for (j=0; j<ly; j++) { */
 /*          proj[i*ly + j].x = i + 0.5; */
 /*          proj[i*ly + j].y = j + 0.5; */
 /*        } */
      integration++;
      if (options[0]>0) Rprintf("Starting iteration %d ", integration);
      if (options[0]>1) Rprintf("\n");
      if (!diff)
        ffb_integrate(options, errorloc);
      else
        diff_integrate(options, errorloc);
      if (errorloc>0) {
        /* ------- free memory on error ------- */
        FREEG1;
        free(proj2);
        UNPROTECT(12);
        UNPROTECT(3); /* rygeom + rcentroidx2 + rcentroidy2*/
        error("error in ffb_integrate/diff_integrate");
        return rans;
      }
      project(centroidx, centroidy, TRUE, options, errorloc,
        n_polycorn, gridexport);
       /* TRUE because we need to project the graticule too. */
      if (errorloc>0) {
        /* ------- free memory on error ------- */
        FREEG1;
        free(proj2);
        UNPROTECT(12);
        UNPROTECT(3); /* rygeom + rcentroidx2 + rcentroidy2*/
        error("error in project");
        return rans;
      }
      /* Overwrite proj with proj2. */
      for (i=0; i<lx; i++)
        for (j=0; j<ly; j++) {
          proj[i*ly + j].x = proj2[i*ly + j].x;
          proj[i*ly + j].y = proj2[i*ly + j].y;
        }
    projtmp = proj;
    proj = proj2;
    proj2 = projtmp;
        if (absrel) {
        curcrit = max_area_err(area_err, cart_area, n_polycorn,
                               cartcorn, &cart_tot_area);
        if (options[0]>0) Rprintf("max. rel. area error: %f\n", curcrit);
      } else {
        curcrit = max_absarea_err(area_err, cart_area, n_polycorn,
                                  cartcorn, &cart_tot_area);
        if (options[0]>0) Rprintf("max. abs. area error: %f\n", curcrit);
      }
    }
    /* Rescale all areas to perfectly match the total area before the          */
    /* integrations started.                                                   */
    correction_factor = sqrt(init_tot_area / cart_tot_area);
    if (options[0]>1) Rprintf("Rescale all areas to perfectly match the total area before the integrations started:\n   correction_factor = %f\n",sqrt(init_tot_area / cart_tot_area));
    for (i=0; i<n_poly; i++)
      for (j=0; j<n_polycorn[i]; j++) {
        cartcorn[i][j].x =
          correction_factor * (cartcorn[i][j].x - 0.5*lx) + 0.5*lx;
        cartcorn[i][j].y =
          correction_factor * (cartcorn[i][j].y - 0.5*ly) + 0.5*ly;
      }
    /* -------**** Rescale centroid coordinate *******------- */
    for (i=0; i<n_reg; i++) {
      centroidx[i] = correction_factor * (centroidx[i] - 0.5*lx) + 0.5*lx;
      centroidy[i] = correction_factor * (centroidy[i] - 0.5*ly) + 0.5*ly;
    }
    /* Run max_area_err() once more so that we print correct absolute areas in */
    /* output_error().                                                      */
      if (absrel) {
        curcrit = max_area_err(area_err, cart_area, n_polycorn,
                               cartcorn, &cart_tot_area);
        if (options[0]>1) Rprintf("after final correction, max. rel. area error: %f\n", curcrit);
       } else {
        curcrit = max_absarea_err(area_err, cart_area, n_polycorn,
                                cartcorn, &cart_tot_area);
        if (options[0]>1) Rprintf("after final correction, max. abs. area error: %f\n", curcrit);
      }
    if (curcrit > MAX_PERMITTED_AREA_ERROR) Rprintf("WARNING criterion: %f > Objective: %f\n Increase maxit or decrease Objective\n", curcrit, MAX_PERMITTED_AREA_ERROR);
    /* -------------------------------------------------------- */
    /* Export the grid if needed */
    /* -------------------------------------------------------- */
    if (gridexport) {
      /* export grid */
      SEXP ransx = PROTECT(allocMatrix(REALSXP, ly, lx));
      SEXP ransy = PROTECT(allocMatrix(REALSXP, ly, lx));
      double *ansx = REAL(ransx);
      double *ansy = REAL(ransy);
      for (i=0; i<lx; i++) {
        for (j=0; j<ly; j++) {
          /* -------**** Rescale grid coordinate *******------- */
          ansx[i*ly + j] = correction_factor * (proj3[i*ly + j].x - 0.5*lx) + 0.5*lx;
          ansy[i*ly + j] = correction_factor * (proj3[i*ly + j].y - 0.5*ly) + 0.5*ly;
          /* please stay on grid (TODO maybe a warning if outside of the grid)*/
          if (ansx[i*ly + j] <0)  ansx[i*ly + j]=0;
          if (ansy[i*ly + j] <0)  ansy[i*ly + j]=0;
          if (ansx[i*ly + j] >lx) ansx[i*ly + j]=lx;
          if (ansy[i*ly + j] >ly) ansy[i*ly + j]=ly;
        }
      }
      SET_VECTOR_ELT(rans, 5, ransx);
      SET_VECTOR_ELT(rans, 6, ransy);
    }
    /* -------------------------------------------------------- */
    /* export cartcorn (ugly global variable) to rygeom an R sf list */
    /* four  steps                                              */
    /* -------------------------------------------------------- */
    /* first return to the original scale */
    inv_rescale_map (centroidx, centroidy, n_polycorn, options);
    /* second: export centroid */
    SET_VECTOR_ELT(rans, 3, rcentroidx2);
    SET_VECTOR_ELT(rans, 4, rcentroidy2);

    /* third: export the cartogram in the R object sf geometry */
    ctrpoly=0;
    double minx, maxx, miny, maxy;
    for (i=0; i<n_rows; i++) {
      rlistcoord = PROTECT(VECTOR_ELT(rygeom, i));
      nbinpoly = length(rlistcoord);
      for (j=0; j<nbinpoly; j++) {
        if (multipoly[i]==0) {
          /* simple sf polygons */
          /* R setup */
          rcoordvert = PROTECT(VECTOR_ELT(rlistcoord, j));
          coordvert = REAL(rcoordvert);
          nbvert = n_polycorn[ctrpoly];
          /* setting values */
          for (k=0; k<nbvert; k++) {
            coordvert[k]=cartcorn[ctrpoly][k].x;
            coordvert[nbvert+k]=cartcorn[ctrpoly][k].y;
            if ((k==0)&&(i==0)) {
              minx=cartcorn[ctrpoly][k].x;
              miny=cartcorn[ctrpoly][k].y;
              maxx=cartcorn[ctrpoly][k].x;
              maxy=cartcorn[ctrpoly][k].y;
            } else {
              minx=fmin2(minx, cartcorn[ctrpoly][k].x);
              miny=fmin2(miny, cartcorn[ctrpoly][k].y);
              maxx=fmax2(maxx, cartcorn[ctrpoly][k].x);
              maxy=fmax2(maxy, cartcorn[ctrpoly][k].y);
            }
          }
          UNPROTECT(1); /* rcoordvert */
          ctrpoly++;
        } else {
          /* multi sf polygons */
          rlistmulti = PROTECT(VECTOR_ELT(rlistcoord, j));
          nblistmulti = length(rlistmulti);
          for (jj=0; jj<nblistmulti; jj++) {
            /* R setup  */
            rcoordvert = PROTECT(VECTOR_ELT(rlistmulti, jj));
            coordvert = REAL(rcoordvert);
            nbvert = n_polycorn[ctrpoly];
            /* setting values */
            for (k=0; k<nbvert; k++) {
              coordvert[k]=cartcorn[ctrpoly][k].x;
              coordvert[nbvert+k]=cartcorn[ctrpoly][k].y;
              if ((k==0)&&(i==0)) {
                minx=cartcorn[ctrpoly][k].x;
                miny=cartcorn[ctrpoly][k].y;
                maxx=cartcorn[ctrpoly][k].x;
                maxy=cartcorn[ctrpoly][k].y;
              } else {
                minx=fmin2(minx, cartcorn[ctrpoly][k].x);
                miny=fmin2(miny, cartcorn[ctrpoly][k].y);
                maxx=fmax2(maxx, cartcorn[ctrpoly][k].x);
                maxy=fmax2(maxy, cartcorn[ctrpoly][k].y);
              }
            }
            UNPROTECT(1); /* rcoordvert */
            ctrpoly++;
          }
          UNPROTECT(1); /* rlistmulti */
        }
      }
      UNPROTECT(1); /* rlistcoord */
    }
    bbox[0]=minx;
    bbox[1]=miny;
    bbox[2]=maxx;
    bbox[3]=maxy;
    /* names of components of the vector rbbox */
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
    /* ygeom is the first component */
    SET_VECTOR_ELT(rans, 0, rygeom);

    /* four:  Areas (on initial scale) */
    for (i=0; i<n_reg; i++) {
      final_area[i]=0;
      for (j=0; j<n_polyinreg[i]; j++)
        final_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
                                      cartcorn[polyinreg[i][j]]);
    }
    SET_VECTOR_ELT(rans, 2, rfinal_area);
    /* -------------------------- end export------------------------- */
    /* ------------------------- Free memory. ----------------------- */
    FREEG1;
    free(proj2);
    if (gridexport) {
      UNPROTECT(16);
      UNPROTECT(3); /* rygeom + rcentroidx2 + rcentroidy2 */
      return rans ;
    } else {
      UNPROTECT(14);
      UNPROTECT(3); /* rygeom + rcentroidx2 + rcentroidy2 */
      return rans ;
    }
  }
  return rans ;
}
