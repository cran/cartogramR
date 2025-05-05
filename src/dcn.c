/** Calculate the grid used in cartogramR for a given L and padding
 * return a list of sfg POINT (the grid)
 *
 *******************************************************************/

/******************************** Inclusions. **********************/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "myomp.h"
#if defined(R_SIGACTION) ||  defined(R_SIGWIN)
#include <signal.h>
#endif
#include "dcn.h"
/**************************** Function prototypes. ******************/
void areacdg(double *x,double *y, int *l1, int *l2,
             double *cdgpx, double *cdgpy, double *areasp,
             int nmp, int *debuts, int *debutsp, int nthreads);
void initpolygons(double *x,double *y, int *l1, int *l2,
          double *cdgpx, double *cdgpy, double *areasp, double *countp, int *regofpol,
          double *count, int nmp, int *debuts, int *debutsp, int nthreads);
void exitpolygons(double *cdgpx, double *cdgpy, double *areasp, double *areasinit,
             int nmp, int np, int *debutsp);
void maxcritbyregion(double *areasp, double *areasObj, double *areaspOld,
                       int *regofpol, int np, int absrel, double *relerr,
                       double *reltol);
/********************************** Functions ***********************/
#if defined(R_SIGACTION) || defined(R_SIGWIN)
static  volatile sig_atomic_t keep_runningg = 1;
static void intHandlerg(int _) {
   (void)_;
    keep_runningg = 0;
}
#endif

/** \fn dcn
 *  \brief Dougenik et al improved algorithm
 *
 * \param  rygeomd: The R list of polygons (each component i is a
 *                     multipolygon) of length nmp
 * \param  rx: vector double [1:n] X-axis coordinates of polygons vertices
 * \param  ry: vector double [1:n] Y-axis coordinates of polygons vertices
 * \param  rcount: vector double [1:nmp] of variable count or density
 *                 in each region
 * \param  rcount: vector double [1:nmp] of variable count or density
 * \param  rl1: vector int [1:n] L1 column of sf::st_coordinates gives
 *         the integer 1,2, 3 got in the following list
 *         POL=(EXT1 HOLE2 HOLE3...)
 * \param  rml2: vector int  [1:n]; the ith coordinate gives the
 *          polygon number  (starting at 1) of vertice i.
 * \param rdup: vector [1:n] if the coordinate k is a vertice with
 *        duplicates (k1, k2 ...) it gives the next indice k1
 *        if no duplicate exists it is 0
 * \param rlast: vector int [1:n] if the coordinate k is a vertice with
 *        duplicates (k1, k2 ..., kk) it gives the last indice kk
 *        of the list of duplicates
 * \param  rparamint vector int
 *           - [0] = maxit
 *           - [1] = number of polygons (not multipolygon which define the region)
 * \param  rparamdouble vector double
 *           - [0] = abs/relerror criterion to stop iterations:
 *                   max_{polygons}abs(area(i) - objective_area)
 *                   or
 *                   max_{polygons}abs((area(i) - objective_area)/obj)
 *           - [1] = abs/reltol criterion to stop iterations
 *                   max_{polygons} abs(area(i) - area(i-1))
 *                   or
 *                   max_{polygons} abs((area(i) - area(i-1))/area(i-1))
 *           - [2] = mp: If a region contains exactly zero population, it will be
 *                   replaced by mp times the smallest positive population
 *                   in any region
 * \param  roptions vector int
 *           - [0] = verbose
 *           - [1] = absolute or relative Error (if 0 absolute else relative)
 *           - [2] = number of threads in cdg calculus (default to 2)
 *           - [3] = number of times algorithm is allowed to increase (default to 3)
 * \param  debuts vector of int (length nmp +1): indicate the first coordinate
 *         (using C enumeration style) of each region (in the x, y vectors).
 *         The last coordinate is n.
 * \param  debutsp vector of int (length nmp): the jth coordinate of this vector
 *         indicate the first coordinate (using C enumeration style)
 *         of POL for the jth region.
 * \return rans : SEXP, The R list
 * [[1]] rygeom the sf object (aka the cartogram)
 * [[2]] original area
 * [[3]] final area
 * [[4]] final x-coordinates of centroids
 * [[5]] final y-coordinates of centroids
* ******************************************************************/

SEXP dcn (SEXP rygeomd, SEXP rx, SEXP ry, SEXP rcount,
          SEXP rl1, SEXP rml2, SEXP rdup, SEXP rlast,
          SEXP rparamint, SEXP rparamdouble, SEXP roptions, SEXP rdebuts,
          SEXP rdebutsp)
{
 /*****************************************************************************/
  /* processing input  from R */
  /*****************************************************************************/
  /* list */
  rygeomd = PROTECT(rygeomd);
  /* double */
  rx = PROTECT(rx);
  ry = PROTECT(ry);
  rcount = PROTECT(rcount);
  rparamdouble = PROTECT(rparamdouble);
  double *x, *y, *count, *paramdouble, relerror, reltol,  MIN_POP_FAC;
  x=REAL(rx);
  y=REAL(ry);
  count=REAL(rcount);
  paramdouble=REAL(rparamdouble);
  relerror = paramdouble[0];
  reltol = paramdouble[1];
  MIN_POP_FAC = paramdouble[2];
  /* integer : l1, l2, paramint, options, debuts, debutsp  */
  rl1 = PROTECT(rl1);
  rml2 = PROTECT(rml2);
  rdup = PROTECT(rdup);
  rlast = PROTECT(rlast);
  rparamint = PROTECT(rparamint);
  roptions = PROTECT(roptions);
  rdebuts = PROTECT(rdebuts);
  rdebutsp = PROTECT(rdebutsp);
  int *l1, *ml2, *dup, *last, *paramint,
    maxit, np, verbose, absrel, nthreads, *options, *debuts,
    *debutsp, maxinc;
  l1 = INTEGER(rl1);
  ml2 = INTEGER(rml2);
  dup = INTEGER(rdup);
  last = INTEGER(rlast);
  paramint = INTEGER(rparamint);
  options = INTEGER(roptions);
  debuts = INTEGER(rdebuts);
  debutsp = INTEGER(rdebutsp);
  maxit = paramint[0];
  np = paramint[1];
  verbose = options[0];
  absrel = options[1];
  nthreads = options[2];
  maxinc=options[3];
  /* signals */
#if defined(R_SIGACTION)
  keep_runningg=1;
  struct sigaction act;
  act.sa_handler = intHandlerg;
  sigemptyset(&act.sa_mask);
  act.sa_flags = 0;
  sigaction(SIGINT, &act, NULL);
#elif defined(R_SIGWIN)
  keep_runningg=1;
  signal(SIGINT, intHandlerg);
#else
  int keep_runningg=1;
#endif
  /*********************************************************************/
  /* input and output from/to R */
  /*********************************************************************/
  /* list output is an R list object */
  /* [[1]] rygeom the sf object */
  /* [[2]] original area */
  /* [[3]] final area */
  /* [[4]] final x-coordinates of centroids */
  /* [[5]] final y-coordinates of centroids */
  SEXP  rygeom = PROTECT(duplicate(rygeomd)), rans;
  rans= PROTECT(allocVector(VECSXP, 6));
  SEXP rareasinit, rareasp, rcdgpx, rcdgpy ;
  double *areasp, *areasinit, *cdgpx, *cdgpy ;
  int nmp = length(rcount);
  rareasinit = PROTECT(allocVector(REALSXP, np));
  rareasp = PROTECT(allocVector(REALSXP, np));
  rcdgpx = PROTECT(allocVector(REALSXP, np));
  rcdgpy = PROTECT(allocVector(REALSXP, np));
  areasinit = REAL(rareasinit);
  areasp = REAL(rareasp);
  cdgpx = REAL(rcdgpx);
  cdgpy = REAL(rcdgpy);
  /************************************************************************/
  /* local variables */
  /************************************************************************/
  int i, j, k, line, iter, cur, n = length(rx), *done, *regofpol, nbofinc=0;
  double *areasObj, *areaspOld, *radiusp, *radiuspObj, *mass, *sizeError,
    min_area, *countp;
  double sumcount=0, sumareas, meansizeError, maxrelError, maxrelErrorOld,
    maxrelTol, forceReductionFactor,
    Fij, Fijx, Fijy, dGjtoV ;
  /* allocations  */
  areasObj = (double *) R_alloc(np, sizeof(double));
  areaspOld = (double *) R_alloc(np, sizeof(double));
  countp = (double *) R_alloc(np, sizeof(double));
  radiusp = (double *) R_alloc(np, sizeof(double));
  radiuspObj = (double *) R_alloc(np, sizeof(double));
  mass = (double *) R_alloc(np, sizeof(double));
  sizeError = (double *) R_alloc(np, sizeof(double));
  done = (int *) R_alloc(n, sizeof(int));
  for (i=0 ; i<n; i++) {
    done[i] = 0;
  }
  regofpol = (int *) R_alloc(np, sizeof(int));
  sumareas=0;
  /* initial state */
  /* areas and centroid */
  initpolygons(x, y, l1, ml2, cdgpx, cdgpy, areasp, countp, regofpol,
               count, nmp, debuts, debutsp, nthreads);
  /* calculus of objective Areas (sums) */
  for (i=0 ; i<np; i++) {
    areasinit[i]=areasp[i];
    sumareas += areasp[i];
    areasObj[i]=countp[i];
    sumcount += countp[i];
  }
  /* objective areas  */
  for (i=0 ; i<np; i++) {
    areasObj[i] *= (sumareas / sumcount);
  }
  /****** Replace target areas equal to zero by a small positive value. ******/
  min_area = areasObj[0];
  for (i=0 ; i<np; i++) {
    if (areasObj[i] > 0.0)
      min_area = fmin2(min_area, areasObj[i]);
    }
  for (i=0; i<np; i++) {
    if (areasObj[i] == 0.0)
      areasObj[i] = MIN_POP_FAC * min_area;
    /* init areaspOld for valgrind */
    areaspOld[i] = areasObj[i];
    }
  /* Radius/Mass/Error and finalization of objective Areas */
  meansizeError=0;
  maxrelError=0;
  for (i=0 ; i<np; i++) {
    /* radius calculus */
    if (areasp[i]<0)
      error("area negative, please check orientation with check_ring");
    radiuspObj[i] = sqrt(areasObj[i] / M_PI);
    radiusp[i] = sqrt(areasp[i] / M_PI);
    /* mass */
    mass[i] = (radiuspObj[i] - radiusp[i]);
    /* size error and its average */
    sizeError[i] = fmax2(areasp[i], areasObj[i]) /
      fmin2(areasp[i], areasObj[i]);
    meansizeError += sizeError[i];

    /* relative error (max) */
  }
  /* we don't have areas of previous iteration  */
  maxcritbyregion(areasp, areasObj, areasp, regofpol,
                  np, absrel, &maxrelError, &maxrelTol);
  /* mean of size error */
  meansizeError /= np;
  /* forcereduction factor */
  forceReductionFactor = 1/(1+  meansizeError);
/* verbose ? */
  if (verbose>1) Rprintf("Initial state\n");
  if (verbose>1) Rprintf(" * Max of abs or relative error %.6g\n", maxrelError);
  if (verbose>1) Rprintf(" * Mean size error %.6f\n", meansizeError);
  if (verbose==1) Rprintf("1 --- MaxErr %.6g\n", maxrelError);
  /* main loop */
  for (iter=0; iter<maxit; iter++) {
    /* update area centroids and error/radius/mass... */
    maxrelErrorOld=maxrelError;
    if (iter>0){
      areacdg(x, y, l1, ml2, cdgpx, cdgpy, areasp, nmp, debuts, debutsp, nthreads);
      meansizeError=0;
      maxrelError=0;
      maxrelTol=0;
      for (i=0 ; i<np; i++) {
        /* radius calculus */
        if (areasp[i] <0 ) {
          error("Negative area in region %d\n the algorithm may have made a crossed polygon ?\n Change algorithm ? Argument method=\"gsm\" in cartogramR() ?\n Returning last value.", i+1);
        }
        radiusp[i] = sqrt(areasp[i] / M_PI);
        /* mass */
        mass[i] = (radiuspObj[i] - radiusp[i]);
        /* size error and its average */
        sizeError[i] = fmax2(areasp[i], areasObj[i]) /
          fmin2(areasp[i], areasObj[i]);
        meansizeError += sizeError[i];
      }
      /* absolute or relative error (max) and tolerance */
      maxcritbyregion(areasp, areasObj, areaspOld, regofpol,
                  np, absrel, &maxrelError, &maxrelTol);
      /* stop if increment between two iterations of main loop is too small */
      if (maxrelTol<reltol) {
        if (verbose>=1) Rprintf("Normal main loop exit. Objective tol. is met: actual tol=%.6g < objective=%.6g\n",
                               maxrelTol, reltol);
        break;
      }
      /* stop if relError is small */
      if (maxrelError<relerror) {
        if (verbose>=1)
          Rprintf("Normal main loop exit. Objective err. is met: actual error=%.6g < objective=%.6g\n",
                  maxrelError, relerror);
        break;
      }
      /* mean of size error */
      meansizeError /= np;
      /* forcereduction factor */
      forceReductionFactor = 1/(1+ meansizeError);
      /* if (verbose==1) Rprintf("forceReductionFactor %.6f -- meansizeError %.6f\n", forceReductionFactor, meansizeError); */
      if (((iter+1) % 10 == 0) && (verbose>1)) {
          Rprintf("Iteration %d\n", iter+1);
          Rprintf(" * Max of abs or relative tolerance %.6g\n", maxrelTol);
          Rprintf(" * Max of abs or relative error %.6g\n", maxrelError);
          Rprintf(" * Mean size error %.8f\n", meansizeError);
      }
      if (((iter+1) % 100 == 0) && (verbose==1)) Rprintf("%d --- MaxErr %.6g\n", iter+1, maxrelError);
      if (!keep_runningg) error("User interruption");
    }
    /* if we do not decrease objective */
      if (maxrelError>maxrelErrorOld) nbofinc = nbofinc + 1;
      if (nbofinc > maxinc) {
        Rprintf(" The maximum number of increases (%d) in the criterion between 2 stages is exceeded (see option maxinc).\n",
                maxinc) ;
        Rprintf(" Main loop exit too early:\n");
        Rprintf("  Objective err. is not met: actual error=%.5g > objective=%.5g\n",
                maxrelError, relerror);
        Rprintf(" If the result does not satisfy your needs, please\n");
        Rprintf("  - increase verbosity level (to understand the problem),\n");
        Rprintf("  - increase maxit,\n");
        Rprintf("  - increase maxinc (risky),\n");
        Rprintf("  - increase maxrelError and maxrelTol\n");
        Rprintf(" in cartogramR() options.\n");
        break;
      }
    /* if we reach objective stop */
    if (maxrelError<relerror) break;
    /* update for increment between two iterations */
    memcpy(areaspOld, areasp, np);
    /* for (k=0; k<np ; k++) { */
    /*       areaspOld[k]=areasp[k]; */
    /*       } */
    /* loop on each polygon vertice */
    for (k=0; k<n ; k++) {
      /* calculation unless it is already done (ie duplicate and done)*/
      if ((dup[k]==(-1) || done[k]==0)) {
        Fijx = 0;
        Fijy = 0;
        /* for each centroid j  */
        for (j=0; j<np; j++){
          /* distance from centroid j to vertice k */
          dGjtoV= sqrt(pow((x[k] - cdgpx[j]),2) +
                       pow((y[k] - cdgpy[j]),2));
          /* Force */
          if (dGjtoV>=radiusp[j]) {
            Fij = mass[j] * radiusp[j]/dGjtoV;
          } else {
            Fij = mass[j] * pow(dGjtoV, 2) /
              pow(radiusp[j], 2) * (4 - 3*dGjtoV/radiusp[j]);
          }
          /* projection of Force on X and Y axis */
          Fijx += Fij / dGjtoV * (x[k] - cdgpx[j]) ;
          Fijy += Fij / dGjtoV * (y[k] - cdgpy[j]) ;
        }
        Fijx *=  forceReductionFactor;
        Fijy *=  forceReductionFactor;
        /* update coordinate vertice */
        x[k] = x[k] + Fijx;
        y[k] = y[k] + Fijy;
        /* update common vertices and mark them done */
        if ((dup[k]!=(-1) && done[k]==0)) {
          cur=dup[k];
          done[k]=1;
          while (cur!=k) {
            x[cur] = x[k];
            y[cur] = y[k];
            done[cur] = 1;
            cur = dup[cur];
          }
        }
      }
      /* if last coord of duplicates, then mark all not done */
      if ((dup[k]!=(-1) && k==last[k])) {
        done[k]=0;
        cur=dup[k];
        while (cur!=k) {
          done[cur] = 0;
          cur = dup[cur];
        }
      }
#if !defined(R_SIGACTION) && !defined(R_SIGWIN)
      R_CheckUserInterrupt();
#endif
    } /* end loop on vertices */
  } /* end main loop */
  /* warning on maxit */
  if ((iter>(maxit-1))&(verbose>=1)) {
    Rprintf("Maximum of iterations reached (%d): increase maxit ?\n", maxit);
    Rprintf("  Atm, options=list(maxit=%d) in cartogramR()\n", maxit);
  }
  /* update region=multipolygon and centroids */
  exitpolygons(cdgpx, cdgpy, areasp, areasinit, nmp, np, debutsp);
  /************************************************************************/
  /* result  R list of np components   */
  /************************************************************************/
  double  minx=0.0, miny=0.0, maxx=0.0, maxy=0.0, *bbox, *exthole;
  int nmpol, npol, ncoo, firstline;
  SEXP rmpol, rpol, rexthole, rbbox, rnamesbbox, rclassbbox;
  iter=0;
   for (i=0; i<nmp; i++) {
     /* MPOL=(POL1, POL2, ....,nmpol) */
     rmpol =  PROTECT(VECTOR_ELT(rygeom, i));
     nmpol = length(rmpol);
      for (j=0; j<nmpol; j++) {
        /* POL1=(EXT, HOLE1, HOLE2, ...,npol) */
        rpol =  PROTECT(VECTOR_ELT(rmpol, j));
        npol = length(rpol);
        for (k=0; k<npol; k++) {
          /* each coordinate k is a matrix */
          rexthole = PROTECT(VECTOR_ELT(rpol, k));
          ncoo = (int) length(rexthole) /2 ;
          /* matrix of size (ncoo x 2) */
          exthole=REAL(rexthole);
          /* the first ncoo-1 lines */
          firstline=iter;
          for (line=0; line<(ncoo-1); line++) {
            exthole[line] = x[iter];
            exthole[ncoo + line] = y[iter];
            /* bounding box */
            if (iter==0) {
          	  minx=x[iter];
              miny=y[iter];
              maxx=x[iter];
              maxy=y[iter];
            } else {
              minx=fmin2(minx, x[iter]);
              miny=fmin2(miny, y[iter]);
              maxx=fmax2(maxx, x[iter]);
              maxy=fmax2(maxy, y[iter]);
            }
            iter++;
          }
          /* the ncooth line: closing polygon */
          exthole[ncoo-1] = x[firstline];
          exthole[ncoo + ncoo -1] = y[firstline];
          UNPROTECT(1) ; /* rexthole : matrix EXT or HOLE */
        } /* end EXT or HOLE*/
        UNPROTECT(1) ; /* rpol */
      } /* end POL */
      UNPROTECT(1) ; /* rmpol */
   } /* end of MPOL */
   /* bbox : last attribute for ans */
   rbbox = PROTECT(allocVector(REALSXP, 4));
   /* assign values to vector rbbox */
   bbox = REAL(rbbox);
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
  /* original area is the second component */
  SET_VECTOR_ELT(rans, 1, rareasinit);
  /* final area is the third component */
  SET_VECTOR_ELT(rans, 2, rareasp);
  /* centroidx is the fourth component */
  SET_VECTOR_ELT(rans, 3, rcdgpx);
  /* centroidy is the fifth component */
  SET_VECTOR_ELT(rans, 4, rcdgpy);
  /*  is the sixth component all criteria*/
  SEXP rcrits;
  rcrits = PROTECT(allocVector(REALSXP, 2));
  double *crits;
  crits=REAL(rcrits);
  crits[0]=maxrelError;
  crits[1]=maxrelTol;
  SET_VECTOR_ELT(rans, 5, rcrits);
   /* unprotect and return */
   UNPROTECT(13); /* args of dcn */
   UNPROTECT(6); /* creation: rans, rygeom, rareasinit rareasp rcdgpx rcdgpy */
   UNPROTECT(3); /* class and attributes */
   UNPROTECT(1); /* rcrits */
   return rans;
}

/** \fn initpolygons
 *  \brief init for polygons: calculate centroid coordinates (cdgpx cdgpy) and area
 *   (areap) of each polygon and fill the countp vector
 *
 *   Each region is a Multipolygon: for example if regions are countries, for
 *   Japan we have several islands and thus several polygons.
 *   Region = MPOL = (POL1, POL2, ...)
 *   For each POLi we have possibly some holes
 *   POLi = (EXT, HOLE1, HOLE2, ...)
 *   We calculate area of each different POLi, and for each POL we affect
 *   in the vector countp the percentage of count calculated in function of area
 *
 * \param  x vector of double X-axis coordinates of polygons vertices
 * \param  y vector of doubleY-axis coordinates of polygons vertices
 * \param  l1 vector of int L1 columun of sf::st_coordinates gives the integer 1,2, 3
 *         got in the following list POL=(EXT1 HOLE2 HOLE3...)
 * \param  l2 vector of int L2 columun of sf::st_coordinates gives the integer 1,2, 3
 *   got in the following list MPOL=(POL1 POL2 POL3...)
 * \param  cdgpx vector of double X-axis coordinates of polygon centroid
 * \param  cdgpy vector of double Y-axis coordinates of polygon centroid
 * \param  areasp vector areas of polygons
 * \param  countp vector of double of length np that will contain for each POL
 *         percentage of count using proportionnality given by areas:
 *         region=(POL1, POL2,...), the count for region (in vector count)
 *         will be split using proportionnality given by area of each POLi
 * \param  regofpol vector of int of length np that will contain for each POL
 *         the number of region/multipolygon it belongs to
 * \param  count vector percentage of count for each region
 * \param  int, the number of regions (MPOL)
 * \param  debuts vector of int (length nmp +1): indicate the first coordinate
 *         (using C enumeration style) of each region (in the x, y vectors).
 *         The last coordinate is n.
 * \param  debutsp vector of int (length nmp): the jth coordinate of this vector
 *         indicate the first coordinate (using C enumeration style)
 *         of POL for the jth region.
 * \param nthreads int the number of threads for open mp
 * \return void (real results are in cdgpx cdgpy areasp, countp and regofpol)
 *******************************************************************/

void initpolygons(double *x,double *y, int *l1, int *l2,
          double *cdgpx, double *cdgpy, double *areasp, double *countp, int *regofpol,
                  double *count, int nmp, int *debuts, int *debutsp, int nthreads)
{
    double areatr2, apartsum, cdgx, cdgy, cdginx , cdginy, area,
        areain, cx, cy, som=0.0;
    int poly, polyinit, hole;
    int imp, i0, i, j, k;
if (nthreads == -1) nthreads= omp_get_num_procs() ;
#pragma omp parallel for  \
  private (areatr2, apartsum, cdgx, cdgy, cdginx , cdginy, area, \
           areain, cx, cy, som, poly, hole,imp, i0, i, k)        \
   num_threads(nthreads)\
  if (nthreads>1)
    for (j=0; j<nmp; j++){
      /* first coordinate */
      /* set level L1 */
      hole=l1[debuts[j]];
      i0=debuts[j];
      apartsum = 0;
      cdginx = 0;
      cdginy = 0;
      /* set  level L2 */
      poly=l2[debuts[j]];
      polyinit=poly;
      area=0;
      cdgx=0;
      cdgy=0;
      imp=debutsp[j];
      regofpol[imp]=j;
      /* loop on the jth region */
      for (i=debuts[j]; i<debuts[j+1]; i++) {
        /* if we are at the last coordinate of the region j
         * skip and finalize calculus */
        if (i<(debuts[j+1]-1)) {
        /*we are not in the last coordinate of the region j*/
          if (l2[i+1]!=poly) {
              /* change of polygon l2*/
              /* area and centroid L1*/
              cdginx /= 3 * apartsum;
              cdginy /= 3 * apartsum;
              areain = apartsum /2 ;
              /* update area and centroid L2 */
              area = area + areain;
              cdgx += cdginx * areain ;
              cdgy += cdginy * areain ;
              /* finalize centroid L2 */
              cdgx /= area;
              cdgy /= area;
              /* update area and centroid L3 */
              areasp[imp] = area;
              cdgpx[imp] = cdgx;
              cdgpy[imp] = cdgy;
              /* set the region number */
              regofpol[imp]=j;
              /* reset level L1 */
              hole=l1[i+1];
              i0=i+1;
              apartsum = 0;
              cdginx = 0;
              cdginy = 0;
              /* reset  level L2*/
              imp++;
              poly=l2[i+1];
              area=0;
              cdgx=0;
              cdgy=0;
              /* end reset  */
          } else {
            /* change l1 */
              if (l1[i+1]!=hole) {
                  /* area and centroid L1*/
                  cdginx /= 3 * apartsum;
                  cdginy /= 3 * apartsum;
                  areain = apartsum /2 ;
                  /* update area and centroid L2 */
                  area = area + areain;
                  cdgx += cdginx * areain ;
                  cdgy += cdginy * areain ;
                  /* reset level L1 */
                  hole=l1[i+1];
                  i0=i+1;
                  apartsum = 0;
                  cdginx = 0;
                  cdginy = 0;
                  /* end reset  */
              } else {
                /* no change, calculus */
                  cx = x[i0] + x[i] + x[i+1];
                  cy = y[i0] + y[i] + y[i+1];
                  areatr2 = (x[i]-x[i0]) *( y[i+1]-y[i0]) -
                      (x[i+1]-x[i0]) * (y[i]-y[i0]);
                  cdginx += areatr2 * cx;
                  cdginy += areatr2 * cy;
                  apartsum += areatr2;
              }
          }   /* end of L2 */
        } else {
        /* change of multipolygon=region */
        /* area and centroid L1*/
        cdginx /= 3 * apartsum;
        cdginy /= 3 * apartsum;
        areain = apartsum /2 ;
        /* update area and centroid L2 */
        area = area + areain;
        cdgx += cdginx * areain ;
        cdgy += cdginy * areain ;
        /* finalize centroid L2 */
        cdgx /= area;
        cdgy /= area;
        /* update area and centroid L2 */
        areasp[imp] = area ;
        cdgpx[imp] = cdgx ;
        cdgpy[imp] = cdgy ;
        /* countp */
        if ((poly-polyinit+1) > 1) {
          som = 0;
          /* total area for all POL of region */
          for (k=0; k<(poly-polyinit+1); k++) {
            som += areasp[imp - k];
          }
          /* percentage */
          for (k=0; k<(poly-polyinit+1); k++) {
            countp[imp - k] = count[j] * areasp[imp - k] / som;
          }
        } else countp[imp] = count[j]; /* only one POL in region */
        /* set the region number */
        regofpol[imp]=j;
        } /* end if else condition on i */
      } /*end of loop inside the jth region: index i*/
    } /* end of loop on all regions: index j */
   return;
}



/** \fn areacdg
 *  \brief calculate centroid area (areap) of each polygon
 *
  *
 * \param  x vector of double X-axis coordinates of polygons vertices
 * \param  y vector of doubleY-axis coordinates of polygons vertices
 * \param  l1 vector of int L1 columun of sf::st_coordinates gives the integer 1,2, 3
 *         got in the following list POL=(EXT1 HOLE2 HOLE3...)
 * \param  l2 vector of int gives the polygon number
 * \param  cdgpx vector of double X-axis coordinates of polygon centroid
 * \param  cdgpy vector of double Y-axis coordinates of polygon centroid
 * \param  areasp vector areas of polygons
 * \param  int, the number of regions (MPOL)
 * \param  debuts vector of int (length nmp +1): indicate the first coordinate
 *         (using C enumeration style) of each region (in the x, y vectors).
 *         The last coordinate is n.
 * \param  debutsp vector of int (length nmp): the jth coordinate of this vector
 *         indicate the first coordinate (using C enumeration style)
 *         of POL for the jth region.
 * \param nthreads int the number of threads for open mp
 * \return void (real results are in cdgpx cdgpy and areasp)
 *******************************************************************/

void areacdg(double *x,double *y, int *l1, int *l2,
             double *cdgpx, double *cdgpy, double *areasp,
             int nmp, int *debuts, int *debutsp, int nthreads)
{
  double areatr2, apartsum, cdgx, cdgy, cdginx , cdginy, area,
    areain, cx, cy;
  /* number of rows in y_geom */
  int poly, hole;
  int imp, i0, i, j;
  if (nthreads == -1) nthreads=  omp_get_num_procs() ;
#pragma omp parallel for \
  private (areatr2, apartsum, cdgx, cdgy, cdginx , cdginy, \
           area, areain, cx, cy, poly, hole, imp, i0, i)   \
  num_threads(nthreads) \
  if (nthreads>1)
    for (j=0; j<nmp; j++){
      /* first coordinate */
      /* set level L1 */
      hole=l1[debuts[j]];
      i0=debuts[j];
      apartsum = 0;
      cdginx = 0;
      cdginy = 0;
      /* set  level L2 */
      poly=l2[debuts[j]];
      area=0;
      cdgx=0;
      cdgy=0;
      imp=debutsp[j];
      /* loop on the jth region */
      for (i=debuts[j]; i<debuts[j+1]; i++) {
        /* if we are at the last coordinate of the region j
         * skip and finalize calculus */
        if (i<(debuts[j+1]-1)) {
        /*we are not in the last coordinate of the region j*/
          if (l2[i+1]!=poly) {
              /* change of polygon l2*/
              /* area and centroid L1*/
              cdginx /= 3 * apartsum;
              cdginy /= 3 * apartsum;
              areain = apartsum /2 ;
              /* update area and centroid L2 */
              area = area + areain;
              cdgx += cdginx * areain ;
              cdgy += cdginy * areain ;
              /* finalize centroid L2 */
              cdgx /= area;
              cdgy /= area;
              /* update area and centroid L3 */
              areasp[imp] = area;
              cdgpx[imp] = cdgx;
              cdgpy[imp] = cdgy;
              /* reset level L1 */
              hole=l1[i+1];
              i0=i+1;
              apartsum = 0;
              cdginx = 0;
              cdginy = 0;
              /* reset  level L2*/
              imp++;
              poly=l2[i+1];
              area=0;
              cdgx=0;
              cdgy=0;
              /* end reset  */
          } else {
            /* change l1 */
              if (l1[i+1]!=hole) {
                  /* area and centroid L1*/
                  cdginx /= 3 * apartsum;
                  cdginy /= 3 * apartsum;
                  areain = apartsum /2 ;
                  /* update area and centroid L2 */
                  area = area + areain;
                  cdgx += cdginx * areain ;
                  cdgy += cdginy * areain ;
                  /* reset level L1 */
                  hole=l1[i+1];
                  i0=i+1;
                  apartsum = 0;
                  cdginx = 0;
                  cdginy = 0;
                  /* end reset  */
              } else {
                /* no change, calculus */
                  cx = x[i0] + x[i] + x[i+1];
                  cy = y[i0] + y[i] + y[i+1];
                  areatr2 = (x[i]-x[i0]) *( y[i+1]-y[i0]) -
                      (x[i+1]-x[i0]) * (y[i]-y[i0]);
                  cdginx += areatr2 * cx;
                  cdginy += areatr2 * cy;
                  apartsum += areatr2;
              }
          }   /* end of L2 */
        } else {
        /* change of multipolygon=region */
        /* area and centroid L1*/
        cdginx /= 3 * apartsum;
        cdginy /= 3 * apartsum;
        areain = apartsum /2 ;
        /* update area and centroid L2 */
        area = area + areain;
        cdgx += cdginx * areain ;
        cdgy += cdginy * areain ;
        /* finalize centroid L2 */
        cdgx /= area;
        cdgy /= area;
        /* update area and centroid L2 */
        areasp[imp] = area ;
        cdgpx[imp] = cdgx ;
        cdgpy[imp] = cdgy ;
        } /* end if else condition on i */
      } /*end of loop inside the jth region: all polygons, index i*/
    } /* end of loop on all regions: index j */
   return;
}

/** \fn exitpolygons
 *  \brief exit for polygons: calculate centroid coordinates (cdgpx cdgpy) and area
 *   (areap) of each region=multipolygon and adjust obtained and objective areas
 *   to multipolygons.
 *
 *   Each region j is a Multipolygon: for example if regions are countries, for
 *   Japan we have several islands and thus several polygons.
 *   Regionj = MPOLj = (POL1, POL2, ...)
 *   For each POLi we have possibly some holes
 *   POLi = (EXT, HOLE1, HOLE2, ...)
 *   We calculate area of each different MPOLj from the final areas of each POLi,
 *   and calculate also the initial area for each MPOLj
 *
 * \param  cdgpx vector of double X-axis coordinates of polygon centroid
 * \param  cdgpy vector of double Y-axis coordinates of polygon centroid
 * \param  areasp vector areas of polygons i
 * \param  areasinit vector percentage of initial areas of each polygon i
 * \param  int, the number of regions (MPOL)
 * \param  int, the number of polygons for all regions
 * \param  debutsp vector of int (length nmp): the jth coordinate of this vector
 *         indicate the first coordinate (using C enumeration style)
 *         of POL for the jth region.
 * \return void (real results are in cdgpx cdgpy areasp and areasinit)
 *******************************************************************/

void exitpolygons(double *cdgpx, double *cdgpy, double *areasp, double *areasinit,
             int nmp, int np, int *debutsp)
{
    double cdgx, cdgy, somi=0.0, somf=0.0;
    int morethan1=0, i, j;
    for (j=0; j<nmp; j++){
      somi=0;
      somf=0;
      cdgx=0;
      cdgy=0;
      morethan1=0;
      if (j<(nmp - 1)) {
        for (i=debutsp[j]; i<debutsp[j+1]; i++){
          if (debutsp[j]== (debutsp[j+1]-1)) {
            /* one pol in region */
            areasinit[j]=areasinit[i];
            areasp[j]=areasp[i];
            cdgpx[j] = cdgpx[i];
            cdgpy[j] = cdgpy[i];
          } else {
            /* several pol in region */
            morethan1=1;
            somi += areasinit[i];
            somf += areasp[i];
            cdgx += cdgpx[i] * areasp[i];
            cdgy += cdgpy[i] * areasp[i];
          }
        }
      } else {
          for (i=debutsp[j]; i<np; i++){
          if (debutsp[j]== (np-1)) {
            /* one pol in region */
            areasinit[j]=areasinit[i];
            areasp[j]=areasp[i];
            cdgpx[j] = cdgpx[i];
            cdgpy[j] = cdgpy[i];
          } else {
            /* several pol in region */
            somi += areasinit[i];
            somf += areasp[i];
            cdgx += cdgpx[i] * areasp[i];
            cdgy += cdgpy[i] * areasp[i];
          }
        }
      }
      if (morethan1>0) {
       areasinit[j]= somi;
       areasp[j]= somf;
       cdgpx[j] = cdgx/somf;
       cdgpy[j] = cdgy/somf;
      }
      }
/* cleaning last entries */
  for (j=nmp; j<np; j++){
        areasinit[j]= 0;
       areasp[j]= 0;
       cdgpx[j] = 0;
       cdgpy[j] = 0;
  }
   return;
}

/** \fn maxcritbyregion
 *  \brief Calculate criteria by region "cumulating" crit by polygon
 *
 *   Each region j is a Multipolygon: for example if regions are countries, for
 *   Japan we have several islands and thus several polygons.
 *   Regionj = MPOLj = (POL1, POL2, ...)
 *   For each POLi we have possibly some holes
 *   POLi = (EXT, HOLE1, HOLE2, ...)
 *   We have criterion for each different POLi and we return max of
 *   criterion for mpolj
 *
 * \param  areasp vector areas of polygons
 * \param  areasObj vector objective areas of polygons
 * \param  areaspOld vector of areas of polygons (previous iteration)
 * \param  regofpol int vector of giving the region number for each poli
 * \param  np int, total number of polygons
 * \param  absrel int, if absrel==0 then absolute error else relative error
 * \param  relerr (pointer on) double, relative error
 * \param  reltol (pointer on) double, relative tolerance
 * \param  absrel int, if absrel==0 then absolute error else relative error
 * \return void (real results are in )
 *******************************************************************/

void maxcritbyregion(double *areasp, double *areasObj, double *areaspOld,
                       int *regofpol, int np, int absrel, double *relerr,
                       double *reltol)
{
  double maxrelError=0, maxrelTol=0, areareg=0, areaobjreg=0, areaoldreg=0;
  int i, poly=0;
  for (i=0 ; i<(np-1); i++) {
    if (regofpol[i+1]!=poly) {
      /* cumulating */
      areareg += areasp[i];
      areaoldreg += areaspOld[i];
      areaobjreg += areasObj[i];
      /* change of region */
      if (absrel==0) {
        maxrelTol=fmax2(maxrelTol,
                        fabs(areareg - areaoldreg));
        maxrelError=fmax2(maxrelError,
                          fabs(areareg - areaobjreg));
      } else {
        maxrelTol=fmax2(maxrelTol,
                        fabs(areareg - areaoldreg)/areaoldreg);
        maxrelError=fmax2(maxrelError,
                          fabs((areareg - areaobjreg)/areaobjreg));
      }
      /* reset */
      areareg=0;
      areaoldreg=0;
      areaobjreg=0;
      poly=regofpol[i+1];
    } else {
      /* no change, cumulating */
      areareg += areasp[i];
      areaoldreg += areaspOld[i];
      areaobjreg += areasObj[i];
    }
  }
  /* last polygon: cumulating */
  areareg += areasp[i];
  areaoldreg += areaspOld[i];
  areaobjreg += areasObj[i];
  /* last polygon: finalizing calculus */
  if (absrel==0) {
    maxrelTol=fmax2(maxrelTol,
                    fabs(areareg - areaoldreg));
    maxrelError=fmax2(maxrelError,
                      fabs(areareg - areaobjreg));
  } else {
    maxrelTol=fmax2(maxrelTol,
                    fabs(areareg - areaoldreg)/areaoldreg);
    maxrelError=fmax2(maxrelError,
                      fabs((areareg - areaobjreg)/areaobjreg));
  }
  /* results */
  *relerr=maxrelError;
  *reltol=maxrelTol;
  return  ;
}

