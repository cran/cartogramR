/** Calculate the grid used in cartogramR for a given L and padding
 * return a list of sfg POINT (the grid)
 *
 *******************************************************************/

/******************************** Inclusions. **********************/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
/**************************** Function prototypes. ******************/
void areacdg (double *x,double *y, int n, int *l1, int *l2, int *l3,
          double *cdgmpx, double *cdgmpy, double *areasmp);
SEXP dcn (SEXP rygeomd, SEXP rx, SEXP ry, SEXP rcount,
          SEXP rl1, SEXP rl2, SEXP rl3,
          SEXP rdup, SEXP rlast, SEXP rparamint,
          SEXP rparamdouble, SEXP roptions);
/********************************** Functions ***********************/

/** \fn dcn
 *  \brief Dougenik et al algorithm
 *
 * \param  rygeomd: The R list of polygons (each component i is a
 *                     multipolygon) of length nmp
 * \param  rx: vector double [1:n] X-axis coordinates of polygons vertices
 * \param  ry: vector double [1:n] Y-axis coordinates of polygons vertices
 * \param  rcount: vector double [1:nmp] of variable count or density
 * \param  rcount: vector double [1:nmp] of variable count or density
 * \param  rl1: vector int [1:n] L1 column of sf::st_coordinates gives
 *         the integer 1,2, 3 got in the following list
 *         POL=(EXT1 HOLE2 HOLE3...)
 * \param  rl2: vector int  [1:n] L2 column of sf::st_coordinates
 *         gives the integer 1,2, 3 got in the following list
 *         MPOL=(POL1 POL2 POL3...)
 * \param  rl3: vector int [1:n] L3 column of sf::st_coordinates gives
 *         the integer 1,2, 3 got in the following list
 *         (MPOL1 MPOL2 MPOL3...)
 * \param rdup: vector [1:n] if the coordinate k is a vertice with
 *        duplicates (k1, k2 ...) it gives the next indice k1
 *        if no duplicate exists it is 0
 * \param rlast: vector int [1:n] if the coordinate k is a vertice with
 *        duplicates (k1, k2 ..., kk) it gives the last indice kk
 *        of the list of duplicates
 * \param  rparamint vector int
 *           - [0] = maxit
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
 * \return rans : SEXP, The R list
 * [[1]] rygeom the sf object (aka the cartogram)
 * [[2]] original area
 * [[3]] final area
 * [[4]] final x-coordinates of centroids
 * [[5]] final y-coordinates of centroids
* ******************************************************************/

SEXP dcn (SEXP rygeomd, SEXP rx, SEXP ry, SEXP rcount,
          SEXP rl1, SEXP rl2, SEXP rl3, SEXP rdup, SEXP rlast,
          SEXP rparamint, SEXP rparamdouble, SEXP roptions)
{
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
  rans= PROTECT(allocVector(VECSXP, 5));
  SEXP rareasinit, rareasmp, rcdgmpx, rcdgmpy ;
  double *areasmp, *areasinit, *cdgmpx, *cdgmpy ;
  int nmp = length(rcount);
  rareasinit = PROTECT(allocVector(REALSXP, nmp));
  rareasmp = PROTECT(allocVector(REALSXP, nmp));
  rcdgmpx = PROTECT(allocVector(REALSXP, nmp));
  rcdgmpy = PROTECT(allocVector(REALSXP, nmp));
  areasinit = REAL(rareasinit);
  areasmp = REAL(rareasmp);
  cdgmpx = REAL(rcdgmpx);
  cdgmpy = REAL(rcdgmpy);
 /*****************************************************************************/
  /* processing input  from R */
  /*****************************************************************************/
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
  /* integer : l1, l2, l3, paramint   */
  rl1 = PROTECT(rl1);
  rl2 = PROTECT(rl2);
  rl3 = PROTECT(rl3);
  rdup = PROTECT(rdup);
  rlast = PROTECT(rlast);
  rparamint = PROTECT(rparamint);
  int *l1, *l2, *l3, *dup, *last, *paramint, maxit, verbose, absrel, *options;
  l1 = INTEGER(rl1);
  l2 = INTEGER(rl2);
  l3 = INTEGER(rl3);
  dup = INTEGER(rdup);
  last = INTEGER(rlast);
  paramint = INTEGER(rparamint);
  options = INTEGER(roptions);
  maxit = paramint[0];
  verbose = options[0];
  absrel = options[1];
  /************************************************************************/
  /* local variables */
  /************************************************************************/
  int i, j, k, line, iter, cur, n = length(rx), *done;
  double *areasObj, *areasmpOld, *radiusmp, *radiusObj, *mass, *sizeError,
    min_area ;
  double sumcount=0, sumareas, meansizeError, maxrelError,
    maxrelTol, forceReductionFactor,
    Fij, Fijx, Fijy, dGjtoV;
  /* allocations  */
  areasObj = (double *) R_alloc(nmp, sizeof(double));
  areasmpOld = (double *) R_alloc(nmp, sizeof(double));
  radiusmp = (double *) R_alloc(nmp, sizeof(double));
  radiusObj = (double *) R_alloc(nmp, sizeof(double));
  mass = (double *) R_alloc(nmp, sizeof(double));
  sizeError = (double *) R_alloc(nmp, sizeof(double));
  done = (int *) R_alloc(n, sizeof(int));
  for (i=0 ; i<n; i++) {
    done[i] = 0;
  }
  sumareas=0;
  /* initial state */
  /* areas and centroid */
  areacdg(x, y, n, l1, l2, l3, cdgmpx, cdgmpy, areasmp);
  /* calculus of objective Areas (sums) */
  for (i=0 ; i<nmp; i++) {
    areasinit[i]=areasmp[i];
    sumareas += areasmp[i];
    areasObj[i]=count[i];
    sumcount += count[i];
  }
  /* objective areas  */
  for (i=0 ; i<nmp; i++) {
    areasObj[i] *= (sumareas / sumcount);
  }
  /****** Replace target areas equal to zero by a small positive value. ******/
  min_area = areasObj[0];
  for (i=0 ; i<nmp; i++) {
    if (areasObj[i] > 0.0)
      min_area = fmin2(min_area, areasObj[i]);
    }
  for (i=0; i<nmp; i++) {
    if (areasObj[i] == 0.0)
      areasObj[i] = MIN_POP_FAC * min_area;
    /* init areasmpOld for valgrind */
    areasmpOld[i] = areasObj[i];
    }
  /* Radius/Mass/Error and finalization of objective Areas */
  meansizeError=0;
  maxrelError=0;
  for (i=0 ; i<nmp; i++) {
    /* radius calculus */
    if (areasmp[i]<0)
      error("area negative, please check orientation with check_ring");
    radiusObj[i] = sqrt(areasObj[i] / M_PI);
    radiusmp[i] = sqrt(areasmp[i] / M_PI);
    /* mass */
    mass[i] = (radiusObj[i] - radiusmp[i]);
    /* size error and its average */
    sizeError[i] = fmax2(areasmp[i], areasObj[i]) /
      fmin2(areasmp[i], areasObj[i]);
    meansizeError += sizeError[i];
    /* relative error (max) */
    if (absrel==0) {
    maxrelError=fmax2(maxrelError,
                      fabs(areasmp[i] - areasObj[i]));
    } else {
    maxrelError=fmax2(maxrelError,
                      fabs((areasmp[i] - areasObj[i])/areasObj[i]));
    }
  }
  /* mean of size error */
  meansizeError /= nmp;
  /* forcereduction factor */
  forceReductionFactor = 1/(1+  meansizeError);
/* verbose ? */
  if (verbose>1) Rprintf("Initial state\n");
  if (verbose>1) Rprintf(" * Max of abs or relative error %.6f\n", maxrelError);
  if (verbose>1) Rprintf(" * Mean size error %.6f\n", meansizeError);
  if (verbose==1) Rprintf("1 * MaxErr %.6f\n", maxrelError);
  /* main loop */
  for (iter=0; iter<maxit; iter++) {
    /* update area centroids and error/radius/mass... */
    if (iter>0){
      areacdg(x, y, n, l1, l2, l3, cdgmpx, cdgmpy, areasmp);
      meansizeError=0;
      maxrelError=0;
      maxrelTol=0;
/* #pragma omp parallel shared(maxrelError, maxrelTol) */
/* #pragma omp for    nowait reduction(max : maxrelError, maxrelTol) */
      for (i=0 ; i<nmp; i++) {
        /* radius calculus */
        radiusmp[i] = sqrt(areasmp[i] / M_PI);
        /* mass */
        mass[i] = (radiusObj[i] - radiusmp[i]);
        /* size error and its average */
        sizeError[i] = fmax2(areasmp[i], areasObj[i]) /
          fmin2(areasmp[i], areasObj[i]);
        meansizeError += sizeError[i];
        /* absolute or relative error (max) */
        if (absrel==0) {
        maxrelError=fmax2(maxrelError,
                          fabs(areasmp[i] - areasObj[i]));
        maxrelTol=fmax2(maxrelTol,
                        fabs(areasmp[i]-areasmpOld[i]));
        } else {
        maxrelError=fmax2(maxrelError,
                          fabs((areasmp[i] - areasObj[i])/areasObj[i]));
        maxrelTol=fmax2(maxrelTol,
                        fabs((areasmp[i]-areasmpOld[i])/areasmpOld[i]));
          }
      }
      /* stop if increment beteween two iterations of main loop is too small */
      if (maxrelTol<reltol) break;
      /* stop if relError is small */
      if (maxrelError<relerror) break;
      /* mean of size error */
      meansizeError /= nmp;
      /* forcereduction factor */
      forceReductionFactor = 1/(1+ meansizeError);
      /* if (verbose==1) Rprintf("forceReductionFactor %.6f -- meansizeError %.6f\n", forceReductionFactor, meansizeError); */
      if (verbose>1) Rprintf("Iteration %d\n", iter+1);
      if (verbose>1) Rprintf(" * Max of abs or relative error %.6f\n", maxrelError);
      if (verbose>1) Rprintf(" * Mean size error %.8f\n", meansizeError);
      if (verbose==1) Rprintf("%d * MaxErr %.6f\n", iter+1, maxrelError);
      R_CheckUserInterrupt();
    }
    /* if we reach objective stop */
    if (maxrelError<relerror) break;
    /* update for increment between two iterations */
    memcpy(areasmpOld, areasmp, nmp);
    /* for (k=0; k<nmp ; k++) { */
    /*       areasmpOld[k]=areasmp[k]; */
    /*       } */
    /* loop on each polygon vertice */
    for (k=0; k<n ; k++) {
      /* calculation unless it is already done (ie duplicate and done)*/
      if (dup[k]==(-1) || done[k]==0) {
        Fijx = 0;
        Fijy = 0;
        /* for each centroid j  */
        for (j=0; j<nmp; j++){
          /* distance from centroid j to vertice k */
          dGjtoV= sqrt(pow((x[k] - cdgmpx[j]),2) +
                       pow((y[k] - cdgmpy[j]),2));
          /* Force */
          if (dGjtoV>=radiusmp[j]) {
            Fij = mass[j] * radiusmp[j]/dGjtoV;
          } else {
            Fij = mass[j] * pow(dGjtoV, 2) /
              pow(radiusmp[j], 2) * (4 - 3*dGjtoV/radiusmp[j]);
          }
          /* projection of Force on X and Y axis */
          Fijx += Fij / dGjtoV * (x[k] - cdgmpx[j]) ;
          Fijy += Fij / dGjtoV * (y[k] - cdgmpy[j]) ;
        }
        Fijx *=  forceReductionFactor;
        Fijy *=  forceReductionFactor;
        /* update coordinate vertice */
        x[k] = x[k] + Fijx;
        y[k] = y[k] + Fijy;
        /* update common vertices and mark them done */
        if (dup[k]!=(-1) && done[k]==0) {
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
      if (dup[k]!=(-1) && k==last[k]) {
        done[k]=0;
        cur=dup[k];
        while (cur!=k) {
          done[cur] = 0;
          cur = dup[cur];
        }
      }
    } /* end loop on vertices */
  } /* end main loop */
  /************************************************************************/
  /* result  R list of nmp components */
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
  SET_VECTOR_ELT(rans, 2, rareasmp);
  /* centroidx is the fourth component */
  SET_VECTOR_ELT(rans, 3, rcdgmpx);
  /* centroidy is the fifth component */
  SET_VECTOR_ELT(rans, 4, rcdgmpy);
   /* unprotect and return */
   UNPROTECT(16);
   UNPROTECT(3); /* class and attributes */
   return rans;
}
/** \fn areacdg
 *  \brief calculate centroid coordinates (cdgmpx cdgmpy) and area
 *   (areamp) of each polygon
 *
 * \param  x vector of double X-axis coordinates of polygons vertices
 * \param  y vector of doubleY-axis coordinates of polygons vertices
 * \param  l1 vector of int L1 columun of sf::st_coordinates gives the integer 1,2, 3
 *         got in the following list POL=(EXT1 HOLE2 HOLE3...)
 * \param  l2 vector of int L2 columun of sf::st_coordinates gives the integer 1,2, 3
 *   got in the following list MPOL=(POL1 POL2 POL3...)
 * \param  l3 vector of int L3 columun of sf::st_coordinates gives the integer 1,2, 3
 *   got in the following list (MPOL1 MPOL2 MPOL3...)
 * \param  cdgmpx vector of double X-axis coordinates of polygon centroid
 * \param  cdgmpy vector of double Y-axis coordinates of polygon centroid
 * \param  areasmp vector areas of polygons
* \return void (real results are in cdgmpx cdgmpy and areasmp)
 *******************************************************************/

void areacdg (double *x,double *y, int n, int *l1, int *l2, int *l3,
          double *cdgmpx, double *cdgmpy, double *areasmp)
{
  /*****************************************************************************/
  double areatr2, apartsum, cdgx, cdgy, cdginx , cdginy, area,
    areain, cx, cy;
  /* number of rows in y_geom */
  int mpoly= 0;
  int poly=0;
  int hole=0;
  int imp=0, i0=0, i=0;
  /* first coordinate */
  /* set level L1 */
  hole=l1[i];
  i0=i;
  apartsum = 0;
  cdginx = 0;
  cdginy = 0;
  /* set  level L2 */
  poly=l2[i];
  area=0;
  cdgx=0;
  cdgy=0;
  /* set L3 */
  mpoly=l3[i];
  areasmp[imp]=0;
  cdgmpx[imp]=0;
  cdgmpy[imp]=0;
  /* end set  */
   for (i=1; i<(n-1); i++) {
     if (l3[i+1]!=mpoly) {
       /* change of multipolygons */
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
       areasmp[imp] += area;
       cdgmpx[imp] += cdgx * area;
       cdgmpy[imp] += cdgy * area;
       /* finalize centroid L3 */
       cdgmpx[imp] /= areasmp[imp];
       cdgmpy[imp] /= areasmp[imp];
       /* reset level L1 */
       i++;
       hole=l1[i];
       i0=i;
       apartsum = 0;
       cdginx = 0;
       cdginy = 0;
       /* reset  level L2 */
       poly=l2[i];
       area=0;
       cdgx=0;
       cdgy=0;
       /* reset L3 */
       mpoly=l3[i];
       imp++;
       areasmp[imp]=0;
       cdgmpx[imp]=0;
       cdgmpy[imp]=0;
       /* end reset  */
     } else {
       if (l2[i+1]!=poly) {
         /* change of polygon */
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
         areasmp[imp] += area;
         cdgmpx[imp] += cdgx * area;
         cdgmpy[imp] += cdgy * area;
         /* reset level L1 */
         i++;
         hole=l1[i];
         i0=i;
         apartsum = 0;
         cdginx = 0;
         cdginy = 0;
         /* reset  level L2*/
         poly=l2[i];
         area=0;
         cdgx=0;
         cdgy=0;
         /* end reset  */
       } else {
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
           i++;
           hole=l1[i];
           i0=i;
           apartsum = 0;
           cdginx = 0;
           cdginy = 0;
           /* end reset  */
         } else {
           cx = x[i0] + x[i] + x[i+1];
           cy = y[i0] + y[i] + y[i+1];
           areatr2 = (x[i]-x[i0]) *( y[i+1]-y[i0]) -
             (x[i+1]-x[i0]) * (y[i]-y[i0]);
           cdginx += areatr2 * cx;
           cdginy += areatr2 * cy;
           apartsum += areatr2;
         }
       }
     }
   }
   /* last coordinate */
   /* change of multipolygons */
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
   areasmp[imp] += area;
   cdgmpx[imp] += cdgx * area;
   cdgmpy[imp] += cdgy * area;
   /* finalize centroid L3 */
   cdgmpx[imp] /= areasmp[imp];
   cdgmpy[imp] /= areasmp[imp];
   return;
 }
