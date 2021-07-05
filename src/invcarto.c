#include <R.h>
#include <Rinternals.h>
#include "invcarto.h"
#include "interpol2.h"
#include "caracmap.h"
#include "inside_functions.h"

/**************************** Function prototypes. ***************************/
POINT affine_transf(int triid, POINT *tri, double x, double y, int lx, int ly);
void inv_project2 (POINT *proj, int lx, int ly, POINT *invproj2, int* options, int *errorloc);
POINT cdgquad(POINT pa, POINT  pb, POINT  pc, POINT  pd);
double areaquad(POINT pa, POINT  pb, POINT  pc, POINT  pd);
/********************************** Functions *******************************/

/** \fn invcarto
 *  \brief Start from a regular grid and apply the inverse cartogram deformation
 * (the flow based cartogram is a transformation T and here we approximate T^{-1})
 *  to that regular grid and we return that transformed grid
 *  and density at these points
 *
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
 * \param  rdensityct Double vector contains the density constant
* \return rans an R list
 * [[1]] rcoordx the x-axis coordinates of the inverse transform
 * [[2]] rcoordy the y-axis coordinates of the inverse transform
 * [[3]] rdensities the density at the x,y coordinate
 *******************************************************************/

SEXP invcarto (SEXP rgridx, SEXP rgridy, SEXP rpadding,
               SEXP rLL, SEXP rbbox, SEXP roptions, SEXP rdensityct)
{
  /*****************************************************************************/
  /* input and output from/to R */
  /*****************************************************************************/
  /* double */
  rgridx = PROTECT(rgridx);
  rgridy = PROTECT(rgridy);
  rbbox = PROTECT(rbbox);
  rpadding = PROTECT(rpadding);
  rdensityct = PROTECT(rdensityct);
  double *gridx, *gridy, *bbox, densityct, padding, map_minx, map_maxx,
    map_miny, map_maxy, *caracmapd;
  double latt_const, new_minx, new_miny;
  gridx = REAL(rgridx);
  gridy = REAL(rgridy);
  bbox = REAL(rbbox);
  padding = REAL(rpadding)[0];
  densityct = REAL(rdensityct)[0];
  /* bounding box (SF order)*/
  map_minx= bbox[0];
  map_miny = bbox[1];
  map_maxx = bbox[2];
  map_maxy = bbox[3];
  /* integer : option(s)   */
  rLL = PROTECT(rLL);
  roptions = PROTECT(roptions);
  /* MAP dimension */
  int LL, *options, *caracmapi, lx, ly, *errorloc;
  LL = INTEGER(rLL)[0];
  options = INTEGER(roptions);
  errorloc = (int *) R_alloc(1, sizeof(int));
  errorloc[0]=0;
  /* Map  */
  /************************************************************************/
  caracmapd = (double *) R_alloc(3, sizeof(double));
  caracmapi = (int *) R_alloc(2, sizeof(int));
  caract_map(caracmapd, caracmapi, padding, LL, map_maxx,
             map_maxy, map_minx,map_miny);
  lx = caracmapi[0];
  ly = caracmapi[1];
  latt_const = caracmapd[0];
  new_minx = caracmapd[1];
  new_miny = caracmapd[2];
  /*------------------------------------------------------------*/
  /* Result: rans  */
  /* [[1]] rcoordx x-axis coordinates */
  /* [[2]] rcoordy y-axis coordinates */
  /* [[3]] rdensities density at point */
  SEXP rans = PROTECT(allocVector(VECSXP, 3));
  SEXP rcoordx = PROTECT(allocMatrix(REALSXP, ly, lx));
  SEXP rcoordy = PROTECT(allocMatrix(REALSXP, ly, lx));
  SEXP rdensities = PROTECT(allocMatrix(REALSXP, ly, lx));
  double *coordx = REAL(rcoordx);
  double *coordy = REAL(rcoordy);
  double *densities = REAL(rdensities);
  /* components names of list rans  */
  SEXP rnamesans;
  rnamesans = PROTECT(allocVector(STRSXP, 3));
  /*******************************************/
  /* program */
  /*******************************************/
  /* local variables  */
  int i, j ;
  double area;
  POINT cdgA, cdgB, cdgC, cdgD;
  /* final grid in point type -> proj */
  /* result in point type -> invproj2 */
  POINT *proj, *invproj2;
  /* allocationS for result  */
  proj = (POINT*) malloc(lx * ly * sizeof(POINT));
  invproj2 = (POINT*) malloc(lx * ly * sizeof(POINT));
  /* populating the proj pointer on POINT with the final grid  */
  for (i=0; i<lx; i++) {
    for (j=0; j<ly; j++) {
      proj[i*ly + j].x = gridx[i*ly + j] ;
      proj[i*ly + j].y = gridy[i*ly + j] ;
    }
  }
  /* inversion of grid */
  inv_project2 (proj, lx, ly, invproj2, options, errorloc);
  if (errorloc[0]>0) {
    UNPROTECT(7); /* function arguments */
    UNPROTECT(5); /* result + name */
   /* free memory */
    FREEIC2;
    error("error in interpol2");
    SEXP nul=NILSXP;
    SET_VECTOR_ELT(rans, 0, nul);
    return rans;
  }
  /* densities in each  point of the grid Gij is calculated   */
  /* by determining the area of polygon ABCD which are the    */
  /* centroids of the adjacent graticules and approximate the */
  /* density as a constant on ABCD                            */
  /*        -------- Gi+1j+2                    */
  /*   Gij+2         /       \ Gi+2j+2           */
  /*   |     +D    /    +C    |                  */
  /*  /          /            |                  */
  /* Gij+1-----Gi+1,j+1 ---- Gi+2j+1             */
  /* |         /             |                   */
  /* Gij  +A  |     +B      /                    */
  /*   \     /             /                     */
  /*     \  |             Gi+2j                  */
  /*       Gi+1,j -----/                         */
  for (i=0; i<lx; i++) {
    for (j=0; j<ly; j++) {
      if ((i==0)&&(j==0)) {
          cdgA.x = 0;
          cdgA.y = 0;
          cdgB.x = 0.5*(invproj2[(i)*ly].x +invproj2[(i+1)*ly].x);
          cdgB.y = 0;
          cdgC = cdgquad(invproj2[i*ly + j], invproj2[(i+1)*ly + j],
                         invproj2[(i+1)*ly + j+1], invproj2[i*ly + j+1]);
          cdgD.x = 0;
          cdgD.y =  0.5*(invproj2[0].y +invproj2[1].y);
        }
      else if ((j==0)&&(i>0)&&(i<(lx-1))) {
          cdgA.x = 0.5*(invproj2[(i-1)*ly].x +invproj2[(i)*ly].x);
          cdgA.y = 0;
          cdgB.x = 0.5*(invproj2[(i)*ly].x +invproj2[(i+1)*ly].x);
          cdgB.y = 0;
          cdgC = cdgquad(invproj2[i*ly + j], invproj2[(i+1)*ly + j],
                         invproj2[(i+1)*ly + j+1], invproj2[i*ly + j+1]);
          cdgD = cdgquad(invproj2[(i-1)*ly + j], invproj2[(i)*ly + j],
                         invproj2[(i)*ly + j+1], invproj2[(i-1)*ly + j+1]);
        }
      else if ((j==0)&&(i==(lx-1))) {
          cdgA.x = 0.5*(invproj2[(i-1)*ly].x +invproj2[(i)*ly].x);
          cdgA.y = 0;
          cdgB.x = lx;
          cdgB.y = 0;
          cdgC.x = lx;
          cdgC.y = 0.5*(invproj2[i*ly + j].y +invproj2[(i)*ly + j+1].y);
          cdgD = cdgquad(invproj2[(i-1)*ly + j], invproj2[(i)*ly + j],
                         invproj2[(i)*ly + j+1], invproj2[(i-1)*ly + j+1]);
        }
      else if ((i==(lx-1))&&(j>0)&&(j<(ly-1))) {
          cdgA = cdgquad(invproj2[(i-1)*ly + j-1], invproj2[(i)*ly + j-1],
                         invproj2[(i)*ly + j], invproj2[(i-1)*ly + j]);
          cdgB.x = lx;
          cdgB.y = 0.5*(invproj2[i*ly + j-1].y +invproj2[(i)*ly + j].y);
          cdgC.x = lx;
          cdgC.y = 0.5*(invproj2[i*ly + j].y +invproj2[(i)*ly + j+1].y);
          cdgD = cdgquad(invproj2[(i-1)*ly + j], invproj2[(i)*ly + j],
                         invproj2[(i)*ly + j+1], invproj2[(i-1)*ly + j+1]);
        }
      else if ((i==(lx -1))&&(j==(ly-1))) {
          cdgA = cdgquad(invproj2[(i-1)*ly + j-1], invproj2[(i)*ly + j-1],
                         invproj2[(i)*ly + j], invproj2[(i-1)*ly + j]);
          cdgB.x = lx;
          cdgB.y = 0.5*(invproj2[i*ly + j-1].y +invproj2[(i)*ly + j].y);
          cdgC.x = lx;
          cdgC.y = ly;
          cdgD.x = 0.5*(invproj2[(i-1)*ly + j].x +invproj2[i*ly +j].x);
          cdgD.y = ly;
        }
      else if ((j==(ly -1))&&(i>0)&&(i<(lx-1))) {
          cdgA = cdgquad(invproj2[(i-1)*ly + j-1], invproj2[(i)*ly + j-1],
                         invproj2[(i)*ly + j], invproj2[(i-1)*ly + j]);
          cdgB = cdgquad(invproj2[i*ly + j-1], invproj2[(i+1)*ly + j-1],
                         invproj2[(i+1)*ly + j], invproj2[i*ly + j]);
          cdgC.x = 0.5*(invproj2[(i)*ly + j].x +invproj2[(i+1)*ly +j].x);
          cdgC.y = ly;
          cdgD.x = 0.5*(invproj2[(i-1)*ly + j].x +invproj2[i*ly +j].x);
          cdgD.y = ly;
        }
      else if ((i==0)&&(j==(ly-1))) {
          cdgA.x = 0;
          cdgA.y = 0.5*(invproj2[j-1].y +invproj2[j].y);
          cdgB = cdgquad(invproj2[i*ly + j-1], invproj2[(i+1)*ly + j-1],
                         invproj2[(i+1)*ly + j], invproj2[i*ly + j]);
          cdgC.x = 0.5*(invproj2[i*ly + j].x +invproj2[(i+1)*ly + j].x);
          cdgC.y = ly;
          cdgD.x = 0;
          cdgD.y = ly;
        }
      else if ((i==0)&&(j>0)&&(j<(ly-1))) {
          cdgA.x = 0;
          cdgA.y = 0.5*(invproj2[j-1].y +invproj2[j].y);
          cdgB = cdgquad(invproj2[i*ly + j-1], invproj2[(i+1)*ly + j-1],
                         invproj2[(i+1)*ly + j], invproj2[i*ly + j]);
          cdgC = cdgquad(invproj2[i*ly + j], invproj2[(i+1)*ly + j],
                         invproj2[(i+1)*ly + j+1], invproj2[i*ly + j+1]);
          cdgD.x = 0;
          cdgD.y = 0.5*(invproj2[j].y +invproj2[j+1].y);
        }
      else {
        cdgA = cdgquad(invproj2[(i-1)*ly + j-1], invproj2[(i)*ly + j-1],
                       invproj2[(i)*ly + j], invproj2[(i-1)*ly + j]);
        cdgB = cdgquad(invproj2[i*ly + j-1], invproj2[(i+1)*ly + j-1],
                       invproj2[(i+1)*ly + j], invproj2[i*ly + j]);
        cdgC = cdgquad(invproj2[i*ly + j], invproj2[(i+1)*ly + j],
                       invproj2[(i+1)*ly + j+1], invproj2[i*ly + j+1]);
        cdgD = cdgquad(invproj2[(i-1)*ly + j], invproj2[(i)*ly + j],
                       invproj2[(i)*ly + j+1], invproj2[(i-1)*ly + j+1]);
      }
      area = 0.5* areaquad(cdgA, cdgB, cdgC, cdgD);
      densities[i*ly + j] = densityct / area;
      /* rescale and set */
      coordx[i*ly + j] = invproj2[i*ly + j].x * latt_const + new_minx;
      coordy[i*ly + j] = invproj2[i*ly + j].y * latt_const + new_miny;;
    }
  }
  /* -------------------------------------------------------- */
  /* Export the results  */
  /* -------------------------------------------------------- */
  SET_VECTOR_ELT(rans, 0,rcoordx);
  SET_VECTOR_ELT(rans, 1, rcoordy);
  SET_VECTOR_ELT(rans, 2, rdensities);
  /* names of components of the list rans */
  SET_STRING_ELT(rnamesans, 0, mkChar("coordx"));
  SET_STRING_ELT(rnamesans, 1, mkChar("coordy"));
  SET_STRING_ELT(rnamesans, 2, mkChar("density"));
  setAttrib(rans, R_NamesSymbol, rnamesans);
  /* -------------------------------------------------------- */
  /* free memory, unprotect and return  */
  /* -------------------------------------------------------- */
  UNPROTECT(7); /* function arguments */
  UNPROTECT(5); /* result + name */
  /* free memory */
  FREEIC2;
  return rans;
}

/*****************************************************************************/
/* How do we compute the inverse projection? We partition the untransformed  */
/* lattice into right triangles. The following figure shows where the        */
/* vertices are and how we label the triangles.                              */
/*                                                                           */
/* ly   -------------------------------------------  ----------------------  */
/*  |  |\       4ly-1      /|\       8ly-1      /|    |\     4lx*ly-1     /| */
/*  |  |   \            /   |   \            /   |    |   \            /   | */
/*  |  |      \      /      |      \      /      |    |      \      /      | */
/* ly- | 4ly-3   \/   4ly-2 | 8ly-3   \/   8ly-2 |... | 4lx*ly  \/  4lx*ly | */
/* 0.5 |         /\         |         /\         |    | -3      /\  -2     | */
/*  |  |      /      \      |      /      \      |    |      /      \      | */
/*  |  |   /            \   |   /            \   |    |   /            \   | */
/*  |  |/       4ly-4      \|/       8ly-4      \|    |/     4lx*ly-4     \| */
/* ly  |--------------------|--------------------|-  -|--------------------| */
/* -1  |          .         |          .         |    |          .         | */
/*  .             .                    .                         .           */
/*  .  |          .         |          .         |    |          .         | */
/*  2  |--------------------|--------------------|-  -|--------------------| */
/*  |  |\         7        /|\       4ly+7      /|    |\    4(lx-1)ly+7   /| */
/*  |  |   \            /   |   \            /   |    |   \            /   | */
/*  |  |      \      /      |      \      /      |    |      \      /      | */
/* 1.5 |  5      \/      6  | 4ly+5   \/   4ly+6 |... | 4(lx-1) \/ 4(lx-1) | */
/*  |  |         /\         |         /\         |    | *ly + 5 /\ *ly + 6 | */
/*  |  |      /      \      |      /      \      |    |      /      \      | */
/*  |  |   /            \   |   /            \   |    |   /            \   | */
/*  |  |/         4        \|/       4ly+4      \|    |/    4(lx-1)ly+4   \| */
/*  1  |--------------------|--------------------|-  -|--------------------| */
/*  |  |\         3        /|\       4ly+3      /|    |\    4(lx-1)ly+3   /| */
/*  |  |   \            /   |   \            /   |    |   \            /   | */
/*  |  |      \      /      |      \      /      |    |      \      /      | */
/* 0.5 |  1      \/      2  | 4ly+1   \/   4ly+2 |... | 4(lx-1) \/ 4(lx-1) | */
/*  |  |         /\         |         /\         |    | *ly + 1 /\ *ly + 2 | */
/*  |  |      /      \      |      /      \      |    |      /      \      | */
/*  |  |   /            \   |   /            \   |    |   /            \   | */
/*  |  |/         0        \|/        4ly       \|    |/     4(lx-1)ly    \| */
/*  0   -------------------------------------------  ----------------------  */
/*     0 ------- 0.5 ------ 1 ------- 1.5 ------ 2 . lx-1 --- lx-0.5 ---- lx */
/*                                                                           */
/* We project all the vertices on this lattice, bearing in mind that we      */
/* already have saved almost half of these coordinates in proj[][].          */
/* Suppose we find that, after the cartogram transformation, a point (x, y)  */
/* is in the projected triangle (a, b, c). We want to find its position in   */
/* the original triangle (p, q, r). We locally approximate the cartogram     */
/* transformation by an affine transformation T such that T(a) = p,          */
/* T(b) = q and T(c) = r. We can think of T as a 3x3 matrix                  */
/*  /t11 t12 t13\                                                            */
/* | t21 t22 t23 |  such that                                                */
/*  \ 0   0   1 /                                                            */
/*  /t11 t12 t13\   /a1 b1 c1\     /p1 q1 r1\                                */
/* | t21 t22 t23 | | a2 b2 c2 | = | p2 q2 r2 | or TA = P. Hence T = PA^{-1}. */
/*  \ 0   0   1 /   \ 1  1  1/     \ 1  1  1/                                */
/*                              /b2-c2 c1-b1 b1*c2-b2*c1\                    */
/* We have A^{-1} = (1/det(A)) | c2-a2 a1-c1 a2*c1-a1*c2 |. By multiplying   */
/*                              \a2-b2 b1-a1 a1*b2-a2*b1/                    */
/* PA^{-1} we obtain t11, t12, t13, t21, t22, t23. The preimage of (x, y) in */
/* the unprojected map is then "pre" with coordinates                        */
/* pre.x = t11*x + t12*y + t13, pre.y = t21*x + t22*y + t23.                 */

/** \fn affine_transf
 *  \brief Function to perform the affine transform
 *  \param triid (int) ID of the triangle, see figure above.
 *  \param tri (pointer on POINT) - the vertices of the
 *         cartogram-transformed triangle
 *  \param x  (double) x-axis coordinates on the cartogram.
 *  \param y  (double) y-axis coordinates on the cartogram.
 *  \param lx (int) grid size on x-axis
 *  \param ly (int) grid size on y-axis
 *******************************************************************/
POINT affine_transf(int triid, POINT *tri, double x, double y, int lx, int ly)
{
  double ainv11, ainv12, ainv13, ainv21, ainv22, ainv23, ainv31, ainv32,
    ainv33, t11, t12, t13, t21, t22, t23, det;
  POINT p, pre, q, r;

  /* Determine the vertices p, q, r of the unprojected triangle from the ID  */
  /* of the triangle. Note that the order of the three points must match the */
  /* order of the vertices in tri[].                                         */

  switch (triid % 4) {
  case 0:
    p.x = triid / (4 * ly);
    p.y = (triid / 4) % ly;
    q.x = p.x + 0.5;
    q.y = p.y + 0.5;
    r.x = p.x + 1;
    r.y = p.y;
    break;
  case 1:
    p.x = triid / (4 * ly);
    p.y = (triid / 4) % ly;
    q.x = p.x;
    q.y = p.y + 1;
    r.x = p.x + 0.5;
    r.y = p.y + 0.5;
    break;
  case 2:
    p.x = triid / (4 * ly) + 0.5;
    p.y = (triid / 4) % ly + 0.5;
    q.x = p.x + 0.5;
    q.y = p.y + 0.5;
    r.x = q.x;
    r.y = q.y - 1;
    break;
  default:
    p.x = triid / (4 * ly);
    p.y = (triid / 4) % ly + 1;
    q.x = p.x + 1;
    q.y = p.y;
    r.x = p.x + 0.5;
    r.y = p.y - 0.5;
  }

  /**************************** Determinant of A. ****************************/

  det = tri[0].x * tri[1].y + tri[1].x * tri[2].y + tri[2].x * tri[0].y
    - tri[1].x * tri[0].y - tri[2].x * tri[1].y - tri[0].x * tri[2].y;

  /*********** Compute det(A) * A^{-1}. We divide by det(A) later. ***********/

  ainv11 = tri[1].y - tri[2].y;
  ainv12 = tri[2].x - tri[1].x;
  ainv13 = tri[1].x * tri[2].y - tri[1].y * tri[2].x;
  ainv21 = tri[2].y - tri[0].y;
  ainv22 = tri[0].x - tri[2].x;
  ainv23 = tri[0].y * tri[2].x - tri[0].x * tri[2].y;
  ainv31 = tri[0].y - tri[1].y;
  ainv32 = tri[1].x - tri[0].x;
  ainv33 = tri[0].x * tri[1].y - tri[0].y * tri[1].x;

  /******************************** Compute T. *******************************/

  t11 = p.x * ainv11 + q.x * ainv21 + r.x * ainv31;
  t12 = p.x * ainv12 + q.x * ainv22 + r.x * ainv32;
  t13 = p.x * ainv13 + q.x * ainv23 + r.x * ainv33;
  t21 = p.y * ainv11 + q.y * ainv21 + r.y * ainv31;
  t22 = p.y * ainv12 + q.y * ainv22 + r.y * ainv32;
  t23 = p.y * ainv13 + q.y * ainv23 + r.y * ainv33;

  /********************* Transform the input coordinates. ********************/

  pre.x = (t11*x + t12*y + t13) / det;
  pre.y = (t21*x + t22*y + t23) / det;

  return pre;
}


/** \fn inv_project2
 *  \brief calculate the inverse projection invproj2
 *
 *  \param proj (pointer on double) which gives the final grid
 *  \param lx (int) grid size on x-axis
 *  \param ly (int) grid size on y-axis
 *  \param invproj2 (pointer on double) will gives the coordinates
 *        of the "inverse projection"
 * \param  options : pointer on int (for verbose mode)
 * \param  error : pointer on int (for error)
 *  \return void
 *******************************************************************/
void inv_project2 (POINT *proj, int lx, int ly, POINT *invproj2,
                   int* options, int *error)
{
  double *xdisp, *ydisp;
  int i, j, k, **xyhalfshift2tri;
  POINT *invproj,  **projgrid, **tri;

  /**************************** Memory allocation. ***************************/

  xdisp = (double*) malloc(lx * ly * sizeof(double));
  ydisp = (double*) malloc(lx * ly * sizeof(double));
  invproj = (POINT*) malloc(lx * ly * sizeof(POINT));
  projgrid = (POINT**) malloc((lx+1) * sizeof(POINT*));
  for (i=0; i<=lx; i++)
    projgrid[i] = (POINT*) malloc((ly+1) * sizeof(POINT));
  tri = (POINT**) malloc(4 * lx * ly * sizeof(POINT*));
  for (i=0; i<4*lx*ly; i++)
    tri[i] = (POINT*) malloc(3 * sizeof(POINT));
  xyhalfshift2tri = (int**) malloc(lx * sizeof(int*));
  for (i=0; i<lx; i++)
    xyhalfshift2tri[i] = (int*) malloc(ly * sizeof(int));

  /* The displacement vector (xdisp[i*ly+j], ydisp[i*ly+j]) is the point     */
  /* that was initially at (i+0.5, j+0.5). We work with (xdisp, ydisp)       */
  /* instead of (proj.x, proj.y) so that we can use the function interpol()  */
  /* defined in integrate.c.                                                 */

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      xdisp[i*ly + j] = proj[i*ly + j].x - i - 0.5;
      ydisp[i*ly + j] = proj[i*ly + j].y - j - 0.5;
    }

  /* projgrid[i][j] is the projected position of (i, j) without half-shift.  */

  for (i=0; i<=lx; i++)
    for (j=0; j<=ly; j++) {
      projgrid[i][j].x = interpol2(i, j, xdisp, 'x', options, error, lx, ly)  + i;
      projgrid[i][j].y = interpol2(i, j, ydisp, 'y', options, error, lx, ly) + j;
      if (error[0]>0) {
        /* ------- free memory on error ------- */
        FREEIC1 ;
        return ;
      }
    }

  /************ Project the triangles shown in the lattice above. ************/

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      tri[4*(i*ly + j)][0].x =                      /* Lower left of square. */
        tri[4*(i*ly + j) + 1][0].x = projgrid[i][j].x;
      tri[4*(i*ly + j)][0].y =
        tri[4*(i*ly + j) + 1][0].y = projgrid[i][j].y;
      tri[4*(i*ly + j) + 1][1].x =                            /* Upper left. */
        tri[4*(i*ly + j) + 3][0].x = projgrid[i][j+1].x;
      tri[4*(i*ly + j) + 1][1].y =
        tri[4*(i*ly + j) + 3][0].y = projgrid[i][j+1].y;
      tri[4*(i*ly + j)][2].x =                               /* Lower right. */
        tri[4*(i*ly + j) + 2][2].x = projgrid[i+1][j].x;
      tri[4*(i*ly + j)][2].y =
        tri[4*(i*ly + j) + 2][2].y = projgrid[i+1][j].y;
      tri[4*(i*ly + j) + 2][1].x =                           /* Upper right. */
        tri[4*(i*ly + j) + 3][1].x = projgrid[i+1][j+1].x;
      tri[4*(i*ly + j) + 2][1].y =
        tri[4*(i*ly + j) + 3][1].y = projgrid[i+1][j+1].y;
      tri[4*(i*ly + j)][1].x =                                  /* Midpoint. */
        tri[4*(i*ly + j) + 1][2].x =
        tri[4*(i*ly + j) + 2][0].x =
        tri[4*(i*ly + j) + 3][2].x = proj[i*ly + j].x;
      tri[4*(i*ly + j)][1].y =
        tri[4*(i*ly + j) + 1][2].y =
        tri[4*(i*ly + j) + 2][0].y =
        tri[4*(i*ly + j) + 3][2].y = proj[i*ly + j].y;
    }

  /***** xyhalfshift2tri[i][j]=k means that (i+0.5, j+0.5) is in tri[k]. *****/

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++)
      xyhalfshift2tri[i][j] = -1;
  for (i=0; i<4*lx*ly; i++)
    set_inside_values_for_polygon(i, 3, tri[i], xyhalfshift2tri);

  /**** Inverse projection for a point at (i+0.5, j+0.5) on the cartogram. ***/

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      k = xyhalfshift2tri[i][j];
      invproj[i*ly + j] = affine_transf(k, tri[k], i+0.5, j+0.5, lx, ly);
    }

  /* Near sharp density gradients, there can be numerical artifacts. We      */
  /* polish them here. If the preimage of (i+0.5, j+0.5) has an x-coordinate */
  /* less than                                                               */
  /* min(invproj[i*ly + j - 1].x, invproj[i*ly + j + 1].x,                   */
  /*     invproj[(i-1)*ly + j].x, invproj[(i+1)*ly + j].x) - 1,              */
  /* greater than                                                            */
  /* max(invproj[i*ly + j - 1].x, invproj[i*ly + j + 1].x,                   */
  /*     invproj[(i-1)*ly + j].x, invproj[(i+1)*ly + j].x) + 1,              */
  /* a y-coordinate less than                                                */
  /* min(invproj[i*ly + j - 1].y, invproj[i*ly + j + 1].y,                   */
  /*     invproj[(i-1)*ly + j].y, invproj[(i+1)*ly + j].y) - 1               */
  /* or greater than                                                         */
  /* max(invproj[i*ly + j - 1].y, invproj[i*ly + j + 1].y,                   */
  /*     invproj[(i-1)*ly + j].y, invproj[(i+1)*ly + j].y) + 1,              */
  /* we replace invproj[i*ly + j] by the centroid of the four neighbouring   */
  /* preimages.                                                              */

  for (j=0; j<ly-1; j++) {
    invproj2[j].x = invproj[j].x;
    invproj2[j].y = invproj[j].y;
  }
  for (i=0; i<lx-1; i++) {
    invproj2[i*ly + ly - 1].x = invproj[i*ly + ly - 1].x;
    invproj2[i*ly + ly - 1].y = invproj[i*ly + ly - 1].y;
  }
  for (j=1; j<ly; j++) {
    invproj2[(lx-1)*ly + j].x = invproj[(lx-1)*ly + j].x;
    invproj2[(lx-1)*ly + j].y = invproj[(lx-1)*ly + j].y;
  }
  for (i=1; i<lx; i++) {
    invproj2[i*ly].x = invproj[i*ly].x;
    invproj2[i*ly].y = invproj[i*ly].y;
  }
  for (i=1; i<lx-1; i++)
    for (j=1; j<ly-1; j++) {
      if (invproj[i*ly + j].x < min4(invproj[i*ly + j - 1].x,
				     invproj[i*ly + j + 1].x,
				     invproj[(i-1)*ly + j].x,
				     invproj[(i+1)*ly + j].x) - 1 ||
	  invproj[i*ly + j].x > max4(invproj[i*ly + j - 1].x,
				     invproj[i*ly + j + 1].x,
				     invproj[(i-1)*ly + j].x,
				     invproj[(i+1)*ly + j].x) + 1 ||
	  invproj[i*ly + j].y < min4(invproj[i*ly + j - 1].y,
				     invproj[i*ly + j + 1].y,
				     invproj[(i-1)*ly + j].y,
				     invproj[(i+1)*ly + j].y) - 1 ||
	  invproj[i*ly + j].y > max4(invproj[i*ly + j - 1].y,
				     invproj[i*ly + j + 1].y,
				     invproj[(i-1)*ly + j].y,
				     invproj[(i+1)*ly + j].y) + 1) {
	invproj2[i*ly + j].x =
	  0.25 * (invproj[i*ly + j - 1].x + invproj[i*ly + j + 1].x +
		  invproj[(i-1)*ly + j].x + invproj[(i+1)*ly + j].x);
	invproj2[i*ly + j].y =
	  0.25 * (invproj[i*ly + j - 1].y + invproj[i*ly + j + 1].y +
		  invproj[(i-1)*ly + j].y + invproj[(i+1)*ly + j].y);
      }
      else {
	invproj2[i*ly + j].x = invproj[i*ly + j].x;
	invproj2[i*ly + j].y = invproj[i*ly + j].y;
      }
    }


  /******************************* Free memory. ******************************/
  FREEIC1 ;
  return;
}


/** \fn cdgquad
 *  \brief calculate the centroid of a quadrilateral
 *
 * \param  pa : POINT, point A
 * \param  pb : POINT, point B
 * \param  pc : POINT, point C
 * \param  pd : POINT, point D
 *            + C
 *           / \
 *          /   \
 *      D +      \
 *        |       \
 *        |        \
 *      A +---------+ B
 *
 * \return rans a POINT : the centroid
 * see Subject 2.02: How can the centroid of a polygon be computed?
 * http://www.faqs.org/faqs/graphics/algorithms-faq/
 **************************************************************************/
POINT cdgquad(POINT pa, POINT pb, POINT pc, POINT pd)
{
  double areatr1,  areatr2, xcdgtr1, xcdgtr2, ycdgtr1, ycdgtr2, sumareatr;
  POINT cdg;
  xcdgtr1 = (pa.x + pb.x + pc.x)/3;
  ycdgtr1 = (pa.y + pb.y + pc.y)/3;
  xcdgtr2 = (pa.x + pc.x + pd.x)/3;
  ycdgtr2 = (pa.y + pc.y + pd.y)/3;
  areatr1 = ((pb.x - pa.x) * (pc.y - pa.y) - (pc.x - pa.x) * (pb.y - pa.y))/2;
  areatr2 = ((pc.x - pa.x) * (pd.y - pa.y) - (pd.x - pa.x) * (pc.y - pa.y))/2;
  sumareatr = areatr1 + areatr2 ;
  cdg.x = (areatr1 * xcdgtr1 + areatr2 * xcdgtr2)/sumareatr ;
  cdg.y = (areatr1 * ycdgtr1 + areatr2 * ycdgtr2)/sumareatr ;
  return cdg;
}
/** \fn areaquad
 *  \brief calculate the centroid of a quadrilateral
 *
 * \param  pa : POINT, point A
 * \param  pb : POINT, point B
 * \param  pc : POINT, point C
 * \param  pd : POINT, point D
 *            + C
 *           / \
 *          /   \
 *      D +      \
 *        |       \
 *        |        \
 *      A +---------+ B
 *
 * \return rans a double : the area of (ABCD)
 *     using the formula below (with x_0=x_n and y_0=y_n)
 *     2 A( P ) = sum_{i=0}^{n-1} (x_i y_{i+1} - y_i x_{i+1})
 **************************************************************************/
double areaquad(POINT pa,POINT  pb,POINT  pc,POINT  pd)
{
  double area;
  area = ((pa.x*pb.y - pa.y* pb.x) + (pb.x*pc.y - pb.y*pc.x) +
          (pc.x*pd.y - pc.y*pd.x) + (pd.x*pa.y - pd.y*pa.x))/2;
  area = fabs(area);
  return area ;
}
