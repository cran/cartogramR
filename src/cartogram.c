/* this file is part of go-cart                     */
/* https://github.com/Flow-Based-Cartograms/go_cart */
/* and is realeased with the following Licence */
/* Adapted from the MIT License */

/* Copyright (c) 2017 Flow-Based-Cartograms */

/* Permission is hereby granted, free of charge, to any person obtaining a copy */
/* of this software and associated documentation files (the "Software"), to deal */
/* in the Software without restriction, including without limitation the rights */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell */
/* copies of the Software, and to permit persons to whom the Software is */
/* furnished to do so, subject to the following conditions: */

/* The above copyright notice and this permission notice shall be included in all */
/* copies or substantial portions of the Software. */

/* Any images generated with the help of the Software shall be referenced to: */
/* Gastner, M., Seguy, V., & More, P. (2018). Fast flow-based algorithm for */
/* creating density-equalizing map projections. Proceedings of the National */
/* Academy of Sciences USA, 115:E2156-E2164. */

/******************************** Inclusions. ********************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cartogram.h"
#include "myomp.h"

/**************************** Function prototypes. ***************************/

double min4 (double a, double b, double c, double d);
double max4 (double a, double b, double c, double d);
POINT affine_transf(int triid, POINT *tri, double x, double y);

/*****************************************************************************/
/********** Function to project the polygons in the input .gen file. *********/

void project (double* centroidx, double* centroidy, Rboolean proj_graticule,
	      int* options, int* error_ptr, int* n_polycorn, Rboolean gridexport)
{
  double *xdisp, x2, *ydisp, y2, ctx, cty;
  int i, j, nthreads;
  int errorloc=0;

  /* The displacement vector (xdisp[i*ly+j], ydisp[i*ly+j]) is the point     */
  /* that was initially at (i+0.5, j+0.5). We work with (xdisp, ydisp)       */
  /* instead of (proj.x, proj.y) so that we can use the function interpol()  */
  /* defined in integrate.c.                                                 */
  nthreads=options[6];
  xdisp = (double*) malloc(lx * ly * sizeof(double));
  ydisp = (double*) malloc(lx * ly * sizeof(double));
  if (nthreads == -1) nthreads= omp_get_num_procs() ;
#pragma omp parallel for                        \
  private(j)                                    \
  num_threads(nthreads)                         \
  if (nthreads>1)
  for (i=0; i<lx; i++) {
    for (j=0; j<ly; j++) {
      xdisp[i*ly + j] = proj[i*ly + j].x - i - 0.5;
      ydisp[i*ly + j] = proj[i*ly + j].y - j - 0.5;
    }
  }
  /********************* Project the polygon coordinates. ********************/

    if (nthreads == -1) nthreads= omp_get_num_procs() ;
#pragma omp parallel for                        \
  private(j, ctx, cty)                          \
  reduction(max: errorloc)                      \
  num_threads(nthreads)                         \
  if (nthreads>1)
  for (i=0; i<n_poly; i++) {
    for (j=0; j<n_polycorn[i]; j++) {
      ctx = polycorn[i][j].x;
      cty = polycorn[i][j].y;
      cartcorn[i][j].x =
        interpol(polycorn[i][j].x, polycorn[i][j].y, xdisp, 'x', options, &errorloc)
        + ctx;
      cartcorn[i][j].y =
        interpol(polycorn[i][j].x, polycorn[i][j].y, ydisp, 'y', options, &errorloc)
        + cty;
    }
  }
  if (errorloc>0) {
    *error_ptr=errorloc;
      free(xdisp);
      free(ydisp);
      return ;
  }
   if (gridexport) {
    for (i=0; i<lx; i++)
      for (j=0; j<ly; j++) {
        ctx = proj3[i*ly + j].x ;
        cty = proj3[i*ly + j].y ;
        /* update grid */
        proj3[i*ly + j].x = interpol(ctx, cty, xdisp, 'x', options, error_ptr) + ctx;
        proj3[i*ly + j].y = interpol(ctx, cty, ydisp, 'y', options, error_ptr) + cty;
      }
  }
   /****************** Project proj2 on the basis of proj. ******************/
  if (proj_graticule) {
    for (i=0; i<lx*ly; i++) {
      x2 = proj2[i].x;
      y2 = proj2[i].y;
      proj2[i].x = interpol(x2, y2, xdisp, 'x', options, error_ptr) + x2;
      if (*error_ptr>0)  {
				free(xdisp);
				free(ydisp);
				return ;
      }
      proj2[i].y = interpol(x2, y2, ydisp, 'y', options, error_ptr) + y2;
      if (*error_ptr>0)  {
				free(xdisp);
				free(ydisp);
				return ;
      	}
    }
	}	
  /****************** Project centroid coordinate ******************/

  for (i=0; i<n_reg; i++) {
    ctx = centroidx[i];
    cty = centroidy[i];
    centroidx[i] = interpol(ctx, cty, xdisp, 'x', options, error_ptr) + ctx;
    if (*error_ptr>0)  {
      free(xdisp);
      free(ydisp);
      return ;
    }
    centroidy[i] = interpol(ctx, cty, ydisp, 'y', options, error_ptr) + cty;
    if (*error_ptr>0)  {
      free(xdisp);
      free(ydisp);
      return ;
    }
  }

  /******************************* Free memory. ******************************/

  free(xdisp);
  free(ydisp);

  return;
}

/*****************************************************************************/
/* Function to return the maximum absolute relative area error. The relative */
/* area error is defined by:                                                 */
/* area_on_cartogram / target_area - 1.                                      */
/* The function also updates the arrays cart_area[] and area_err[] that are  */
/* passed by reference.                                                      */

double max_area_err (double *area_err, double *cart_area,  int* n_polycorn,
		     POINT **corn, double *sum_cart_area)
{
  double max, obj_area, sum_target_area;
  int i, j;

  for (i=0; i<n_reg; i++) {
    cart_area[i] = 0.0;
    for (j=0; j<n_polyinreg[i]; j++)
      cart_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
				   corn[polyinreg[i][j]]);
  }
  for (i=0, sum_target_area=0.0; i<n_reg; i++)
    sum_target_area += target_area[i];
  for (i=0, *sum_cart_area=0.0; i<n_reg; i++)
    *sum_cart_area += cart_area[i];
  for (i=0; i<n_reg; i++) {
    obj_area =                         /* Objective area in cartogram units. */
      target_area[i] * (*sum_cart_area) / sum_target_area;
    area_err[i] = cart_area[i] / obj_area - 1.0;
  }
  max = 0.0;                   /* Determine the maximum absolute area error. */
  for (i=0; i<n_reg; i++)
    max = MAX(max, fabs(area_err[i]));

  return max;
}
/*****************************************************************************/
/* Function to return the maximum absolute area error. The area error is */
/* defined by:      abs(area_on_cartogram - target_area).  */
/* The function also updates the arrays cart_area[] and area_err[] that are  */
/* passed by reference.                                                      */

double max_absarea_err (double *area_err, double *cart_area,  int* n_polycorn,
		     POINT **corn, double *sum_cart_area)
{
  double max, obj_area, sum_target_area;
  int i, j;

  for (i=0; i<n_reg; i++) {
    cart_area[i] = 0.0;
    for (j=0; j<n_polyinreg[i]; j++)
      cart_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
				   corn[polyinreg[i][j]]);
  }
  for (i=0, sum_target_area=0.0; i<n_reg; i++)
    sum_target_area += target_area[i];
  for (i=0, *sum_cart_area=0.0; i<n_reg; i++)
    *sum_cart_area += cart_area[i];
  for (i=0; i<n_reg; i++) {
    obj_area =                         /* Objective area in cartogram units. */
      target_area[i] * (*sum_cart_area) / sum_target_area;
    area_err[i] = cart_area[i] - obj_area ;
  }
  max = 0.0;                   /* Determine the maximum absolute area error. */
  for (i=0; i<n_reg; i++)
    max = MAX(max, fabs(area_err[i]));

  return max;
}

