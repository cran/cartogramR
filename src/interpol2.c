/* function are adapted from these in go-cart       */
/* https://github.com/Flow-Based-Cartograms/go_cart */
/* and is realeased with the following Licence */
/* Adapted from the MIT License */

/* Copyright (c) 2017 Flow-Based-Cartograms */
/* Copyright (c) 2021 cartogramR */

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

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h> /* for fmin2 */
#include "interpol2.h"




/** \fn interpol2
 *  \brief Calculates bilinear interpolation
 *
 * Function to bilinearly interpolate a numerical array grid[0..lx*ly-1]
 * whose entries are numbers for the positions:
 * x = (0.5, 1.5, ..., lx-0.5), y = (0.5, 1.5, ..., ly-0.5).
 * The final argument "zero" can take two possible values: "x" or "y". If
 * zero==x, the interpolated function is forced to return 0 if x=0 or x=lx.
 * This option is suitable fo interpolating from gridvx because there can be
 * no flow through the boundary. If zero==y, the interpolation returns 0 if
 * y=0 or y=ly, suitable for gridvy. The unconstrained boundary will be
 * determined by continuing the function value at 0.5 (or lx-0.5 or ly-0.5)
 * all the way to the edge (i.e. the slope is 0 consistent with a cosine
 * transform). As it is mainly used to interpolate a displacement, forcing
 * it to be 0 on grid boundary is useful to stay on the map...
 * \param  x : coordx to interpolate
 * \param  y : coordy to interpolate
 * \param  grid : pointer on double, the grid
 * \param  zero : char, two possible values: "x" or "y"
 * \param  options : pointer on int (for verbose mode)
 * \param  errorloc : int (errorloc >0 if problem occurs)
 * \param  lx : int, number of grid points on x-axis
 * \param  ly : int, number of grid points on y-axis
 * \return value interpolated
 *******************************************************************/

double interpol2 (double x, double y, double *grid, char zero, int* options,
		 int *errorloc, int lx, int ly)
{
  double delta_x, delta_y, fx0y0, fx0y1, fx1y0, fx1y1, x0, x1, y0, y1;

  if (x<0 || x>lx || y<0 || y>ly) {
    if (options[0]>0) {
      Rprintf("ERROR: coordinate outside bounding box in interpol2().\n");
      Rprintf("x=%f, y=%f\n", x, y);
    }
    *errorloc=1;
    return -1;
  }
  if (zero != 'x' && zero != 'y' ) {
    if (options[0]>0) Rprintf("ERROR: unknown argument zero in interpol2().\n");
    *errorloc=2;
    return -1;
  }

  x0 =                             /* Nearest grid point smaller than x.     */
    MAX(0.0, floor(x+0.5) - 0.5);  /* Exception: if x<0.5, x0 becomes 0.0.   */
  x1 =                             /* Nearest grid point larger than x.      */
    MIN(lx, floor(x+0.5) + 0.5);   /* Exception: if x>lx-0.5, x1 becomes lx. */
  y0 = MAX(0.0, floor(y+0.5) - 0.5);                     /* Similarly for y. */
  y1 = MIN(ly, floor(y+0.5) + 0.5);
  delta_x = (x-x0) / (x1-x0);   /* On a scale from 0 to 1, how far is x (or */
  delta_y = (y-y0) / (y1-y0);   /* y) away from x0 (or y0)? 1 means x=x1.   */

  /* Function value at (x0, y0). */

  if ((x<0.5 && y<0.5) || (x<0.5 && zero == 'x') ||
      (y<0.5 && zero == 'y'))
    fx0y0 = 0.0;
  else
    fx0y0 = grid[(int)x0*ly + (int)y0];

  /* Function value at (x0, y1). */

  if ((x<0.5 && y>=ly-0.5) || (x<0.5 && zero == 'x') ||
      (y>=ly-0.5 && zero == 'y'))
    fx0y1 = 0.0;
  else if (x>=0.5 && y>=ly-0.5 && zero == 'x')
    fx0y1 = grid[(int)x0*ly + ly -1];
  else
    fx0y1 = grid[(int)x0*ly + (int)y1];

  /* Function value at (x1, y0). */

  if ((x>=lx-0.5 && y<0.5) || (x>=lx-0.5 && zero == 'x') ||
      (y<0.5 && zero == 'y'))
    fx1y0 = 0.0;
  else if (x>=lx-0.5 && y>=0.5 && zero == 'y')
    fx1y0 = grid[(lx-1)*ly + (int)y0];
  else
    fx1y0 = grid[(int)x1*ly + (int)y0];

  /* Function value at (x1, y1). */

  if ((x>=lx-0.5 && y>=ly-0.5) || (x>=lx-0.5 && zero == 'x') ||
      (y>=ly-0.5 && zero == 'y'))
    fx1y1 = 0.0;
  else if (x>=lx-0.5 && y<ly-0.5 && zero == 'y')
    fx1y1 = grid[(lx-1)*ly + (int)y1];
  else if (x<lx-0.5 && y>=ly-0.5 && zero == 'x')
    fx1y1 = grid[(int)x1*ly + ly - 1];
  else
    fx1y1 = grid[(int)x1*ly + (int)y1];

  return (1-delta_x)*(1-delta_y)*fx0y0 + (1-delta_x)*delta_y*fx0y1
    + delta_x*(1-delta_y)*fx1y0 + delta_x*delta_y*fx1y1;
}
