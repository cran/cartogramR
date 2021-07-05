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

#include "inside_functions.h"
#include <math.h>



/*****************************************************************************/
/* Function to set values of inside[][], used in                             */
/* set_inside_values_for_polygon() below. It sets the value in inside[][]    */
/* for all x-values between poly_minx and the x-value (of the point on the   */
/* line connecting the given two coordinates) that corresponds to the        */
/* current y-value l.                                                        */

void set_inside_value_at_y (int region, POINT pk, POINT pn, int l,
          double poly_minx, int **inside)
{
  double intersection;
  int m;

  /* x-value of the intersection between y = l and the line formed by the    */
  /* coordinates (pkx, pky) and (pnx, pny).                                  */

  intersection = (pn.x-0.5 - (pk.x-0.5)) * (l - (pk.y-0.5)) /
          (pn.y-0.5 - (pk.y-0.5)) + (pk.x-0.5);
  for (m = (int) poly_minx; m < intersection; m++)
    inside[m][l] = region - inside[m][l] - 1;

  return;
}

/*****************************************************************************/
/* Function that takes two polygon coordinates and loops over the y-values   */
/* between the two input y-coordinates. It updates the value of inside[][]   */
/* for all points between polyminx (for this polygon) and the x-value at all */
/* coordinates on the horizontal line to the left of the line segment        */
/* connecting the input coordinates.                                         */

void set_inside_values_between_points (int region, POINT pk, POINT pn,
               double poly_minx, int **inside)
{
  int l;

  /* Loop over all integer y-values between the two input y-coordinates.     */

  for (l = ceil(MIN(pn.y, pk.y) - 0.5); l < MAX(pn.y - 0.5, pk.y - 0.5); l++)
    set_inside_value_at_y(region, pk, pn, l, poly_minx, inside);

  return;
}

/*****************************************************************************/
/**** Function to set values in inside[][] for a particular polygon in a  ****/
/**** region.                                                             ****/

void set_inside_values_for_polygon (int region, int n_polycorn,
            POINT *polycorn, int **inside)
{
  double poly_minx = polycorn[0].x;
  int k, n;

  /************ Determine the minimum x-coordinate of the polygon. ***********/

  for (k=0; k<n_polycorn; k++)
    poly_minx = MIN(poly_minx, polycorn[k].x);

  /* Loop over all pairs of consecutive coordinates of polygon.              */

  for (k=0, n=n_polycorn-1; k<n_polycorn; n=k++)
    set_inside_values_between_points(region, polycorn[k], polycorn[n],
             poly_minx, inside);
  return;
}
/*************** Helper function: return min/max of 4 numbers. ***************/
double min4 (double a, double b, double c, double d)
{
  if (a <= b && a <= c && a <= d)
    return a;
  if (b <= a && b <= c && b <= d)
    return b;
  if (c <= a && c <= b && c <= d)
    return c;
  return d;
}
double max4 (double a, double b, double c, double d)
{
  if (a >= b && a >= c && a >= d)
    return a;
  if (b >= a && b >= c && b >= d)
    return b ;
  if (c >= a && c >= b && c >= d)
    return c;
  return d;
}
