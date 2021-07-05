/* this file is part of go-cart                     */
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


#include <stdio.h>
#include <stdlib.h>
#include "cartogram.h"

/**************************** Function prototypes. ***************************/

void read_poly (int* n_polycorn, double* coordvertices);
void make_region (int* cumsum_polyinreg);

/*****************************************************************************/
/******************* Function to count the number of polygons. ***************/
void read_poly (int* n_polycorn, double* coordvertices)
{
  /*****************************************************************************/
  /**** Function to read polygon corners. The first and last vertex of each ****/
  /**** polygon must be identical.                                          ****/
  int i,j,coordctr=0;
  /* allocation step 1*/
  polycorn = (POINT**) malloc(n_poly * sizeof(POINT*));
  for (i=0; i<n_poly; i++) { 
    /* allocation */
    polycorn[i] = (POINT*) malloc(n_polycorn[i] * sizeof(POINT));
    /* setting values */
    for (j=0; j<n_polycorn[i]; j++) { 
      polycorn[i][j].x=coordvertices[coordctr];
      coordctr++;
      polycorn[i][j].y=coordvertices[coordctr];
      coordctr++;
    }
  }
  return;
}

/*****************************************************************************/
/* Function to make regions from polygons. Region IDs in the .gen file must  */
/* be nonnegative.                                                           */

void make_region (int* nb_polyinreg)
{
  
  int i, j, numligne;
  n_polyinreg = (int*) calloc(n_reg, sizeof(int));  /* Which polygons con-   */
  polyinreg = (int**) malloc(n_reg * sizeof(int*));
  numligne = 0;
  for (i=0; i<n_reg; i++) {
      n_polyinreg[i]=nb_polyinreg[i];
      polyinreg[i] = (int*) malloc(n_polyinreg[i] * sizeof(int));
      for (j=0; j<n_polyinreg[i]; j++) {
	polyinreg[i][j] = numligne;
	numligne++;
      }
  }
  return;
}

/*****************************************************************************/
/******************* Function to process map information. ********************/

void read_gen (int* nb_polyinreg, int* options)
{
  /* read_polyRsf(n_polycorn, rygeom, n_rows, nbpoly, multipoly); */
  /* read_poly(n_polycorn, coordvertices); */
  make_region(nb_polyinreg);
  if (options[0]>1) Rprintf("%i polygon(s), %i region(s)\n", n_poly, n_reg);
  return;
}
