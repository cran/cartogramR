/* this file is part of go-cart                     */
/* https://github.com/Flow-Based-Cartograms/go_cart */
/* with adaptation to GNU-R  */
/* and are realeased with th following Licence */
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

#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <string.h>
#include "cartogram.h"
#include "inside_functions.h"
#include "myomp.h"

/***************************** Global variables. *****************************/

int **xyhalfshift2reg;

/**************************** Function prototypes. ***************************/

void rescale_map (double* centroidx, double* centroidy,
		  int* n_polycorn, int* options);
void interior (int* n_polycorn);
void gaussian_blur (double tot_init_area, double avg_dens);

/*****************************************************************************/
/* Function to change coordinates from (minx, miny, maxx, maxy) to           */
/* (0, 0, LX, LY).                                                           */

void rescale_map (double* centroidx, double* centroidy,
		  int* n_polycorn, int* options)
{
  double latt_const, new_maxx, new_maxy, new_minx, new_miny;
  int i, j;
  
  /* Minimum dimensions that leave enough space between map and rectangular  */
  /* boundaries.                                                             */

  new_maxx = 0.5 * ((1.0+PADDING)*map_maxx + (1.0-PADDING)*map_minx);
  new_minx = 0.5 * ((1.0-PADDING)*map_maxx + (1.0+PADDING)*map_minx);
  new_maxy = 0.5 * ((1.0+PADDING)*map_maxy + (1.0-PADDING)*map_miny);
  new_miny = 0.5 * ((1.0-PADDING)*map_maxy + (1.0+PADDING)*map_miny);  
  if (map_maxx-map_minx > map_maxy-map_miny) {    
    lx = L;
    latt_const = (new_maxx-new_minx) / L;
    ly = 1 << ((int)ceil(log2((new_maxy-new_miny)/latt_const)));
    new_maxy = 0.5*(map_maxy+map_miny) + 0.5*ly*latt_const;
    new_miny = 0.5*(map_maxy+map_miny) - 0.5*ly*latt_const;
  }
  else {
    ly = L;
    latt_const = (new_maxy-new_miny) / L;
    lx = 1 << ((int) ceil(log2((new_maxx-new_minx) / latt_const)));
    new_maxx = 0.5*(map_maxx+map_minx) + 0.5*lx*latt_const;
    new_minx = 0.5*(map_maxx+map_minx) - 0.5*lx*latt_const;
  }
  if (options[0]>1) Rprintf("Using a %d x %d lattice with bounding box\n\t(%f %f %f %f).\n",
	 lx, ly, new_minx, new_miny, new_maxx, new_maxy);

  /********************* Rescale all polygon coordinates. ********************/

  for (i=0; i<n_poly; i++)
    for (j=0; j<n_polycorn[i]; j++) {
      polycorn[i][j].x = (polycorn[i][j].x - new_minx) / latt_const;
      polycorn[i][j].y = (polycorn[i][j].y - new_miny) / latt_const;
    }

  /********************* Rescale all centroid coordinates. ********************/
  for (i=0; i<n_reg; i++) {
    centroidx[i] = (centroidx[i] - new_minx) / latt_const;
    centroidy[i] = (centroidy[i] - new_miny) / latt_const;
  }
  
  return;
}

/* Function to determine if a grid point at (x+0.5, y+0.5) is inside one of  */
/* the polygons and, if yes, in which region. The result is stored in the    */
/* array xyhalfshift2reg[x][y]. If (x+0.5, y+0.5) is outside all polygons,   */
/* then xyhalfshift2reg[x][y] = -1. If it is inside region i, then           */
/* xyhalfshift2reg[x][y] = i.                                                */

void interior (int* n_polycorn)
{
  int i, j, poly;

  /********************* Initialize xyhalfshift2reg[][]. *********************/

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++)
      xyhalfshift2reg[i][j] = -1;

  /**** Use the concept of the "crossing number" to determine the correct ****/
  /**** values of xyhalfshift2reg[][].                                    ****/

  for (i=0; i<n_reg; i++)
    for (j=0; j<n_polyinreg[i]; j++) {
      poly = polyinreg[i][j];
      set_inside_values_for_polygon(i, n_polycorn[poly], polycorn[poly],
				    xyhalfshift2reg);
    }
  return;
}

/*****************************************************************************/
/* Function to determine polygon area. This is needed to determine the       */
/* average population.                                                       */
/* The problem in short is to find the area of a polygon whose vertices are  */
/* given. Recall Stokes' theorem in 3d for a vector field v:                 */
/* integral[around closed curve dA]v(x,y,z).ds =                             */
/*                                          integral[over area A]curl(v).dA. */
/* Now let v(x,y,z) = (0,Q(x,y),0) and dA = (0,0,dx*dy). Then                */
/* integral[around closed curve dA]Q(x,y)dy =                                */
/*                                         integral[over area A]dQ/dx*dx*dy. */
/* If Q = x:                                                                 */
/* A = integral[over area A]dx*dy = integral[around closed curve dA]x dy.    */
/* For every edge from (x[i],y[i]) to (x[i+1],y[i+1]) there is a             */
/* parametrization                                                           */
/* (x(t),y(t)) = ((1-t)x[i]+t*x[i+1],(1-t)y[i]+t*y[i+1]), 0<t<1              */
/* so that the path integral along this edge is                              */
/* int[from 0 to 1]{(1-t)x[i]+t*x[i+1]}(y[i+1]-y[i])dt =                     */
/*                                          0.5*(y[i+1]-y[i])*(x[i]+x[i+1]). */
/* Summing over all edges yields:                                            */
/* Area = 0.5*[(x[0]+x[1])(y[1]-y[0]) + (x[1]+x[2])(y[2]-y[1]) + ...         */
/*               ... + (x[n-1]+x[n])(y[n]-y[n-1]) + (x[n]+x[0])(y[0]-y[n])]. */
/* ArcGIS treats a clockwise direction as positive, so that there is an      */
/* additional minus sign.                                                    */

double polygon_area (int ncrns, POINT *polygon)
{
  double area = 0.0;
  int i;

  for (i=0; i<ncrns-1; i++)
    /* no minus sign as SF polygon are in the right counterclockwise direction */
    /* area -=  */
    area +=
      0.5 * (polygon[i].x+polygon[i+1].x) * (polygon[i+1].y-polygon[i].y);
  /* no minus sign as SF polygon are in the right counterclockwise direction */
  /*  return area -= 0.5 * (polygon[ncrns-1].x+polygon[0].x) *
      (polygon[0].y-polygon[ncrns-1].y);*/
  return area += 0.5 * (polygon[ncrns-1].x+polygon[0].x) *
          (polygon[0].y-polygon[ncrns-1].y);
}

/*****************************************************************************/
/************ Function to smoothen the density near polygon edges. ***********/

void gaussian_blur (double tot_init_area, double avg_dens)
{
  double prefactor, scale_i, scale_j;
  fftw_plan plan_bwd;
  int i, j;

  /* Prepare the backward transform (i.e. from Fourier to real space). */

  plan_bwd = fftw_plan_r2r_2d(lx, ly, rho_ft, rho_init,
  			  FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);

  /* Upon forward and backward transform, rho_init is multiplied by the      */
  /* factor 4*lx*ly.                                                         */

  for (i=0; i<lx*ly; i++)
    rho_init[i] /= 4*lx*ly;

  /* At the start of this function, rho_ft is the unsmoothed density. We     */
  /* first have to transform it to Fourier space.                            */

  fftw_execute(plan_fwd);

  /* Now perform the Gaussian blur.                                          */

  prefactor = -0.5 * BLUR_WIDTH * BLUR_WIDTH * M_PI * M_PI;
  for (i=0; i<lx; i++) {
    scale_i = (double) i / lx;
    for (j=0; j<ly; j++) {
      scale_j = (double) j /ly;
      rho_ft[i*ly+j] *=
	exp(prefactor * (scale_i*scale_i + scale_j*scale_j));
    }
  }
  fftw_execute(plan_bwd);      /* Transform back from Fourier to real space. */

  /* Free memory. */

  fftw_destroy_plan(plan_bwd);

  return;
}

/*****************************************************************************/
/* Function to fill the lx-times-ly-grid with the initial density. It reads  */
/* the input target areas and produces a .eps image of the input polygons.   */
/* It also performs a Gaussian blur on the input density. This function      */
/* should only be called once, namely before the first round of integration. */
/* Afterwards use fill_with_density2() below.                                */

void fill_with_density1 (double* centroidx, double* centroidy,
			 int* n_polycorn, double* varregion, int* nb_polyinreg,
			 int* options, double* original_area)
{
  double  avg_dens, *dens, *init_area, min_area, tot_init_area,
          tot_target_area;
  int i,  j;

  /************************** Make region from data. **************************/

  read_gen(nb_polyinreg, options);

  /************************** Areas (on initial scale) **************************/
  for (i=0; i<n_reg; i++) {
    original_area[i]=0;
    for (j=0; j<n_polyinreg[i]; j++)
      original_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
				   polycorn[polyinreg[i][j]]);
  }
  /***************** Fit the map on an (lx)*(ly)-square grid. ****************/

  rescale_map(centroidx, centroidy, n_polycorn, options);

  /***************************** Allocate memory. ****************************/

  xyhalfshift2reg = (int**) malloc(lx * sizeof(int*));
  for (i=0; i<lx; i++)
    xyhalfshift2reg[i] = (int*) malloc(ly * sizeof(int));
  dens = (double*) malloc(n_reg * sizeof(double));
  target_area = (double*) malloc(n_reg * sizeof(double));
  init_area = (double*) calloc(n_reg, sizeof(double));

  /******* Determine inside which regions the grid points are located. *******/

  interior(n_polycorn);

  /*********************** Read target areas from file. **********************/
  for (i=0; i<n_reg; i++) target_area[i] = varregion[i];

  /****** Replace target areas equal to zero by a small positive value. ******/
  
  min_area = target_area[0];
  for (i=1; i<n_reg; i++) 
    if (target_area[i] > 0.0)
      min_area = MIN(min_area, target_area[i]);
  for (i=0; i<n_reg; i++)
    if (target_area[i] == 0.0)
      target_area[i] = MIN_POP_FAC * min_area;

  /**************** Calculate all initial areas and densities. ***************/

  for (i=0; i<n_reg; i++)
    for (j=0; j<n_polyinreg[i]; j++)
      init_area[i] += polygon_area(n_polycorn[polyinreg[i][j]],
				   polycorn[polyinreg[i][j]]);
  /* no division by area if data is already a density */
  for (i=0; i<n_reg; i++) {
       dens[i] = target_area[i] / init_area[i];
  }
  /******************** Calculate the "average density" = ********************/
  /**************** (total target area)/(total initial area). ****************/

  for (i=0, tot_init_area=0.0; i<n_reg; i++)
    tot_init_area += init_area[i];
  for (i=0, tot_target_area=0.0; i<n_reg; i++)
    tot_target_area += target_area[i];
  avg_dens = tot_target_area / tot_init_area;

  /***** Allocate memory for the Fourier transform of the input density. *****/

  rho_ft = (double*) fftw_malloc(lx * ly * sizeof(double));
  rho_init = (double*) fftw_malloc(lx * ly * sizeof(double));

  /************************** Digitize the density. **************************/

  for (i=0; i<lx; i++)
    for (j=0; j<ly; j++) {
      if (xyhalfshift2reg[i][j]==-1)
	rho_init[i*ly+j] = avg_dens;
      else
	rho_init[i*ly+j] = dens[xyhalfshift2reg[i][j]];
    }

  /*** Plan the Fourier transform already now so that it can be shared with **/
  /*** gaussian_blur().                                                     **/

  plan_fwd = fftw_plan_r2r_2d(lx, ly, rho_init, rho_ft,
  			  FFTW_REDFT10, FFTW_REDFT10, FFTW_ESTIMATE);

  /** Smoothen the density profile to avoid uncontrolled distortions around **/
  /** the edges of the polygons.                                            **/

  gaussian_blur(tot_init_area, avg_dens);

  /* Compute rho_ft[], the two dimensional cosine forward Fourier transform  */
  /* of rho_init[]. After the transform we have                              */
  /* rho_ft[i*ly+j] = 4 \sum_{m=0}^{lx-1} \sum_{n=0}^{ly-1}                  */
  /*           rho_init[m][n] cos(pi i (m+1/2) / lx) cos(pi j (n+1/2) / ly). */
  /* We must bear in mind that rho_0[m][n] is the density associated with    */
  /* the coordinates (m+1/2, n+1/2).                                         */

  fftw_execute(plan_fwd);

  /******************************* Free memory. ******************************/

  for (i=0; i<lx; i++)
    free(xyhalfshift2reg[i]);
  free(xyhalfshift2reg);
  free(dens);
  free(init_area);

  return;
}

/*****************************************************************************/
/* Function to fill the lx-times-ly grid with density *after* the first      */
/* round of integration. The main differences compared to                    */
/* fill_with_density1() are that fill_with_density2()                        */
/* - does not assign the target areas,                                       */
/* - does not produce a map,                                                 */
/* - does not perform a Gaussian blur.                                       */

void fill_with_density2 (int* n_polycorn, int *options)
{
  double avg_dens, *dens, *tmp_area, tot_target_area, tot_tmp_area;
  int i, j, nthreads;

  nthreads=options[6];
  /* Copy cartcorn[][] to polycorn[][]. */

   for (i=0; i<n_poly; i++)
     memcpy(polycorn[i], cartcorn[i], n_polycorn[i] * sizeof(POINT));
  /*   for (j=0; j<n_polycorn[i]; j++) */
  /*     polycorn[i][j] = cartcorn[i][j]; */

  /***************************** Allocate memory. ****************************/

  xyhalfshift2reg = (int**) malloc(lx * sizeof(int*));
  for (i=0; i<lx; i++)
    xyhalfshift2reg[i] = (int*) malloc(ly * sizeof(int));
  dens = (double*) malloc(n_reg * sizeof(double));
  tmp_area = (double*) calloc(n_reg, sizeof(double));

  /******* Determine inside which regions the grid points are located. *******/

  interior(n_polycorn);

  /** Calculate all region areas and densities up to this point in the      **/
  /** algorithm.                                                            **/

  for (i=0; i<n_reg; i++)
    for (j=0; j<n_polyinreg[i]; j++)
      tmp_area[i] +=
	polygon_area(n_polycorn[polyinreg[i][j]], polycorn[polyinreg[i][j]]);
  for (i=0; i<n_reg; i++) dens[i] = target_area[i] / tmp_area[i];

  /******************** Calculate the "average density" = ********************/
  /**************** (total target area)/(total initial area). ****************/

  for (i=0, tot_tmp_area=0.0; i<n_reg; i++)
    tot_tmp_area += tmp_area[i];
  for (i=0, tot_target_area=0.0; i<n_reg; i++)
    tot_target_area += target_area[i];
  avg_dens = tot_target_area / tot_tmp_area;

  /************************** Digitize the density. **************************/

  if (nthreads == -1) nthreads= omp_get_num_procs() ;
#pragma omp parallel for \
  private(j)  \
  num_threads(nthreads) \
  if (nthreads>1)
  for (i=0; i<lx; i++) {
    for (j=0; j<ly; j++) {
      if (xyhalfshift2reg[i][j]==-1)
        rho_init[i*ly+j] = avg_dens;
      else
        rho_init[i*ly+j] = dens[xyhalfshift2reg[i][j]];
    }
  }
  fftw_execute(plan_fwd);

  /******************************* Free memory. ******************************/

  for (i=0; i<lx; i++)
    free(xyhalfshift2reg[i]);
  free(xyhalfshift2reg);
  free(dens);
  free(tmp_area);

  return;
}
