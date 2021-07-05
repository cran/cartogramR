#include "cartogram.h"
/********************************** Functions *******************************/
/** \fn inv_rescale_map
*  \brief recalculate coordinates on initial scale
 *
 * \param  centroidx : pointer on double, vector of x-coordinates of centroids
 * \param  centroidx : pointer on double, vector of y-coordinates of centroids
 * \param  n_polycorn : pointer on int, vector of numbers of vertices for each 
 *                     polygon
 * \param options : pointer on int, vector of options
 * \return void
 *******************************************************************/

void inv_rescale_map (double* centroidx, double* centroidy, int* n_polycorn, int* options)
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
      cartcorn[i][j].x = cartcorn[i][j].x*latt_const + new_minx;
      cartcorn[i][j].y = cartcorn[i][j].y*latt_const + new_miny;
    }

  /********************* Rescale all centroid coordinates. ********************/
  for (i=0; i<n_reg; i++) {
    centroidx[i] = centroidx[i]*latt_const + new_minx;
    centroidy[i] = centroidy[i]*latt_const + new_miny;
  }
 
  return;
}

double scale_map_factor (void)
{
  double latt_const, new_maxx, new_maxy, new_minx, new_miny;

  /* Minimum dimensions that leave enough space between map and rectangular  */
  /* boundaries.                                                             */

  new_maxx = 0.5 * ((1.0+PADDING)*map_maxx + (1.0-PADDING)*map_minx);
  new_minx = 0.5 * ((1.0-PADDING)*map_maxx + (1.0+PADDING)*map_minx);
  new_maxy = 0.5 * ((1.0+PADDING)*map_maxy + (1.0-PADDING)*map_miny);
  new_miny = 0.5 * ((1.0-PADDING)*map_maxy + (1.0+PADDING)*map_miny);
  if (map_maxx-map_minx > map_maxy-map_miny) {
    lx = L;
    latt_const = (new_maxx-new_minx) / L;
  }
  else {
    ly = L;
    latt_const = (new_maxy-new_miny) / L;
  }

  return latt_const;
}
