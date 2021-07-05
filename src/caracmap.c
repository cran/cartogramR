/******************************** Inclusions. **********************/
#include <math.h>

/********************************** Functions *******************************/
/** \fn caract_map
 *  \brief Calculates map caracteristics
 *
 * \param  caracmapd : pointer on double, will contains results
 *                   (scale, centerx, centery)
 * \param  caracmapi : pointer on int,  will contains results
 *                    (lx. ly) number of points
 * \param  padding : double, the padding used in cartogram 
 *         Determines space between map and boundary (default to 1.5)
 * \param  LL : the value of L in cartogram  (default is 512), 
 *         must be a power of two (for fftw) 
 * \param  map_maxx : double, Bounding box
 * \param  map_maxy : double, Bounding box
 * \param  map_minx : double, Bounding box
 * \param  map_miny : double, Bounding box
 * \return void
 *******************************************************************/

 void caract_map (double *caracmapd, int *caracmapi, double padding,
		  int L, double map_maxx, double map_maxy,
		  double map_minx, double map_miny)
{
  double latt_const, new_maxx, new_maxy, new_minx, new_miny;
  int lx, ly;

  
  /* Minimum dimensions that leave enough space between map and rectangular  */
  /* boundaries.                                                             */

  new_maxx = 0.5 * ((1.0+padding)*map_maxx + (1.0-padding)*map_minx);
  new_minx = 0.5 * ((1.0-padding)*map_maxx + (1.0+padding)*map_minx);
  new_maxy = 0.5 * ((1.0+padding)*map_maxy + (1.0-padding)*map_miny);
  new_miny = 0.5 * ((1.0-padding)*map_maxy + (1.0+padding)*map_miny);  
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
  caracmapd[0] = latt_const;
  caracmapd[1] = new_minx;
  caracmapd[2] = new_miny;
  caracmapi[0] = lx;
  caracmapi[1] = ly;
  return ;
}
