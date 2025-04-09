#ifndef __ICARTOGRAM_H_
#define __ICARTOGRAM_H_
/******************************** Inclusions. ********************************/

#include <R.h>
#include <Rinternals.h>
#include <fftw3.h>
#include "typepoint.h"

/******************************** Definitions. *******************************/

/* Areas on cartogram differ at most by an absolute relative error of        */
/* MAX_PERMITTED_AREA_ERROR. That is,                                        */
/* |area_on_cartogram / target_area - 1| <= MAX_PERMITTED_AREA_ERROR.        */
extern double MAX_PERMITTED_AREA_ERROR;
extern int L;
/* Maximum dimension of the FFT lattice is L x L. */

/* If a region contains exactly zero population, it will be replaced by      */
/* MIN_POP_FAC times the smallest positive population in any region.         */
extern double MIN_POP_FAC, PADDING, BLUR_WIDTH;
  /* MIN_POP_FAC  (0.2)  Replace area 0 by the minimum times this. */
  /* PADDING (1.5)     Determines space between map and boundary. */
  /* BLUR_WIDTH (5e0)  Width of Gaussian blur to smoothen the density. */
#define MAX_STRING_LENGTH (1000)
/* #define PI (3.14159265358979323846264338327950288419716939937510) */

/********************************** Macros. **********************************/

#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)>(b)) ? (b) : (a))


/*********************************** Types. **********************************/

/* POINT see typepoint.h */

/* Already in R */
/* typedef enum {          /\* Declares an enumeration data type called BOOLEAN. *\/ */
/*   FALSE,                /\* FALSE = 0, TRUE = 1 *\/ */
/*   TRUE */
/* } BOOLEAN;  */

/***************************** Global variables. *****************************/

/* Variables for map. */

extern double *area_err, *cart_area, map_maxx, map_maxy, map_minx, map_miny,
  *target_area, *bbox, *coordvertices;
extern int max_id, min_id, n_poly, *n_polyinreg, n_reg, **polyinreg;
extern POINT **cartcorn, **origcorn, **polycorn, *proj, *proj2, *proj3;

/* Variables for digitizing the density. */

extern double *rho_ft, *rho_init;
extern fftw_plan plan_fwd;
extern int lx, ly;

/**************************** Function prototypes. ***************************/
SEXP cartogramR (SEXP rcentroidx, SEXP rcentroidy, SEXP rygeomd, SEXP rvarregion, 
		 SEXP rnb_polyinreg, SEXP rn_polycorn,
		 SEXP rdimpoly, SEXP rbbox, SEXP rparamsdouble, SEXP rparamsint,
		 SEXP roptions, SEXP rmultipoly) ;
void set_inside_values_for_polygon (int region, int n_polycorn,
				    POINT *polycorn, int **inside);
double polygon_area (int ncrns, POINT *polygon);
void fill_with_density1 (double* centroidx, double* centroidy,
			 int* n_polycorn, double* varregion,
			 int* nb_polyinreg, int* options,double*  original_area);
void fill_with_density2 (int* n_polycorn, int *options);
double interpol (double x, double y, double *grid, char zero, int* options, int* error_ptr) ;
void read_gen (int* nb_polyinreg, int* options);
void ffb_integrate (int* options, int* error_ptr);
void diff_integrate (int* options, int* error_ptr);
void project (double* centroidx, double* centroidy, Rboolean proj_graticule,
	      int* options, int* error_ptr, int* n_polycorn, Rboolean gridexport);
double max_area_err (double *area_err, double *cart_area,  int* n_polycorn,
		     POINT **corn, double *sum_cart_area);
double max_absarea_err (double *area_err, double *cart_area,  int* n_polycorn,
		     POINT **corn, double *sum_cart_area);
/* IO */
void inv_rescale_map (double* centroidx, double* centroidy, int* n_polycorn, int* options);
double scale_map_factor (void);
#endif
