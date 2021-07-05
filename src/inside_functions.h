#ifndef __INSIDE_FUNCTIONS_H_
#define __INSIDE_FUNCTIONS_H_
#include "typepoint.h"
/********************************** Macros. **********************************/
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)>(b)) ? (b) : (a))
/**************************** Function prototypes. ***************************/
double min4 (double a, double b, double c, double d);
double max4 (double a, double b, double c, double d);
void set_inside_value_at_y (int region, POINT pk, POINT pn, int l, double poly_minx, int **inside);
void set_inside_values_between_points (int region, POINT pk, POINT pn, double poly_minx, int **inside);
void set_inside_values_for_polygon (int region, int n_polycorn, POINT *polycorn, int **inside);
#endif // __INSIDE_FUNCTIONS_H_
