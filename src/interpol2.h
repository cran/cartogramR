#ifndef __INTERPOL2_H_
#define __INTERPOL2_H_
/********************************** Macros. **********************************/
#define MAX(a,b) (((a)>(b)) ? (a) : (b))
#define MIN(a,b) (((a)>(b)) ? (b) : (a))
/**************************** Function prototypes. ***************************/
double interpol2 (double x, double y, double *grid, char zero, int* options,
		  int * errorloc, int lx, int ly) ;
#endif // __INTERPOL2_H_
