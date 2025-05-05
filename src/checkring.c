/** \brief Checks  if polygons are counter clockwise or not and
 *         corrects if asked
 *
 * For each polygon (EXT, HOLE1, HOLE2) area of EXT must be positive
 * and areas of HOLE1, HOLE2 must be negative. The result is either
 * the corrected polygons or a vector of 0/1 that indicates for each
 * line i if the polygons are in the right order (=1) or not (=0)
 *******************************************************************/

/******************************** Inclusions. **********************/

#include <R.h>
#include <Rinternals.h>


/** \fn checkringcorrect
 *  \brief Checks  if polygons are counter clockwise or not and
 *         corrects if asked
 *
 * \param  rygeom  The R list of polygons (each component i is either a
 *                     multipolygon or a polygon) of length n_rows
 * \param  rmultipoly Integer vector  of length n_rows
 *                        Indicates if the component i is a multipolygon (=1)
 *                        or a polygon (=0)
 * \return rygeom2 The corrected R list of polygons
 *******************************************************************/

SEXP checkringcorrect (SEXP rygeom, SEXP rmultipoly)
{
  /*****************************************************************************/
  /* input and output from/to R */
  /*****************************************************************************/
  /* list output is an R object */
  SEXP  rygeom2;
  /*****************************************************************************/
  /* processing input  from R */
  /*****************************************************************************/
  /* rmultipoly is a vector of 0/1 which indicates if row i is */
  /* - or a simple polygon (EXT, HOLE1, HOLE2, ...) ->  rmultipoly[i]=0*/
  /* - a multipolygon ((EXT, HOLE1, HOLE2, ...), (EXT, HOLE1, ...), ...) */
  /* (integer) */
  int *multipoly;
  rmultipoly = PROTECT(rmultipoly);
  multipoly  = INTEGER(rmultipoly);
  rygeom = PROTECT(rygeom);

  /************************************************************************/
  /* local variables */
  /************************************************************************/
  int i, j, jj, k; /* check will be initialized */
  /* process options */
  /* number of rows in y_geom */
  int n_rows;
  n_rows = length(rygeom);

  /************************************************************************/
  /* Result: rygeom2  */
  /************************************************************************/
    rygeom2 = PROTECT(duplicate(rygeom));
  /************************************************************************/
  /* Read polygon from R list to n_polycorn (and assign memory)*/
  /************************************************************************/
  int nbinpoly, nbvert, two_nbvert_minus1, nblistmulti;


  SEXP rlistcoord,  rcoordvert, rlistmulti;
  SEXP rlistcoord2,  rcoordvert2, rlistmulti2;/* rlistmulti2 will be initialized */
  double *coordvert, *coordvert2, area;
  for (i=0; i<n_rows; i++) {
    rlistcoord = PROTECT(VECTOR_ELT(rygeom, i));
    rlistcoord2 = PROTECT(VECTOR_ELT(rygeom2, i));
    nbinpoly = length(rlistcoord);
    for (j=0; j<nbinpoly; j++) {
      if (multipoly[i]==0) {
	/********************************************/
	/* simple sf polygons */
	/* get values directly it is an sf POLYGON*/
	/********************************************/
	rcoordvert = PROTECT(VECTOR_ELT(rlistcoord, j));
	rcoordvert2 = PROTECT(VECTOR_ELT(rlistcoord2, j));
	coordvert = REAL(rcoordvert);
	coordvert2 = REAL(rcoordvert2);
	nbvert = (int) length(rcoordvert)/2;
	two_nbvert_minus1= nbvert+nbvert-1;
	/********************************************/
	/* area calculus */
	/********************************************/
	area = 0;
	for (k=0; k<nbvert -1; k++) {
	  area += 0.5 * (coordvert[k+1]+coordvert[k]) *
	    (coordvert[nbvert+k+1] - coordvert[nbvert+k]);
	}
	area += 0.5 * (coordvert[0]+coordvert[nbvert-1]) *
	  (coordvert[nbvert] - coordvert[two_nbvert_minus1]);
	/********************************************/
	/* change */
	/* outer polygon => area must be positive */
	/********************************************/
	if (j==0) {
	  if (area<0) {
      for (k=0; k<nbvert; k++) {
        coordvert2[k]=coordvert[nbvert-1-k];
        coordvert2[nbvert+k]=coordvert[two_nbvert_minus1-k];
      }
	  }
	} else {
	/********************************************/
	/* holes polygon => area must be negative */
	/********************************************/
	  if (area>0) {
      for (k=0; k<nbvert; k++) {
        coordvert2[k]=coordvert[nbvert-1-k];
        coordvert2[nbvert+k]=coordvert[two_nbvert_minus1-k];
      }
	  }
	}
	UNPROTECT(1); /* rcoordvert */
	UNPROTECT(1); /* rcoordvert2 */
      } else {
	/********************************************/
	/* multi sf polygons */
	/********************************************/
	rlistmulti = PROTECT(VECTOR_ELT(rlistcoord, j));
	rlistmulti2 = PROTECT(VECTOR_ELT(rlistcoord2, j));
	nblistmulti = length(rlistmulti);
	for (jj=0; jj<nblistmulti; jj++) {
	  /********************************************/
	  /* get values it is an sf POLYGON*/
	  /********************************************/
	  rcoordvert = PROTECT(VECTOR_ELT(rlistmulti, jj));
	  rcoordvert2 = PROTECT(VECTOR_ELT(rlistmulti2, jj));
	  coordvert = REAL(rcoordvert);
	  coordvert2 = REAL(rcoordvert2);
	  nbvert = (int) length(rcoordvert)/2;
	  two_nbvert_minus1= nbvert+nbvert-1;
	  /********************************************/
	  /* area calculus */
	  /********************************************/
	  area = 0;
	  for (k=0; k<nbvert -1; k++) {
	    area += 0.5 * (coordvert[k+1]+coordvert[k]) *
	      (coordvert[nbvert+k+1] - coordvert[nbvert+k]);
	  }
	  area += 0.5 * (coordvert[0]+coordvert[nbvert-1]) *
	    (coordvert[nbvert] - coordvert[two_nbvert_minus1]);
	  /********************************************/
	  /* change */
	  /* outer polygon => area positive */
	  /********************************************/
	  if (jj==0) {
	    if (area<0) {
        for (k=0; k<nbvert; k++) {
          coordvert2[k]=coordvert[nbvert-1-k];
          coordvert2[nbvert+k]=coordvert[two_nbvert_minus1-k];
        }
	    }
	  } else {
	    /********************************************/
	    /* holes polygon => area negative */
	    /********************************************/
	    if (area>0) {
        for (k=0; k<nbvert; k++) {
          coordvert2[k]=coordvert[nbvert-1-k];
          coordvert2[nbvert+k]=coordvert[two_nbvert_minus1-k];
        }
	    }
	  }
	  UNPROTECT(1); /* rcoordvert */
	  UNPROTECT(1); /* rcoordvert2 */
	}
	UNPROTECT(1); /* rlistmulti */
	UNPROTECT(1); /* rlistmulti2 */
      }
    }
    UNPROTECT(1); /* rlistcoord */
    UNPROTECT(1); /* rlistcoord2 */
  }
  UNPROTECT(3); /* rygeom + rygeom2 + rmultipoly*/
  return rygeom2;
}


/** \fn checkringnocorrect
 *  \brief Checks  if polygons are counter clockwise or not
 *
 * \param  rygeom  The R list of polygons (each component i is either a
 *                     multipolygon or a polygon) of length n_rows
 * \param  rmultipoly Integer vector  of length n_rows
 *                        Indicates if the component i is a multipolygon (=1)
 *                        or a polygon (=0)
 * \return rcheck an integer vector  of length n_rows which indicates
 *                     if line i is in the right direction (1) or not (0)
 *******************************************************************/

SEXP checkringnocorrect (SEXP rygeom, SEXP rmultipoly)
{
  /*****************************************************************************/
  /* input and output from/to R */
  /*****************************************************************************/
  /* list output is an R object */
  SEXP  rcheck;
  /*****************************************************************************/
  /* processing input  from R */
  /*****************************************************************************/
  /* rmultipoly is a vector of 0/1 which indicates if row i is */
  /* - or a simple polygon (EXT, HOLE1, HOLE2, ...) ->  rmultipoly[i]=0*/
  /* - a multipolygon ((EXT, HOLE1, HOLE2, ...), (EXT, HOLE1, ...), ...) */
  /* (integer) */
  int *multipoly;
  rmultipoly = PROTECT(rmultipoly);
  multipoly  = INTEGER(rmultipoly);
  rygeom = PROTECT(rygeom);

  /************************************************************************/
  /* local variables */
  /************************************************************************/
  int i, j, jj, k, *check; /* check will be initialized */
  /* process options */
  /* number of rows in y_geom */
  int n_rows;
  n_rows = length(rygeom);

  /************************************************************************/
  /* Result: rygeom2  */
  /************************************************************************/
    rcheck = PROTECT(allocVector(INTSXP, n_rows));
    check = INTEGER(rcheck);
  /************************************************************************/
  /* Read polygon from R list to n_polycorn (and assign memory)*/
  /************************************************************************/
  int nbinpoly, nbvert, two_nbvert_minus1, nblistmulti;


  SEXP rlistcoord,  rcoordvert, rlistmulti;
  double *coordvert,  area;
  for (i=0; i<n_rows; i++) {
    rlistcoord = PROTECT(VECTOR_ELT(rygeom, i));
    check[i]=1;
    nbinpoly = length(rlistcoord);
    for (j=0; j<nbinpoly; j++) {
      if (multipoly[i]==0) {
	/********************************************/
	/* simple sf polygons */
	/* get values directly it is an sf POLYGON*/
	/********************************************/
	rcoordvert = PROTECT(VECTOR_ELT(rlistcoord, j));
	coordvert = REAL(rcoordvert);
	nbvert = (int) length(rcoordvert)/2;
	two_nbvert_minus1= nbvert+nbvert-1;
	/********************************************/
	/* area calculus */
	/********************************************/
	area = 0;
	for (k=0; k<nbvert -1; k++) {
	  area += 0.5 * (coordvert[k+1]+coordvert[k]) *
	    (coordvert[nbvert+k+1] - coordvert[nbvert+k]);
	}
	area += 0.5 * (coordvert[0]+coordvert[nbvert-1]) *
	  (coordvert[nbvert] - coordvert[two_nbvert_minus1]);
	/********************************************/
	/* outer polygon => area must be positive */
	/********************************************/
	if (j==0) {
	  if (area<0) {
	    check[i]=0;
	  }
	} else {
	/********************************************/
	/* holes polygon => area must be negative */
	/********************************************/
	  if (area>0) {
	    check[i]=0;
	  }
	}
	UNPROTECT(1); /* rcoordvert */
      } else {
	/********************************************/
	/* multi sf polygons */
	/********************************************/
	rlistmulti = PROTECT(VECTOR_ELT(rlistcoord, j));
	nblistmulti = length(rlistmulti);
	for (jj=0; jj<nblistmulti; jj++) {
	  /********************************************/
	  /* get values it is an sf POLYGON*/
	  /********************************************/
	  rcoordvert = PROTECT(VECTOR_ELT(rlistmulti, jj));
	  coordvert = REAL(rcoordvert);
	  nbvert = (int) length(rcoordvert)/2;
	  two_nbvert_minus1= nbvert+nbvert-1;
	  /********************************************/
	  /* area calculus */
	  /********************************************/
	  area = 0;
	  for (k=0; k<nbvert -1; k++) {
	    area += 0.5 * (coordvert[k+1]+coordvert[k]) *
	      (coordvert[nbvert+k+1] - coordvert[nbvert+k]);
	  }
	  area += 0.5 * (coordvert[0]+coordvert[nbvert-1]) *
	    (coordvert[nbvert] - coordvert[two_nbvert_minus1]);
	  /********************************************/
	  /* outer polygon => area positive */
	  /********************************************/
	  if (jj==0) {
	    if (area<0) {
	      check[i]=0;
	    }
	  } else {
	    /********************************************/
	    /* holes polygon => area negative */
	    /********************************************/
	    if (area>0) {
	      check[i]=0;
	    }
	  }
	  UNPROTECT(1); /* rcoordvert */
	}
	UNPROTECT(1); /* rlistmulti */
      }
    }
    UNPROTECT(1); /* rlistcoord */
  }
      UNPROTECT(3); /* rygeom + rcheck + rmultipoly */
    return rcheck;
}
