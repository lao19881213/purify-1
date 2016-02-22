
#ifndef PURIFY_CHIRP
#define PURIFY_CHIRP
#include "purify_config.h"
#include <stdio.h>


/*!  
 * Definition of chirp operator.
 */
typedef struct  {
   /*! Number of image pixels in first dimension. */
 int nx;
  /*! Number of image pixels in second dimension. */
 int ny;
  /*! Number of measurements */
 int M;
   /*! complex values of chirpoperator in form M * nx * ny */
 complex double *cval;
} purify_chirp;

// second structure for chirp in sparse matrix format

void purify_chirp_init(purify_chirp *chirp, 				 
				 void **data, int pc, double usfac);

void purify_chirp_convolution(purify_sparsemat_row *gmat, 				 
				 void **data);

void purify_chirp_write_complex(complex double *chirp,
				 int nx, 
				 int ny, 
				 int part,
				const char*filename);

void purify_chirp_write( double *matrix, 
					int nx, 
					int ny, 
					const char*filename);

void purify_chirp_write_matrix_ascii(complex double *matrix,
				 int nx, 
				 int ny, 
				 int part, 
				 const char*filename);


#endif
