/*! 
 * \file purify_chirp.c
 * Functionality to set up chirp matrix.
 */
#include "purify_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h> 
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>
#include PURIFY_BLAS_H
#include "purify_image.h"
#include "purify_sparsemat.h"
#include "purify_visibility.h"
#include "purify_error.h"
#include "purify_types.h"
#include "purify_utils.h" 
#include "purify_measurement.h" 
#include "purify_ran.h"  
#include "purify_chirp.h"


/*!
 * 
 * 
 * \param[out] (purify_chirp*) Chirp matrix in Fourier space.
 * \param[in] data 
 * - data[0] (purify_visibility): visibiliy coverage
 * - data[2] (purify_image): input image
 * - data[3] (purify_measurment_cparams): measurement operator parameters

 * \authors Laura Wolz
 */
void purify_chirp_init(purify_chirp *chirpop, 				 
				 void **data){

	double fovx, fovy, dlx, dly, Lx, Ly, usfac, a0;
	int  M, ox, oy, nx, ny, nxo, nyo, nxou, nyou, Npix;

	//temporary storage for each line of the Chirpoperator of size 2*Nx * 2*Ny;
	complex double  *chirpfft, *chirpfftshift; 
	complex double *chirpmatrix, *temp;
	//Chirpoperator for full chirp matrix (Nx*Ny) x M

	purify_image *img;
	purify_visibility *vis;
	purify_measurement_cparam *params;
  	fftw_plan planfft;


	vis= (purify_visibility*)data[0];
	img = (purify_image*)data[1];
	params= (purify_measurement_cparam*)data[2];

	ox=params->ofx;
	oy=params->ofy;

	nx=img->nx; //dimension of original input image
	ny=img->ny;

	nxo=nx*ox; //multiplied by oversampling of DDE
	nyo=ny*oy;

	usfac=2.0; //for now times 2 ->needs to be dynamically adapted to maxmimum w
	a0=1.0/usfac;

	nxou=nxo*usfac; //upsampled and oversampled
	nyou=nyo*usfac;

	fovx=img->fov_x;
	fovy=img->fov_y;

	Lx=2.0*sin(fovx/2.0);
	Ly=2.0*sin(fovy/2.0);

//Divide by oversampling factor
	dlx=Lx/(double)(nx);// 
	dly=Ly/(double)(ny);// 

	
	M=params->nmeas;

	Npix=2*nxou*2*nyou; //total number of pixels needed for chirpmatrix

	printf("\n**********************\n");
  	printf("Chirp initialisation\n");
  	printf("**********************\n");

  	printf("Number of measurements %i \n", M);
  	printf("Field-of-view in x-dir = %f  and y-dir = %f \n", fovx, fovy);
  	printf("Length in x-dir = %f  and y-dir = %f \n", Lx, Ly);
  	printf("dl in x-dir = %f  and y-dir = %f \n", dlx, dly);


  	printf("Number of pixels for chirp %i \n", Npix);


  	chirpmatrix = (complex double*)malloc((Npix) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(chirpmatrix);

  	chirpfft = (complex double*)malloc((Npix) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(chirpfft);

  	chirpfftshift = (complex double*)malloc((Npix) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(chirpfftshift);

  	temp = (complex double*)malloc((nyou*nxou) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(temp);

  	for(int i=0; i<Npix; i++){
  		chirpfft[i]=0.0+I*0.0;
	  	chirpfftshift[i]=0.0+I*0.0;
	  	chirpmatrix[i]=0.0+I*0.0;

  	}
  	for(int i=0; i<nxou*nyou; i++){
  		temp[i]=1e5+I*1e5;
  	}

	chirpop->nx=nxou;
  	chirpop->ny=nyou;
  	chirpop->M=M;
  	chirpop->cval = (complex double*)malloc((M*nyou*nxou) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(chirpop->cval);
  
	
  	planfft = fftw_plan_dft_2d(2*nxou, 2*nyou, chirpmatrix, chirpfft, FFT_FORWARD, FFTW_MEASURE);


	double l, m , lmsq, max;
	int ii, jj;
 
	for(int p=0; p<M; p++){

	  	for(int i=0; i<2*nxou; i++){

	  		for(int j=0; j<2*nyou; j++){

	  			l =(double) ( i - nxou + 1 );
	  			l= l*dlx;
	  			m = (double) ( j - nyou +1 );
	  			m = m*dly;
	  			lmsq = pow(l,2.0)+pow(m, 2.0);


	  			// double x=-2.0*PURIFY_PI * vis->w[p] *(sqrt (1.0 - lmsq) -1.0 );
	  			double x=PURIFY_PI * vis->w[p] *lmsq;

	  			chirpmatrix[i* (2*nyou) + j ]=cexp(I* x);//do we need normalizing here?
	  		}
	  	}
	  	if(p==0){
	  		printf("uvw %f %f %f\n", vis->u[p], vis->v[p], vis->w[p]);
	  		purify_chirp_write_complex(chirpmatrix, 2*nxou, 2*nyou, 0, "data/test/chirp/chirp_0_abs.fits");
	  		purify_chirp_write_complex(chirpmatrix, 2*nxou, 2*nyou, 1, "data/test/chirp/chirp_0_real.fits");
			purify_chirp_write_complex(chirpmatrix, 2*nxou, 2*nyou, 2, "data/test/chirp/chirp_0_imag.fits");

	  	}
	  	
		fftw_execute(planfft);

		//fft gives unnormalized result
		for(int i=0; i<2*nxou; i++){
	  		for(int j=0; j<2*nyou; j++){			
	  			chirpfft[i* (2*nyou) + j]=chirpfft[i* (2*nyou) + j]/sqrt(Npix);
	  		}
		}
		
		if(p==0){
	  		purify_chirp_write_complex(chirpfft, 2*nxou, 2*nyou, 0, "data/test/chirp/chirp_fft_0_abs.fits");
	  		purify_chirp_write_complex(chirpfft, 2*nxou, 2*nyou, 1, "data/test/chirp/chirp_fft_0_real.fits");
			purify_chirp_write_complex(chirpfft, 2*nxou, 2*nyou, 2, "data/test/chirp/chirp_fft_0_imag.fits");

	  	}
		//not sure if shift works correctly
		purify_utils_fftshift_2d(chirpfftshift, chirpfft, 2*nxou, 2*nxou);
		
		if(p==0){
	  		purify_chirp_write_complex(chirpfftshift, 2*nxou, 2*nyou, 0, "data/test/chirp/chirp_fftshift_0_abs.fits");
	  		purify_chirp_write_complex(chirpfftshift, 2*nxou, 2*nyou, 1, "data/test/chirp/chirp_fftshift_0_real.fits");
			purify_chirp_write_complex(chirpfftshift, 2*nxou, 2*nyou, 2, "data/test/chirp/chirp_fftshift_0_imag.fits");

	  	}
		for(int i=0; i<2*nxou; i++){
	  		for(int j=0; j<2*nyou; j++){

	  			//masking 
	  			if(i>(nxou- nxo/2 -1 ) && i<(nxou + nxo/2 ) && j>(nyou - nyo/2 - 1 ) && j<(nyou + nyo/2 )){
	  				
	  			}else chirpfftshift[ i* (2*nyou) + j ] =0.0 + I*0.0;
	  			
	  			ii=i-( nxou/2 );
	  			jj=j-( nyou/2 );

	  			//cutting the fftchirp with dimension nxou * nyou
	  			//finding maximum of resulting chirp
	  			if(i>( nxou/2 -1 ) && i<(nxou + nxou/2 ) && j>( nyou/2 - 1) && j<(nyou + nyou/2 )){

					temp[ii* (nyou) + jj]=chirpfftshift[ i* (2*nyou) + j];

					if(ii==0 && jj==0) {
						max=cabs(temp[ii* (nyou) + jj ]);
					}else if (cabs(temp[ii* (nyou) + jj ])>max)	
						max=cabs(temp[ii* (nyou) + jj]);
				}
	  		}
	  	}
		if(p==0){
	  		purify_chirp_write_complex(temp, nxou, nyou, 0, "data/test/chirp/chirp_cut_0_abs.fits");
	  		purify_chirp_write_complex(temp, nxou, nyou, 1, "data/test/chirp/chirp_cut_0_real.fits");
			purify_chirp_write_complex(temp, nxou, nyou, 2, "data/test/chirp/chirp_cut_0_imag.fits");
			

	  	}

 		for(int i=0; i<nxou*nyou; i++){
 			//different cutting out!
 			if(cabs(temp[i])< max*1e-10) temp[i]=0.0 + I*0.0;

 			chirpop->cval[p*(nxou*nyou) + i]=temp[i];
 		}


	}


//set up sparse matrix

	fftw_destroy_plan(planfft);
	free(chirpmatrix);
	free(chirpfftshift);
	free(chirpfft);
	free(temp);
}



/*!
 * 
 * 
 * \param[out] (purify_sparsemat_row) convoloution kernel
 * \param[in] data 
 * - data[0] (purify_visibility): visibiliy coverage
 * - data[2] (purify_chirp): chirp
 * - data[3] (purify_sparsemat_row): interpolation kernel for continuous visibilities
 * \authors Laura Wolz
 */
void purify_chirp_convolution(purify_sparsemat_row *gmatout, 				 
				 void **data){

	printf("\n**********************\n");
  	printf("Chirp convolution\n");
  	printf("**********************\n");

	

	purify_visibility *vis;
	purify_chirp *chirp;
	purify_sparsemat_row *gmat;

	vis= (purify_visibility*)data[0];
	chirp= (purify_chirp*)data[1];
	gmat= (purify_sparsemat_row*)data[2];
	
	printf("rows= %i \n", gmat->nrows); 
	printf("cols= %i \n", gmat->ncols); 
	printf("number of non-zero elements= %i \n",gmat->nvals); 
	printf("real 1 or im 0  %i \n", gmat->real); 
	

	int nx=chirp->nx;
	int ny=chirp->ny;

	
	int Npix=nx*ny;

	int nxu=nx/2;
	int nyu=ny/2;

	int zerofreq=(chirp->nx)*chirp->ny/2;
	double *testmat;
	testmat= (double*)malloc((Npix) * sizeof(double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(testmat);

  	complex double *newG;
	newG = (complex double*)malloc((Npix) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(newG);

	purify_chirp_write_complex(chirp->cval, chirp->nx*chirp->ny,chirp->M, 1, "data/test/compl_chirp_op_real.fits");

	int c, rr;
	double *A;
    A = (double*)calloc(gmat->nrows * gmat->ncols, sizeof(double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(A);

/*
  // Construct explicit matrix.
  	for (c = 0; c < gmat->nrows; c++){
    	for (rr = gmat->rowptr[c]; rr < gmat->rowptr[c+1]; rr++){
		    A[c * gmat->ncols + gmat->colind[rr]] = gmat->vals[rr];
		    if(gmat->vals[rr]!=0.0) printf(" G matrix %i %i %i %e \n",c, rr,gmat->colind[rr] , gmat->vals[rr]);
		}
	}
	purify_chirp_write(A, gmat->nrows, gmat->ncols, "data/test/chirp/Gmat.fits");
*/

	for(int i=0; i<nx; i++){
		for(int j=0; j<ny; j++){

			//store the index of the elements of the non-zero elements of the shifted chirp
			//double check this loop with output etc.
			for(int ii=0; ii<nx; ii++){
				for(int jj=0; jj<ny; jj++){

					int sign, oldpixi, oldpixj;
					if(i>=nx/2) sign=-1;
					if(i<nx/2) sign=+1;

					//if inside shifted area, store OLD pixel number.
					if( ii<(i+nx/2) && ii>(i-nx/2) &&  jj<(j+ny/2) && jj>(j-ny/2) ){
						oldpixi=sign*nx/2+ii;
						oldpixj=sign*ny/2+jj;
						testmat[ii*ny+jj]= (double) oldpixi*ny+oldpixj;
					}else{//if outside new area put to -1;
						testmat[ii*ny+jj]=-1.0;
					}
					//printf("sign %i oldpix %i,%i  index in array %f \n", sign, oldpixi, oldpixj, testmat[ii*ny+jj]);
				}
			}
			if(i==0 &&j==10) 	purify_chirp_write(testmat, nx,ny, "data/test/chirp/index_shiftedchirp_matrix.fits");
			if(j==0) printf("\n in pixel row %i \n\n", i);
			//loop over all measurements
			for(int m=0; m<chirp->M; m++){

				complex double Gtemp=0.0+I*0.0;
				complex double G0temp=0.0+I*0.0, chirptemp=0.0+I*0.0;

				int from, to;
				from=gmat->rowptr[m];
				if(m==(chirp->M-1)){
					to=gmat->nvals;
				}else
					to=gmat->rowptr[m+1];

				//only loop over the non-zero gmat elements
				for(int pix=from; pix<to; pix++ ){





					 //which pixel in shifted chirp does this pixel index correspond to
					int chirppixindex=(int) testmat[gmat->colind[pix]];

					//if shifted chirp is NOT zero at that position
					if(chirppixindex>=0){
						chirptemp=chirp->cval[m*Npix + chirppixindex];

						G0temp=gmat->vals[pix]+I*0.0;

						Gtemp=Gtemp + ( G0temp * chirptemp );
					}
				}


				newG[m*Npix+i*ny+j]=Gtemp;

			}
		}
	} 
 
		//make newG into a sparse matrix again



}





/*!
 * 
 * 
 * \param[in] (complex double*) chirp matrix in 1d array .
 * \param[in] (int) nx
 * \param[in] (int) ny
 * \param[in] (int) part: 0 absolute value, 1 real part, 2 imaginary part
 * \retval error Zero return indicates no errors.
 *
 * \authors Laura Wolz
 */
void purify_chirp_write_complex(complex double *chirp, int nx, int ny, int part,  const char*filename){

	printf("writing chirp to file \n");
	purify_image img;
	purify_image_filetype filetype=1;
	int x;
	img.nx=nx;
	img.ny=ny;
	img.pix = (double*)malloc(img.nx * img.ny * sizeof(double));
	if(part==0){ //absolute value
		for(int i=0; i<nx*ny; i++) img.pix[i]=cabs(chirp[i]);
	}else if(part==1){ //real part
		for(int i=0; i<nx*ny; i++) img.pix[i]=creal(chirp[i]);

	}else if(part==2){ //imaginary part
		for(int i=0; i<nx*ny; i++) img.pix[i]=cimag(chirp[i]);

	}else{
		PURIFY_ERROR_GENERIC("Wrong part in purify_chirp_write, can only be 0,1 or 2 \n");
	}

	x=purify_image_writefile(&img, filename,filetype);
	if(x!=0) printf("Possible error when printing chirp matrix \n");
	
}

/*!
 * 
 * 
 * \param[in] (complex double*) chirp matrix in 1d array .
 * \param[in] (int) nx
 * \param[in] (int) ny
 * \param[in] (int) part: 0 absolute value, 1 real part, 2 imaginary part
 * \retval error Zero return indicates no errors.
 *
 * \authors Laura Wolz
 */
void purify_chirp_write( double *matrix, int nx, int ny, const char*filename){

	printf("writing matrix to file \n");
	purify_image img;
	purify_image_filetype filetype=1;
	int x;
	img.nx=nx;
	img.ny=ny;
	img.pix = (double*)malloc(img.nx * img.ny * sizeof(double));

	for(int i=0; i<nx*ny; i++) img.pix[i]=(matrix[i]);
	
	x=purify_image_writefile(&img, filename,filetype);

	if(x!=0) printf("Possible error when printing chirp matrix \n");
	
}


