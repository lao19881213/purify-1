/*! 
 * \file purify_chirp.c
 * Functionality to set up chirp matrix.
 */
#include "purify_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>
#include <math.h>
#include <assert.h>
#include <time.h> 
#ifdef _OPENMP 
  #include <omp.h>
#endif 

#include PURIFY_BLAS_H
#include "purify_visibility.h"
#include "purify_sparsemat.h"
#include "purify_image.h"
#include "purify_measurement.h"
#include "purify_types.h"
#include "purify_error.h"
#include "purify_ran.h"
#include "sopt_image.h"
#include "sopt_utility.h"
#include "sopt_prox.h" 
#include "sopt_tv.h"
#include "sopt_l1.h"
#include "sopt_wavelet.h"
#include "sopt_sara.h" 
#include "sopt_ran.h" 
#include "purify_utils.h" 
#include "purify_chirp.h"

//#define VERBOSE 0
/*!
 * 
 * 
 * \param[out] (purify_chirp*) Chirp matrix in Fourier space.
 * \param[in] data 
 * - data[0] (purify_visibility): visibiliy coverage
 * - data[2] (purify_image): input image
 * - data[3] (purify_measurment_cparams): measurement operator parameters
 * pc : percentage of maximum w (constant)
 * \authors Laura Wolz
 */
void purify_chirp_init(purify_chirp *chirpop, 				 
				 void **data, int pc, double usfac){

	double fovx, fovy, dlx, dly, Lx, Ly, Bx, By, a0;
	int  M, ox, oy, nx, ny, nxou, nyou, nxu, nyu, Npix;

	purify_image *img;
	purify_visibility *vis;
	purify_measurement_cparam *params;

	vis= (purify_visibility*)data[0];
	img = (purify_image*)data[1];
	params= (purify_measurement_cparam*)data[2];
	ox=params->ofx;
	oy=params->ofy;

	nx=img->nx; //dimension of original input image
	ny=img->ny;
	printf(" Original (upsampled) image dim %i %i \n", nx, ny);

	a0=1.0/usfac;

	nxu=nx;
	nyu=ny;

	

	nxou=nx*ox; //upsampled and oversampled
	nyou=ny*oy;
	printf(" After over sampling image dim %i %i \n", nxou, nyou);

	fovx=img->fov_x;
	fovy=img->fov_y;

	Lx=2.0*sin(fovx/2.0);
	Ly=2.0*sin(fovy/2.0);

	Bx=nx*a0/(2.0*Lx); By=Bx;

	dlx=Lx/(double)(nxu);
	dly=Ly/(double)(nyu);
	
	M=params->nmeas;

	Npix=nxou*nyou; //total number of pixels needed for chirpmatrix

	printf("\n**********************\n");
  	printf("Chirp initialisation\n");
  	printf("**********************\n");

	FILE *fp;

	fp=fopen("./data/test/chirp5pc/chirp_params.txt", "w");
  	fprintf(fp,"Number of measurements %i \n\n", M);
	fprintf(fp,"Original image dim %i X %i \n\n", nx, ny);
	fprintf(fp,"After over sampling image dim %i X %i \n\n", nxou, nyou);
  	fprintf(fp,"Field-of-view in x-dir = %f  and y-dir = %f \n\n", fovx, fovy);
  	fprintf(fp,"Length in x-dir = %f  and y-dir = %f \n\n", Lx, Ly);
  	fprintf(fp,"dl in x-dir = %f  and y-dir = %f \n\n", dlx, dly);
  	fprintf(fp,"Bandwidth B in x-dir = %f  and y-dir = %f \n\n", Bx, By);
	fprintf(fp," w_max= Nx_orig /Lx^2 = %f \n\n",(double) nx / pow(Lx, 2.0) );

	fclose(fp);

	chirpop->cval = (complex double*)malloc((M*nyou*nxou) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(chirpop->cval);


	chirpop->nx=nxou;
  	chirpop->ny=nyou;
  	chirpop->M=M;
 	

	//SET TO ZERO FOR TESTs
	 for(int p=0; p<M; p++) vis->w[p]=0.0; 

	//temporary storage for each line of the Chirpoperator of size 2*Nx * 2*Ny;
	complex double  *chirpfft, *chirpfftshift; 
	complex double *chirpmatrix, *temp;
	//Chirpoperator for full chirp matrix (Nx*Ny) x M

  	chirpmatrix = (complex double*)malloc((Npix) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(chirpmatrix);

  	chirpfft = (complex double*)malloc((Npix) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(chirpfft);

  	chirpfftshift = (complex double*)malloc((Npix) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(chirpfftshift);

  	temp = (complex double*)malloc((nyou*nxou) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(temp);


  	
  	for(int i=0; i<nxou*nyou; i++){
  		temp[i]=1e5+I*1e5;
  	}

 	fftw_plan planfft;

  	planfft = fftw_plan_dft_2d(nxou, nyou, chirpmatrix, chirpfft, FFTW_FORWARD, FFTW_MEASURE);
	
	for(int i=0; i<Npix; i++){
  		chirpfft[i]=0.0+I*0.0;
	  	chirpfftshift[i]=0.0+I*0.0;
	  	chirpmatrix[i]=0.0+I*0.0;
	  	if(pc==0) chirpmatrix[i]=1.0+I*0.0;
  	}

	for(int p=0; p<M; p++){

		int ii, jj;
		double l, m , lmsq, max;

		if(p%100==0) printf(" Set up w kernel for m=%i \n", p);
	  	for(int i=0; i<nxu; i++){

	  		for(int j=0; j<nyu; j++){

	  			l =(double) ( i - nxu/2 + 1 );
	  			l= l*dlx;
	  			m = (double) ( j - nyu/2 +1 );
	  			m = m*dly;
	  			lmsq = pow(l,2.0)+pow(m, 2.0);

				vis->w[p]=(double) nx / pow(Lx, 2.0) * (double)pc/100.0;

	  			//double x=-2.0*PURIFY_PI * vis->w[p] *(sqrt (1.0 - lmsq) -1.0 );
	  			double x=-PURIFY_PI * vis->w[p] *lmsq;
	  			int ii=i+ (nxou-nxu)/2;
	  			int jj=j+ (nyou-nyu)/2;
	  			chirpmatrix[ii* nyou + jj ]=cexp(I* x);//do we need normalizing here?
	  		}
	  	}
	  		  	

		fftw_execute(planfft);
		double sumchirp=0.0, testsum=0.0;
		//fft gives unnormalized result
		for(int i=0; i<nxou; i++){
	  		for(int j=0; j<nyou; j++){			
	  			chirpfft[i* nyou + j]=chirpfft[i* nyou + j]/(Npix); //CHECK IF CORRECT WHEN w!=0
	  			sumchirp+=cabs(chirpfft[i* nyou + j]);
	  		}
		}
		if(p==0) printf("sumchirp = %e \n", sumchirp);
		purify_utils_fftshift_2d(chirpfftshift, chirpfft, nxou, nxou);
		
		//maybe sparsifying thing here


 		for(int i=0; i<nxou*nyou; i++){
 			//different cutting out!
 			chirpop->cval[p*(nxou*nyou) + i]=chirpfftshift[i]/sumchirp;
 			testsum+=cabs(chirpop->cval[p*(nxou*nyou) + i]);
 		}
 		if(p==0) printf("test sumchirp = %e \n", testsum);



 		if(p==0)purify_chirp_write_matrix_ascii(chirpfftshift,  nxou,  nyou, 1, "data/test/chirp.txt");
	}
	
	

//set up sparse matrix ?

	printf("\n**********************\n");
  	printf("Chirp init complete \n");
  	printf("**********************\n");

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
	
	printf("Gmat rows= %i \n", gmat->nrows); 
	printf("Gmat cols= %i \n", gmat->ncols); 
	printf("Gmat number of non-zero elements= %i \n",gmat->nvals); 
	printf("real 1 or im 0  %i \n", gmat->real); 
	

	int nx=chirp->nx;
	int ny=chirp->ny;

	printf("Chirp rows= %i \n", nx); 
	printf("Chirp cols= %i \n", ny); 

	int Npix=nx*ny;

	
	int zerofreq=(chirp->nx)*chirp->ny/2;
	
  	complex double *newG;
	newG = (complex double*)malloc((Npix*chirp->M) * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(newG);

	complex double	*Gtemp_mat, *G_shifted;
	Gtemp_mat = (complex double*)malloc(Npix * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(Gtemp_mat);

	G_shifted = (complex double*)malloc(Npix * sizeof(complex double));
  	PURIFY_ERROR_MEM_ALLOC_CHECK(G_shifted);
	//purify_chirp_write_complex(chirp->cval, chirp->nx*chirp->ny,chirp->M, 1, "data/test/compl_chirp_op_real.fits");

	int c, rr;

	int sparsecount=0;

  	for(int i=0; i<Npix*chirp->M; i++) newG[i]=0.0+I*0.0;

	//loop over every pixel in every column
	//where the column indices would be ny*i+j
	//#pragma omp parallel for
	#ifdef _OPENMP 
	#pragma omp parallel for
	#endif

	//loop over all measurements m (i.e. every row of the G matrix )
	for(int m=0; m<chirp->M; m++){//chirp->M

		if(m%100==0) printf("In M %i \n", m);

		//loop over every pixel
		for(int i=0; i<nx; i++){//nx
			for(int j=0; j<ny; j++){ //ny

	
				complex double Gtemp=0.0+I*0.0;
				complex double G0temp=0.0+I*0.0;
				complex double chirptemp=0.0+I*0.0;

				int from, to;
				from=gmat->rowptr[m];
				if(m==(chirp->M-1)){
					to=gmat->nvals;
				}else
					to=gmat->rowptr[m+1];

				//only loop over the non-zero gmat elements
				for(int pix=from; pix<to; pix++ ){

					//express the column index as two-dimensional indices in image plane
					int ii, jj, i_fftshift, j_fftshift;

					jj= gmat->colind[pix] % ny; 
					ii= (gmat->colind[pix] - jj)/ny;

					if(ii<nx/2) i_fftshift=ii+nx/2;
					if(ii>=nx/2) i_fftshift=ii-nx/2;
					if(jj<ny/2) j_fftshift=jj+ny/2;
					if(jj>=ny/2) j_fftshift=jj-ny/2;
					
					
					int oldpixi, oldpixj;
					
					//translate the chirp matrix for m to center on the pixel (i,j)
					//store old pixel indices of Chirp 
					oldpixi=nx/2-i+i_fftshift;//sign*nx/2+ii;
					oldpixj=ny/2-j+j_fftshift;//sign*ny/2+jj;

					//index of the chirp which translates to (ii,jj)
					int chirppixindex= (double) oldpixi*ny+oldpixj;

					chirptemp=chirp->cval[m*Npix + chirppixindex];

					//only add if within the overlap between chirp and Gmat
					//no circular convolution

					if(oldpixi>=0 && oldpixi<nx){
							if(oldpixj>=0 && oldpixj<ny){

							double temp = gmat->vals[pix];//(int) (gmat->vals[pix] * pow(10,6));
							
							//temp /= pow(10,6);

							G0temp= temp +I*0.0;

							Gtemp=Gtemp + ( G0temp * chirptemp );

						}
					}		
				}
				Gtemp_mat[i*ny+j]=Gtemp;
				
				}
		}

		purify_utils_fftshift_2d(G_shifted, Gtemp_mat, nx, ny);
		for(int i=0; i<Npix; i++){
			newG[m*Npix+i]=G_shifted[i];

			if(cabs(newG[m*Npix+i])!=0.0){
				sparsecount++;
			}
		} 

	}

	printf("elements of gmat after sparse %i\n", sparsecount);

	gmatout->nrows = chirp->M;
	gmatout->ncols = Npix;
	gmatout->nvals = sparsecount;
	gmatout->real = 0;
	gmatout->cvals = (complex double*)malloc(gmatout->nvals * sizeof(complex double));
	gmatout->vals = NULL;
	PURIFY_ERROR_MEM_ALLOC_CHECK(gmatout->cvals);
	gmatout->colind = (int*)malloc(gmatout->nvals * sizeof(int));
	PURIFY_ERROR_MEM_ALLOC_CHECK(gmatout->colind);
	gmatout->rowptr = (int*)malloc((gmatout->nrows + 1) * sizeof(int));
	PURIFY_ERROR_MEM_ALLOC_CHECK(gmatout->rowptr);
	
	printf("after Alloc %i\n", sparsecount);
	

	int valcount=0, colcount=0;
	gmatout->rowptr[0]=0;

	for(int m=0; m<chirp->M; m++){ //chirp->M //ROW


		
		for(int i=0; i<nx; i++){//nx  COLUMN
			for(int j=0; j<ny; j++){ //ny 
				
				if(cabs( newG[m*Npix+i*ny+j] )!=0.0){

					gmatout->cvals[valcount]=newG[m*Npix+i*ny+j];	
					gmatout->colind[valcount]=i*ny+j; 
					valcount++;
					
				}
			}
		}		
		if(m < (chirp->M - 1) ){
		 gmatout->rowptr[m+1]=valcount;
		}

	}


	printf("finished sparsification! \n");
	
		
	//copy convolved newG in old sparse gmat
	// only for testing purposes;
 	for(int m=0; m<chirp->M; m++){ //
 		int from, to;
		from=gmat->rowptr[m];

		if(m==(chirp->M-1)){
			to=gmat->nvals;
		}else
			to=gmat->rowptr[m+1];

		for(int pix=from; pix<to; pix++ ){
	 		
	 		int ii, jj;
			jj= gmat->colind[pix] % ny;
			ii= (gmat->colind[pix] - jj)/ny;

	 		//gmat->vals[pix]=creal(newG[m*Npix+ii*ny+jj]);

	 		
		}
	}
		
	//purify_chirp_write_matrix_ascii(newG,chirp->M, Npix, 1, "./data/test/chirp5pc/Gmatconv.txt");


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
 * \param[in] (complex double*) matrix in 1d array .
 * \param[in] (int) nx
 * \param[in] (int) ny
 *
 *
 * \authors Laura Wolz
 */
void purify_chirp_write_matrix_ascii(complex double *matrix, int nx, int ny, int part, const char*filename){


	FILE *fp;
	fp=fopen(filename, "w");

	for(int i =0; i<nx; i++){
		for(int j=0; j<ny; j++){
			if(part==0){ //complex value
				fprintf(fp, "%e+i%e  ", creal(matrix[i*ny+j]), cimag(matrix[i*ny+j]));
			}else if(part==1){ //real part
				fprintf(fp, "%e  ", creal(matrix[i*ny+j]));
			}else if(part==2){ //imaginary part
				fprintf(fp, "%e  ", cimag(matrix[i*ny+j]));
			}

		}
		fprintf(fp, " \n");
	}
	fclose(fp);

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


