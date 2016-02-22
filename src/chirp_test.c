/*! 
 * \file chirp_tests.c
 * Image recontruction example from continuous visibilities using SARA
 * Test image: M31 128x128.
 * Coverage: random variable density with M=0.2N visibilities.
 *  author: Laura Wolz
 */
#include "purify_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>
#include <math.h>
#include <assert.h>
#include <time.h> 
#include <string.h>
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

#define VERBOSE 1


int main(int argc, char *argv[]) {

    printf("\n\n starting \n\n");


  int i, j, Nx, Ny, Nr, Nb;
  int seedn=54;
  double sigma;
  double a;
  double mse;
  double snr;
  double snr_out;
  double gamma=0.001;
  double aux1, aux2, aux3, aux4;
  complex double alpha;
  

  purify_image img, img_copy;
  purify_visibility_filetype filetype_vis;
  purify_image_filetype filetype_img;
  complex double *xinc;
  complex double *y0;
  complex double *y;
  complex double *noise;
  double *xout;
  double *w;
  double *error;
  complex double *xoutc;
  double *wdx;
  double *wdy;
  double *dummyr;
  complex double *dummyc;


  //parameters for the continuos Fourier Transform
  double *deconv;
  purify_sparsemat_row gmat;
  purify_visibility vis_test;
  purify_measurement_cparam param_m1;
  purify_measurement_cparam param_m2;
  complex double *fft_temp1;
  complex double *fft_temp2;
  void *datafwd[6];
  void *dataadj[6];
  fftw_plan planfwd;
  fftw_plan planadj;
  complex double *shifts;

  //Structures for sparsity operator
  sopt_wavelet_type *dict_types;
 
  sopt_sara_param param1;

  void *datas[1];
 

  //Structures for the opmization problems
  sopt_l1_sdmmparam param4;
  sopt_l1_rwparam param5;
  
clock_t start, stop;
  double t = 0.0;
  #ifdef _OPENMP 
  double start1, stop1;
  #endif
  int dimy, dimx;

  filetype_img = PURIFY_IMAGE_FILETYPE_FITS;

   //Read input image
  purify_image_readfile(&img, "data/images/M31_64pix.fits", filetype_img);
  printf("Image dimension: %i, %i \n\n", img.nx, img.ny); 



  //Image dimension of the zero padded image 
  //zero padded in real space -> over-sampling!!!
  //Dimensions should be power of 2
  dimx =img.nx;
  dimy =img.ny;
  
  double upsampling_factor=2.0;
  // upsampling
  //zero pad in Fourier
  fftw_plan planfft, planfftup;

  complex double *tempimg, *imgfft, *imgfftup, *imgup, *imgfftshift, *imgfftupshift;
  tempimg = (complex double*)malloc((dimy*dimx) * sizeof(complex double));
  imgfft = (complex double*)malloc((dimy*dimx) * sizeof(complex double));
  imgfftup = (complex double*)malloc((upsampling_factor*dimy*upsampling_factor*dimx) * sizeof(complex double));
  imgup = (complex double*)malloc((upsampling_factor*dimy*upsampling_factor*dimx) * sizeof(complex double));
  imgfftshift = (complex double*)malloc((dimy*dimx) * sizeof(complex double));
  imgfftupshift = (complex double*)malloc((upsampling_factor*dimy*upsampling_factor*dimx) * sizeof(complex double));

  planfft = fftw_plan_dft_2d(dimx, dimy, tempimg, imgfft, FFTW_FORWARD, FFTW_MEASURE);
  planfftup = fftw_plan_dft_2d(upsampling_factor*dimx, upsampling_factor*dimy, imgfftup, imgup, FFTW_BACKWARD, FFTW_MEASURE);
  
  
  for(int i =0; i<dimx*dimy; i++){
    tempimg[i]=img.pix[i]+0.0*I;
  }

  fftw_execute(planfft);

  purify_utils_fftshift_2d(imgfftshift, imgfft, dimx, dimy);

//ONLY WORKS FOR upsampling_fatcor=2 AT THE MOMENT

  for(int i=0; i<dimx; i++){
    for(int j=0; j<dimy; j++){
      int ii=i+dimx/2;
      int jj=j+dimy/2;
      imgfftupshift[ii*2*dimy+jj]=imgfftshift[i*dimy+j];

    }
  }
  purify_utils_fftshift_2d(imgfftup, imgfftupshift, upsampling_factor*dimx, upsampling_factor*dimy);

  fftw_execute(planfftup);

  img.nx=upsampling_factor*dimx;
  img.ny=upsampling_factor*dimy;
  dimx =img.nx;
  dimy =img.ny;
  free(img.pix);
  img.pix=( double*)malloc((img.ny*img.nx) * sizeof( double));

  for(int i=0; i<img.nx*img.ny;i++){
    img.pix[i]=creal(imgup[i]);
  }


//set field of view on degrees
  double FOVx, FOVy, Lx, Ly, Bx, By;
  FOVx=1.0/ 180.0 * PURIFY_PI;  //in radians!!
  FOVy=FOVx;

  Lx=2.0*sin(FOVx/2.0); Ly=Lx;

  Bx=img.nx/(2.0*Lx); By=Bx;  //Is only half of the FULL bandwidth -B to B

  // Input image.
  img.fov_x = FOVx ;
  img.fov_y = FOVy ;
 

  printf("\n Field-of-view in x-dir = %f  and y-dir = %f \n", FOVx, FOVy);
  printf("Length in x-dir = %f  and y-dir = %f \n\n", Lx, Ly);

//Define parameters

  filetype_vis = PURIFY_VISIBILITY_FILETYPE_PROFILE_VIS_NODUMMY;

  //Read coverage
  purify_visibility_readfile(&vis_test,
             "./data/images/Coverages/uv_p01.vis",
             filetype_vis); 
  
  for(int i=0; i<20; i++){
    printf("Vis = ( %e , %e )\n",vis_test.u[i], vis_test.v[i]);

  }

  //Read coverage
/*
  filetype_vis = PURIFY_VISIBILITY_FILETYPE_PROFILE_WIS;
  filetype_img = PURIFY_IMAGE_FILETYPE_FITS;

  purify_visibility_readfile(&vis_test,
             "./data/images/Coverages/uvw.txt",
             filetype_vis); 
  
  for(int i=0; i<vis_test.nmeas; i++){
    vis_test.u[i]=vis_test.u[i]*2.0*PURIFY_PI/Bx;
    vis_test.v[i]=vis_test.v[i]*2.0*PURIFY_PI/By;

  }
  for(int i=0; i<20; i++){
    printf("Vis = ( %e , %e )\n",vis_test.u[i], vis_test.v[i]);

  }
  */
  printf("Number of visibilities: %i \n\n", vis_test.nmeas);  

 

  param_m1.nmeas = vis_test.nmeas;
  param_m1.ny1 = dimy;
  param_m1.nx1 = dimx;
  param_m1.ofy = 2;
  param_m1.ofx = 2;
  param_m1.ky = 24;
  param_m1.kx = 24;

  param_m2.nmeas = vis_test.nmeas;
  param_m2.ny1 = dimy;
  param_m2.nx1 = dimx;
  param_m2.ofy = 2;
  param_m2.ofx = 2;
  param_m2.ky = 24;
  param_m2.kx = 24;

  Nb = 9;
  Nx=param_m2.ny1*param_m2.nx1;
  printf("sqrt(Nx) %i \n ", Nx);
  Nr=Nb*Nx;
  Ny=param_m2.nmeas;

  //Memory allocation for the different variables
  deconv = (double*)malloc((Nx) * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(deconv);
  xinc = (complex double*)malloc((Nx) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(xinc);
  xout = (double*)malloc((Nx) * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(xout);
  y = (complex double*)malloc((vis_test.nmeas) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(y);
  y0 = (complex double*)malloc((vis_test.nmeas) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(y0);
  noise = (complex double*)malloc((vis_test.nmeas) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(noise);
  shifts = (complex double*)malloc((vis_test.nmeas) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(shifts);
  w = (double*)malloc((Nr) * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(w);
  error = (double*)malloc((Nx) * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(error);
  xoutc = (complex double*)malloc((Nx) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(xoutc);

  wdx = (double*)malloc((Nx) * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(wdx);
  wdy = (double*)malloc((Nx) * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(wdy);

  dummyr = malloc(Nr * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(dummyr);
  dummyc = malloc(Nr * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(dummyc);

  for (i=0; i < Nx; i++){
    xinc[i] = 0.0 + 0.0*I;
  }

 
//read image into complex array (xinc) with zero imaginary part
  for (i=0; i < img.ny; i++){
    for (j=0; j < img.nx; j++){
      xinc[i+j*param_m1.ny1] = img.pix[i+j*img.ny] + 0.0*I;
    }
  }
  
  //Initialize griding matrix gmat
  start = clock();
  assert(start != -1);
  purify_measurement_init_cft(&gmat, deconv, shifts,
      vis_test.u, vis_test.v, &param_m1);
  stop = clock();
  t = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("Time initalization: %f \n\n", t);
  
  ///*****************************************************************************///
  ///*********************************** CHIRP ***********************************///
  ///*****************************************************************************///


 
  purify_sparsemat_row gmatout;

  printf(" G matrix rows cols %i %i  \n",gmat.nrows , gmat.ncols);

 //Set up the chirp operator
  void *datachirp[3];
  datachirp[0] = (void*)&vis_test;
  datachirp[1] = (void*)&img;
  datachirp[2] = (void*)&param_m1;


  purify_chirp chirpop;
  start = clock();
  assert(start != -1);
  purify_chirp_init(&chirpop, datachirp, 5, upsampling_factor);

  stop = clock();
  t = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("\n\n Time Chirp operator : %f \n\n", t);

  //Set up Gmatrix through convolution of gmat with chirp
  //save complex outcome in gmatout
  void *datachirp2[3];
  datachirp2[0] = (void*)&vis_test;
  datachirp2[1] = (void*)&chirpop;
  datachirp2[2] = (void*)&gmat;

  start = clock();
  assert(start != -1);
  
  purify_chirp_convolution(&gmatout, datachirp2);

  stop = clock();
  t = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("\n\n Time Chrip convolution : %f \n\n", t);


  ///*****************************************************************************///

  //Memory allocation for the fft
  i = Nx*param_m1.ofy*param_m1.ofx;
  fft_temp1 = (complex double*)malloc((i) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(fft_temp1);
  fft_temp2 = (complex double*)malloc((i) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(fft_temp2);

  //FFT plan  
  planfwd = fftw_plan_dft_2d(param_m1.nx1*param_m1.ofx, param_m1.ny1*param_m1.ofy, 
              fft_temp1, fft_temp1, 
              FFTW_FORWARD, FFTW_MEASURE);

  planadj = fftw_plan_dft_2d(param_m1.nx1*param_m1.ofx, param_m1.ny1*param_m1.ofy, 
              fft_temp2, fft_temp2, 
              FFTW_BACKWARD, FFTW_MEASURE);


  datafwd[0] = (void*)&param_m1;
  datafwd[1] = (void*)deconv;
  datafwd[2] = (void*)&gmatout;
  datafwd[3] = (void*)&planfwd;
  datafwd[4] = (void*)fft_temp1;
  datafwd[5] = (void*)shifts;

  dataadj[0] = (void*)&param_m2;
  dataadj[1] = (void*)deconv;
  dataadj[2] = (void*)&gmatout;
  dataadj[3] = (void*)&planadj;
  dataadj[4] = (void*)fft_temp2;
  dataadj[5] = (void*)shifts;


  printf("FFT plan done \n\n");
  
  
  start = clock();
  assert(start != -1);
  purify_measurement_cftfwd((void*)y0, (void*)xinc, datafwd);
  stop = clock();
  t = (double) (stop-start)/CLOCKS_PER_SEC;
  printf("Time forward operator: %f \n\n", t);

  
 
  ///*****************************************************************************///


//Noise realization
  //Input snr
  snr = 30.0;
  a = cblas_dznrm2(Ny, (void*)y0, 1);
  sigma = a*pow(10.0,-(snr/20.0))/sqrt(Ny);

    
  for (i=0; i < Ny; i++) {
      noise[i] = (sopt_ran_gasdev2(seedn) + sopt_ran_gasdev2(seedn)*I)*(sigma/sqrt(2));
      y[i] = y0[i] + noise[i];
  }

  printf("After noise init \n");

  //Rescaling the measurements

  aux4 = (double)Ny/(double)Nx;

  for (i=0; i < Ny; i++) {
    y[i] = y[i]/sqrt(aux4);
  }

  for (i=0; i < Nx; i++) {
    deconv[i] = deconv[i]/sqrt(aux4);
  }
  
  
  // Output image.
  img_copy.fov_x = FOVx ;
  img_copy.fov_y = FOVy ;
  img_copy.nx = param_m1.nx1;
  img_copy.ny = param_m1.ny1;
  printf("\n fovx %f \n\n ", img_copy.fov_x);
  for (i=0; i < Nx; i++){
    xoutc[i] = 0.0 + 0.0*I;
  }
  
  //Dirty image
  purify_measurement_cftadj((void*)xoutc, (void*)y, dataadj);

  printf("after adjacent operator \n");

  for (i=0; i < Nx; i++) {
    xout[i] = creal(xoutc[i]);
  }
    printf("after xout operator \n");

  aux1 = purify_utils_maxarray(xout, Nx);



   img_copy.pix = (double*)malloc((Nx) * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(img_copy.pix);

  for (i=0; i < Nx; i++){
    img_copy.pix[i] = creal(xoutc[i]);
  }
    printf("after img copy operator \n");




  purify_image_writefile(&img_copy, "./data/test/chirp5pc_up/m31dirty.fits", filetype_img);
  
  printf("Written dirty image to file \n");


  ///*****************************************************************************///

  //SARA structure initialization
    printf("sara struc \n");
    param1.ndict = Nb;
    param1.real = 0;

    dict_types = malloc(param1.ndict * sizeof(sopt_wavelet_type));
    PURIFY_ERROR_MEM_ALLOC_CHECK(dict_types);


    dict_types[0] = SOPT_WAVELET_DB1;
    dict_types[1] = SOPT_WAVELET_DB2;
    dict_types[2] = SOPT_WAVELET_DB3;
    dict_types[3] = SOPT_WAVELET_DB4;
    dict_types[4] = SOPT_WAVELET_DB5;
    dict_types[5] = SOPT_WAVELET_DB6;
    dict_types[6] = SOPT_WAVELET_DB7;
    dict_types[7] = SOPT_WAVELET_DB8;
    dict_types[8] = SOPT_WAVELET_Dirac;
    

    sopt_sara_initop(&param1, param_m1.ny1, param_m1.nx1, 4, dict_types);

    datas[0] = (void*)&param1;


//Scaling constants in the diferent representation domains
    //only for SARA

  sopt_sara_analysisop((void*)dummyc, (void*)xoutc, datas);

  for (i=0; i < Nr; i++) {
    dummyr[i] = creal(dummyc[i]);
  }

  aux2 = purify_utils_maxarray(dummyr, Nr);


  // Output image.
  img_copy.fov_x = FOVx;
  img_copy.fov_y = FOVy;
  img_copy.nx = param_m1.nx1;
  img_copy.ny = param_m1.ny1;

    printf("img copy \n");


  //Initial solution and weights
  for (i=0; i < Nx; i++) {
      xoutc[i] = 0.0 + 0.0*I;
      wdx[i] = 1.0;
      wdy[i] = 1.0;
  }
  for (i=0; i < Nr; i++){
    w[i] = 1.0;
  }
  //Copy true image in xout
  for (i=0; i < Nx; i++) {
    xout[i] = creal(xinc[i]);
  }
  printf("**********************\n");
  printf("SARA reconstruction\n");
  printf("**********************\n");


  //Structure for the L1 solver      
  param4.verbose = 2;
  param4.max_iter = 300;
  param4.gamma = gamma*aux2;
  param4.rel_obj = 0.001;
  param4.epsilon = sqrt(Ny + 2*sqrt(Ny))*sigma/sqrt(aux4);
  param4.epsilon_tol = 0.01;
  param4.real_data = 0;
  param4.cg_max_iter = 100;
  param4.cg_tol = 0.000001;

  
  //Initial solution
  for (i=0; i < Nx; i++) {
      xoutc[i] = 0.0 + 0.0*I;
  }
 

  //Structure for the RWL1 solver    
  param5.verbose = 2;
  param5.max_iter = 5;
  param5.rel_var = 0.001;
  param5.sigma = sigma*sqrt((double)Ny/(double)Nr);
  param5.init_sol = 1;

  
  #ifdef _OPENMP 
    start1 = omp_get_wtime();
  #else
    start = clock();
    assert(start != -1);
  #endif
  sopt_l1_rwsdmm((void*)xoutc, Nx,
                   &purify_measurement_cftfwd,
                   datafwd,
                   &purify_measurement_cftadj,
                   dataadj,
                   &sopt_sara_synthesisop,
                   datas,
                   &sopt_sara_analysisop,
                   datas,
                   Nr,
                   (void*)y, Ny, param4, param5);

  #ifdef _OPENMP 
    stop1 = omp_get_wtime();
    t = stop1 - start1;
  #else
    stop = clock();
    t = (double) (stop-start)/CLOCKS_PER_SEC;
  #endif

  printf("Time SARA: %f \n\n", t); 

  
    //SNR
    for (i=0; i < Nx; i++) {
        error[i] = creal(xoutc[i])-xout[i];
    }
    mse = cblas_dnrm2(Nx, error, 1);
    a = cblas_dnrm2(Nx, xout, 1);
    snr_out = 20.0*log10(a/mse);
    printf("SNR: %f dB\n\n", snr_out);

  for (i=0; i < Nx; i++){
    img_copy.pix[i] = creal(xoutc[i]);
  }
  
  purify_image_writefile(&img_copy,"./data/test/chirp5pc_up/m31sara.fits", filetype_img);

  

   //Residual image

  purify_measurement_cftfwd((void*)y0, (void*)xoutc, datafwd);
  alpha = -1.0 +0.0*I;
  cblas_zaxpy(Ny, (void*)&alpha, y, 1, y0, 1);
  purify_measurement_cftadj((void*)xinc, (void*)y0, dataadj);

  for (i=0; i < Nx; i++){
    img_copy.pix[i] = creal(xinc[i]);
  }

  purify_image_writefile(&img_copy, "./data/test/chirp5pc_up/m31sarares.fits", filetype_img);

  //Error image
  for (i=0; i < Nx; i++){
    img_copy.pix[i] = error[i];
  }

  purify_image_writefile(&img_copy,"./data/test/chirp5pc_up/m31saraerror.fits", filetype_img);


 
  //Free all memory
  purify_image_free(&img);
  purify_image_free(&img_copy);
  free(deconv);
  purify_visibility_free(&vis_test);
  free(y);
  free(xinc);
  free(xout);
  free(w);
  free(noise);
  free(y0);
  free(error);
  free(xoutc);
  free(wdx);
  free(wdy);
  free(shifts);

  sopt_sara_free(&param1);
  
  free(dict_types);

  //free(fft_temp1);
  //free(fft_temp2);
  fftw_destroy_plan(planfwd);
  fftw_destroy_plan(planadj);
  purify_sparsemat_freer(&gmat);

  free(dummyr);
  free(dummyc);


  return 0;

}















