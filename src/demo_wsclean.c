//
// Created by Vijay Kartik on 27/03/15.
//


#include <stdio.h>
#include <stdlib.h>
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <purify_measurement.h>
#include <purify_interfaces.h>
//#include <CoreGraphics/CoreGraphics.h>
//#include <Foundation/Foundation.h>
#include <sopt_utility.h>
#include <time.h>
#include <CoreGraphics/CoreGraphics.h>
#include <GLKit/GLKit.h>

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
#include "purify_measurement.h"
#include "purify_interfaces.h"
#include "wscleaninterface.h"
#define VERBOSE 1
int main(int argc, char *argv[]) {

  int i;
  double gamma=0.01;
  double aux1, aux2;
  complex double alpha;
  complex double *y;
  //purify_visibility vis_test;
  purify_measurement_cparam param_m1;
  purify_measurement_cparam param_m2;

  //WSCLEAN INTERFACE
  void *userdata;
  purify_domain_info wscleanParams;
  purify_domain_data_format wscleanFormat;
  //WSCLEAN INTERFACE


  //Structures for sparsity operator
  sopt_wavelet_type *dict_types;
  sopt_sara_param param1;
  void *datas[1];

  //Structures for the optimization problems
  sopt_l1_sdmmparam param4;
  //sopt_l1_rwparam param5;

  int dimy = 128;
  int dimx = 128;

  // Rescale the visibilities.
  double dx = 120; // arcsec
  dx = dx / 60.0 / 60.0 / 180.0 * PURIFY_PI; // convert to radians
  double umax = 1.0 / (2.0 * dx);
  double uscale = 2.0*PURIFY_PI / umax;


  //WSCLEAN OPERATOR
  if(argc < 2) {
    printf("Syntax: %s <ms>\n", argv[0]);
    exit(-1);
  }
  wscleanParams.msPath = argv[1];
  printf("Opening: %s\n", argv[1]);
  wscleanParams.imageWidth = dimx;
  wscleanParams.imageHeight = dimy;
  wscleanParams.pixelScaleX = uscale;
  wscleanParams.pixelScaleY = uscale;
  wscleanParams.extraParameters = "-weight natural";
  purify_interface_initialiseOperator(&userdata, &wscleanParams, &wscleanFormat);
  int Ny = wscleanFormat.data_size;
  y = malloc(Ny*sizeof(*y));
  PURIFY_ERROR_MEM_ALLOC_CHECK(y);
  double *w = malloc(Ny*sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(w);
  double* weights = malloc(Ny*sizeof(*weights));
  PURIFY_ERROR_MEM_ALLOC_CHECK(weights);
  purify_interface_readData(userdata, y, weights);
  for (i=0;i<Ny;i++) y[i] *= weights[i];
  //WSCLEAN OPERATOR

  // Define dimensional parameters.
  int Nb = 9;
  int Nx = dimx * dimy;
  int Nr = Nb * Nx;

  // Memory allocation for the different variables.
  double* xout = (double*)malloc((Nx) * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(xout);
  complex double* xoutc = (complex double*)malloc((Nx) * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(xoutc);
  double* dummyr = malloc(Nr * sizeof(double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(dummyr);
  complex double* dummyc = malloc(Nr * sizeof(complex double));
  PURIFY_ERROR_MEM_ALLOC_CHECK(dummyc);

  for (i=0; i < Nx; i++){
    xoutc[i] = 0.0 + 0.0*I;
  }

  // Define parameters for the convolutional gridding.
  param_m1.nmeas = Ny;
  param_m1.ny1 = dimy;
  param_m1.nx1 = dimx;
  param_m1.ofy = 2;
  param_m1.ofx = 2;
  param_m1.ky = 24;
  param_m1.kx = 24;


  //SARA structure initialization

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


  //Scaling constants in the different representation domains

  sopt_sara_analysisop((void*)dummyc, (void*)xoutc, datas);

  for (i=0; i < Nr; i++) {
    dummyr[i] = creal(dummyc[i]);
  }

  aux2 = purify_utils_maxarray(dummyr, Nr);

  printf("Scale of gamma = %f \n\n", aux2);

  //Initial solution and weights
  for (i=0; i < Nx; i++) {
    xoutc[i] = 0.0 + 0.0*I;
  }
  for (i=0; i < Nr; i++){
    w[i] = 1.0;
  }

  printf("**********************\n");
  printf("BPSA reconstruction\n");
  printf("**********************\n");

  //Structure for the L1 solver
  param4.verbose = 2;
  param4.max_iter = 5;
  param4.gamma = gamma*aux2;//*sqrt(aux4);
  param4.rel_obj = 0.0001;
  param4.epsilon = 0.01*aux1; //sqrt(Ny + 2*sqrt(Ny))*sigma/sqrt(aux4);
  param4.epsilon_tol = 0.01;
  param4.real_data = 0;
  param4.cg_max_iter = 100;
  param4.cg_tol = 0.000001;

  sopt_l1_sdmm((void*)xoutc, Nx,
          &purify_interface_fwdOp,
          &userdata,
          &purify_interface_adjOp,
          &userdata,
          &sopt_sara_synthesisop,
          datas,
          &sopt_sara_analysisop,
          datas,
          Nr,
          (void*)y, Ny, w, param4);

  printf("Time BPSA: %f \n\n", t);


  free(y0);
  free(y);
  free(xout);
  free(xoutc);
  free(w);
  free(dummyr);
  free(dummyc);
  free(dict_types);
  sopt_sara_free(&param1);

return 0;
}
