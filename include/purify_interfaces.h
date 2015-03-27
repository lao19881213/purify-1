//
// Created by Vijay Kartik on 25/03/15.
//

#ifndef PURIFY_INTERFACES_H_
#define PURIFY_INTERFACES_H_
#define WSCLEAN

#include "complex.h"
#include "purify_measurement.h"
#include "purify_sparsemat.h"

#include <wscleaninterface.h>

//just supporting WSClean for now
/*!
 * Structure storing parameters for other measurement operator.
 *
 */
//typedef struct {
    /*! Measurement Set filename. */
//    const char *msFilename;
    /*! Number of columns in the discrete image. */
//    unsigned int nX;
    /*! Number of rows in the discrete image. */
//    unsigned int nY;
    /*! Pixel scale in x. */
//    double pixelScaleX;
    /*! Pixel scale in y. */
//    double pixelScaleY;
    /*! String storing other flags that are read by WSClean. */
//    const char *flags;

//} purify_domain_info;

//typedef struct {
//    int dataSize;
//} purify_domain_data_format;

void purify_interface_initialiseOperator(void **userdata, const purify_domain_info *wscleanParams, purify_domain_data_format *wscleanformat);
void purify_interface_finaliseOperator(void *userdata);

void purify_interface_fwdOp(void *out, void *in, void **data);

void purify_interface_adjOp(void *out, void *in, void **data);

void purify_interface_readData(void* userData, complex double *data, double* weights);

#endif //PURIFY_INTERFACES_H_
