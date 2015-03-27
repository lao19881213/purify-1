//
// Created by Vijay Kartik on 25/03/15.
//

#include "purify_interfaces.h"
#include "wscleaninterface.h"


void purify_interface_initialiseOperator(void *userdata, const purify_domain_info *interfaceParams, purify_domain_data_format *interfaceFormat) {
#ifdef WSCLEAN
    wsclean_initialize(&userdata, interfaceParams, interfaceFormat);

#endif
}


void purify_interface_finaliseOperator(void *userdata) {
#ifdef WSCLEAN
    wsclean_deinitialize(userdata);
#endif
}

void purify_interface_readData(void* userData, complex double *data, double* weights) {
#ifdef WSCLEAN
    wsclean_read(userData, data, weights);
#endif
}

void purify_interface_fwdOp(void *out, void *in, void **data)
{
#ifdef WSCLEAN
    wsclean_operator_A(in, out, *data);
#endif
}

void purify_interface_adjOp(void *out, void *in, void **data)
{
#ifdef WSCLEAN
    wsclean_operator_At(in, out, *data);
#endif
}
