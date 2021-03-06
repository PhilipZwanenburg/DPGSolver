// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__array_free_h__INCLUDED
#define DPG__array_free_h__INCLUDED

#include <stddef.h>
#include <complex.h>

#include "S_OpCSR.h"
#include "matrix_structs.h"

extern void array_free2_c     (size_t const iMax, char           **A);
extern void array_free2_ui    (size_t const iMax, unsigned int   **A);
extern void array_free2_i     (size_t const iMax, int            **A);
extern void array_free2_d     (size_t const iMax, double         **A);
extern void array_free2_cmplx (size_t const iMax, double complex **A);
extern void array_free3_c     (size_t const iMax, size_t const jMax, char         ***A);
extern void array_free3_ui    (size_t const iMax, size_t const jMax, unsigned int ***A);
extern void array_free3_i     (size_t const iMax, size_t const jMax, int          ***A);
extern void array_free3_d     (size_t const iMax, size_t const jMax, double       ***A);
extern void array_free4_ui    (size_t const iMax, size_t const jMax, size_t const kMax, unsigned int ****A);
extern void array_free4_d     (size_t const iMax, size_t const jMax, size_t const kMax, double       ****A);
extern void array_free5_d     (size_t const iMax, size_t const jMax, size_t const kMax, size_t const lMax, double *****A);

extern void array_free1_CSR_d (struct S_OpCSR *A);
extern void array_free4_CSR_d (size_t const iMax, size_t const jMax, size_t const kMax, struct S_OpCSR ****A);
extern void array_free5_CSR_d (size_t const iMax, size_t const jMax, size_t const kMax, size_t const lMax, struct S_OpCSR *****A);

#endif // DPG__array_free_h__INCLUDED
