// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_unit_jacobian_boundary.h"

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"

#include "test_code_fluxes_inviscid.h"
#include "test_support.h"
#include "array_norm.h"
#include "jacobian_boundary_conditions.h"
#include "boundary_conditions_c.h"

/*
 *	Purpose:
 *		Test correctness of implementation of functions relating to jacobians of inviscid fluxes.
 *
 *	Comments:
 *		Correctness is assessed based on comparison with derivatives computed using the complex step method.
 *
 *	Notation:
 *
 *	References:
 *		Squire(1998)-Using_Complex_Variables_to_Estimate_Derivatives_of_Real_Functions
 */

static void compute_dWdW_cs(const unsigned int Neq, const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                            double *W, double *dWdW_cs, double *nL, const char *BType)
{
	unsigned int   i, iMax, n, var, var2, dim, Nvar, NnTotal, IndW, IndWB, InddWdW;
	double         h;
	double complex Wp[Nn*Nel*Neq], *WB;

	h = EPS*EPS;

	Nvar    = Neq;
	NnTotal = Nn*Nel;

	WB = malloc(NnTotal*Nvar * sizeof *WB); // free

	for (i = 0, iMax = NnTotal*Nvar; i < iMax; i++)
		WB[i] = W[i];

	for (var = 0; var < Nvar; var++) {
		IndW = NnTotal*var;
		for (n = 0; n < NnTotal; n++) {
			Wp[IndW+n] += h*I;

		if (strstr(BType,"SlipWall"))
			boundary_SlipWall_c(Nn,Nel,Wp,WB,nL,d);
		else
			printf("Error: Unsupported BType.\n"), EXIT_MSG;

		for (var2 = 0; var2 < Nvar; var2++) {
			IndWB = NnTotal*var2;
			InddWdW = (var*Nvar+var2)*NnTotal;
			for (n = 0; n < NnTotal; n++)
				dWdW_cs[InddWdW+n] = cimag(WB[IndWB+n])/h;
		}

		for (n = 0; n < NnTotal; n++) {
			Wp[IndW+n] -= h*I;
	}
	free(WB);
}

static unsigned int compare_jacobian_boundary(const unsigned int Nn, const unsigned int Nel, const unsigned int d,
                                              const unsigned int Neq, double *W, double *nL, const char *BType)
{
	unsigned int pass = 0;

	unsigned int NnTotal, Nvar;
	double       *dWdW, *dWdW_cs;

	NnTotal = Nn*Nel;
	Nvar    = Neq;

	dWdW    = malloc(NnTotal*Nvar*Neq * sizeof *dWdW);    // free
	dWdW_cs = malloc(NnTotal*Nvar*Neq * sizeof *dWdW_cs); // free

	if (strstr(BType,"SlipWall"))
		jacobian_boundary_SlipWall(Nn,Nel,W,dWdW,nL,d);
	else
		printf("Error: Unsupported BType.\n"), EXIT_MSG;

	compute_dWdW_cs(Neq,Nn,Nel,d,W,dWdW_cs,nL,BType);

	if (array_norm_diff_d(NnTotal*Nvar*Neq,dWdW,dWdW_cs,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	free(dWdW);
	free(dWdW_cs);

	return pass;
}

void test_unit_jacobian_boundary(void)
{
	unsigned int pass;

	/*
	 *	Input:
	 *		Jacobians computed using the linearized code and the complex step method.
	 *
	 *	Expected Output:
	 *		No difference between the results.
	 *
	 */

printf("Add comment to jacobian_boundary that dWdW outputs are transposes of the standard jacobian matrix notation (Wikipedia).\n");

	unsigned int NBTypes = 2;

	char         *BType[NBTypes];
	unsigned int i, Nn, Nel, d, Neq;
	double       *W;

	for (i = 0; i < NBTypes; i++)
		BType[i] = malloc(STRLEN_MIN * sizeof *BType[i]); // free

	strcpy(BType[0],"SlipWall");
	strcpy(BType[1],"Riemann ");

	for (i = 0; i < NBTypes; i++) {
	for (d = 1; d <= 3; d++) {
		Neq = d+2;

		W    = initialize_W(&Nn,&Nel,d); // free
		nL   = initialize_n(Nn,Nel,d);   // free
		pass = compare_jacobian_boundary(Nn,Nel,d,Neq,W,nL,BType[i]);

		if (d == 1) {
			//     0         10        20        30        40        50
			printf("jacobian_boundary_%s (d = %d):              ",BType[i],d);
		} else {
			printf("                           (d = %d):              ",d);
		}
		test_print(pass);

		free(W);
		free(nL);
	}}

	for (i = 0; i < NBTypes; i++)
		free(BType[i]);
}
