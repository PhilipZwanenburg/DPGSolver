// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "test_integration_poisson.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "Parameters.h"
#include "Test.h"
#include "S_DB.h"

#include "petscmat.h"

/*
#include "test_code_integration.h"
#include "compute_errors.h"
*/
/*
#include "S_VOLUME.h"

#include "test_support.h"
#include "compute_errors.h"
#include "adaptation.h"
#include "array_norm.h"
*/

/*
 *	Purpose:
 *		Test various aspects of the Poisson solver implementation:
 *			1) Linearization
 *			2) Optimal convergence orders
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void compute_A_cs(Mat *A, Vec *b, Vec *x, const unsigned int assemble_type)
{
	if (!assemble_type) {
		compute_A_cs(A,b,x,1);
		compute_A_cs(A,b,x,2);
		compute_A_cs(A,b,x,3);

		finalize_Mat(A,1);
		return;
	}

	// Initialize DB Paramters
	unsigned int Nvar = DB.Nvar;

	// Standard datatypes
	unsigned int   i, j, iMax, jMax, side, NvnS[2], IndA[2], nnz_d;
	double         h;
	double complex *RHS_c;

	struct S_VOLUME *VOLUME, *VOLUME2;
	struct S_FACET  *FACET;

	PetscInt    *m, *n;
	PetscScalar *vv;

	h = EPS*EPS;

	if (*A == NULL)
		initialize_KSP(A,b,x);

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		NvnS[0] = VOLUME->NvnS;

		if (VOLUME->uhat_c)
			free(VOLUME->uhat_c);

		VOLUME->uhat_c = calloc(NvnS[0]*Nvar , sizeof *(VOLUME->uhat_c));
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			VOLUME->uhat_c[i] += h*I;

			compute_qhat_VOLUME_c();
			compute_qhat_FACET_c();
			finalize_qhat_c();
			switch (assemble_type) {
			default: // 0
				compute_uhat_VOLUME_c();
				compute_uhat_FACET_c();
				break;
			case 1:
				compute_uhat_VOLUME_c();
				break;
			case 2:
			case 3:
				compute_uhat_FACET_c();
				break;
			}

			if (assemble_type == 1) {
				for (VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
					if (VOLUME->indexg != VOLUME2->indexg)
						continue;

					IndA[1] = VOLUME2->IndA;
					NvnS[1] = VOLUME2->NvnS;
					nnz_d   = VOLUME2->nnz_d;

					m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
					n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
					vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

					RHS_c = VOLUME2->RHS_c;
					for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
						m[j]  = IndA[1]+j;
						n[j]  = IndA[0]+i;
						vv[j] = cimag(RHS_c[j])/h;
					}

					MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
					free(m); free(n); free(vv);
				}
			}

			if (assemble_type == 2 || assemble_type == 3) {
				for (FACET = DB.FACET; FACET; FACET = FACET->next) {
				for (side = 0; side < 2; side++) {
					if (assemble_type == 2) {
						if (side == 0) {
							VOLUME2 = FACET->VIn;
							RHS_c = FACET->RHSIn_c;
						} else {
							if (FACET->Boundary)
								continue;
							VOLUME2 = FACET->VOut;
							RHS_c = FACET->RHSOut_c;
						}
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						nnz_d   = VOLUME2->nnz_d;

						m  = malloc(NvnS[0]*Nvar * sizeof *m);  // free
						n  = malloc(NvnS[0]*Nvar * sizeof *n);  // free
						vv = malloc(NvnS[0]*Nvar * sizeof *vv); // free

						for (j = 0, jMax = NvnS[0]*Nvar; j < jMax; j++) {
							m[j]  = IndA[0]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					} else if (assemble_type == 3) {
						if (FACET->Boundary)
							continue;

						if (side == 0) {
							if (VOLUME->indexg != FACET->VIn->indexg)
								continue;
							VOLUME2 = FACET->VOut;
							RHS_c = FACET->RHSOut_c;
						} else {
							if (VOLUME->indexg != FACET->VOut->indexg)
								continue;
							VOLUME2 = FACET->VIn;
							RHS_c = FACET->RHSIn_c;
						}

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						nnz_d   = VOLUME2->nnz_d;

						m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
						n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
						vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

						for (j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
							m[j]  = IndA[1]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					}
				}}
			}
			VOLUME->uhat_c[i] -= h*I;
		}
	}
}

void test_integration_poisson(int nargc, char **argv)
{
	unsigned int pass;
	char         **argvNew;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free

	strcpy(argvNew[0],argv[0]);

	TestDB.TestCase = malloc(STRLEN_MAX * sizeof *(TestDB.TestCase)); // free

	/*
	 *	Input:
	 *
	 *		ToBeModified
	 *
	 *	Expected Output:
	 *
	 *		ToBeModified
	 *
	 */

	unsigned int P, ML, PMin, PMax, MLMin, MLMax;

	Mat A = NULL, A_cs = NULL, A_csc = NULL;
	Vec b = NULL, b_cs = NULL, b_csc = NULL,
	    x = NULL, x_cs = NULL, x_csc = NULL;

	// **************************************************************************************************** //
	// TRIs (change to Mixed) ToBeModified
	strcpy(argvNew[1],"test/Test_poisson_TRI");
	strcpy(TestDB.TestCase,"Poisson");

	TestDB.PG_add = 0;
	TestDB.IntOrder_mult = 2;

	// Linearization
	TestDB.PGlobal = DB.PGlobal;
	TestDB.ML      = DB.ML;

	code_startup(nargc,argvNew,0,1);

	implicit_info_Poisson();

	finalize_LHS(&A,&b,1);
//	finalize_LHS(&A,&b,2);
//	finalize_LHS(&A,&b,3);
	finalize_Mat(&A,1);

	compute_A_cs(&A_cs,&b_cs,&x_cs,1);
//	compute_A_cs(&A_cs,&b_cs,&x_cs,2);
//	compute_A_cs(&A_cs,&b_cs,&x_cs,3);
	finalize_Mat(&A_cs,1);

	MatView(A,PETSC_VIEWER_STDOUT_SELF);
	MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);
	EXIT_MSG;

//	finalize_LHS(&A,&b,&x,0);
//	compute_A_cs(&A_cs,&b_cs,&x_cs,0);
//	compute_A_cs_complete(&A_csc,&b_csc,&x_csc);

//	MatView(A_csc,PETSC_VIEWER_STDOUT_SELF);

	pass = 0;
	if (PetscMatAIJ_norm_diff_d(DB.dof,A,A_cs,"Inf")  < EPS &&
	    PetscMatAIJ_norm_diff_d(DB.dof,A,A_csc,"Inf") < EPS)
		pass = 1, TestDB.Npass++;

	//     0         10        20        30        40        50
	printf("Linearization Poisson (2D - TRI  ):              ");
	test_print(pass);

	finalize_ksp(&A,&b,&x,2);
	finalize_ksp(&A_cs,&b_cs,&x_cs,2);
	finalize_ksp(&A_csc,&b_csc,&x_csc,2);
	code_cleanup();

	// Convergence orders
	PMin = 0;  PMax = 4;
	MLMin = 0; MLMax = 4;

	for (P = PMin; P <= PMax; P++) {
	for (ML = MLMin; ML <= MLMax; ML++) {
		TestDB.PGlobal = P;
		TestDB.ML = ML;

		code_startup(nargc,argvNew,0,1);

		solver_poisson();
		compute_errors_global();

		code_cleanup();
	}}
	// test with various boundary conditions and fluxes (ToBeDeleted)


	free(argvNew[0]); free(argvNew[1]); free(argvNew);
	free(TestDB.TestCase);
}
