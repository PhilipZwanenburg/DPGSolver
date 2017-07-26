// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_integration_linearization.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <complex.h>
#include <time.h>

#include "petscmat.h"

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"
#include "Test.h"

#include "test_code_integration.h"
#include "test_support.h"

#include "solver_symmetric_functions.h"
#include "finalize_LHS.h"
#include "compute_VOLUME_RLHS_DG.h"
#include "compute_FACE_RLHS_DG.h"
#include "compute_GradW_DG.h"

#include "solver.h"
#include "explicit_VOLUME_info_c.h"
#include "explicit_FACE_info_c.h"
#include "explicit_GradW_c.h"
#include "finalize_RHS_c.h"

#include "array_norm.h"

/*
 *	Purpose:
 *		Provide functions for linearization integration testing.
 *
 *	Comments:
 *		The linearization is verified by comparing with the output when using the complex step method.
 *
 *		For second order equations, the verification of the linearization of the weak gradients is also performed.
 *
 *		When collocation is enabled, symmetry of the global system matrix is lost for the diffusion operator due to the
 *		premultiplication of operators by inverse cubature weights. To recover the symmetry, the RHS/LHS terms are
 *		corrected before calling the linear solver.
 *
 *		For second order equations, the VOLUME->RHS_c contributes off-diagonal terms to the global system matrix due to
 *		the use of the fully corrected gradient. A flag is provided to avoid the computation of these terms when
 *		checking only the diagonal VOLUME/FACE contributions to the global system using compute_A_cs with assembly_type
 *		equal to 1 or 2.
 *
 *		For 2nd order equations, modifying What in a given VOLUME results in a modification of Qhat in all adjacent
 *		VOLUMEs, which in turn results in explicit_VOLUME_info generating off-diagonal terms in the global system. These
 *		terms are verified with CheckOffDiagonal is set to true.
 *
 *	Notation:
 *
 *	References:
 *		Martins(2003)-The Complex-Step Derivative Approximation
 *		Squire(1998)-Using Complex Variables to Estimate Derivatives of Real Functions
 */

static void update_VOLUME_FACEs        (void);
static void compute_A_cs               (Mat *A, unsigned int const assemble_type, bool const AllowOffDiag);
static void compute_A_cs_complete      (Mat *A);
static void finalize_LHS_Qhat          (Mat*const A, unsigned int const assemble_type, unsigned int const dim);
static void compute_A_Qhat_cs          (Mat*const A, unsigned int const assemble_type, unsigned int const dim);
static void compute_A_Qhat_cs_complete (Mat*const A, unsigned int const dim);

static void set_test_linearization_data(struct S_linearization *const data, char const *const TestName)
{
	// default values
	TestDB.CheckOffDiagonal = 0; // Should not be modified.

	data->PrintEnabled           = false; // Here used to output the matrix to a file when enabled.
	data->PrintTimings           = false;
	data->CheckFullLinearization = true;
	data->CheckWeakGradients     = false;

	data->PG_add        = 1;
	data->IntOrder_mult = 2;
	data->IntOrder_add  = 0;

	data->PGlobal = 3;
	data->ML      = 0;

	data->Nref        = 2;
	data->update_argv = true;

	strcpy(data->argvNew[1],TestName);
	if (strstr(TestName,"Advection")) {
		; // Do nothing
	} else if (strstr(TestName,"Poisson")) {
		; // Do nothing
	} else if (strstr(TestName,"Euler")) {
//		data->PrintTimings = true;
		data->update_argv = false;
	} else if (strstr(TestName,"NavierStokes")) {
//		data->PrintTimings = true;
		data->CheckWeakGradients = true;
	} else {
		printf("%s\n",TestName); EXIT_UNSUPPORTED;
	}
}

static void check_passing(struct S_linearization const *const data, unsigned int *pass)
{
	double diff_cs = 0.0;
	if (data->CheckFullLinearization)
		diff_cs = PetscMatAIJ_norm_diff_d(DB.dof,data->A_cs,data->A_csc,"Inf",false);

	if (diff_cs > 2e1*EPS)
		*pass = 0;

	double const diff_LHS = PetscMatAIJ_norm_diff_d(DB.dof,data->A_cs,data->A,"Inf",true);
	if (diff_LHS > 2e1*EPS)
		*pass = 0;

	if (DB.Symmetric) {
		PetscBool Symmetric = 0;
		MatIsSymmetric(data->A,1e1*EPS,&Symmetric);
		if (!Symmetric) {
			*pass = 0;
			printf("Failed symmetric.\n");
		}
	}

	if (!(*pass)) {
		if (data->CheckFullLinearization)
			printf("diff_cs:  %e\n",diff_cs);
		printf("diff_LHS: %e\n",diff_LHS);
	}
}

static void update_times (clock_t *const tc, clock_t *const tp, double *const times, unsigned int const counter,
                          bool const PrintTimings)
{
	if (!PrintTimings)
		return;

	*tc = clock();
	times[counter] = (*tc-*tp)/(1.0*CLOCKS_PER_SEC);
	*tp = *tc;
}

static void print_times (double const *const times, bool const PrintTimings)
{
	if (!PrintTimings)
		return;

	printf("Timing ratios: ");
	for (size_t i = 0; i < 2; i++)
		printf("% .2f ",times[i+1]/times[0]);
	printf("\n");
}

static Mat allocate_A (void)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	solver_info.create_RHS = false;
	initialize_petsc_structs(&solver_info);

	return solver_info.A;
}

static void compute_A (Mat A, const unsigned int CheckLevel)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	solver_info.create_RHS = false;
	solver_info.A = A;

	if (CheckLevel == 1)
		solver_info.compute_FACE = false;

	compute_RLHS(&solver_info);
}

static void assemble_A (Mat A)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	solver_info.create_RHS = false;
	solver_info.A = A;
	assemble_petsc_structs(&solver_info);
}

static void destroy_A (Mat A)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	solver_info.create_RHS = false;
	solver_info.A = A;
	destroy_petsc_structs(&solver_info);
}

void test_linearization(struct S_linearization *const data, char const *const TestName)
{
	/*
	 *	Expected Output:
	 *		Correspondence of LHS matrices computed using complex step and exact linearization.
	 *
	 *	Comments:
	 *
	 *		DG:
	 *
	 *		By default, the complete linearization is checked here. However, linearizations of the individual
	 *		contributions listed below may be checked separately:
	 *			0) complete check (default)
	 *			1) diagonal VOLUME contributions
	 *			2) diagonal FACE contributions
	 *			3) off-diagonal FACE (and VOLUME for 2nd order equations) contributions
	 *		Further, the complete linearization can be checked either through the assembly of the individual
	 *		contributions or more simply by using the assembled RHS directly.
	 *
	 *
	 *		HDG:
	 *
	 *		Similar functionality to that of the DG scheme is provided.
	 */

	set_test_linearization_data(data,TestName);

	int  const               nargc   = data->nargc;
	char const *const *const argvNew = (char const *const *const) data->argvNew;

	bool const CheckFullLinearization = data->CheckFullLinearization,
	           CheckWeakGradients     = data->CheckWeakGradients,
	           PrintEnabled           = data->PrintEnabled;

	unsigned int const Nref        = data->Nref,
	                   update_argv = data->update_argv;

	TestDB.PGlobal       = data->PGlobal;
	TestDB.ML            = data->ML;
	TestDB.PG_add        = data->PG_add;
	TestDB.IntOrder_add  = data->IntOrder_add;
	TestDB.IntOrder_mult = data->IntOrder_mult;

	code_startup(nargc,argvNew,Nref,update_argv);
	update_VOLUME_FACEs();
	compute_dof();

	// Perturb the solution for nonlinear equations
	if (strstr(TestName,"Euler") || strstr(TestName,"NavierStokes")) {
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int const Nvar = DB.Nvar,
			                   NvnS = VOLUME->NvnS;
			for (size_t i = 0, iMax = NvnS*Nvar; i < iMax; i++)
				VOLUME->What[i] += 1e3*EPS*((double) rand() / ((double) RAND_MAX+1));
		}
	} else if (strstr(TestName,"Advection") || strstr(TestName,"Poisson")) {
		; // Do nothing
	} else {
		EXIT_UNSUPPORTED;
	}

	for (size_t nTest = 0; nTest < 2; nTest++) {
		if (nTest == 0) { // Check weak gradients
			if (!CheckWeakGradients)
				continue;

			if (!strstr(TestName,"NavierStokes"))
				EXIT_UNSUPPORTED;

			set_PrintName("linearization (weak gradient)",data->PrintName,&data->TestTRI);

			struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
			compute_GradW_DG(&solver_info);

			Mat A[DMAX]     = { NULL },
			    A_cs[DMAX]  = { NULL },
			    A_csc[DMAX] = { NULL };

			unsigned int const d = DB.d;
			if (!CheckFullLinearization) {
				for (size_t dim = 0; dim < d; dim++) {
					const unsigned int CheckLevel = 3;

					A[dim] = allocate_A();
					for (size_t i = 1; i <= CheckLevel; i++)
						finalize_LHS_Qhat(&A[dim],i,dim);
					assemble_A(A[dim]);

					A_cs[dim] = allocate_A();
					for (size_t i = 1; i <= CheckLevel; i++)
						compute_A_Qhat_cs(&A_cs[dim],i,dim);
					assemble_A(A_cs[dim]);
				}
			} else {
				for (size_t dim = 0; dim < d; dim++) {
					clock_t tp = clock(), tc;
					double times[3];
					unsigned int counter = 0;

					A[dim] = allocate_A();
					finalize_LHS_Qhat(&A[dim],0,dim);
					update_times(&tc,&tp,times,counter++,data->PrintTimings);

					A_cs[dim] = allocate_A();
					compute_A_Qhat_cs(&A_cs[dim],0,dim);
					update_times(&tc,&tp,times,counter++,data->PrintTimings);

					A_csc[dim] = allocate_A();
					compute_A_Qhat_cs_complete(&A_csc[dim],dim);
					update_times(&tc,&tp,times,counter++,data->PrintTimings);

					print_times(times,data->PrintTimings);
				}
			}

			unsigned int pass = 1;
			for (size_t dim = 0; dim < d; dim++) {
				data->A     = A[dim];
				data->A_cs  = A_cs[dim];
				data->A_csc = A_csc[dim];

				check_passing(data,&pass);
			}
			test_print2(pass,data->PrintName);

			for (size_t dim = 0; dim < d; dim++) {
				destroy_A(A[dim]);
				destroy_A(A_cs[dim]);
				destroy_A(A_csc[dim]);
			}
		} else { // Standard linearization
			set_PrintName("linearization",data->PrintName,&data->TestTRI);

			Mat A     = allocate_A(),
			    A_cs  = allocate_A(),
			    A_csc = allocate_A();

			if (!CheckFullLinearization) {
				correct_collocated_for_symmetry();

				// Note: Poisson fails symmetric for CheckLevel = 1 and this is as expected. There is a FACE
				//       contribution from Qhat which has not yet been balanced by the FACE contribution from the 2nd
				//       equation.
				const unsigned int CheckLevel = 3;
				compute_A(A,CheckLevel);
				assemble_A(A);

				bool AllowOffDiag = 0;
				if (CheckLevel == 3)
					AllowOffDiag = 1;

				for (size_t i = 1; i <= CheckLevel; i++)
					compute_A_cs(&A_cs,i,AllowOffDiag);
				assemble_A(A_cs);
			} else {
				clock_t tp = clock(), tc;
				double times[3];
				unsigned int counter = 0;

				compute_A(A,3);
				assemble_A(A);
				update_times(&tc,&tp,times,counter++,data->PrintTimings);

				compute_A_cs(&A_cs,0,true);
				update_times(&tc,&tp,times,counter++,data->PrintTimings);

				compute_A_cs_complete(&A_csc);
				update_times(&tc,&tp,times,counter++,data->PrintTimings);

				print_times(times,data->PrintTimings);
			}

			if (PrintEnabled) {
				// See commented MatView options in solver_implicit to output in more readable format.
				MatView(A,PETSC_VIEWER_STDOUT_SELF);
				MatView(A_cs,PETSC_VIEWER_STDOUT_SELF);
			}

			data->A = A; data->A_cs = A_cs; data->A_csc = A_csc;

			unsigned int pass = 1;
			check_passing(data,&pass);
			test_print2(pass,data->PrintName);

			destroy_A(A);
			destroy_A(A_cs);
			destroy_A(A_csc);
		}
	}
	code_cleanup();
}

static void flag_nonzero_LHS(unsigned int *A)
{
	unsigned int const Nvar = DB.Nvar,
	                   dof  = DB.dof;
	unsigned int       NvnS[2], IndA[2];
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {

			// Diagonal entries
			for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				if (VOLUME->indexg != VOLUME2->indexg)
					continue;

				IndA[1] = VOLUME2->IndA;
				NvnS[1] = VOLUME2->NvnS;

				size_t const Indn = (IndA[0]+i)*dof;
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					size_t const m = IndA[1]+j;
					A[Indn+m] = 1;
				}

			}

			// Off-diagonal entries
			for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
			for (size_t side = 0; side < 2; side++) {
				struct S_VOLUME *VOLUME2;

				if (FACE->Boundary)
					continue;

				if (side == 0) {
					if (VOLUME->indexg != FACE->VL->indexg)
						continue;
					VOLUME2 = FACE->VR;
				} else {
					if (VOLUME->indexg != FACE->VR->indexg)
						continue;
					VOLUME2 = FACE->VL;
				}

				IndA[1] = VOLUME2->IndA;
				NvnS[1] = VOLUME2->NvnS;

				size_t const Indn = (IndA[0]+i)*dof;
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					size_t const m = IndA[1]+j;
					A[Indn+m] = 1;
				}
			}}
		}
	}
}

static void update_VOLUME_FACEs (void)
{
	/*
	 *	Purpose:
	 *		For each VOLUME, set the list of pointers to all adjacent FACEs.
	 */

	// Ensure that the lists are reset (in case the mesh has been updated)
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		for (size_t i = 0; i < NFMAX*NSUBFMAX; i++)
			VOLUME->FACE[i] = NULL;
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		update_data_FACE(FACE);

		struct S_SIDE_DATA data = FACE->data[0];
		data.VOLUME->FACE[data.Indsfh] = FACE;

		if (!FACE->Boundary) {
			data = FACE->data[1];
			data.VOLUME->FACE[data.Indsfh] = FACE;
		}
	}
}

static void compute_A_cs(Mat *const A, unsigned int const assemble_type, bool const AllowOffDiag)
{
	if (AllowOffDiag) {
		if (DB.Viscous)
			TestDB.CheckOffDiagonal = 1;
	}

	if (assemble_type == 0) {
		compute_A_cs(A,1,AllowOffDiag);
		compute_A_cs(A,2,AllowOffDiag);
		compute_A_cs(A,3,AllowOffDiag);

		assemble_A(*A);
		return;
	}

	unsigned int const Nvar = DB.Nvar;
	unsigned int       NvnS[2], IndA[2];
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		NvnS[0] = VOLUME->NvnS;

		// Note: Initialize with zeros for linear cases.
		if (strstr(DB.TestCase,"Advection") || strstr(DB.TestCase,"Poisson")) {
			if (VOLUME->What_c)
				free(VOLUME->What_c);

			VOLUME->What_c = calloc(NvnS[0]*Nvar , sizeof *(VOLUME->What_c));
		} else if (strstr(DB.TestCase,"Euler") || strstr(DB.TestCase,"NavierStokes")) {
			if (VOLUME->What_c)
				free(VOLUME->What_c);
			VOLUME->What_c = malloc(NvnS[0]*Nvar * sizeof *(VOLUME->What_c));

			for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
				VOLUME->What_c[i] = VOLUME->What[i];
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			double const h = EPS*EPS;
			VOLUME->What_c[i] += h*I;

			explicit_GradW_c(VOLUME,false);
			if (assemble_type == 1) {
				explicit_VOLUME_info_c(VOLUME,0);
				correct_collocated_for_symmetry_c(VOLUME,0,true,false);
			} else if (assemble_type == 2 || assemble_type == 3) {
				explicit_FACE_info_c(VOLUME,0);
				correct_collocated_for_symmetry_c(VOLUME,false,false,true);
			} else {
				EXIT_UNSUPPORTED;
			}

			if (assemble_type == 1) {
				unsigned int const dof = DB.dof;

				bool CheckOffDiagonal = TestDB.CheckOffDiagonal;
				unsigned int *A_nz  = NULL;
				if (CheckOffDiagonal) {
					A_nz = calloc(dof*dof , sizeof *A_nz); // free
					flag_nonzero_LHS(A_nz);
				}

				for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
					if (!CheckOffDiagonal) {
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;
					}

					IndA[1] = VOLUME2->IndA;
					if (CheckOffDiagonal && !A_nz[(IndA[0]+i)*dof+IndA[1]])
						continue;

					NvnS[1] = VOLUME2->NvnS;
					unsigned int const nnz_d = VOLUME2->nnz_d;

					PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
					PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
					PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

					double complex const *const RHS_c = VOLUME2->RHS_c;
					for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
						m[j]  = IndA[1]+j;
						n[j]  = IndA[0]+i;
						vv[j] = cimag(RHS_c[j])/h;
					}

					MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
					free(m); free(n); free(vv);
				}

				if (CheckOffDiagonal)
					free(A_nz);
			} else if (assemble_type == 2 || assemble_type == 3) {
				for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
				for (size_t side = 0; side < 2; side++) {
					double complex const *RHS_c;

					struct S_VOLUME *VOLUME2;
					if (assemble_type == 2) {
						if (side == 0) {
							VOLUME2 = FACE->VL;
							RHS_c = FACE->RHSL_c;
						} else {
							if (FACE->Boundary)
								continue;
							VOLUME2 = FACE->VR;
							RHS_c = FACE->RHSR_c;
						}
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						unsigned int const nnz_d = VOLUME2->nnz_d;

						PetscInt    *const m  = malloc(NvnS[0]*Nvar * sizeof *m);  // free
						PetscInt    *const n  = malloc(NvnS[0]*Nvar * sizeof *n);  // free
						PetscScalar *const vv = malloc(NvnS[0]*Nvar * sizeof *vv); // free

						for (size_t j = 0, jMax = NvnS[0]*Nvar; j < jMax; j++) {
							m[j]  = IndA[0]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					} else if (assemble_type == 3) {
						if (FACE->Boundary)
							continue;

						if (side == 0 && VOLUME == FACE->VL) {
							VOLUME2 = FACE->VR;
							RHS_c = FACE->RHSR_c;
						} else if (side == 1 && VOLUME == FACE->VR) {
							VOLUME2 = FACE->VL;
							RHS_c = FACE->RHSL_c;
						} else {
							continue;
						}

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						unsigned int const nnz_d = VOLUME2->nnz_d;

						PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
						PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
						PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

						for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
							m[j]  = IndA[1]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(RHS_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					}
				}}
			} else {
				EXIT_UNSUPPORTED;
			}
			VOLUME->What_c[i] -= h*I;
		}
	}
}

static void compute_A_cs_complete(Mat *A)
{
	unsigned int const dof = DB.dof;

	unsigned int *A_nz = calloc(dof*dof , sizeof *A_nz); // free
	flag_nonzero_LHS(A_nz);

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const Nvar = DB.Nvar;
		unsigned int       NvnS[2], IndA[2];

		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			double const h = EPS*EPS;

			VOLUME->What_c[i] += h*I;

			explicit_GradW_c(VOLUME,0);
			explicit_VOLUME_info_c(VOLUME,0);
			explicit_FACE_info_c(VOLUME,0);
			finalize_RHS_c(VOLUME,0);

			// correct_FACE = false as FACE terms were added to VOLUME RHS in finalize_RHS_c
			correct_collocated_for_symmetry_c(VOLUME,0,true,false);

			for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				IndA[1] = VOLUME2->IndA;
				if (!A_nz[(IndA[0]+i)*dof+IndA[1]])
					continue;

				NvnS[1] = VOLUME2->NvnS;
				unsigned int const nnz_d = VOLUME2->nnz_d;

				PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
				PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
				PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

				double complex const *const RHS_c = VOLUME2->RHS_c;
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					m[j]  = IndA[1]+j;
					n[j]  = IndA[0]+i;
					vv[j] = cimag(RHS_c[j])/h;
				}

				MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
				free(m); free(n); free(vv);
			}

			VOLUME->What_c[i] -= h*I;
		}
	}
	free(A_nz);

	assemble_A(*A);
}

static void finalize_LHS_Qhat(Mat *const A, unsigned int const assemble_type, unsigned int const dim)
{
	/*
	 *	Purpose:
	 *		Assembles the various contributions to the weak gradient in A.
	 *
	 *	Comments:
	 *		assemble_type is a flag for which parts of the global matrix are assembled.
	 */

	unsigned int Nvar = DB.Nvar,
	             Neq  = DB.Neq;

	switch (assemble_type) {
	default: // 0
		finalize_LHS_Qhat(A,1,dim);
		finalize_LHS_Qhat(A,2,dim);
		finalize_LHS_Qhat(A,3,dim);

		assemble_A(*A);
		break;
	case 1: // diagonal VOLUME contributions
		for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			unsigned int const IndA = VOLUME->IndA,
			                   NvnS = VOLUME->NvnS;

			PetscInt *const m = malloc(NvnS * sizeof *m), // free
			         *const n = malloc(NvnS * sizeof *n); // free

			double *zeros = calloc(NvnS*NvnS , sizeof *zeros); // free

			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const Indm = IndA + eq*NvnS;
				for (size_t i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (size_t var = 0; var < Nvar; var++) {
					size_t const Indn = IndA + var*NvnS;
					for (size_t i = 0; i < NvnS; i++)
						n[i] = Indn+i;

					PetscScalar const *vv;
					if (eq == var)
						vv = VOLUME->QhatV_What[dim];
					else
						vv = zeros;

					MatSetValues(*A,NvnS,m,NvnS,n,vv,ADD_VALUES);
				}
			}
			free(m);
			free(n);
			free(zeros);
		}
		break;
	case 2: // diagonal FACE contributions
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		for (size_t side = 0; side < 2; side++) {
			double const *QhatF_What = NULL;

			struct S_VOLUME const *VOLUME;
			if (side == 0) {
				VOLUME = FACE->VL;
				QhatF_What = FACE->QhatL_WhatL[dim];
			} else {
				if (FACE->Boundary)
					continue;
				VOLUME = FACE->VR;
				QhatF_What = FACE->QhatR_WhatR[dim];
			}

			unsigned int IndA = VOLUME->IndA,
			             NvnS = VOLUME->NvnS;

			PetscInt *const m = malloc(NvnS * sizeof *m), // free
			         *const n = malloc(NvnS * sizeof *n); // free

			double *zeros = calloc(NvnS*NvnS , sizeof *zeros); // free

			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const Indm = IndA + eq*NvnS;
				for (size_t i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (size_t var = 0; var < Nvar; var++) {
					size_t const Indn = IndA + var*NvnS;
					for (size_t i = 0; i < NvnS; i++)
						n[i] = Indn+i;

					PetscScalar const *vv;
					if (FACE->Boundary) {
						vv = &QhatF_What[(eq*Nvar+var)*NvnS*NvnS];
					} else {
						if (eq == var)
							vv = &QhatF_What[0];
						else
							vv = zeros;
					}

					MatSetValues(*A,NvnS,m,NvnS,n,vv,ADD_VALUES);
				}
			}
			free(m);
			free(n);
			free(zeros);
		}}
		break;
	case 3: // off-diagonal contributions
		for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		for (size_t side = 0; side < 2; side++) {
			if (FACE->Boundary)
				continue;

			double const *QhatF_What = NULL;

			struct S_VOLUME const *VOLUME, *VOLUME2;
			if (side == 0) {
				VOLUME  = FACE->VR;
				VOLUME2 = FACE->VL;
				QhatF_What = FACE->QhatR_WhatL[dim];
			} else {
				VOLUME  = FACE->VL;
				VOLUME2 = FACE->VR;
				QhatF_What = FACE->QhatL_WhatR[dim];
			}

			unsigned int const IndA  = VOLUME->IndA,
			                   IndA2 = VOLUME2->IndA,
			                   NvnS  = VOLUME->NvnS,
			                   NvnS2 = VOLUME2->NvnS;

			PetscInt *const m = malloc(NvnS  * sizeof *m), // free
			         *const n = malloc(NvnS2 * sizeof *n); // free

			double *zeros = calloc(NvnS*NvnS2 , sizeof *zeros); // free

			for (size_t eq = 0; eq < Neq; eq++) {
				size_t const Indm = IndA + eq*NvnS;
				for (size_t i = 0; i < NvnS; i++)
					m[i] = Indm+i;

				for (size_t var = 0; var < Nvar; var++) {
					size_t const Indn = IndA2 + var*NvnS2;
					for (size_t i = 0; i < NvnS2; i++)
						n[i] = Indn+i;

					PetscScalar const *vv;
					if (eq == var)
						vv = &QhatF_What[0];
					else
						vv = zeros;

					MatSetValues(*A,NvnS,m,NvnS2,n,vv,ADD_VALUES);
				}
			}
			free(m);
			free(n);
			free(zeros);
		}}
		break;
	}
}

static void compute_A_Qhat_cs(Mat *const A, unsigned int const assemble_type, unsigned int const dim)
{
	if (!strstr(DB.TestCase,"NavierStokes"))
		EXIT_UNSUPPORTED;

	if (assemble_type == 0) {
		compute_A_Qhat_cs(A,1,dim);
		compute_A_Qhat_cs(A,2,dim);
		compute_A_Qhat_cs(A,3,dim);

		assemble_A(*A);
		return;
	}

	unsigned int const Nvar = DB.Nvar;
	unsigned int       NvnS[2], IndA[2];
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		NvnS[0] = VOLUME->NvnS;

		// Note: Initialize with zeros for linear cases.
		if (strstr(DB.TestCase,"NavierStokes")) {
			if (VOLUME->What_c)
				free(VOLUME->What_c);
			VOLUME->What_c = malloc(NvnS[0]*Nvar * sizeof *(VOLUME->What_c));

			for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++)
				VOLUME->What_c[i] = VOLUME->What[i];
		} else {
			EXIT_UNSUPPORTED;
		}
	}

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			double const h = EPS*EPS;

			VOLUME->What_c[i] += h*I;

			explicit_GradW_c(VOLUME,0);

			if (assemble_type == 1) {
				for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
					if (VOLUME->indexg != VOLUME2->indexg)
						continue;

					IndA[1] = VOLUME2->IndA;
					NvnS[1] = VOLUME2->NvnS;
					unsigned int const nnz_d = VOLUME2->nnz_d;

					PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
					PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
					PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

					double complex const *const QhatV_c = VOLUME2->QhatV_c[dim];
					for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
						m[j]  = IndA[1]+j;
						n[j]  = IndA[0]+i;
						vv[j] = cimag(QhatV_c[j])/h;
					}

					MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
					free(m); free(n); free(vv);
				}
			}

			if (assemble_type == 2 || assemble_type == 3) {
				for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
				for (size_t side = 0; side < 2; side++) {
					double complex const *QhatF_c;

					struct S_VOLUME *VOLUME2;
					if (assemble_type == 2) {
						if (side == 0) {
							VOLUME2 = FACE->VL;
							QhatF_c = FACE->QhatL_c[dim];
						} else {
							if (FACE->Boundary)
								continue;
							VOLUME2 = FACE->VR;
							QhatF_c = FACE->QhatR_c[dim];
						}
						if (VOLUME->indexg != VOLUME2->indexg)
							continue;

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						unsigned int const nnz_d = VOLUME2->nnz_d;

						PetscInt    *const m  = malloc(NvnS[0]*Nvar * sizeof *m);  // free
						PetscInt    *const n  = malloc(NvnS[0]*Nvar * sizeof *n);  // free
						PetscScalar *const vv = malloc(NvnS[0]*Nvar * sizeof *vv); // free

						for (size_t j = 0, jMax = NvnS[0]*Nvar; j < jMax; j++) {
							m[j]  = IndA[0]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(QhatF_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					} else if (assemble_type == 3) {
						if (FACE->Boundary)
							continue;

						if (side == 0 && VOLUME == FACE->VL) {
							VOLUME2 = FACE->VR;
							QhatF_c = FACE->QhatR_c[dim];
						} else if (side == 1 && VOLUME == FACE->VR) {
							VOLUME2 = FACE->VL;
							QhatF_c = FACE->QhatL_c[dim];
						} else {
							continue;
						}

						IndA[1] = VOLUME2->IndA;
						NvnS[1] = VOLUME2->NvnS;
						unsigned int const nnz_d = VOLUME2->nnz_d;

						PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
						PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
						PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

						for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
							m[j]  = IndA[1]+j;
							n[j]  = IndA[0]+i;
							vv[j] = cimag(QhatF_c[j])/h;
						}

						MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
						free(m); free(n); free(vv);
					}
				}}
			}
			VOLUME->What_c[i] -= h*I;
		}
	}
}

static void compute_A_Qhat_cs_complete(Mat *const A, unsigned int const dim)
{
	if (!strstr(DB.TestCase,"NavierStokes"))
		EXIT_UNSUPPORTED;

	unsigned int const dof = DB.dof;

	unsigned int *A_nz = calloc(dof*dof , sizeof *A_nz); // free
	flag_nonzero_LHS(A_nz);

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const Nvar = DB.Nvar;
		unsigned int       NvnS[2], IndA[2];

		IndA[0] = VOLUME->IndA;
		NvnS[0] = VOLUME->NvnS;
		for (size_t i = 0, iMax = NvnS[0]*Nvar; i < iMax; i++) {
			double const h = EPS*EPS;

			VOLUME->What_c[i] += h*I;

			explicit_GradW_c(VOLUME,0);

			for (struct S_VOLUME *VOLUME2 = DB.VOLUME; VOLUME2; VOLUME2 = VOLUME2->next) {
				IndA[1] = VOLUME2->IndA;
				if (!A_nz[(IndA[0]+i)*dof+IndA[1]])
					continue;

				NvnS[1] = VOLUME2->NvnS;
				unsigned int const nnz_d = VOLUME2->nnz_d;

				PetscInt    *const m  = malloc(NvnS[1]*Nvar * sizeof *m);  // free
				PetscInt    *const n  = malloc(NvnS[1]*Nvar * sizeof *n);  // free
				PetscScalar *const vv = malloc(NvnS[1]*Nvar * sizeof *vv); // free

				double complex const *const Qhat_c = VOLUME2->Qhat_c[dim];
				for (size_t j = 0, jMax = NvnS[1]*Nvar; j < jMax; j++) {
					m[j]  = IndA[1]+j;
					n[j]  = IndA[0]+i;
					vv[j] = cimag(Qhat_c[j])/h;
				}

				MatSetValues(*A,nnz_d,m,1,n,vv,ADD_VALUES);
				free(m); free(n); free(vv);
			}

			VOLUME->What_c[i] -= h*I;
		}
	}
	free(A_nz);

	assemble_A(*A);
}
