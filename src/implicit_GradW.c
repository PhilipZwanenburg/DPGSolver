// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "implicit_GradW.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>

#include "Macros.h"
#include "Parameters.h"
#include "S_DB.h"
#include "S_VOLUME.h"
#include "S_FACE.h"

#include "solver_functions.h"
#include "matrix_functions.h"
#include "array_free.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Compute weak gradient contributions required for the computation of viscous fluxes for implicit runs.
 *
 *	Comments:
 *		See comments of explicit_GradW.c.
 *
 *	Notation:
 *
 *	References:
 */

static void implicit_GradW_VOLUME   (void);
static void implicit_GradW_FACE     (void);
static void implicit_GradW_finalize (void);

void implicit_GradW(void)
{
	if (!DB.Viscous)
		return;

	implicit_GradW_VOLUME();
	implicit_GradW_FACE();
	implicit_GradW_finalize();
}

static void implicit_GradW_VOLUME(void)
{
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	struct S_OPERATORS_V *OPS[2];

	struct S_VDATA *const VDATA = malloc(sizeof *VDATA); // free
	VDATA->OPS = (struct S_OPERATORS_V const *const *) OPS;

	for (size_t i = 0; i < 2; i++)
		OPS[i] = malloc(sizeof *OPS[i]); // free

	struct S_Dxyz *const DxyzInfo = malloc(sizeof *DxyzInfo); // free

	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		init_VDATA(VDATA,VOLUME);

		unsigned int const NvnS = VDATA->OPS[0]->NvnS,
		                   NvnI = VDATA->OPS[0]->NvnI;

		DxyzInfo->Nbf = VDATA->OPS[0]->NvnS;
		DxyzInfo->Nn  = VDATA->OPS[0]->NvnI;
		DxyzInfo->D   = (double const *const *const) VDATA->OPS[0]->D_Weak;
		DxyzInfo->C   = VOLUME->C_vI;

		double const *const ChiS_vI = VDATA->OPS[0]->ChiS_vI;

		double **const QhatV_What = VOLUME->QhatV_What;
		for (size_t dim = 0; dim < d; dim++) {
			DxyzInfo->dim = dim;
			double *const Dxyz = compute_Dxyz(DxyzInfo,d); // free

			// Note: The detJ_vI term cancels with the gradient operator (Zwanenburg(2016), eq. (B.2))
			if (DB.Collocated) { // ChiS_vI == I
				for (size_t i = 0; i < NvnS*NvnS; i++)
					QhatV_What[dim][i] = Dxyz[i];
			} else {
				mm_d(CBRM,CBNT,CBNT,NvnS,NvnS,NvnI,1.0,0.0,Dxyz,ChiS_vI,QhatV_What[dim]);
			}
			free(Dxyz);

			// Compute intermediate Qhat contribution
			mkl_dimatcopy('R','T',NvnS,NvnS,1.0,QhatV_What[dim],NvnS,NvnS);
			mm_CTN_d(NvnS,Nvar,NvnS,QhatV_What[dim],VOLUME->What,VOLUME->QhatV[dim]);

			for (size_t i = 0; i < NvnS*Nvar; i++)
				VOLUME->Qhat[dim][i] = VOLUME->QhatV[dim][i];
		}
	}

	free(VDATA);
	for (size_t i = 0; i < 2; i++)
		free(OPS[i]);
	free(DxyzInfo);
}

static void implicit_GradW_FACE(void)
{
	/*
	 *	Comments:
	 *		dnSolNumdW(L/R) may only have dependence on all solution variables through the boundary condition.
	 */

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar,
	                   Neq  = DB.Neq;

	// Standard datatypes
	struct S_OPERATORS_F *OPSL[2], *OPSR[2];
	struct S_FDATA       *const FDATAL = malloc(sizeof *FDATAL), // free
	                     *const FDATAR = malloc(sizeof *FDATAR); // free
	FDATAL->OPS = (struct S_OPERATORS_F const *const *) OPSL;
	FDATAR->OPS = (struct S_OPERATORS_F const *const *) OPSR;

	struct S_NumericalFlux *const NFluxData = malloc(sizeof *NFluxData); // free
	FDATAL->NFluxData = NFluxData;
	FDATAR->NFluxData = NFluxData;

	for (size_t i = 0; i < 2; i++) {
		OPSL[i]  = malloc(sizeof *OPSL[i]);  // free
		OPSR[i]  = malloc(sizeof *OPSR[i]);  // free
	}

	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		init_FDATA(FDATAL,FACE,'L');
		init_FDATA(FDATAR,FACE,'R');

		// Compute WL_fIL and WR_fIL (i.e. as seen from the (L)eft VOLUME)
		unsigned int const IndFType = FDATAL->IndFType,
		                   NfnI     = OPSL[IndFType]->NfnI;

		FDATAL->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAL->W_fIL)), // free
		FDATAR->W_fIL = malloc(NfnI*Nvar * sizeof *(FDATAR->W_fIL)); // free

		coef_to_values_fI(FDATAL,'W','I');
		compute_WR_fIL(FDATAR,FDATAL->W_fIL,FDATAR->W_fIL);

		// Compute numerical flux as seen from the left VOLUME
		NFluxData->WL_fIL         = FDATAL->W_fIL;
		NFluxData->WR_fIL         = FDATAR->W_fIL;
		NFluxData->nSolNum_fI     = malloc(d * sizeof *(NFluxData->nSolNum_fI));     // free
		NFluxData->dnSolNumdWL_fI = malloc(d * sizeof *(NFluxData->dnSolNumdWL_fI)); // free
		NFluxData->dnSolNumdWR_fI = malloc(d * sizeof *(NFluxData->dnSolNumdWR_fI)); // free
		for (size_t dim = 0; dim < d; dim++) {
			NFluxData->nSolNum_fI[dim]     = malloc(NfnI*Neq      * sizeof *(NFluxData->nSolNum_fI[dim]));     // free
			if (FACE->Boundary) {
				NFluxData->dnSolNumdWL_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnSolNumdWL_fI[dim])); // free
				NFluxData->dnSolNumdWR_fI[dim] = malloc(NfnI*Neq*Nvar * sizeof *(NFluxData->dnSolNumdWR_fI[dim])); // free
			} else {
				NFluxData->dnSolNumdWL_fI[dim] = malloc(NfnI * sizeof *(NFluxData->dnSolNumdWL_fI[dim])); // free
				NFluxData->dnSolNumdWR_fI[dim] = malloc(NfnI * sizeof *(NFluxData->dnSolNumdWR_fI[dim])); // free
			}
		}

		compute_numerical_solution(FDATAL,'I');
		add_Jacobian_scaling_FACE(FDATAL,'I','Q');

		finalize_QhatF_Weak(FDATAL,FDATAR,'L','E');
		finalize_QhatF_Weak(FDATAL,FDATAR,'L','I');
		if (!FACE->Boundary) {
			finalize_QhatF_Weak(FDATAL,FDATAR,'R','E');
			finalize_QhatF_Weak(FDATAL,FDATAR,'R','I');
		}

		free(FDATAL->W_fIL);
		free(FDATAR->W_fIL);
		array_free2_d(d,NFluxData->nSolNum_fI);
		array_free2_d(d,NFluxData->dnSolNumdWL_fI);
		array_free2_d(d,NFluxData->dnSolNumdWR_fI);
	}

	for (size_t i = 0; i < 2; i++) {
		free(OPSL[i]);
		free(OPSR[i]);
	}
	free(NFluxData);
	free(FDATAL);
	free(FDATAR);
}

static void finalize_Qhat(struct S_VOLUME const *const VOLUME, unsigned int const NvnS, double *const *const Qhat)
{
	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	if (DB.Collocated) {
		double const *const detJV_vI = VOLUME->detJV_vI;
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NvnS; n++) {
				Qhat[dim][var*NvnS+n] /= detJV_vI[n];
			}}
		}
	} else {
		double *Qhat_tmp = malloc(NvnS*Nvar * sizeof *Qhat_tmp); // free
		for (size_t dim = 0; dim < d; dim++) {
			mm_CTN_d(NvnS,Nvar,NvnS,VOLUME->MInv,Qhat[dim],Qhat_tmp);
			for (size_t var = 0; var < Nvar; var++) {
			for (size_t n = 0; n < NvnS; n++) {
				Qhat[dim][var*NvnS+n] = Qhat_tmp[var*NvnS+n];
			}}
		}
		free(Qhat_tmp);
	}
}

static void finalize_Qhat_What(struct S_VOLUME const *const VOLUME, unsigned int const NRows, unsigned int const NCols,
                               bool const variable_eqvar, double *const *const Qhat_What)
{
	/*
	 *	Comments:
	 *		The variable_eqvar flag is provided as QhatV_What is constant for all equations and variables. The same is
	 *		true for Qhat_What for FACEs which are not on the boundary. In the case of boundary FACEs, while the
	 *		contribution may vary depending on the boundary condition employed.
	 */

	unsigned int const d = DB.d;

	unsigned int eqMax, varMax;

	if (variable_eqvar) {
		eqMax  = DB.Neq;
		varMax = DB.Nvar;
	} else {
		eqMax  = 1;
		varMax = 1;
	}

	if (DB.Collocated) {
		double const *const detJV_vI = VOLUME->detJV_vI;
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < eqMax; eq++) {
			for (size_t var = 0; var < varMax; var++) {
				size_t const Indeqvar = (eq*varMax+var)*NRows*NCols;
				for (size_t i = 0; i < NRows; i++) {
				for (size_t j = 0; j < NCols; j++) {
					Qhat_What[dim][Indeqvar+i*NCols+j] /= detJV_vI[i];
				}}
			}}
		}
	} else {
		double *Qhat_tmp = malloc(NRows*NCols * sizeof *Qhat_tmp); // free
		for (size_t dim = 0; dim < d; dim++) {
			for (size_t eq = 0; eq < eqMax; eq++) {
			for (size_t var = 0; var < varMax; var++) {
				size_t const Indeqvar = (eq*varMax+var)*NRows*NCols;
				mm_d(CBRM,CBNT,CBNT,NRows,NCols,NRows,1.0,0.0,VOLUME->MInv,&Qhat_What[dim][Indeqvar],Qhat_tmp);
				for (size_t i = 0; i < NRows; i++) {
				for (size_t j = 0; j < NCols; j++) {
					Qhat_What[dim][Indeqvar+i*NCols+j] = Qhat_tmp[i*NCols+j];
				}}
			}}
		}
		free(Qhat_tmp);
	}
}

static void implicit_GradW_finalize(void)
{
	/*
	 *	Comments:
	 *		The contributions to Qhat_What from VOLUMEs and FACEs are not summed as they are treated in the FACE info
	 *		function.
	 */

	unsigned int const d    = DB.d,
	                   Nvar = DB.Nvar;

	// Add FACE contributions to VOLUME->Qhat then multiply by MInv
	for (struct S_FACE *FACE = DB.FACE; FACE; FACE = FACE->next) {
		struct S_VOLUME const *const VL = FACE->VIn,
		                      *const VR = FACE->VOut;

		unsigned int const NvnSL = VL->NvnS,
		                   NvnSR = VR->NvnS;

		for (size_t dim = 0; dim < d; dim++) {
			for (size_t i = 0; i < NvnSL*Nvar; i++)
				VL->Qhat[dim][i] += FACE->QhatL[dim][i];

			if (!FACE->Boundary) {
				for (size_t i = 0; i < NvnSR*Nvar; i++)
					VR->Qhat[dim][i] += FACE->QhatR[dim][i];
			}
		}

		finalize_Qhat(VL,NvnSL,FACE->QhatL);
		if (FACE->Boundary) {
			finalize_Qhat_What(VL,NvnSL,NvnSL,1,FACE->Qhat_WhatLL);
		} else {
			finalize_Qhat(VR,NvnSR,FACE->QhatR);
			finalize_Qhat_What(VL,NvnSL,NvnSL,0,FACE->Qhat_WhatLL);
			finalize_Qhat_What(VL,NvnSL,NvnSR,0,FACE->Qhat_WhatRL);
			finalize_Qhat_What(VR,NvnSR,NvnSL,0,FACE->Qhat_WhatLR);
			finalize_Qhat_What(VR,NvnSR,NvnSR,0,FACE->Qhat_WhatRR);
		}
	}

	// Multiply VOLUME Qhat terms by MInv
	for (struct S_VOLUME *VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		unsigned int const NvnS = VOLUME->NvnS;

		finalize_Qhat(VOLUME,NvnS,VOLUME->Qhat);
		finalize_Qhat(VOLUME,NvnS,VOLUME->QhatV);

//		finalize_Qhat_What(VOLUME,NvnS,NvnS,0,VOLUME->QhatV_What);
	}
printf("iGf: disabled functionality.\n");
}