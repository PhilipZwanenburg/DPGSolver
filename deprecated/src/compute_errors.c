// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "compute_errors.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "Parameters.h"
#include "Macros.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"
#include "S_OpCSR.h"
#include "Test.h"

#include "element_functions.h"
#include "matrix_functions.h"
#include "variable_functions.h"
#include "exact_solutions.h"
#include "update_VOLUMEs.h"
#include "initialize_test_case.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Compute L2 errors for supported test cases.
 *
 *	Comments:
 *		The OUTPUT_GRADIENTS flag is provided if it is desired to also output gradient errors for the TaylorCouette
 *		case. It was found that gradients converged optimally whenever the solution converged optimally and this
 *		functionality was disabled.
 *
 *	Notation:
 *
 *	References:
 */

#define OUTPUT_GRADIENTS 0

static void output_errors (const double *L2Error, const unsigned int NvarError, const unsigned int DOF, const double Vol);
static void collect_errors (const unsigned int NvarError);
static char *set_fname(const unsigned int collect);

struct S_OPERATORS {
	unsigned int NvnS, NvnI;
	double       *I_vG_vI, *w_vI, *ChiS_vI;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	// Standard datatypes
	unsigned int P, type, curved;
	struct S_ELEMENT *ELEMENT, *ELEMENT_OPS;

	P = VOLUME->P;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);
	if (1 || type == TRI || type == TET || type == PYR)
		ELEMENT_OPS = ELEMENT;
	else if (type == LINE || type == QUAD || type == HEX || type == WEDGE)
		ELEMENT_OPS = ELEMENT->ELEMENTclass[IndClass];

	OPS->NvnS = ELEMENT_OPS->NvnS[P];
	if (!curved) {
		OPS->NvnI = ELEMENT_OPS->NvnIs[P];

		OPS->I_vG_vI = ELEMENT_OPS->I_vGs_vIs[1][P][0];
		OPS->w_vI    = ELEMENT_OPS->w_vIs[P];
		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIs[P][P][0];
	} else {
		OPS->NvnI = ELEMENT_OPS->NvnIc[P];

		OPS->I_vG_vI = ELEMENT_OPS->I_vGc_vIc[P][P][0];
		OPS->w_vI    = ELEMENT_OPS->w_vIc[P];
		OPS->ChiS_vI = ELEMENT_OPS->ChiS_vIc[P][P][0];
	}
}

void compute_errors(struct S_VOLUME *VOLUME, double *L2Error2, double *Vol, unsigned int *DOF,
                    const unsigned int solved)
{
	// Initialize DB Parameters
	char         *TestCase        = DB.TestCase;
	unsigned int d                = DB.d,
	             Nvar             = DB.Nvar,
	             TestL2projection = DB.TestL2projection;

	// Standard datatypes
	unsigned int i, iMax, j, NvnS, NvnI, IndU, dim, Indq;
	double       *XYZ_vI,
	             *rho, *p, *s, *rhoEx, *pEx, *sEx, *U, *UEx, *What, *W, *u, *q, *uEx, *qEx,
	             *detJV_vI, *w_vI, *ChiS_vI, *wdetJV_vI, err,
	             *ChiSwdetJV_vI, *MInvChiSwdetJV_vI, *L2_Ex_vI, *diag_wdetJV_vI;

	struct S_OPERATORS *OPS;

	OPS = malloc(sizeof *OPS); // free

	init_ops(OPS,VOLUME,0);

	NvnS    = OPS->NvnS;
	NvnI    = OPS->NvnI;
	w_vI    = OPS->w_vI;
	ChiS_vI = OPS->ChiS_vI;

//printf("%d %d\n",VOLUME->indexg,VOLUME->curved);
	detJV_vI = VOLUME->detJV_vI;

	wdetJV_vI = malloc(NvnI * sizeof *wdetJV_vI); // free
	for (i = 0; i < NvnI; i++)
		wdetJV_vI[i] = w_vI[i]*detJV_vI[i];

	XYZ_vI = malloc(NvnI*d * sizeof *XYZ_vI); // free
	mm_CTN_d(NvnI,d,VOLUME->NvnG,OPS->I_vG_vI,VOLUME->XYZ,XYZ_vI);

	if (strstr(TestCase,"Advection")) {
		for (i = 0, iMax = 1; i < iMax; i++)
			L2Error2[i] = 0.0;

		u = malloc(NvnI * sizeof *u); // free
		mm_CTN_d(NvnI,1,NvnS,ChiS_vI,VOLUME->What,u);

		uEx = malloc(NvnI * sizeof *uEx); // free
		compute_exact_solution(NvnI,XYZ_vI,uEx,solved);

		for (i = 0, iMax = 1; i < iMax; i++) {
			for (j = 0; j < NvnI; j++) {
				err = u[j]-uEx[j];
				L2Error2[i] += err*err*wdetJV_vI[j];
			}
		}
		free(u);
		free(uEx);
	} else if (strstr(TestCase,"Poisson")) {
		for (i = 0, iMax = DMAX+1; i < iMax; i++)
			L2Error2[i] = 0.0;

		u = malloc(NvnI      * sizeof *u); // free
		q = calloc(NvnI*DMAX , sizeof *q); // free
		mm_CTN_d(NvnI,1,NvnS,ChiS_vI,VOLUME->What,u);
		for (dim = 0; dim < d; dim++)
			mm_CTN_d(NvnI,1,NvnS,ChiS_vI,VOLUME->Qhat[dim],&q[dim*NvnI]);

		uEx = malloc(NvnI      * sizeof *uEx); // free
		qEx = calloc(NvnI*DMAX , sizeof *qEx); // free

		compute_exact_solution(NvnI,XYZ_vI,uEx,solved);
		compute_exact_gradient(NvnI,XYZ_vI,qEx);

		for (i = 0, iMax = DMAX+1; i < iMax; i++) {
			if (i == 0) { // u
				for (j = 0; j < NvnI; j++) {
					err = u[j]-uEx[j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			} else { // q
				Indq = (i-1)*NvnI;
				for (j = 0; j < NvnI; j++) {
					err = q[Indq+j]-qEx[Indq+j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			}
		}
		free(u);
		free(q);
		free(uEx);
		free(qEx);
	} else if (strstr(TestCase,"PeriodicVortex") ||
	           strstr(TestCase,"EllipticPipe") ||
	           strstr(TestCase,"ParabolicPipe") ||
	           strstr(TestCase,"SinusoidalPipe") ||
	           strstr(TestCase,"SupersonicVortex")) {
		for (i = 0, iMax = NVAR3D+1; i < iMax; i++)
			L2Error2[i] = 0.0;

		UEx = malloc(NvnI*NVAR3D * sizeof *UEx); // free

		compute_solution(NvnI,XYZ_vI,UEx,solved);

		W = malloc(NvnI*Nvar   * sizeof *W); // free
		U = malloc(NvnI*NVAR3D * sizeof *U); // free
		if (!TestL2projection) {
			What = VOLUME->What;
			mm_CTN_d(NvnI,Nvar,NvnS,ChiS_vI,What,W);

			convert_variables(W,U,d,3,NvnI,1,'c','p');
		} else {
			compute_inverse_mass(VOLUME);

			ChiSwdetJV_vI     = malloc(NvnS*NvnI * sizeof *ChiSwdetJV_vI);     // free
			MInvChiSwdetJV_vI = malloc(NvnS*NvnI * sizeof *MInvChiSwdetJV_vI); // free
			L2_Ex_vI          = malloc(NvnI*NvnI * sizeof *L2_Ex_vI);          // free

			diag_wdetJV_vI = diag_d(wdetJV_vI,NvnI); // free

			mm_d(CBRM,CBT,CBNT,NvnS,NvnI,NvnI,1.0,0.0,ChiS_vI,diag_wdetJV_vI,ChiSwdetJV_vI);
			mm_d(CBRM,CBNT,CBNT,NvnS,NvnI,NvnS,1.0,0.0,VOLUME->MInv,ChiSwdetJV_vI,MInvChiSwdetJV_vI);
			mm_d(CBRM,CBNT,CBNT,NvnI,NvnI,NvnS,1.0,0.0,ChiS_vI,MInvChiSwdetJV_vI,L2_Ex_vI);

			mm_CTN_d(NvnI,NVAR3D,NvnI,L2_Ex_vI,UEx,U);
//			mm_d(CBCM,CBT,CBNT,NvnI,NVAR3D,NvnI,1.0,0.0,L2_Ex_vI,UEx,U); // ToBeDeleted

			free(ChiSwdetJV_vI);
			free(MInvChiSwdetJV_vI);
			free(L2_Ex_vI);
			free(diag_wdetJV_vI);
		}

		rhoEx = &UEx[NvnI*0];
		pEx   = &UEx[NvnI*(NVAR3D-1)];
		rho   = &U[NvnI*0];
		p     = &U[NvnI*(NVAR3D-1)];

		sEx = malloc(NvnI * sizeof *sEx); // free
		s   = malloc(NvnI * sizeof *s);   // free
		for (i = 0; i < NvnI; i++) {
			sEx[i] = pEx[i]/pow(rhoEx[i],GAMMA);
			s[i]   = p[i]/pow(rho[i],GAMMA);
		}

		for (i = 0; i <= NVAR3D; i++) {
			IndU = i*NvnI;
			if (i == 0 || i == NVAR3D-1) { // rho, p
				for (j = 0; j < NvnI; j++) {
					err = (U[IndU+j]-UEx[IndU+j])/UEx[IndU+j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			} else if (i > 0 && i < NVAR3D-1) { // u, v, w (Not normalized as variables may be negative)
				for (j = 0; j < NvnI; j++) {
					err = U[IndU+j]-UEx[IndU+j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			} else if (i == NVAR3D) { // s
				for (j = 0; j < NvnI; j++) {
					err = (s[j]-sEx[j])/sEx[j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			}
		}
		free(W);
		free(U);
		free(s);
		free(UEx);
		free(sEx);
	} else if (strstr(TestCase,"GaussianBump") ||
	           strstr(TestCase,"SubsonicNozzle")) {
		L2Error2[0] = 0.0;

		W = malloc(NvnI*Nvar   * sizeof *W); // free
		U = malloc(NvnI*NVAR3D * sizeof *U); // free

		What = VOLUME->What;
		mm_CTN_d(NvnI,Nvar,NvnS,ChiS_vI,What,W);

		convert_variables(W,U,d,3,NvnI,1,'c','p');

		rho   = &U[NvnI*0];
		p     = &U[NvnI*(NVAR3D-1)];

		sEx = malloc(NvnI * sizeof *sEx); // free
		s   = malloc(NvnI * sizeof *s);   // free
		for (i = 0; i < NvnI; i++) {
			// Should s = log(P/rho^GAMMA)? ToBeModified
			sEx[i] = DB.pInf/pow(DB.rhoInf,GAMMA);
			s[i]   = p[i]/pow(rho[i],GAMMA);
		}

		for (j = 0; j < NvnI; j++) {
			err = (s[j]-sEx[j])/sEx[j];
			L2Error2[0] += err*err*wdetJV_vI[j];
		}

		free(W);
		free(U);
		free(sEx);
		free(s);
	} else if (strstr(TestCase,"TaylorCouette")) {
		if (d != 2)
			EXIT_UNSUPPORTED;

		if (!OUTPUT_GRADIENTS) {
			for (i = 0, iMax = 3; i < iMax; i++)
				L2Error2[i] = 0.0;
		} else {
			for (i = 0, iMax = 3*DMAX; i < iMax; i++)
				L2Error2[i] = 0.0;
		}

		UEx = malloc(NvnI*NVAR3D * sizeof *UEx); // free
		compute_solution(NvnI,XYZ_vI,UEx,solved);

		W = malloc(NvnI*Nvar   * sizeof *W); // free
		U = malloc(NvnI*NVAR2D * sizeof *U); // free

		mm_CTN_d(NvnI,Nvar,NvnS,ChiS_vI,VOLUME->What,W);

		convert_variables(W,U,d,2,NvnI,1,'c','p');

		rhoEx = &UEx[NvnI*0];
		pEx   = &UEx[NvnI*(NVAR3D-1)];
		rho   = &U[NvnI*0];
		p     = &U[NvnI*(NVAR2D-1)];

		double *const TEx = malloc(NvnI * sizeof *TEx), // free
		       *const T   = malloc(NvnI * sizeof *T);   // free

		double const Rg = DB.Rg;
		for (j = 0; j < NvnI; j++) {
			TEx[j] = pEx[j]/(rhoEx[j]*Rg);
			T[j]   = p[j]/(rho[j]*Rg);
		}

		// Solution errors
		for (i = 0; i < 3; i++) {
			IndU = (i+1)*NvnI;
			if (i < 2) { // u, v
				for (j = 0; j < NvnI; j++) {
					err = U[IndU+j]-UEx[IndU+j]; // u, v (Not normalized as variables may be negative)
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			} else if (i == 2) { // T
				for (j = 0; j < NvnI; j++) {
					err = (T[j]-TEx[j])/TEx[j];
					L2Error2[i] += err*err*wdetJV_vI[j];
				}
			}
		}

		if (OUTPUT_GRADIENTS) {
			// Compute exact gradients for the u, v velocity components and temperature
			double *const QuvTEx = malloc(NvnI*3*d * sizeof *QuvTEx); // free
			compute_exact_gradient(NvnI,XYZ_vI,QuvTEx);

			// Compute approximate gradients
			double *const Q    = malloc(NvnI*Nvar * sizeof *Q),    // free
			       *const QuvT = malloc(NvnI*3*d  * sizeof *QuvT); // free

			for (size_t dim = 0; dim < d; dim++) {
				mm_CTN_d(NvnI,Nvar,NvnS,ChiS_vI,VOLUME->Qhat[dim],Q);
				size_t const IndQ = dim*3*NvnI;
				for (size_t n = 0; n < NvnI; n++) {
					double const drho  = Q[0*NvnI+n],
					             drhou = Q[1*NvnI+n],
					             drhov = Q[2*NvnI+n],
					             dE    = Q[(d+1)*NvnI+n];

					double const rho_inv = 1.0/rho[n],
					             u       = U[1*NvnI+n],
					             v       = U[2*NvnI+n],
					             E       = W[(d+1)*NvnI+n],
					             du      = rho_inv*(drhou-drho*u),
					             dv      = rho_inv*(drhov-drho*v);

					double const dEoRho  = rho_inv*rho_inv*(dE*rho[n]-E*drho),
					             dV2     = 2.0*(u*du+v*dv),
					             dT      = GM1/Rg*(dEoRho-0.5*dV2);

					QuvT[IndQ+0*NvnI+n] = du;
					QuvT[IndQ+1*NvnI+n] = dv;
					QuvT[IndQ+2*NvnI+n] = dT;
				}
			}

			// Gradient errors
			for (size_t i = 0; i < 2*3; i++) {
				size_t const IndL2 = i+3;
				for (size_t j = 0; j < NvnI; j++) {
					err = QuvT[i*NvnI+j]-QuvTEx[i*NvnI+j];
					L2Error2[IndL2] += err*err*wdetJV_vI[j];
				}
			}

			free(Q);
			free(QuvT);
			free(QuvTEx);
		}
		free(W);
		free(U);
		free(UEx);
		free(TEx);
		free(T);
	} else {
		EXIT_UNSUPPORTED;
	}

	*DOF = NvnS;
	*Vol = 0.0;
	for (i = 0; i < NvnI; i++)
		*Vol += wdetJV_vI[i];

	free(wdetJV_vI);
	free(XYZ_vI);

	free(OPS);
}

void compute_errors_global(void)
{
	/*
	 *	Purpose:
	 *		Compute and output global L2 solution errors.
	 *
	 *	Comments:
	 *		The VOLUME loop is placed outside of the compute_errors function as the same function is used on a VOLUME
	 *		basis in adapt_initial.
	 */

	// Initialize DB Parameters
	char *TestCase = DB.TestCase,
	     *Geometry = DB.Geometry;
	int  MPIrank   = DB.MPIrank;

	// Standard datatypes
	unsigned int i, DOF, DOF_l, NvarError, CurvedOnly;
	double       Vol, Vol_l, *L2Error2, *L2Error2_l;

	struct S_VOLUME *VOLUME;

	// silence
	DOF = 0; Vol = 0.0;

	CurvedOnly = 0;
	if (strstr(TestCase,"Advection")) {
		NvarError = 1;
	} else if (strstr(TestCase,"Poisson")) {
		NvarError = DMAX+1;
	} else if (strstr(TestCase,"PeriodicVortex") ||
	           strstr(TestCase,"SupersonicVortex") ||
	           strstr(TestCase,"ParabolicPipe") ||
	           strstr(TestCase,"EllipticPipe")) {
		NvarError = NVAR3D+1;
	} else if (strstr(TestCase,"GaussianBump") ||
	           strstr(TestCase,"SubsonicNozzle")) {
		if (strstr(Geometry,"NacaSymmetric"))
			CurvedOnly = 1; // Avoid trailing edge singularity when computing entropy error.
		NvarError = 1;
	} else if (strstr(TestCase,"TaylorCouette")) {
		if (!OUTPUT_GRADIENTS)
			NvarError = 3; // u, v, T
		else
			NvarError = 3*DMAX; // u, v, T and Gradients in x, y directions
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}

	L2Error2   = calloc(NvarError , sizeof *L2Error2);   // free
	L2Error2_l = malloc(NvarError * sizeof *L2Error2_l); // free

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		if (CurvedOnly && !VOLUME->curved)
			continue;
		compute_errors(VOLUME,L2Error2_l,&Vol_l,&DOF_l,1);

		Vol += Vol_l;
		DOF += DOF_l;

		for (i = 0; i < NvarError; i++)
			L2Error2[i] += L2Error2_l[i];
	}
	free(L2Error2_l);

	// Write to files and collect
	output_errors(L2Error2,NvarError,DOF,Vol);
	free(L2Error2);

	MPI_Barrier(MPI_COMM_WORLD);
	if (!MPIrank)
		collect_errors(NvarError);
}

static void output_errors(const double *L2Error2, const unsigned int NvarError, const unsigned int DOF, const double Vol)
{
	/*
	 * Comments:
	 *		Squared L2 errors are output here.
	 */

	// Initialize DB Parameters
	char *TestCase = DB.TestCase;

	// standard datatypes
	unsigned int i;
	char         *f_name;

	FILE *fID;

	f_name = set_fname(0); // free
	if ((fID = fopen(f_name,"w")) == NULL)
		printf("Error: File: %s, did not open.\n",f_name), EXIT_MSG;

	if (strstr(TestCase,"Advection")) {
		fprintf(fID,"DOF         Vol         L2u2\n");
	} else if (strstr(TestCase,"Poisson")) {
		fprintf(fID,"DOF         Vol         L2u2        L2q1_2      L2q2_2      L2q3_2\n");
	} else if (strstr(TestCase,"PeriodicVortex") ||
	           strstr(TestCase,"SupersonicVortex") ||
	           strstr(TestCase,"EllipticPipe") ||
	           strstr(TestCase,"ParabolicPipe")) {
		fprintf(fID,"DOF         Vol         L2rho2      L2u2        L2v2        L2w2        L2p2        L2s2\n");
	} else if (strstr(TestCase,"GaussianBump") ||
	           strstr(TestCase,"SubsonicNozzle")) {
		fprintf(fID,"DOF         Vol         L2s2\n");
	} else if (strstr(TestCase,"TaylorCouette")) {
		if (!OUTPUT_GRADIENTS)
			fprintf(fID,"DOF         Vol         L2u2        L2v2        L2T2\n");
		else
			fprintf(fID,"DOF         Vol         L2u2        L2v2        L2T2        L2Gradu2                L2Gradv2                L2GradT2 \n");
	} else {
		EXIT_UNSUPPORTED;
	}
	fprintf(fID,"%-10d  %.4e  ",DOF,Vol);
	for (i = 0; i < NvarError; i++)
		fprintf(fID,"%.4e  ",L2Error2[i]);

	fclose(fID);
	free(f_name);
}

static void collect_errors(const unsigned int NvarError)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase;
	int  MPIsize   = DB.MPIsize;

	// Standard datatypes
	char         *f_name, StringRead[STRLEN_MAX], *data;
	int          rank, offset;
	unsigned int i, DOF;
	double       tmp_d, *L2Error2, Vol, *L2Error;

	FILE *fID;

	L2Error2 = calloc(NvarError , sizeof *L2Error2); // free
	Vol = 0;
	for (rank = 0; rank < MPIsize; rank++) {
		f_name = set_fname(0); // free
		if ((fID = fopen(f_name,"r")) == NULL)
			printf("Error: File: %s, did not open.\n",f_name), EXIT_MSG;

		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			i = 0;
			data = StringRead;
			if (sscanf(data," %d%n",&DOF,&offset) == 1)
				data += offset;
			if (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
				Vol += tmp_d;
				data += offset;
			}
			while (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
				L2Error2[i++] += tmp_d;
				data += offset;
			}
		}

		fclose(fID);
		free(f_name);
	}

	L2Error = malloc(NvarError * sizeof *L2Error); // free
	for (i = 0; i < NvarError; i++)
		L2Error[i] = sqrt(L2Error2[i]/Vol);
	free(L2Error2);

	f_name = set_fname(1); // free
	if ((fID = fopen(f_name,"w")) == NULL)
		printf("Error: File: %s, did not open.\n",f_name), EXIT_MSG;

	if (strstr(TestCase,"Advection")) {
		fprintf(fID,"DOF         L2u\n");
	} else if (strstr(TestCase,"Poisson")) {
		fprintf(fID,"DOF         L2u         L2q1        L2q2        L2q3\n");
	} else if (strstr(TestCase,"PeriodicVortex") ||
	           strstr(TestCase,"SupersonicVortex") ||
	           strstr(TestCase,"EllipticPipe") ||
	           strstr(TestCase,"ParabolicPipe")) {
		fprintf(fID,"DOF         L2rho       L2u         L2v         L2w         L2p         L2s\n");
	} else if (strstr(TestCase,"GaussianBump") ||
	           strstr(TestCase,"SubsonicNozzle")) {
		fprintf(fID,"DOF         L2s\n");
	} else if (strstr(TestCase,"TaylorCouette")) {
		if (!OUTPUT_GRADIENTS) {
			fprintf(fID,"DOF         L2u         L2v         L2T\n");
		} else {
			fprintf(fID,"DOF         L2u         L2v         L2T         L2ux        L2uy        L2vx        L2vy        L2Tx        L2Ty        \n");
		}
	} else {
		EXIT_UNSUPPORTED;
	}
	fprintf(fID,"%-10d  ",DOF);
	for (i = 0; i < NvarError; i++)
		fprintf(fID,"%.4e  ",L2Error[i]);
	free(L2Error);

	fclose(fID);
	free(f_name);
}

static char *set_fname(const unsigned int collect)
{
	// Initialize DB Parameters
	char *TestCase = DB.TestCase,
	     *MeshType = DB.MeshType;
	int  MPIrank = DB.MPIrank;

	// standard datatypes
	char *f_name, string[STRLEN_MAX];

	f_name = malloc(STRLEN_MAX * sizeof *f_name); // keep (requires external free)

	if (!collect)
		strcpy(f_name,"errors/");
	else
		strcpy(f_name,"results/");

	strcat(f_name,TestCase); strcat(f_name,"/");
	strcat(f_name,MeshType); strcat(f_name,"/");
	if (!collect)
		strcat(f_name,"L2errors2_");
	else
		strcat(f_name,"L2errors_");
	sprintf(string,"%dD_",DB.d);   strcat(f_name,string);
	                               strcat(f_name,MeshType);
	if (DB.Adapt == ADAPT_0) {
		sprintf(string,"_ML%d",DB.ML); strcat(f_name,string);
		sprintf(string,"P%d",DB.PGlobal), strcat(f_name,string);
	} else {
		sprintf(string,"_ML%d",TestDB.ML); strcat(f_name,string);
		sprintf(string,"P%d",TestDB.PGlobal), strcat(f_name,string);
	}
	if (!collect)
		sprintf(string,"_%d",MPIrank), strcat(f_name,string);
	strcat(f_name,".txt");

	return f_name;
}
