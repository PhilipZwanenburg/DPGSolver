// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_integration.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

#include "petscsys.h"
#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"
#include "Test.h"
#include "S_DB.h"
#include "S_ELEMENT.h"
#include "S_VOLUME.h"

#include "test_support.h"
#include "solver.h"
#include "initialization.h"
#include "setup_parameters.h"
#include "setup_mesh.h"
#include "setup_operators.h"
#include "setup_operators_HDG.h"
#include "setup_structures.h"
#include "setup_geometry.h"
#include "initialize_test_case.h"
#include "adaptation.h"
#include "setup_Curved.h"
#include "element_functions.h"
#include "memory_free.h"
#include "array_norm.h"
#include "array_free.h"

#include "array_print.h"

/*
 *	Purpose:
 *		Provide functions for integration testing.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

static void update_TestCase(void)
{
	/*
	 *	Comments:
	 *		Appending '_Test' for certain integration test cases is necessary to allow for the mesh refinement in
	 *		adapt_initial to proceed for all VOLUMEs.
	 */

	if (strstr(DB.TestCase,"Advection")) {
		if (strstr(DB.TestCase,"Default")) {
			strcpy(DB.TestCase,"Advection_Default");
		} else if (strstr(DB.TestCase,"Peterson")) {
			strcpy(DB.TestCase,"Advection_Peterson");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else if (strstr(DB.TestCase,"Poisson_Ringleb") ||
	           strstr(DB.TestCase,"Poisson_n-Cube") ||
	           strstr(DB.TestCase,"Poisson_n-Ball") ||
	           strstr(DB.TestCase,"Poisson_n-Ellipsoid") ||
	           strstr(DB.TestCase,"Poisson_HoldenRamp") ||
	           strstr(DB.TestCase,"Poisson_GaussianBump")) {
		strcpy(DB.TestCase,"Poisson");
	} else if (strstr(DB.TestCase,"update_h")) {
		strcpy(DB.TestCase,"Poisson_Test");
	} else if (strstr(DB.TestCase,"L2_proj")) {
		strcpy(DB.TestCase,"Euler_PeriodicVortex_Test");
	} else if (strstr(DB.TestCase,"linearization")) {
		strcpy(DB.TestCase,"Euler_SupersonicVortex_Test");
	} else if (strstr(DB.TestCase,"Euler")) {
		if (strstr(DB.TestCase,"SupersonicVortex")) {
			strcpy(DB.TestCase,"Euler_SupersonicVortex");
		} else if (strstr(DB.TestCase,"PeriodicVortex")) {
			if (strstr(DB.TestCase,"Stationary"))
				strcpy(DB.TestCase,"Euler_PeriodicVortex_Stationary");
			else
				strcpy(DB.TestCase,"Euler_PeriodicVortex");
		} else if (strstr(DB.TestCase,"EllipticPipe")) {
			strcpy(DB.TestCase,"Euler_EllipticPipe");
		} else if (strstr(DB.TestCase,"ParabolicPipe")) {
			strcpy(DB.TestCase,"Euler_ParabolicPipe");
		} else if (strstr(DB.TestCase,"SinusoidalPipe")) {
			strcpy(DB.TestCase,"Euler_SinusoidalPipe");
		} else if (strstr(DB.TestCase,"GaussianBump")) {
			strcpy(DB.TestCase,"Euler_InternalSubsonic_GaussianBump");
		} else {
			EXIT_UNSUPPORTED;
		}
/*		if (strstr(DB.Geometry,"Ellipsoidal_Section") ||
		    strstr(DB.Geometry,"GaussianBump") ||
		    strstr(DB.Geometry,"EllipsoidalBump"))
			strcpy(DB.TestCase,"SubsonicNozzle");
		else if (strstr(DB.Geometry,"ExpansionCorner"))
			strcpy(DB.TestCase,"PrandtlMeyer");
*/
	} else if (strstr(DB.TestCase,"NavierStokes")) {
		if (strstr(DB.TestCase,"TaylorCouette")) {
			strcpy(DB.TestCase,"NavierStokes_TaylorCouette");
		} else if (strstr(DB.TestCase,"PlaneCouette")) {
			strcpy(DB.TestCase,"NavierStokes_PlaneCouette");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		printf("%s\n",DB.TestCase);
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

void code_startup
	(const int nargc, const char*const*const argv, const unsigned int n_ref, const unsigned int modify_params)
{
	int  MPIrank, MPIsize;

	// Start MPI and PETSC
	PetscInitialize((int *) &nargc,(char ***) &argv,PETSC_NULL,PETSC_NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
	MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

	// Test memory leaks only from Petsc and MPI using valgrind
	//PetscFinalize(), exit(1);

	DB.MPIsize = MPIsize;
	DB.MPIrank = MPIrank;

	// Initialization
	// - Here, we will first read the control file and set the parameters
	// 		from the control file (the mesh level, P, ...)
	//		then, overwrite certain values using the test parameters (stored
	//		in TestDB). Then, create the name of the mesh file again using this.

	initialization(nargc,argv);
	update_TestCase();
	if (modify_params) {
		DB.PGlobal = TestDB.PGlobal;
		if (modify_params == 1) {
			DB.ML = TestDB.ML;
			set_MeshFile(-1);
		}
	}

	setup_parameters();
	if (modify_params)
		setup_parameters_L2proj();

	// Generate the new mesh file name using the test parameters
	// that were copied into the DB structure and based on the parameters
	// that were set. Also, in the case we are using a Bezier mesh,
	// input the order of the mesh to be used (in PGc)
	if (DB.BezierMesh)
		set_MeshFile(DB.PGc[DB.PGlobal]);

	initialize_test_case_parameters();
	setup_mesh();
	if (strstr(DB.MeshType,"Curved"))
		strcat(DB.Geometry,"Curved");
	setup_operators();
	setup_operators_HDG();
	setup_structures();
	setup_geometry();

	initialize_test_case(n_ref);

	if (modify_params == 2)
		mesh_to_level(TestDB.ML);
}

void code_cleanup(void)
{
	memory_free();
}

void code_startup_mod_prmtrs(int const nargc, char const *const *const argv, unsigned int const Nref,
                             unsigned int const update_argv, unsigned int const phase)
{
	int  MPIrank, MPIsize;

	if (phase == 1) {
		// Start MPI and PETSC
		PetscInitialize((int *) &nargc,(char ***) &argv,PETSC_NULL,PETSC_NULL);
		MPI_Comm_size(MPI_COMM_WORLD,&MPIsize);
		MPI_Comm_rank(MPI_COMM_WORLD,&MPIrank);

		// Test memory leaks only from Petsc and MPI using valgrind
		//PetscFinalize(), exit(1);

		DB.MPIsize = MPIsize;
		DB.MPIrank = MPIrank;

		// Initialization
		initialization(nargc,argv);
		update_TestCase();
		if (update_argv) {
			DB.PGlobal = TestDB.PGlobal;
			if (update_argv == 1)
				DB.ML      = TestDB.ML;
		}

		initialize_test_case_parameters();
		setup_parameters();
		if (update_argv)
			setup_parameters_L2proj();
	} else if (phase == 2) {
		setup_mesh();
		if (strstr(DB.MeshType,"Curved"))
			strcat(DB.Geometry,"Curved");
		setup_operators();
		setup_structures();
		setup_geometry();

		initialize_test_case(Nref);

		if (update_argv == 2)
			mesh_to_level(TestDB.ML);
	} else {
		printf("Error: Unsupported.\n"), EXIT_MSG;
	}
}

void code_startup_mod_ctrl(int const nargc, char const *const *const argv, unsigned int const Nref,
                           unsigned int const update_argv, unsigned int const phase)
{
	if (phase == 1) {
		// Start MPI and PETSC
		PetscInitialize((int *) &nargc,(char ***) &argv,PETSC_NULL,PETSC_NULL);
		MPI_Comm_size(MPI_COMM_WORLD,&DB.MPIsize);
		MPI_Comm_rank(MPI_COMM_WORLD,&DB.MPIrank);

		// Initialization
		initialization(nargc,argv);
	} else if (phase == 2) {
		initialize_test_case_parameters();
		setup_parameters();
		if (update_argv)
			setup_parameters_L2proj();
	} else if (phase == 3) {
		setup_mesh();
		setup_operators();
		setup_structures();
		setup_geometry();

		initialize_test_case(Nref);
	} else {
		EXIT_UNSUPPORTED;
	}
}

void check_convergence_orders(const unsigned int MLMin, const unsigned int MLMax, const unsigned int PMin,
                              const unsigned int PMax, unsigned int *pass, const bool PrintEnabled)
{
	// Initialize DB Parameters
	char         *TestCase = DB.TestCase,
	             *MeshType = DB.MeshType;
	unsigned int d         = DB.d;

	// Standard datatypes
	char         f_name[STRLEN_MAX], string[STRLEN_MAX], StringRead[STRLEN_MAX], *data;
	unsigned int i, ML, P, NVars, NML, NP, Indh, *VarsToCheck;
	int          offset;
	double       **L2Errors, **ConvOrders, *h, tmp_d;

	FILE *fID;

	// silence
	NVars = 0;

	if (strstr(TestCase,"Advection")) {
		NVars = 1;
	} else if (strstr(TestCase,"Poisson")) {
		NVars = DMAX+1;
	} else if (strstr(TestCase,"SupersonicVortex") ||
	           strstr(TestCase,"PeriodicVortex") ||
	           strstr(TestCase,"EllipticPipe") ||
	           strstr(TestCase,"ParabolicPipe")) {
		NVars = DMAX+2+1;
	} else if (strstr(TestCase,"GaussianBump") ||
	           strstr(TestCase,"SubsonicNozzle")) {
		NVars = 1;
	} else if (strstr(TestCase,"TaylorCouette")) {
		NVars = 3;
	} else {
		printf("Error: Unsupported TestCase.\n"), EXIT_MSG;
	}

	VarsToCheck    = malloc(NVars * sizeof *VarsToCheck);    // free
	double *const OrderIncrement = malloc(NVars * sizeof *OrderIncrement); // free
	if (strstr(TestCase,"Advection")) {
		for (size_t i = 0; i < NVars; i++) {
			OrderIncrement[i] = 1;
			VarsToCheck[i]    = 1;
		}
	} else if (strstr(TestCase,"Poisson")) {
		for (i = 0; i < NVars; i++) {
			OrderIncrement[i] = 0;
			if (i == 0) {
				OrderIncrement[i] = 1;
				VarsToCheck[i]    = 1;
			} else if (i <= d) {
				VarsToCheck[i]    = 1;
			} else {
				VarsToCheck[i]    = 0;
			}
		}
	} else if (strstr(TestCase,"SupersonicVortex") ||
	           strstr(TestCase,"PeriodicVortex")  ||
	           strstr(TestCase,"ParabolicPipe")  ||
	           strstr(TestCase,"EllipticPipe")  ||
	           strstr(TestCase,"InviscidChannel")  ||
	           strstr(TestCase,"GaussianBump")  ||
	           strstr(TestCase,"SubsonicNozzle")) {
		for (i = 0; i < NVars; i++) {
			OrderIncrement[i] = 1;
			if (i <= d || i > DMAX) {
				VarsToCheck[i] = 1;
			} else {
				VarsToCheck[i] = 0;
			}
		}
	} else if (strstr(TestCase,"TaylorCouette")) {
		for (i = 0; i < NVars; i++) {
			OrderIncrement[i] = 1;
			VarsToCheck[i] = 1;
		}
	} else {
		EXIT_UNSUPPORTED;
	}


	NML = MLMax-MLMin+1;
	NP  = PMax-PMin+1;

	L2Errors   = malloc(NVars * sizeof *L2Errors);   // free
	ConvOrders = malloc(NVars * sizeof *ConvOrders); // free
	for (i = 0; i < NVars; i++) {
		L2Errors[i]   = calloc(NML*NP , sizeof *L2Errors[i]);   // free
		ConvOrders[i] = calloc(NML*NP , sizeof *ConvOrders[i]); // free
	}
	h = calloc(NML*NP , sizeof *h); // free

	// Read in data and compute convergence orders
	for (ML = MLMin; ML <= MLMax; ML++) {
	for (P = PMin; P <= PMax; P++) {
		Indh = (ML-MLMin)*NP+(P-PMin);

		strcpy(f_name,"../cases/results/");
		strcat(f_name,TestCase); strcat(f_name,"/");
		strcat(f_name,MeshType); strcat(f_name,"/");
		strcat(f_name,"L2errors_");
		sprintf(string,"%dD_",d);   strcat(f_name,string);
		                            strcat(f_name,MeshType);
		sprintf(string,"_ML%d",ML); strcat(f_name,string);
		sprintf(string,"P%d",P);    strcat(f_name,string);
		strcat(f_name,".txt");

		if ((fID = fopen(f_name,"r")) == NULL)
			printf("Error: File: %s, did not open.\n",f_name), EXIT_MSG;

		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) { ; }
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			i = 0;
			data = StringRead;
			if (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
				data += offset;
				h[Indh] = 1.0/pow(tmp_d,1.0/d);
			}
			while (sscanf(data," %lf%n",&tmp_d,&offset) == 1) {
				L2Errors[i++][Indh] = tmp_d;
				data += offset;
			}
		}
		fclose(fID);
	}}

	for (ML = MLMin+1; ML <= MLMax; ML++) {
	for (P = PMin; P <= PMax; P++) {
		Indh = (ML-MLMin)*NP+(P-PMin);
		for (i = 0; i < NVars; i++) {
			if (fabs(L2Errors[i][Indh]) > EPS && fabs(L2Errors[i][Indh-NP]) > EPS)
				ConvOrders[i][Indh] = log10(L2Errors[i][Indh]/L2Errors[i][Indh-NP])/log10(h[Indh]/h[Indh-NP]);
		}
	}}

	*pass = 1;

	bool Peterson_warning = 0;
	for (i = 0; i < NVars; i++) {
		if (!VarsToCheck[i])
			continue;

		for (P = PMin; P <= PMax; P++) {
			Indh = (MLMax-MLMin)*NP+(P-PMin);
			if ((ConvOrders[i][Indh]-(P+OrderIncrement[i])) < -0.125) {
				if (strstr(DB.PDESpecifier,"Peterson") && (ConvOrders[i][Indh]-(P+OrderIncrement[i]-0.5)) > -0.125) {
					Peterson_warning = 1;
					continue;
				}
				*pass = 0;

				printf("i = %d, P = %d, ConvOrder = (% .3e, % .3e)\n",i,P,ConvOrders[i][Indh],1.0*(P+OrderIncrement[i]));
				break;
			}
		}
	}

	if (Peterson_warning)
		test_print_warning("Convergence orders suboptimal by up to half of an order on Peterson meshes");

	free(VarsToCheck);
	free(OrderIncrement);

	if (!(*pass)) {
		test_print_warning("Suboptimal Convergence");
		for (i = 0; i < NVars; i++)
			array_print_d(NML,NP,ConvOrders[i],'R');
		printf("\n\n\n");
	}

	if (PrintEnabled) {
//		printf("ViscousFlux, Blending, Parametrization: %d, %d, %d\n",DB.ViscousFluxType,DB.Blending,DB.Parametrization);
		printf("h:\n");
		array_print_d(NML,NP,h,'R');
		printf("L2Errors: \n");
		for (i = 0; i < NVars; i++)
		array_print_d(NML,NP,L2Errors[i],'R');
		printf("Conv Orders (h): \n");
		for (i = 0; i < NVars; i++)
		array_print_d(NML,NP,ConvOrders[i],'R');

		unsigned int u1 = 1;
		double **L2ErrorsP;

		L2ErrorsP = malloc(NVars * sizeof *L2ErrorsP); // free

		//printf("Conv (p): \n");
		for (i = 0; i < NVars; i++) {
			L2ErrorsP[i] = calloc(NML*NP , sizeof *L2ErrorsP[i]);
			for (ML = MLMin; ML <= MLMax; ML++) {
			for (P = max(PMin,u1); P <= PMax; P++) {
				Indh = (ML-MLMin)*NP+(P-PMin);
				if (L2Errors[i][Indh] > EPS)
					L2ErrorsP[i][Indh] = log(L2Errors[i][Indh])/((double) P);
			}}
		//	array_print_d(NML,NP,L2ErrorsP[i],'R');
		}
		array_free2_d(NVars,L2ErrorsP);
	}

	for (i = 0; i < NVars; i++) {
		free(L2Errors[i]);
		free(ConvOrders[i]);
	}
	free(L2Errors);
	free(ConvOrders);
	free(h);
}

void evaluate_mesh_regularity(double *mesh_quality)
{
	/*
	 *	Purpose:
	 *		Compute the ratio of enclosing to enclosed spheres by each TET ELEMENT and return the maximum.
	 *
	 *	Comments:
	 *		For straight meshes, the mesh regularity check is only necessary for TET refinement as all other ELEMENTs
	 *		refine into pieces having the same shape but with all edge lengths scaled by 1/2 when isotropic h-refinement
	 *		is performed. However, on curved meshes, the movement of the vertices to the exact boundary can cause other
	 *		ELEMENT types to become less regular with mesh refinement.
	 *
	 *		NOTE: IT IS ONLY BE NECESSARY TO COMPUTE THE RATIO OF THE CIRCUMSPHERE TO THE IN-SPHERE FOR THE ESTIMATE.
	 *		THIS IS CONCLUDED BASED ON DISCUSSIONS IN CIARLET(1972)-General_Lagrange_and_Hermite... eq. (2.11) and
	 *		(2.12).
	 *		      -> DELETE CHECKS FOR OTHER TYPES (ToBeModified)
	 *
	 *		TET algorithm:
	 *
	 *		While testing this function, several meshes generated by gmsh returned elements having extremely small
	 *		internal radius (nearly coplanar TET vertices). These elements, being of the worst quality, resulted in a
	 *		growing ratio of the internal to external radii with the mesh refinement despite obtaining optimal
	 *		convergence orders (L2 norm). As the L2 norm is global, it was thus decided to take the average of ratios
	 *		from each element. Correct results were then obtained (Suboptimal using refine by split, optimal using
	 *		non-nested unstructured meshes).
	 *			-> Potentially look into why gmsh gives such poor quality elements. (ToBeModified)
	 *
	 *		The exact radii of the maximal radius enclosed and minimal radius enclosing spheres are calculated as
	 *		follows:
	 *
	 *		Enclosed:
	 *			Define the planes of each face by a*x+b*y+c*z = d:
	 *				1) Compute [a b c] as the normal vector to the plane and find d by substituting one coordinate.
	 *				2) Normalize [a b c] (and d) and invert the normal if not pointing outwards:
	 *					Given distance from point to plane = |a*x+b*y+c*z-d|/norm([a b c],2), ensure that when
	 *					substituting the coordinate of the vertex not in the plane that the result is negative.
	 *				3) Using the distance formula again:
	 *					r = n1*x_c+n2*y_c+n3*z_c-d (4 eqns and 4 unknowns)
	 *					Solve: [ones(d+1,1) nr]*[r x_c y_c z_c]' = d
	 *
	 *		Enclosing:
	 *			It may be sufficient for the mesh quality measure to simply use the circumsphere radius. (ToBeModified)
	 *
	 *			There are three possibilities to consider for this case.
	 *
	 *			1) Two   corners of the TET touch the sphere.
	 *			2) Three corners of the TET touch the sphere.
	 *			3) Four  corners of the TET touch the sphere (circumsphere).
	 *
	 *			The respective radii for the cases above are determined based on the vertices of the longest line, the
	 *			triangle formed from the three longest edges, and the circumsphere.
	 *
	 *			Algorithm:
	 *				Check if all four corners fit within the sphere formed from (in the following order):
	 *				- Sphere  of radius (lmax/2) centered at the midpoint of the longest line;
	 *				- Spheres of radius determined from the triangles formed from each face in asending order of radius;
	 *					In this case, the sphere equation is determined from the three coordinates and from the fact
	 *					that the sphere center must lie in the plane of the three vertices.
	 *
	 *				As soon as a condition is met, the appropriate sphere has been found, if neither of the above
	 *				two conditions are met, the circumsphere provides the correct radius.
	 *
	 *		2D algorithm (TRI/QUADs):
	 *
	 *			Similar to TET but using enclosing and enclosed circles.
	 *
	 */

	// Initialize DB Parameters
	char         *MeshType = DB.MeshType;
	unsigned int d         = DB.d;

	// Standard datatypes
	unsigned int i, j, k, l, f, e, c, Nf, Ne, Nve, iMax, jMax, kMax, lMax, Ecount,
	             *IndsE, *IndsF, Found, TETcount, NormType, *VeEcon, *VeFcon;
	double       r_ratio, r, rIn, rOut, rTmp, d1, d2,
	             *XYZ, *XYZdiff, *n, *nNorm, **LHS, **RHS, *lenE, *XYZc, *abcF, *rF, *XYZcE, *d_p, *XYZ_vV;

	struct S_ELEMENT *ELEMENT;
	struct S_VOLUME  *VOLUME;

	NormType = 0; // Options: 0 (Inf), 2 (L2)

	XYZ     = malloc(NVEMAX*d * sizeof *XYZ);     // free
	XYZdiff = malloc(d        * sizeof *XYZdiff); // free

	n     = malloc(NFMAX*d * sizeof *n);     // free
	nNorm = malloc(NFMAX   * sizeof *nNorm); // free

	LHS = malloc(3 * sizeof *LHS); // free
	RHS = malloc(3 * sizeof *RHS); // free
	lapack_int **piv = malloc(3 * sizeof *piv); // free

	LHS[0] = malloc(DMAX*DMAX             * sizeof *LHS[0]); // free
	LHS[1] = malloc(NFMAX*DMAX*NFMAX*DMAX * sizeof *LHS[1]); // free
	LHS[2] = malloc((DMAX+1)*(DMAX+1)     * sizeof *LHS[2]); // free
	RHS[0] = malloc(DMAX*1                * sizeof *RHS[0]); // free
	RHS[1] = malloc(NFMAX*DMAX*1          * sizeof *RHS[1]); // free
	RHS[2] = malloc((DMAX+1)*1            * sizeof *RHS[2]); // free
	piv[0] = malloc(DMAX*1                * sizeof *piv[0]); // free
	piv[1] = malloc(NFMAX*DMAX*1          * sizeof *piv[1]); // free
	piv[2] = malloc((DMAX+1)*1            * sizeof *piv[2]); // free

	IndsE = calloc(28     , sizeof *IndsE); // free (HEX determines size = d*(7+6+5+...) = d*28
	IndsF = malloc(NFMAX  * sizeof *IndsF); // free
	lenE  = calloc(28     , sizeof *lenE);  // free
	XYZc  = malloc(DMAX   * sizeof *XYZc);  // free
	XYZcE = malloc(DMAX   * sizeof *XYZcE); // free
	d_p   = malloc(NFMAX  * sizeof *d_p);   // free

	abcF = malloc(NFMAX*DMAX * sizeof *abcF); // free
	rF   = malloc(NFMAX*1    * sizeof *rF);   // free

	r_ratio  = 1.0*d; // Minimum value for regular TRI/TET.
	if (d == 3) {
		// Currently only being used for TETs.
		TETcount = 0;
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			if (VOLUME->type != TET)
				continue;

			ELEMENT = get_ELEMENT_type(VOLUME->type);

			Nf  = ELEMENT->Nf;
			Ne  = ELEMENT->Ne;
			Nve = ELEMENT->Nve;

			VeEcon = ELEMENT->VeEcon;
			VeFcon = ELEMENT->VeFcon;

			// Obtain vertex coordinates (Row-major)
			if (strstr(MeshType,"ToBeCurved"))
				XYZ_vV = VOLUME->XYZ_vVc;
			else
				XYZ_vV = VOLUME->XYZ_vV;

			for (i = 0; i < Nve; i++) {
			for (j = 0; j < d; j++) {
				XYZ[i*d+j] = XYZ_vV[j*Nve+i];
			}}

			// XYZ coordinates of TET center
			for (j = 0; j < d; j++) {
				XYZcE[j] = 0.0;
				for (i = 0; i < Nve; i++)
					XYZcE[j] += XYZ[i*d+j];
				XYZcE[j] /= 1.0*Nve;
			}

			// Evaluate radius of enclosed sphere

			// 1) Normal vector (normalized) computation
			for (f = 0; f < Nf; f++) {
				// Cross-product
				compute_plane(&XYZ[VeFcon[f*NFVEMAX+0]*d],&XYZ[VeFcon[f*NFVEMAX+1]*d],&XYZ[VeFcon[f*NFVEMAX+2]*d],
				              RHS[0],&d_p[f]);

				// Normalize
				nNorm[f] = array_norm_d(d,RHS[0],"L2");
				for (i = 0; i < d; i++)
					n[f*d+i] = RHS[0][i]/nNorm[f];
				d_p[f] /= nNorm[f];

				// Ensure that normal points outwards
				r = -d_p[f];
				for (i = 0; i < d; i++)
					r += n[f*d+i]*XYZ[f*d+i];

				if (r > 0.0) {
					for (i = 0; i < d; i++)
						n[f*d+i] *= -1.0;
					d_p[f] *= -1.0;
				}
			}

			// 2) Assemble LHS, RHS and solve system to obtain rIn
			for (i = 0, iMax = d+1; i < iMax; i++) {
				for (j = 0, jMax = d+1; j < jMax; j++) {
					if (j == 0)
						LHS[2][i*jMax+j] = 1.0;
					else
						LHS[2][i*jMax+j] = n[i*d+(j-1)];
				}
				RHS[2][i] = d_p[i];
			}

			if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
				printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

			rIn = RHS[2][0];


			// Evaluate radius of enclosing sphere

			// 1) Compute length of all edges
			for (e = 0; e < Ne; e++) {
				for (i = 0; i < d; i++)
					XYZdiff[i] = XYZ[VeEcon[e*NEVEMAX+0]*d+i]-XYZ[VeEcon[e*NEVEMAX+1]*d+i];
				lenE[e]  = array_norm_d(d,XYZdiff,"L2");
				IndsE[e] = e;
			}
			PetscSortRealWithPermutation(Ne,lenE,(int *) IndsE);

			// Sphere radius for case 1.
			rOut = 0.5*lenE[IndsE[Ne-1]];

			Found = 0;
			for (k = 0; k < d && !Found; k++) {
				if (k == 0) {
					// Case 1: 2 Corners touch the sphere

					// Compute sphere center
					for (j = 0; j < d; j++)
						XYZc[j] = 0.5*(XYZ[VeEcon[IndsE[Ne-1]*NEVEMAX+0]*d+j]+XYZ[VeEcon[IndsE[Ne-1]*NEVEMAX+1]*d+j]);

					// Check if all corners are enclosed by the sphere
					Found = 1;
					for (c = 0; c < Nve; c++) {
						for (i = 0; i < d; i++)
							XYZdiff[i] = XYZ[c*d+i]-XYZc[i];
						r = array_norm_d(d,XYZdiff,"L2");
						if (r > rOut+EPS)
							Found = 0;
					}
//					printf("rOut (%d) % .3e\n",k,rOut);
					if (Found)
						break;
				} else if (k == 1) {
					// Case 2: 3 Corners touch the sphere

					// Find the spheres generated by each face and check in order of ascending radius
					jMax = d+1;
					for (f = 0; f < Nf; f++) {
						// Equation of the plane for each face already determined above (coefs == n*nNorm, d == d_p)

						// Solve for sphere center and radius
						for (i = 0; i < d; i++) {
							for (j = 0; j < d; j++)
								LHS[2][i*jMax+j] = 2.0*XYZ[VeFcon[f*NFVEMAX+i]*d+j];
							LHS[2][i*jMax+d] = 1.0;
							RHS[2][i] = pow(array_norm_d(d,&XYZ[VeFcon[f*NFVEMAX+i]*d],"L2"),2.0);
						}
						i = d;
						for (j = 0; j < d; j++)
							LHS[2][i*jMax+j] = n[f*d+j];
						LHS[2][i*jMax+d] = 0.0;
						RHS[2][i] = d_p[f];

						if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
							printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

						for (i = 0; i < d; i++)
							abcF[f*d+i] = RHS[2][i];
						rF[f] = sqrt(RHS[2][3]+pow(array_norm_d(d,RHS[2],"L2"),2.0));

						IndsF[f] = f;
					}
					PetscSortRealWithPermutation(Nf,rF,(int *) IndsF);

					for (f = 0; f < Nf; f++) {
						for (i = 0; i < d; i++)
							XYZc[i] = abcF[IndsF[f]*d+i];
						rOut = rF[IndsF[f]];

						// Check if all corners are enclosed by the sphere
						Found = 1;
						for (c = 0; c < Nve; c++) {
							for (i = 0; i < d; i++)
								XYZdiff[i] = XYZ[c*d+i]-XYZc[i];
							r = array_norm_d(d,XYZdiff,"L2");
							if (r-rOut > 1e2*EPS)
								Found = 0;
						}
//						printf("rOut (%d, %d) % .3e\n",k,f,rOut);
						if (Found)
							break;
					}
				} else if (k == 2) {
					// Case 3: 4 Corners touch the sphere

					// Compute sphere radius
					for (i = 0; i < d+1; i++) {
						for (j = 0; j < d; j++)
							LHS[2][i*jMax+j] = 2.0*XYZ[i*d+j];
						LHS[2][i*jMax+d] = 1.0;
						RHS[2][i] = pow(array_norm_d(d,&XYZ[i*d],"L2"),2.0);
					}

					if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
						printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

					rOut = sqrt(RHS[2][3]+pow(array_norm_d(d,RHS[2],"L2"),2.0));
//					printf("rOut (%d) % .3e\n",k,rOut);
				}
			}

			if (NormType == 0) {
				if (rOut/rIn > r_ratio) {
					r_ratio = rOut/rIn;
				}
			} else if (NormType == 2) {
				r_ratio += rOut/rIn;
			}
			TETcount += 1;

//printf("%d % .3e % .3e % .3e % .3e % .3e\n",
//       k,rOut/rIn,rOut,rIn,sqrt(0.375)*lenE[IndsE[Ne-1]],rOut-sqrt(0.375)*lenE[IndsE[Ne-1]]);

			if (rOut - sqrt(0.375)*lenE[IndsE[Ne-1]] > 1e2*EPS)
				printf("Error: rOut larger than upper bound on sphere radius.\n"), EXIT_MSG;
		}
		if (TETcount) {
			if (NormType == 0)
				*mesh_quality = r_ratio;
			else if (NormType == 2)
				*mesh_quality = r_ratio/TETcount;
			else
				printf("Error: Unsupported.\n"), EXIT_MSG;
		} else {
			printf("Update mesh regularity check for element types other than TET.\n");
			*mesh_quality = 3.0; // ratio for regular TET
		}
	} else if (d == 2) {
		unsigned int IndcF[4*3]   = { 0, 1, 2,  0, 1, 3,  0, 2, 3,  1, 2, 3 },   // (Ind)ices of (c)ircle (F)aces
		             IndcVe2[6*2] = { 0, 1,  0, 2,  1, 2,  0, 3,  1, 3,  2, 3 }, // (Ind)ices of (c)ircle (Ve)rtices
		             IndcVe3[4*3] = { 0, 1, 2,  0, 1, 3,  0, 2, 3,  1, 2, 3 };   // (Ind)ices of (c)ircle (Ve)rtices

		Ecount = 0;
		for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
			ELEMENT = get_ELEMENT_type(VOLUME->type);

			Nf  = ELEMENT->Nf;
			Ne  = ELEMENT->Ne;
			Nve = ELEMENT->Nve;

			VeEcon = ELEMENT->VeEcon;
			VeFcon = ELEMENT->VeFcon;

			// Obtain vertex coordinates (Row-major)
			if (strstr(MeshType,"ToBeCurved"))
				XYZ_vV = VOLUME->XYZ_vVc;
			else
				XYZ_vV = VOLUME->XYZ_vV;

			for (i = 0; i < Nve; i++) {
			for (j = 0; j < d; j++) {
				XYZ[i*d+j] = XYZ_vV[i+Nve*j];
			}}

			// XYZ coordinates of ELEMENT center
			for (j = 0; j < d; j++) {
				XYZcE[j] = 0.0;
				for (i = 0; i < Nve; i++)
					XYZcE[j] += XYZ[i*d+j];
				XYZcE[j] /= 1.0*Nve;
			}

			// Evaluate radius of enclosed circle

			// 1) Normal vector (normalized) computation
			for (f = 0; f < Nf; f++) {
				n[f*d+0] =   XYZ[VeFcon[f*NFVEMAX+0]*d+1]-XYZ[VeFcon[f*NFVEMAX+1]*d+1];
				n[f*d+1] = -(XYZ[VeFcon[f*NFVEMAX+0]*d+0]-XYZ[VeFcon[f*NFVEMAX+1]*d+0]);

				nNorm[f] = array_norm_d(d,&n[f*d],"L2");

				// Normalize
				for (i = 0; i < d; i++)
					n[f*d+i] /= nNorm[f];

				// Ensure that normal points outwards
				d1 = 0.0;
				for (i = 0; i < d; i++)
					d1 += pow(XYZ[VeFcon[f*NFVEMAX+0]*d+i]-n[f*d+i]*nNorm[f]-XYZcE[i],2.0);
				d1 = sqrt(d1);

				d2 = 0.0;
				for (i = 0; i < d; i++)
					d2 += pow(XYZ[VeFcon[f*NFVEMAX+0]*d+i]+n[f*d+i]*nNorm[f]-XYZcE[i],2.0);
				d2 = sqrt(d2);

				if (d1 > d2) {
					for (i = 0; i < d; i++)
						n[f*d+i] *= -1.0;
				}

				// Compute d for the equation of the line: a*x+b*y = d
				d_p[f] = 0.0;
				for (i = 0; i < d; i++)
					d_p[f] += n[f*d+i]*XYZ[VeFcon[f*NFVEMAX+0]*d+i];
			}

			// 2) Assemble LHS, RHS and solve system to obtain rIn
			kMax = 0;
			if (VOLUME->type == TRI)
				kMax = 1;
			else if (VOLUME->type == QUAD)
				kMax = 4;

			rIn = 0.0;
			for (k = 0; k < kMax; k++) {
				for (i = 0, iMax = d+1; i < iMax; i++) {
					for (j = 0, jMax = d+1; j < jMax; j++) {
						if (j == 0)
							LHS[2][i*jMax+j] = 1.0;
						else
							LHS[2][i*jMax+j] = n[IndcF[k*(d+1)+i]*d+(j-1)];
					}
					RHS[2][i] = d_p[IndcF[k*(d+1)+i]];
				}

				if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
					printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

				if (k == 0 || (RHS[2][0] < rIn))
					rIn = RHS[2][0];
			}

			// Evaluate radius of enclosing sphere

			// 1) Compute length of all possible 2 vertex circles
			kMax = 0;
			if (VOLUME->type == TRI)
				kMax = 3;
			else if (VOLUME->type == QUAD)
				kMax = 6;

			for (k = 0; k < kMax; k++) {
				for (i = 0; i < d; i++)
					XYZdiff[i] = XYZ[IndcVe2[k*2+0]*d+i]-XYZ[IndcVe2[k*2+1]*d+i];
				lenE[k]  = array_norm_d(d,XYZdiff,"L2");
				IndsE[k] = k;
			}
			PetscSortRealWithPermutation(kMax,lenE,(int *) IndsE);

			// Circle radius for case 1.
			rOut = 0.5*lenE[IndsE[kMax-1]];

			Found = 0;
			for (k = 0; k < d && !Found; k++) {
				if (k == 0) {
					// Case 1: 2 Corners touch the circle

					// Compute circle center
					for (j = 0; j < d; j++)
						XYZc[j] = 0.5*(XYZ[IndcVe2[IndsE[Ne-1]*2+0]*d+j]+XYZ[IndcVe2[IndsE[Ne-1]*2+1]*d+j]);

					// Check if all corners are enclosed by the circle
					Found = 1;
					for (c = 0; c < Nve; c++) {
						for (i = 0; i < d; i++)
							XYZdiff[i] = XYZ[c*d+i]-XYZc[i];
						r = array_norm_d(d,XYZdiff,"L2");
						if (r > rOut+EPS) {
							Found = 0;
							break;
						}
					}
//					printf("rOut (%d) % .3e\n",k,rOut);
					if (Found)
						break;
				} else if (k == 1) {
					// Case 2: 3 Corners touch the circle
					lMax = 0;
					if (VOLUME->type == TRI)
						lMax = 1;
					else if (VOLUME->type == QUAD)
						lMax = 4;

					rOut = 0.0;
					for (l = 0; l < lMax; l++) {
						// Compute circle radius
						for (i = 0; i < d+1; i++) {
							for (j = 0; j < d; j++)
								LHS[2][i*jMax+j] = 2.0*XYZ[IndcVe3[l*3+i]*d+j];
							LHS[2][i*jMax+d] = 1.0;
							RHS[2][i] = pow(array_norm_d(d,&XYZ[IndcVe3[l*3+i]*d],"L2"),2.0);
						}

						if (LAPACKE_dgesv(LAPACK_ROW_MAJOR,d+1,1,LHS[2],d+1,piv[2],RHS[2],1) > 0)
							printf("Error: mkl LAPACKE_dgesv failed.\n"), EXIT_MSG;

						rTmp = sqrt(RHS[2][d]+pow(array_norm_d(d,RHS[2],"L2"),2.0));
						if (rTmp > rOut)
							rOut = rTmp;
					}
//					printf("rOut (%d) % .3e\n",k,rOut);
				}
			}

			if (rOut/rIn < 1.0)
				printf("Error: Invalid value (% .3e, % .3e % .3e).\n",rIn,rOut,rOut/rIn), EXIT_MSG;

			if (NormType == 0) {
				if (rOut/rIn > r_ratio) {
					r_ratio = rOut/rIn;
				}
			} else if (NormType == 2) {
				r_ratio += rOut/rIn;
			}
			Ecount++;
		}
		if (NormType == 0)
			*mesh_quality = r_ratio;
		else if (NormType == 2)
			*mesh_quality = r_ratio/Ecount;
		else
			printf("Error: Unsupported.\n"), EXIT_MSG;
	} else if (d == 1) {
		*mesh_quality = 1.0;
	}

//	printf("%d % .3e\n",DB.ML,*mesh_quality);

	free(XYZ);
	free(XYZdiff);
	free(n);
	free(nNorm);

	array_free2_d(3,LHS);
	array_free2_d(3,RHS);
	array_free2_i(3,piv);

	free(IndsE);
	free(IndsF);
	free(lenE);
	free(XYZc);
	free(XYZcE);
	free(d_p);

	free(abcF);
	free(rF);
}

void check_mesh_regularity(const double *mesh_quality, const unsigned int NML, unsigned int *pass,
                           const bool PrintEnabled)
{
	/*
	 *	Purpose:
	 *		Checks that the sequence of refined meshes used satisfy the mesh regularity condition required for optimal
	 *		convergence to be obtained.
	 *
	 *	Comments:
	 *		For the moment, the slopes of the last two intervals are used to assess this.
	 */

	unsigned int i;
	double       slope_quality[2];

	if (NML > 2) {
		for (i = 0; i < 2; i++)
			slope_quality[i] = mesh_quality[NML-2+i]-mesh_quality[NML-3+i];

		if (slope_quality[1]-slope_quality[0] > 1e-3 && PrintEnabled) {
			test_print_warning("Potential mesh regularity issue");
		}

		if (slope_quality[1] > 0.0) {
			if (slope_quality[1]/slope_quality[0] > 5e0)
				*pass = 0;
		}

		if (PrintEnabled) {
			printf("\nMesh quality:");
			array_print_d(1,NML,mesh_quality,'R');
		}
	}
}

void set_PrintName(char *name_type, char *PrintName, bool *omit_root)
{
	if (!(*omit_root)) {
		*omit_root = 1;
		if (strstr(name_type,"conv_orders")) {
			strcpy(PrintName,"Convergence Orders         (");
		} else if (strstr(name_type,"equiv_rc")) {
			strcpy(PrintName,"Equivalence Real/Complex   (");
		} else if (strstr(name_type,"equiv_alg")) {
			strcpy(PrintName,"Equivalence Algorithms     (");
		} else if (strstr(name_type,"linearization")) {
			if (strstr(name_type,"weak gradient"))
				strcpy(PrintName,"Linearization (weak grad.) (");
			else
				strcpy(PrintName,"Linearization              (");
		} else {
			EXIT_UNSUPPORTED;
		}
	} else {
		strcpy(PrintName,"                           (");
	}

	char *method = malloc(STRLEN_MIN * sizeof *method); // free
	sprintf(method,"Method: %d",DB.Method);
	strcat(PrintName,method); strcat(PrintName,", ");
	free(method);

	strcat(PrintName,DB.PDE); strcat(PrintName,", ");
	if (!strstr(DB.PDESpecifier,"NONE")) {
		strcat(PrintName,DB.PDESpecifier); strcat(PrintName,", ");
	}
	strcat(PrintName,DB.MeshType); strcat(PrintName,") : ");
}

void compute_dof (void)
{
	struct S_solver_info solver_info = constructor_solver_info(false,false,false,'I',DB.Method);
	set_global_indices(&solver_info);

	DB.dof = solver_info.dof;
}
