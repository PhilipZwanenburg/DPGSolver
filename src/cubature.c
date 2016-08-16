// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include "cubature.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "petscsys.h"
#include "mkl.h"
 
#include "Parameters.h"
#include "Macros.h"

#include "array_norm.h"
#include "array_sort.h"
#include "matrix_functions.h"

void cubature_TP(double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                 const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType)
{
	/*
	 *	Purpose:
	 *		Return nodes, weights and symmetries for tensor-product cubature depending on the nodetype.
	 *
	 *	Comments:
	 *		Note that GL nodes are NOT the optimal integration nodes for d > 1. This was discussed in Hughes' FEM book
	 *		and the pyFR code has some nodal sets available for low orders. Possibly implement this later (ToBeDeleted).
	 *		The order of rst, w is important for minimizing memory stride while computing the length 2 Discrete Fourier
	 *		Transform (ToBeModified).
	 *		Ordering convention:
	 *			GL/GLL: 2-blocks of -/+ node with abs(node) going from 1 to 0 followed by optional 0 node.
	 *		rst is stored in memory as r, s, then t. This was motivated by the desire to be consistent with how XYZ
	 *		coordinates are stored in each element, with this ordering being required so that vectorized operations can
	 *		be performed on XYZ arrays.
	 *
	 *	Notation:
	 *		rst   : Nodes array of dimension Nn*d (column-major storage)
	 *		w     : Weights array of dimension Nn*1
	 *		symms : Symmetries array
	 *		        This always returns 1d symmetries in the current implementation (ToBeModified)
	 *		Nn    : (N)umber of (n)odes
	 *		Ns    : (N)umber of (s)ymmetries
	 *
	 *	References:
	 *		GL  : http://www.mathworks.com/matlabcentral/fileexchange/26737-legendre-laguerre-and-hermite-gauss-quadrature/
	 *		      content/GaussLegendre.m
	 *		GLL : http://www.mathworks.com/matlabcentral/fileexchange/4775-legende-gauss-lobatto-nodes-and-weights/content/
	 *		      lglnodes.m
	 */

	// Standard datatypes
	unsigned int i, j, k, iMax, jMax, kMax, dim, u1,
	             N, rInd, row, Nrows, NsOut,
	             *symmsOut, *Indices;
	int          sd, sN;
	double       norm, swapd,
	             *r, *rold, *rdiff, *wOut, *r_d, *wOut_d, *r_std, *wOut_std,
	             *V, *a, *CM, *eigs;

	// Arbitrary initializations for variables defined in conditionals (to eliminate compiler warnings)
	Nrows = 0;

	N = P+1;

	u1 = 1;
	sd = d;
	sN = N;

	r    = malloc(N * sizeof *r); // free
	wOut = malloc(N * sizeof *wOut); // free

	// Note: GLL must be first as "GL" is in "GLL"
	if (strstr(NodeType,"GLL")) {
		if (P == 0)
			printf("Error: Cannot use GLL nodes of order P0.\n"), exit(1);

		rold  = malloc(N * sizeof *rold); // free
		rdiff = malloc(N * sizeof *rdiff); // free

		// Use the Chebyshve-Guass-Lobatto nodes as the first guess
		for (i = 0; i < N; i++)
			r[i] = -cos(PI*(i)/P);
// array_print_d(1,N,r,'R');

		// Legendre Vandermonde Matrix
		V = malloc(N*N * sizeof *V); // free
		for (i = 0, iMax = N*N; i < iMax; i++)
			V[i] = 0.;

		/* Compute P_(N) using the recursion relation. Compute its first and second derivatives and update r using the
		Newton-Raphson method */
		for (i = 0; i < N; i++) rold[i] = 2.;
		for (i = 0; i < N; i++) rdiff[i] = r[i]-rold[i];
		norm = array_norm_d(N,rdiff,"Inf");

		while (norm > EPS) {
			for (i = 0; i < N; i++)
				rold[i] = r[i];
			for (j = 0; j < N; j++) {
				V[0*N+j] = 1.;
				V[1*N+j] = r[j];
			}

		for (i = 1; i < P; i++) {
		for (j = 0; j < N; j++) {
			V[(i+1)*N+j] = ((2*i+1)*r[j]*V[i*N+j] - i*V[(i-1)*N+j])/(i+1);
		}}

		for (j = 0; j < N; j++)
			r[j] = rold[j] - (r[j]*V[P*N+j]-V[(P-1)*N+j])/(N*V[P*N+j]);

		for (i = 0; i < N; i++)
			rdiff[i] = r[i]-rold[i];
		norm = array_norm_d(N,rdiff,"Inf");
		}

		for (j = 0; j < N; j++)
			wOut[j] = 2./(P*N*pow(V[P*N+j],2));

		free(rold);
		free(rdiff);
		free(V);

// array_print_d(1,N,r,'R');
// array_print_d(1,N,wOut,'R');
	} else if (strstr(NodeType,"GL")) {
		// Build the companion matrix CM
		/* CM is defined such that det(rI-CM)=P_n(r), with P_n(r) being the Legendre poynomial under consideration.
		 * Moreover, CM is constructed in such a way so as to be symmetrical.
		 */
		a = malloc(P * sizeof *a); // free
		for (i = 1; i < N; i++)
			a[i-1] = (1.*i)/sqrt(4.*pow(i,2)-1);

		CM = malloc(N*N * sizeof *CM); // free
		for (i = 0, iMax = N*N; i < iMax; i++)
			CM[i] = 0.;

		for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			if (i == j+1) CM[i*N+j] = a[j];
			if (j == i+1) CM[i*N+j] = a[i];
		}}
		free(a);

		// Determine the abscissas (r) and weights (w)
		/* Because det(rI-CM) = P_n(r), the abscissas are the roots of the characteristic polynomial, the eigenvalues of
		 * CM. The weights can then be derived from the corresponding eigenvectors.
		 */

		eigs    = malloc(N *sizeof *eigs); // free
		Indices = malloc(N *sizeof *Indices); // free
		for (i = 0; i < N; i++)
			Indices[i] = i;

		if (LAPACKE_dsyev(LAPACK_ROW_MAJOR,'V','U',(MKL_INT) N,CM,(MKL_INT) N,eigs) > 0)
			printf("Error: mkl LAPACKE_sysev failed to compute eigenvalues.\n"), exit(1);

// array_print_d(1,N,eigs,'R');
// array_print_d(N,N,CM,'R');

		array_sort_d(1,N,eigs,Indices,'R','N');

		for (i = 0; i < N; i++) {
		for (j = 0; j < N; j++) {
			swapd = CM[i*N+j];
			CM[i*N+j]          = CM[Indices[i]*N+j];
			CM[Indices[i]*N+j] = swapd;
		}}
		free(Indices);

		for (j = 0; j < N; j++) {
			r[j]    = eigs[j];
			wOut[j] = 2*pow(CM[0*i+j],2);
		}
		free(CM);
		free(eigs);

// array_print_d(1,N,r,'R');
// array_print_d(1,N,wOut,'R');
	}

	// Re-arrange r and w for GL/GLL nodes
	r_std    = malloc(N * sizeof *r_std); // free
	wOut_std = malloc(N * sizeof *wOut_std); // free

	for (i = 0; i < N; i++) {
		r_std[i]    = r[i];
		wOut_std[i] = wOut[i];
	}

	if (N % 2 == 1) {
		k = (unsigned int) floor(N/2);

		r[N-1]    = r_std[k];
		wOut[N-1] = wOut_std[k];

		j = k;
		for (i = 0, iMax = N-1; i < iMax; ) {
			r[i]    = r_std[k-j];
			wOut[i] = wOut_std[k-j];
			i++;

			r[i]    = r_std[k+j];
			wOut[i] = wOut_std[k+j];
			i++;

			j--;
		}
	} else {
		k = N/2;
		j = k;
		for (i = 0, iMax = N; i < iMax; ) {
			r[i]    = r_std[k-j];
			wOut[i] = wOut_std[k-j];
			i++;

			r[i]    = r_std[k+j-1];
			wOut[i] = wOut_std[k+j-1];
			i++;

			j--;
		}
	}

	free(r_std);
	free(wOut_std);


	r_d    = malloc(pow(N,d)*d * sizeof *r_d);    // keep (requires external free)
	wOut_d = malloc(pow(N,d)   * sizeof *wOut_d); // free/keep (Conditional return_w)

	row = 0; Nrows = pow(N,d);
	for (k = 0, kMax = (unsigned int) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
	for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
	for (i = 0, iMax = min(max((d-0)*N,u1),N); i < iMax; i++) {
		wOut_d[row] = wOut[i];
		if (d == 2) wOut_d[row] *= wOut[j];
		if (d == 3) wOut_d[row] *= wOut[j]*wOut[k];
		for (dim = 0; dim < d; dim++) {
			if (dim == 0) rInd = i;
			if (dim == 1) rInd = j;
			if (dim == 2) rInd = k;
			r_d[dim*Nrows+row] = r[rInd];
		}
		row++;
	}}}
	free(r);
	free(wOut);

	// Compute symmetries
	NsOut = (unsigned int) ceil(N/2.0);
	symmsOut = malloc(NsOut * sizeof *symmsOut); // keep (requires external free)
	if (N % 2 == 1) {
		for (i = 0; i < NsOut-1; i++)
			symmsOut[i] = 2;
		symmsOut[NsOut-1] = 1;
	} else {
		for (i = 0; i < NsOut; i++)
			symmsOut[i] = 2;
	}

	*rst   = r_d;
	*symms = symmsOut;

	*Nn = pow(N,d);
	*Ns = NsOut;

	if (return_w != 0) *w = wOut_d;
	else               *w = NULL, free(wOut_d);

// array_print_d(pow(N,d),d,r_d,'R');
// array_print_d(pow(N,d),1,wOut_d,'R');
}

void cubature_TRI(double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                  const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType)
{
	/*
	 *	Purpose:
	 *		Return nodes, weights and symmetries for triangular cubature depending on the nodetype.
	 *
	 *	Comments:
	 *		Check that AO nodes correspond to those in pyfr after writing the converter script. (ToBeDeleted)
	 *		The WSH and WV nodes were determined based off of those from the pyfr code (pyfr/quadrules/tri) after being
	 *		transfered to the equilateral reference TRI used in this code.
	 *		The order of rst, w is important for minimizing memory stride while computing the length 3 Discrete Fourier
	 *		Transform (ToBeModified).
	 *		Ordering convention:
	 *			3-blocks of symmetric nodes going from farthest from center towards the center, followed by 1-block of
	 *			center node if present.
	 *		rst is stored in memory as r, s, then t (See cubature_TP for the motivation).
	 *
	 *		WSH nodes have the following [order, cubature strength]:
	 *			[0,1], [1,2], [2,4], [3,5], [4,7], [5,8], [6,10], [7,12], [8,14]
	 *
	 *		Input P for WV nodes is the cubature strength desired.
	 *
	 *	Notation:
	 *		rst   : Nodes array of dimension Nn*d (column-major storage)
	 *		w     : Weights array of dimension Nn*1
	 *		symms : Symmetries array
	 *		Nn    : (N)umber of (n)odes
	 *		Ns    : (N)umber of (s)ymmetries
	 *
	 *	References:
	 *		pyfr code : http://www.pyfr.org
	 *
	 *		AO  : Hesthaven(2008)-Nodal_Discontinuous_Galerkin_Methods (ToBeModified: Likely use same conversion)
	 *		WSH : pyfr/quadrules/tri + conversion to barycentric coordinates (ToBeModified: See python script)
	 *		WV  : pyfr/quadrules/tri + conversion to barycentric coordinates (ToBeModified: See python script)
	 */

	// Standard datatypes
	static unsigned int perms[18] = { 0, 1, 2,
	                                  2, 0, 1,
	                                  1, 2, 0,
	                                  0, 2, 1,
	                                  1, 0, 2,
	                                  2, 1, 0};
	unsigned int symms31[2] = { 0, 0};
	double rst_c[6] = { -1.0,            1.0,           0.0,
	                    -1.0/sqrt(3.0), -1.0/sqrt(3.0), 2.0/sqrt(3.0)};
	unsigned int i, iMax, j, jMax, k,
	             IndB, IndBC, IndGroup, GroupCount, Nc,
	             PMax, NnOut, NsOut, Ngroups, Nsymms;
	unsigned int *symmsOut, *symms_Nperms, *symms_count;
	char         *StringRead, *strings, *stringe, *CubFile, *Pc;
	double       *rstOut, *wOut, *BCoords, *BCoords_complete, *w_read;

	FILE *fID;

	// Silence compiler warnings
	PMax = 0; Nsymms = 0; Ngroups = 0;
	w_read = malloc(0 * sizeof *w_read); // silence (ToBeModified: change this to NULL and remove free below)
	wOut = NULL;

	if (strstr(NodeType,"AO")) {
		if (return_w)
			printf("Error: Invalid value for return_w in cubature_TRI.\n"), exit(1);

		PMax = 15;
	} else if (strstr(NodeType,"WSH")) {
		PMax = 8;
	} else if (strstr(NodeType,"WV")) {
		PMax = 20;
	}

	if (P > PMax)
		printf("Error: %s TRI nodes of order %d are not available.\n",NodeType,P), exit(1);

	Nc = 3;

	CubFile = malloc(STRLEN_MAX * sizeof *CubFile); // free
	Pc      = malloc(STRLEN_MIN * sizeof *Pc);      // free
	sprintf(Pc,"%d",P);

	strcpy(CubFile,"../cubature/tri/");
	strcat(CubFile,NodeType);
	strcat(CubFile,Pc);
	strcat(CubFile,".txt");
	free(Pc);

	if ((fID = fopen(CubFile,"r")) == NULL)
		printf("Error: Cubature file %s not found.\n",CubFile), exit(1);
	free(CubFile);

	StringRead = malloc(STRLEN_MAX * sizeof *StringRead); // free

	// Read header
	if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		strings = StringRead;
		Ngroups = strtol(strings,&stringe,10); strings = stringe;
		Nsymms  = strtol(strings,&stringe,10); strings = stringe;
	}

	symms_count  = malloc(Nsymms * sizeof *symms_count);  // free
	symms_Nperms = malloc(Nsymms * sizeof *symms_Nperms); // free

	for (i = 0; i < Nsymms; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			symms_count[i]  = strtol(strings,&stringe,10); strings = stringe;
			symms_Nperms[i] = strtol(strings,&stringe,10); strings = stringe;
		}
	}

	for (iMax = 1; iMax--; ) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			; // skip 1 line(s)
		}
	}

	BCoords = malloc(Ngroups*Nc * sizeof * BCoords); // free

	for (i = 0; i < Ngroups; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			for (j = 0; j < Nc; j++) {
				BCoords[i*Nc+j] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	// Read weights
	if (return_w) {
		for (iMax = 1; iMax--; ) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				; // skip 1 line(s)
			}
		}

		free(w_read);
		w_read = malloc(Ngroups * sizeof *w_read); // free

		for (i = 0; i < Ngroups; i++) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				strings = StringRead;
				w_read[i] = strtod(strings,&stringe); strings = stringe;
			}
		}
//array_print_d(Ngroups,1,w_read,'R');
	}

	fclose(fID);
	free(StringRead);

	// Convert barycentric coordinates to rst nodal coordinates
	// Also set weights in array of size NnOut x 1 if applicable
	NnOut = 0;
	for (i = 0; i < Nsymms; i++)
		NnOut += symms_count[i]*symms_Nperms[i];

	if (return_w)
		wOut = malloc(NnOut * sizeof *wOut); // free/keep (conditional return_w)

	BCoords_complete = malloc(NnOut*Nc * sizeof *BCoords_complete); // free

	IndB = 0; IndBC = 0;
	GroupCount = 0;

	IndGroup = 0;
	for (i = 0; i < Nsymms; i++) {
		if (symms_count[i] == 0)
			IndGroup++;
		else
			break;
	}

	for (i = 0; i < Ngroups; i++) {
		GroupCount++;

		jMax = symms_Nperms[IndGroup];
		// Count 3/1 symmetries
		if      (jMax == 6) symms31[0] += 2;
		else if (jMax == 3) symms31[0] += 1;
		else if (jMax == 1) symms31[1] += 1;

		for (j = 0; j < jMax; j++) {
			if (return_w)
				wOut[IndBC] = w_read[IndB];

			for (k = 0; k < Nc; k++) {
				BCoords_complete[IndBC*Nc+k] = BCoords[IndB*Nc+perms[j*Nc+k]];
			}
			IndBC++;
		}
		IndB++;

		if (symms_count[IndGroup] == GroupCount) {
			GroupCount = 0;
			IndGroup += 1;
		}
	}
	free(symms_count);
	free(symms_Nperms);
	free(BCoords);
	free(w_read);

	rstOut = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,NnOut,d,Nc,1.0,BCoords_complete,rst_c);
	// keep (requires external free)
	free(BCoords_complete);

//array_print_d(NnOut,d,rstOut,'C');

	NsOut = 0;
	for (i = 0; i < 2; i++)
		NsOut += symms31[i];

	symmsOut = malloc(NsOut * sizeof *symmsOut); // keep (requires external free)
	k = 0;
	for (i = 0; i < 2; i++) {
	for (j = 0; jMax = symms31[i], j < jMax; j++) {
		if (i == 0) symmsOut[k] = 3;
		else        symmsOut[k] = 1;
		k++;
	}}

	*rst   = rstOut;
	*symms = symmsOut;

	*Nn = NnOut;
	*Ns = NsOut;

	if (return_w) {
		*w = wOut;
	} else {
		*w = NULL;
	}
}

void cubature_TET(double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                  const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType)
{
	/*
	 *	Purpose:
	 *		Return nodes, weights and symmetries for tetrahedral cubature depending on the nodetype.
	 *
	 *	Comments:
	 *		Check that AO nodes correspond to those in pyfr after writing the converter script. (ToBeDeleted)
	 *		The WSH and WV nodes were determined based off of those from the pyfr code (pyfr/quadrules/tet) after being
	 *		transfered to the regular TET used in this code.
	 *		The order of rst, w is important for minimizing memory stride while computing the length 3 Discrete Fourier
	 *		Transform (Note: as w is a 1d matrix, its ordering is actually not relevant) (ToBeModified).
	 *		Ordering convention:
	 *			3-blocks of symmetric nodes in the same rotational order, followed by 1-block(s) of center node(s) if
	 *			present.
	 *		rst is stored in memory as r, s, then t (See cubature_TP for the motivation).
	 *
	 *		WSH nodes have the following [order, cubature strength]:
	 *			[0,1], [1,2], [2,3], [3,5], [4,6], [5,8], [6,9]
	 *
	 *		Input P for WV nodes is the cubature strength desired.
	 *
	 *	Notation:
	 *		rst   : Nodes array of dimension Nn*d (column-major storage)
	 *		w     : Weights array of dimension Nn*1
	 *		symms : Symmetries array
	 *		Nn    : (N)umber of (n)odes
	 *		Ns    : (N)umber of (s)ymmetries
	 *
	 *	References:
	 *		pyfr code : http://www.pyfr.org
	 *
	 *		AO  : Hesthaven(2008)-Nodal_Discontinuous_Galerkin_Methods (ToBeModified: Likely use same conversion)
	 *		WSH : pyfr/quadrules/tet + conversion to barycentric coordinates (ToBeModified: See python script)
	 *		WV  : pyfr/quadrules/tet + conversion to barycentric coordinates (ToBeModified: See python script)
	 */

	// Standard datatypes
	static unsigned int perms_TRI[18] = { 0, 1, 2,
	                                      2, 0, 1,
	                                      1, 2, 0,
	                                      0, 2, 1,
	                                      1, 0, 2,
	                                      2, 1, 0},
	                    perms_TET[16] = { 0, 1, 2, 3,
	                                      3, 0, 1, 2,
	                                      2, 3, 0, 1,
	                                      1, 2, 3, 0};
	unsigned int symms31[2] = { 0, 0}, NTRIsymms[3], TETperms[4];
	double rst_c[12] = { -1.0,            1.0,            0.0,           0.0,
	                     -1.0/sqrt(3.0), -1.0/sqrt(3.0),  2.0/sqrt(3.0), 0.0,
	                     -1.0/sqrt(6.0), -1.0/sqrt(6.0), -1.0/sqrt(6.0), 3.0/sqrt(6.0)},
	       BCoords_tmp[4];
	unsigned int i, iMax, j, jMax, k, kMax, l, lMax, TRIsymm,
	             IndB, IndBC, IndGroup, Ind1, Indperm, GroupCount, Nc, N1,
	             PMax, NnOut, NsOut, Ngroups, Nsymms, Nperms;
	unsigned int *symmsOut, *symms_Nperms, *symms_count;
	char         *StringRead, *strings, *stringe, *CubFile, *Pc;
	double       *rstOut, *wOut, *BCoords, *BCoords_complete, *w_read;

	FILE *fID;

	// Silence compiler warnings
	PMax = 0; Nsymms = 0; Ngroups = 0;
	w_read = malloc(0 * sizeof *w_read); // silence

	if (strstr(NodeType,"AO")) {
		if (return_w)
			printf("Error: Invalid value for return_w in cubature_TET.\n"), exit(1);

		PMax = 15;
	} else if (strstr(NodeType,"WSH")) {
		PMax = 6;
	} else if (strstr(NodeType,"WV")) {
		PMax = 10;
	}

	if (P > PMax)
		printf("Error: %s TET nodes of order %d are not available.\n",NodeType,P), exit(1);

	Nc = 4;

	CubFile = malloc(STRLEN_MAX * sizeof *CubFile); // free
	Pc      = malloc(STRLEN_MIN * sizeof *Pc);      // free
	sprintf(Pc,"%d",P);

	strcpy(CubFile,"../cubature/tet/");
	strcat(CubFile,NodeType);
	strcat(CubFile,Pc);
	strcat(CubFile,".txt");
	free(Pc);

	if ((fID = fopen(CubFile,"r")) == NULL)
		printf("Error: Cubature file %s not found.\n",CubFile), exit(1);
	free(CubFile);

	StringRead = malloc(STRLEN_MAX * sizeof *StringRead); // free

	// Read header
	if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		strings = StringRead;
		Ngroups = strtol(strings,&stringe,10); strings = stringe;
		Nsymms  = strtol(strings,&stringe,10); strings = stringe;
	}

	symms_count  = malloc(Nsymms * sizeof *symms_count);  // free
	symms_Nperms = malloc(Nsymms * sizeof *symms_Nperms); // free

	for (i = 0; i < Nsymms; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			symms_count[i]  = strtol(strings,&stringe,10); strings = stringe;
			symms_Nperms[i] = strtol(strings,&stringe,10); strings = stringe;
		}
	}

	for (iMax = 1; iMax--; ) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			; // skip 1 line(s)
		}
	}

	BCoords = malloc(Ngroups*Nc * sizeof * BCoords); // free

	for (i = 0; i < Ngroups; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			for (j = 0; j < Nc; j++) {
				BCoords[i*Nc+j] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	// Read weights
	if (return_w) {
		for (iMax = 1; iMax--; ) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				; // skip 1 line(s)
			}
		}

		free(w_read);
		w_read = malloc(Ngroups * sizeof *w_read); // free

		for (i = 0; i < Ngroups; i++) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				strings = StringRead;
				w_read[i] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	fclose(fID);
	free(StringRead);

	// Convert barycentric coordinates to rst nodal coordinates
	// Also set weights in array of size NnOut x 1 if applicable

	// Find number of 1 symmetries and total number of nodes
	for (i = 0, N1 = 0, NnOut = 0; i < Nsymms; i++) {
		if (i >= Nsymms-2)
			N1 += symms_count[i];
		NnOut += symms_count[i]*symms_Nperms[i];
	}

	if (return_w)
		wOut = malloc(NnOut * sizeof *wOut); // free/keep (conditional return_w)

	BCoords_complete = malloc(NnOut*Nc * sizeof *BCoords_complete); // free

	IndB = 0; IndBC = 0;
	GroupCount = 0;

	IndGroup = 0;
	for (i = 0; i < Nsymms; i++) {
		if (symms_count[i] == 0)
			IndGroup++;
		else
			break;
	}

	Ind1 = 0;
	for (i = 0; i < Ngroups; i++) {
		GroupCount++;

		Nperms = symms_Nperms[IndGroup];
		// Count 3/1 symmetries and establish TRI symmetries to loop over
		if (Nperms == 24) {
			symms31[0] += 8;

			NTRIsymms[0] = 0; NTRIsymms[1] = 0; NTRIsymms[2] = 4;
			TETperms[0]  = 0; TETperms[1]  = 1; TETperms[2]  = 2; TETperms[3]  = 3;
		} else if (Nperms == 12) {
			symms31[0] += 4;

			NTRIsymms[0] = 0; NTRIsymms[1] = 2; NTRIsymms[2] = 1;
			TETperms[0]  = 3; TETperms[1]  = 2; TETperms[2]  = 0; TETperms[3]  = 9;
		} else if (Nperms == 6) {
			symms31[0] += 2;

			NTRIsymms[0] = 0; NTRIsymms[1] = 2; NTRIsymms[2] = 0;
			TETperms[0]  = 0; TETperms[1]  = 3; TETperms[2]  = 9; TETperms[3]  = 9;
		} else if (Nperms == 4) {
			symms31[0] += 1, symms31[1] += 1;

			NTRIsymms[0] = 1; NTRIsymms[1] = 1; NTRIsymms[2] = 0;
			TETperms[0]  = 3; TETperms[1]  = 0; TETperms[2]  = 9; TETperms[3]  = 9;
		} else if (Nperms == 1) {
			symms31[1] += 1;

			NTRIsymms[0] = 1; NTRIsymms[1] = 0; NTRIsymms[2] = 0;
			TETperms[0]  = 0; TETperms[1]  = 9; TETperms[2]  = 9; TETperms[3]  = 9;
		}

		Indperm = 0;
		if (NTRIsymms[0]) {
			for (j = 0; j < Nc; j++)
				BCoords_tmp[j] = BCoords[IndB*Nc+perms_TET[TETperms[Indperm]*Nc+j]];
			Indperm++;

			for (k = 0; k < Nc; k++)
				BCoords_complete[(NnOut-N1+Ind1)*Nc+k] = BCoords_tmp[k];

			if (return_w)
				wOut[NnOut-N1+Ind1] = w_read[IndB];

			Ind1++;
		}

		for (TRIsymm = 1; TRIsymm <= 2; TRIsymm++) {
			for (l = 0, lMax = NTRIsymms[TRIsymm]; l < lMax; l++) {

				for (j = 0; j < Nc; j++)
					BCoords_tmp[j] = BCoords[IndB*Nc+perms_TET[TETperms[Indperm]*Nc+j]];
				Indperm++;

				jMax = 0;
				for (k = 1, kMax = TRIsymm+1; k <= kMax; k++)
					jMax += k;

				for (j = 0; j < jMax; j++) {
					if (return_w)
						wOut[IndBC] = w_read[IndB];

					for (k = 0; k < Nc-1; k++)
						BCoords_complete[IndBC*Nc+k] = BCoords_tmp[perms_TRI[j*(Nc-1)+k]];

					k = Nc-1;
					BCoords_complete[IndBC*Nc+k] = BCoords_tmp[3];

					IndBC++;
				}
			}
		}
		IndB++;

		if (symms_count[IndGroup] == GroupCount) {
			GroupCount = 0;
			IndGroup++;
			while (IndGroup < Nsymms && symms_count[IndGroup] == 0)
				IndGroup++;
		}
	}
//array_print_d(NnOut,Nc,BCoords_complete,'R');

	free(symms_count);
	free(symms_Nperms);
	free(BCoords);
	free(w_read);

	rstOut = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,NnOut,d,Nc,1.0,BCoords_complete,rst_c);
	// keep (requires external free)
	free(BCoords_complete);

//array_print_d(NnOut,d,rstOut,'C');

	NsOut = 0;
	for (i = 0; i < 2; i++)
		NsOut += symms31[i];

	symmsOut = malloc(NsOut * sizeof *symmsOut); // keep (requires external free)
	k = 0;
	for (i = 0; i < 2; i++) {
	for (j = 0; jMax = symms31[i], j < jMax; j++) {
		if (i == 0) symmsOut[k] = 3;
		else        symmsOut[k] = 1;
		k++;
	}}

	*rst   = rstOut;
	*symms = symmsOut;

	*Nn = NnOut;
	*Ns = NsOut;

	if (return_w) {
		*w = wOut;
	} else {
		*w = NULL;
	}
}

void cubature_PYR(double **rst, double **w, unsigned int **symms, unsigned int *Nn, unsigned int *Ns,
                  const unsigned int return_w, const unsigned int P, const unsigned int d, const char *NodeType)
{
	/*
	 *	Purpose:
	 *		Return nodes, weights and symmetries for pyramidal cubature depending on the nodetype.
	 *
	 *	Comments:
	 *		Barycentric coordinates for the PYR element may be negative, unlike for the SI elements. This does not pose
	 *		any challenges or result in any limitations, but is simply a result of the last additional condition (chosen
	 *		for good conditioning of the matrix) when computing the barycentric coordinates of the nodes in the pyfr
	 *		code.
	 *
	 *		Note:
	 *
	 *			[ones]_{Nn x 1}  = [BCoords]_{Nn x Nc} * [ones]_{Nc x 1}          => Partition of unity (1 condition)
	 *			[rst]_{Nn x d}   = [BCoords]_{Nn x Nc} * [rst_(c)orners]_{Nc x d} => Linear precision   (d conditions)
	 *			[r*s*t]_{Nn x 1} = [BCoords]_{Nn x Nc} * [{r*s*t}_c]_{Nc x 1}     => Arbitrary          (1 condition)
	 *
	 *			Then
	 *				[BCoords] = [ones(Nn,1) rst r*s*t]*inv([ones(Nc,1) rst_c {r*s*t}_c])
	 *
	 *			where
	 *				cond([ones(Nc,1) rst_c {r*s*t}_c]) = 2.0
	 *
	 *		All nodes were determined based off of those from the pyfr code (pyfr/quadrules/pyr) after being transfered
	 *		to the regular PYR used in this code. (ToBeModified if other nodes are added)
	 *			Options:
	 *				GL     : GL nodes,                       no weights
	 *				GLL    : GLL nodes,                      no weights
	 *				GLW    : GL nodes,                       with weights
	 *				GLLW   : GLL nodes,                      with weights
	 *				WV     : WV PYR nodes,                   with weights
	 *					Exact integration to lower order than expected.
	 *				WVHToP : WV HEX nodes transfered to PYR, with weights
	 *					Exact integration to lower order than expected.
	 *
	 *			After implementing the traditional PYR orthogonal basis
	 *			(Chan(2015)-Orthogonal_Bases_for_Vertex_Mapped_Pyramids, eq. 2.1) and the orthogonal basis from the pyfr
	 *			code (Witherden(2015,Thesis), eq. 3.20), the following conclusions were drawn:
	 *				Despite being lower-order in the "c" term, the traditional basis mass matrix is not integrated
	 *				exactly using the WV nodes, while the basis in the pyfr code is (eventually); this is odd. Further,
	 *				given the maximum strength WV rule (WV10), even the P3 PYR mass matrix is not given exactly with
	 *				constant Jacobian. Further testing revealed that the WV nodes integrate exactly (with the expected
	 *				order) only when there is variation either in a OR in b, but not in both; the nodes cannot exactly
	 *				integrate polynomials on QUAD cross-sections of the pyramid, unlike the GL nodes.
	 *
	 *				Using GL nodes transfered to the PYR element, the cubature strength for exact mass matrix is as
	 *				expected, after accounting for the added contribution to the "c" term from the weight (w =
	 *				w_HEX*pow(1-c,2); GL nodes of order P+1 are required for exact integration of all terms). Using the
	 *				WV HEX nodes transfered to the PYR element does not result in exact integration (similarly to the
	 *				PYR nodes) either for the traditional or pyfr PYR basis.
	 *
	 *		The order of rst, w is important for minimizing memory stride while computing the length 4 Discrete Fourier
	 *		Transform (Note: as w is a 1d matrix, its ordering is actually not relevant) (ToBeModified).
	 *		Ordering convention:
	 *			4-blocks of symmetric nodes in the same rotational order, followed by 1-block(s) of center node(s) if
	 *			present.
	 *		rst is stored in memory as r, s, then t (See cubature_TP for the motivation).
	 *
	 *		Input P for WV nodes is the cubature strength desired.
	 *
	 *	Notation:
	 *		rst   : Nodes array of dimension Nn*d (column-major storage)
	 *		w     : Weights array of dimension Nn*1
	 *		symms : Symmetries array
	 *		Nn    : (N)umber of (n)odes
	 *		Ns    : (N)umber of (s)ymmetries
	 *
	 *	References:
	 *		ToBeModified: Add reference to Witherden(2015)-On_the_Identification_of_Symmetric_...
	 *		pyfr code : http://www.pyfr.org
	 *
	 *		GL  : pyfr/quadrules/pyr + conversion to barycentric coordinates (ToBeModified: See python script)
	 *		GLL : pyfr/quadrules/pyr + conversion to barycentric coordinates (ToBeModified: See python script)
	 *		WV  : pyfr/quadrules/pyr + conversion to barycentric coordinates (ToBeModified: See python script)
	 */

	// Standard datatypes
	static unsigned int perms_QUAD[16] = { 0, 1, 2, 3,
	                                       3, 0, 1, 2,
	                                       2, 3, 0, 1,
	                                       1, 2, 3, 0};
	unsigned int symms41[2] = { 0, 0}, NQUADsymms = 0;
	double rst_c[15] = { -1.0,            1.0,            1.0,           -1.0,           0.0,
	                     -1.0,           -1.0,            1.0,            1.0,           0.0,
	                     -sqrt(2.0)/5.0, -sqrt(2.0)/5.0, -sqrt(2.0)/5.0, -sqrt(2.0)/5.0, 4.0/5.0*sqrt(2.0)},
	       BCoords_tmp[5];
	unsigned int i, iMax, j, jMax, k, QUADsymm,
	             IndB, IndBC, IndGroup, Ind1, GroupCount, Nc, N1,
	             PMax, NnOut, NsOut, Ngroups, Nsymms, Nperms;
	unsigned int *symmsOut, *symms_Nperms, *symms_count;
	char         *StringRead, *strings, *stringe, *CubFile, *Pc;
	double       *rstOut, *wOut, *BCoords, *BCoords_complete, *w_read;

	FILE *fID;

	// Silence compiler warnings
	PMax = 0; Nsymms = 0; Ngroups = 0;
	w_read = malloc(0 * sizeof *w_read); // silence
	wOut = NULL;

	if (strstr(NodeType,"GL")) {
		if (strstr(NodeType,"W") == NULL) {
			if (return_w)
				printf("Error: Invalid value for return_w in cubature_PYR.\n"), exit(1);

			PMax = 6;
		} else {
			PMax = 6;
		}
	} else if (strstr(NodeType,"WV")) {
		if (strstr(NodeType,"HToP") == NULL)
			PMax = 10;
		else
			PMax = 11;
	}

	if (P > PMax)
		printf("Error: %s PYR nodes of order %d are not available.\n",NodeType,P), exit(1);

	Nc = 5;

	CubFile = malloc(STRLEN_MAX * sizeof *CubFile); // free
	Pc      = malloc(STRLEN_MIN * sizeof *Pc);      // free
	sprintf(Pc,"%d",P);

	strcpy(CubFile,"../cubature/pyr/");
	strcat(CubFile,NodeType);
	strcat(CubFile,Pc);
	strcat(CubFile,".txt");
	free(Pc);

	if ((fID = fopen(CubFile,"r")) == NULL)
		printf("Error: Cubature file %s not found.\n",CubFile), exit(1);
	free(CubFile);

	StringRead = malloc(STRLEN_MAX * sizeof *StringRead); // free

	// Read header
	if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
		strings = StringRead;
		Ngroups = strtol(strings,&stringe,10); strings = stringe;
		Nsymms  = strtol(strings,&stringe,10); strings = stringe;
	}

	symms_count  = malloc(Nsymms * sizeof *symms_count);  // free
	symms_Nperms = malloc(Nsymms * sizeof *symms_Nperms); // free

	for (i = 0; i < Nsymms; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			symms_count[i]  = strtol(strings,&stringe,10); strings = stringe;
			symms_Nperms[i] = strtol(strings,&stringe,10); strings = stringe;
		}
	}

	for (iMax = 1; iMax--; ) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			; // skip 1 line(s)
		}
	}

	BCoords = malloc(Ngroups*Nc * sizeof * BCoords); // free

	for (i = 0; i < Ngroups; i++) {
		if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
			strings = StringRead;
			for (j = 0; j < Nc; j++) {
				BCoords[i*Nc+j] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	// Read weights
	if (return_w) {
		for (iMax = 1; iMax--; ) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				; // skip 1 line(s)
			}
		}

		free(w_read);
		w_read = malloc(Ngroups * sizeof *w_read); // free

		for (i = 0; i < Ngroups; i++) {
			if (fscanf(fID,"%[^\n]\n",StringRead) == 1) {
				strings = StringRead;
				w_read[i] = strtod(strings,&stringe); strings = stringe;
			}
		}
	}

	fclose(fID);
	free(StringRead);

	// Convert barycentric coordinates to rst nodal coordinates
	// Also set weights in array of size NnOut x 1 if applicable

	// Find number of 1 symmetries and total number of nodes
	for (i = 0, N1 = 0, NnOut = 0; i < Nsymms; i++) {
		if (i >= Nsymms-1)
			N1 += symms_count[i];
		NnOut += symms_count[i]*symms_Nperms[i];
	}

	if (return_w)
		wOut = malloc(NnOut * sizeof *wOut); // free/keep (conditional return_w)

	BCoords_complete = malloc(NnOut*Nc * sizeof *BCoords_complete); // free

	IndB = 0; IndBC = 0;
	GroupCount = 0;

	IndGroup = 0;
	for (i = 0; i < Nsymms; i++) {
		if (symms_count[i] == 0)
			IndGroup++;
		else
			break;
	}

	Ind1 = 0;
	for (i = 0; i < Ngroups; i++) {
		GroupCount++;

		Nperms = symms_Nperms[IndGroup];
		// Count 4/1 symmetries and establish QUAD symmetries to loop over
		if (Nperms == 8) {
			symms41[0] += 2;
			NQUADsymms = 1;
		} else if (Nperms == 4) {
			symms41[0] += 1;
			NQUADsymms = 0;
		} else if (Nperms == 1) {
			symms41[1] += 1;
			NQUADsymms = 0;
		}

		if (Nperms == 1) {
			for (j = 0; j < Nc; j++)
				BCoords_tmp[j] = BCoords[IndB*Nc+j];

			for (k = 0; k < Nc; k++)
				BCoords_complete[(NnOut-N1+Ind1)*Nc+k] = BCoords_tmp[k];

			if (return_w)
				wOut[NnOut-N1+Ind1] = w_read[IndB];

			Ind1++;
		} else {
			for (QUADsymm = 0; QUADsymm <= NQUADsymms; QUADsymm++) {
				if (QUADsymm == 0) {
					for (j = 0; j < Nc; j++)
						BCoords_tmp[j] = BCoords[IndB*Nc+j];
				} else {
					BCoords_tmp[0] = BCoords[IndB*Nc+1];
					BCoords_tmp[1] = BCoords[IndB*Nc+0];
					BCoords_tmp[2] = BCoords[IndB*Nc+3];
					BCoords_tmp[3] = BCoords[IndB*Nc+2];
					BCoords_tmp[4] = BCoords[IndB*Nc+4];
				}

				for (j = 0; j < 4; j++) {
					if (return_w)
						wOut[IndBC] = w_read[IndB];

					for (k = 0; k < Nc-1; k++)
						BCoords_complete[IndBC*Nc+k] = BCoords_tmp[perms_QUAD[j*(Nc-1)+k]];

					k = Nc-1;
					BCoords_complete[IndBC*Nc+k] = BCoords_tmp[4];

					IndBC++;
				}
			}
		}
		IndB++;

		if (symms_count[IndGroup] == GroupCount) {
			GroupCount = 0;
			IndGroup += 1;
		}
	}
//array_print_d(NnOut,Nc,BCoords_complete,'R');

	free(symms_count);
	free(symms_Nperms);
	free(BCoords);
	free(w_read);

	rstOut = mm_Alloc_d(CblasColMajor,CblasTrans,CblasNoTrans,NnOut,d,Nc,1.0,BCoords_complete,rst_c);
	// keep (requires external free)
	free(BCoords_complete);

//array_print_d(NnOut,d,rstOut,'C');

	NsOut = 0;
	for (i = 0; i < 2; i++)
		NsOut += symms41[i];

	symmsOut = malloc(NsOut * sizeof *symmsOut); // keep (requires external free)
	k = 0;
	for (i = 0; i < 2; i++) {
	for (j = 0; jMax = symms41[i], j < jMax; j++) {
		if (i == 0) symmsOut[k] = 4;
		else        symmsOut[k] = 1;
		k++;
	}}

	*rst   = rstOut;
	*symms = symmsOut;

	*Nn = NnOut;
	*Ns = NsOut;

	if (return_w) {
		*w = wOut;
	} else {
		free(wOut);
		*w = NULL;
	}
}