// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "plotting_element_info.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "mkl.h"

#include "Parameters.h"
#include "Macros.h"

#include "matrix_functions.h"

/*
 *	Purpose:
 *		Return nodes and connectivity information required for plotting with paraview using the VTK format.
 *
 *	Comments:
 *		rst is stored in memory as r, s, then t.
 *		The ordering of 'Con' is determined by the by the cell type ordering in VTK for the following elements:
 *			QUAD, HEX, WEDGE, PYR.
 *
 *	Notation:
 *
 *	References:
 *		VTK User's Guide, File Formats for VTK Version 4.2 (Figure 2. Linear cell types found in VTK)
 */

static unsigned int get_Ne_type(const unsigned int type)
{
	switch (type) {
		case LINE:  return 0; break;
		case TRI:   return 3; break;
		case QUAD:  return 4; break;
		case TET:   return 6; break;
		case HEX:   return 12; break;
		case WEDGE: return 9; break;
		case PYR:   return 8; break;
		default:
			printf("Error: Unsupported.\n"), EXIT_MSG;
			break;
	}
}

void plotting_element_info(double **rst, unsigned int **connect, unsigned int **types, unsigned int **connectE,
                           unsigned int *Nn, unsigned int *NE, const unsigned int P, const unsigned int typeIn)
{
	unsigned int i, j, k, l, m, iMax, jMax, kMax, lMax, row, j0, jS, jC, count, N,
	             d, Nc, NnOut, NEOut, Ne,
	             layer, lBs, lTs, sum, u1 = 1,
	             *connectOut, *typesOut, *connectEOut;
	double di, dj, dk, *rstOut;

	if (P == 0)
		printf("Error: Input P must be greater than 0 for plotting nodes.\n"), EXIT_MSG;

	Ne = get_Ne_type(typeIn);

	if (typeIn == LINE)
		connectEOut = NULL;
	else
		connectEOut = calloc(Ne*(P+1) , sizeof *connectEOut); // keep

	if (typeIn == LINE || typeIn == QUAD || typeIn == HEX) {
		unsigned int dim,
		             Indr, Indc,
		             N2,
		             nLINE[2], nQUAD[4], nHEX[8];
		int          sd, sP, sN;
		double *r;

		// Arbitrary initializations for variables defined in conditionals (to eliminate compiler warnings)
		d = 0;

		if      (typeIn == LINE) d = 1;
		else if (typeIn == QUAD) d = 2;
		else if (typeIn == HEX)  d = 3;

		N = P+1;
		NnOut = pow(N,d);
		NEOut = pow(P,d);

		u1 = 1;
		sd = d;
		sP = P;
		sN = N;

		r      = malloc(N       * sizeof *r);      // free
		rstOut = malloc(NnOut*d * sizeof *rstOut); // keep (requires external free)

		connectOut  = calloc(NEOut*8 , sizeof *connectOut); // keep (requires external free)
		typesOut    = malloc(NEOut   * sizeof *typesOut);   // keep (requires external free)

		if (P == 0) {
			r[0] = 0;
		} else {
			for (i = 0; i < N; i++)
				r[i] = -1.0 + (2.0*i)/P;
		}

		row = 0;
		for (k = 0, kMax = (unsigned int) min(max((sd-2)*sN,1),sN); k < kMax; k++) {
		for (j = 0, jMax = min(max((d-1)*N,u1),N); j < jMax; j++) {
		for (i = 0, iMax = min(max((d-0)*N,u1),N); i < iMax; i++) {
			for (dim = 0; dim < d; dim++) {
				if (dim == 0) Indr = i;
				if (dim == 1) Indr = j;
				if (dim == 2) Indr = k;
				rstOut[dim*NnOut+row] = r[Indr];
			}
			row++;
		}}}
		free(r);

		N2 = pow(N,2);

		Indc = 0;
		for (k = 0, kMax = (unsigned int) max(sP*min(sd-2,1),1); k < kMax; k++) {
			for (j = 0, jMax = max(P*min(d-1,u1),u1); j < jMax; j++) {
				for (i = 0; i < P; i++) {
					nLINE[0] = i;
					nLINE[1] = i+1;

					for (l = 0; l < 2; l++)             nQUAD[l] = nLINE[l]+N*j;
					for (l = 2, m = 1; l < 4; l++, m--) nQUAD[l] = nLINE[m] + N*(j+1);

					for (l = 0; l < 4; l++)             nHEX[l] = nQUAD[l] + N2*k;
					for (l = 4, m = 0; l < 8; l++, m++) nHEX[l] = nQUAD[m] + N2*(k+1);

					for (l = 0, lMax = pow(2,d); l < lMax; l++)
						connectOut[Indc*8+l] = nHEX[l];
					Indc++;
				}
			}
		}

		if (typeIn == LINE) {
			for (i = 0; i < NEOut; i++)
				typesOut[i] = 3;
		} else if (typeIn == QUAD) {
			for (i = 0; i < NEOut; i++)
				typesOut[i] = 9;
			for (i = 0; i < Ne; i++) {
				switch (i) {
					case 0: j0 = 0;   jS = N; break;
					case 1: j0 = P;   jS = N; break;
					case 2: j0 = 0;   jS = 1; break;
					case 3: j0 = P*N; jS = 1; break;
					default:
						printf("Error: Unsupported.\n"), EXIT_MSG; break;
				}
				for (j = j0, count = 0; count < N; j += jS, count++)
					connectEOut[i*N+count] = j;
			}
		} else if (typeIn == HEX) {
			for (i = 0; i < NEOut; i++)
				typesOut[i] = 12;
			for (i = 0; i < Ne; i++) {
				switch (i) {
					case 0:  j0 = 0;         jS = 1;   break;
					case 1:  j0 = P*N;       jS = 1;   break;
					case 2:  j0 = P*N*N;     jS = 1;   break;
					case 3:  j0 = (N*N-1)*N; jS = 1;   break;
					case 4:  j0 = 0;         jS = N;   break;
					case 5:  j0 = P;         jS = N;   break;
					case 6:  j0 = P*N*N;     jS = N;   break;
					case 7:  j0 = P*N*N+P;   jS = N;   break;
					case 8:  j0 = 0;         jS = N*N; break;
					case 9:  j0 = P;         jS = N*N; break;
					case 10: j0 = P*N;       jS = N*N; break;
					case 11: j0 = P*N+P;     jS = N*N; break;
					default:
						printf("Error: Unsupported.\n"), EXIT_MSG; break;
				}
				for (j = j0, count = 0; count < N; j += jS, count++)
					connectEOut[i*N+count] = j;
			}
		}

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connectE = connectEOut;
		*Nn       = NnOut;
		*NE       = NEOut;
	} else if (typeIn == TRI) {
		d = 2;
		Nc = 3;
		NnOut = 1.0/2.0*(P+1)*(P+2);
		NEOut = 0;
		for (i = 0; i < P; i++)
			NEOut += 2*i+1;

		unsigned int iStart;
		double rst_V[Nc*d], BCoords[NnOut*Nc];

		connectOut = calloc(NEOut*8 , sizeof *connectOut); // keep (requires external free)
		typesOut   = malloc(NEOut   * sizeof *typesOut); // keep (requires external free)

		// Determine barycentric coordinates of equally spaced nodes
		row = 0;
		for (j = 0; j <= P; j++) {
		for (i = 0, iMax = P-j; i <= iMax; i++) {
			di = i;
			dj = j;

			BCoords[row*Nc+2] = dj/P;
			BCoords[row*Nc+1] = di/P;
			BCoords[row*Nc+0] = 1.0 - (BCoords[row*Nc+1]+BCoords[row*Nc+2]);

			row++;
		}}

		// TRI vertex nodes
		rst_V[0*d+0] = -1.0; rst_V[0*d+1] = -1.0/sqrt(3.0);
		rst_V[1*d+0] =  1.0; rst_V[1*d+1] = -1.0/sqrt(3.0);
		rst_V[2*d+0] =  0.0; rst_V[2*d+1] =  2.0/sqrt(3.0);

		rstOut = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NnOut,d,Nc,1.0,BCoords,rst_V);
		// keep (requires external free)

		// Convert to column-major ordering
		mkl_dimatcopy('R','T',NnOut,d,1.0,rstOut,d,NnOut);

//array_print_d(NnOut,Nc,BCoords,'R');
//array_print_d(NnOut,d,rstOut,'R');

		row = 0;
		// Regular TRIs
		for (j = P; j; j--) {
			iStart = 0;
			for (l = P+1, m = j; P-m; m++, l--)
				iStart += l;

			for (i = iStart, iMax = iStart+j; i < iMax; i++) {
				connectOut[row*8+0] = i;
				connectOut[row*8+1] = i+1;
				connectOut[row*8+2] = i+1+j;
				row++;
			}
		}

		// Inverted  TRIs
		for (j = P; j; j--) {
			iStart = 1;
			for (l = P+1, m = j; P-m; m++, l--)
				iStart += l;

			for (i = iStart, iMax = iStart+j-1; i < iMax; i++) {
				connectOut[row*8+0] = i;
				connectOut[row*8+1] = i+j;
				connectOut[row*8+2] = i+j+1;
				row++;
			}
		}

//array_print_ui(NEOut,8,connectOut,'R');

		for (i = 0; i < NEOut; i++)
			typesOut[i] = 5;

		N = P+1;
		for (i = 0; i < Ne; i++) {
			switch (i) {
				case 0:  j0 = P; jS = P; break;
				case 1:  j0 = 0; jS = N; break;
				case 2:  j0 = 0; jS = 1; break;
				default:
					printf("Error: Unsupported.\n"), EXIT_MSG; break;
			}
			for (j = j0, count = 0; count < N; j += jS, count++) {
				connectEOut[i*N+count] = j;

				if (j != j0 && jS > 1)
					jS--;
			}
		}

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connectE = connectEOut;
		*Nn       = NnOut;
		*NE       = NEOut;
	} else if (typeIn == TET) {
		d = 3;
		Nc = 4;
		NnOut = 1.0/6.0*(P+1)*(P+2)*(P+3);
		NEOut = 0;
		for (i = 1; i <= P; i++) {
			// Regular TETs
			for (j = 1; j <= i; j++)
				NEOut += j;

			// PYRs
			NEOut += i*(i-1);

			// Inverted TETs
			if (i > 2)
				for (j = 1; j <= (i-2); j++)
					NEOut += j;
		}

		double BCoords[NnOut*Nc], rst_V[Nc*d];

		connectOut = calloc(NEOut*8 , sizeof *connectOut); // keep (requires external free)
		typesOut   = malloc(NEOut   * sizeof *typesOut); // keep (requires external free)

		// Determine barycentric coordinates of equally spaced nodes
		row = 0;
		for (k = 0; k <= P; k++) {
		for (j = 0, jMax = P-k; j <= jMax; j++) {
		for (i = 0, iMax = P-(j+k); i <= iMax; i++) {
			di = i;
			dj = j;
			dk = k;

			BCoords[row*Nc+3] = dk/P;
			BCoords[row*Nc+2] = dj/P;
			BCoords[row*Nc+1] = di/P;
			BCoords[row*Nc+0] = 1.0 - (BCoords[row*Nc+1]+BCoords[row*Nc+2]+BCoords[row*Nc+3]);

			row++;
		}}}

		// TET vertex nodes
		rst_V[0*d+0] = -1.0; rst_V[0*d+1] = -1.0/sqrt(3.0); rst_V[0*d+2] = -1.0/sqrt(6);
		rst_V[1*d+0] =  1.0; rst_V[1*d+1] = -1.0/sqrt(3.0); rst_V[1*d+2] = -1.0/sqrt(6);
		rst_V[2*d+0] =  0.0; rst_V[2*d+1] =  2.0/sqrt(3.0); rst_V[2*d+2] = -1.0/sqrt(6);
		rst_V[3*d+0] =  0.0; rst_V[3*d+1] =  0.0          ; rst_V[3*d+2] =  3.0/sqrt(6);

		rstOut = mm_Alloc_d(CblasRowMajor,CblasNoTrans,CblasNoTrans,NnOut,d,Nc,1.0,BCoords,rst_V);
		// keep (requires external free)

		// Convert to column-major ordering
		mkl_dimatcopy('R','T',NnOut,d,1.0,rstOut,d,NnOut);

//array_print_d(NnOut,Nc,BCoords,'R');
//array_print_d(NnOut,d,rstOut,'R');

		row = 0;
		for (layer = P; layer; layer--) {
			lBs = 0;
			for (i = P; i > layer; i--) {
				sum = 0;
				for (j = 1; j <= i+1; j++)
					sum += j;
				lBs += sum;
			}

			lTs = 0;
			for (i = P; i > (layer-1); i--) {
				sum = 0;
				for (j = 1; j <= i+1; j++)
					sum += j;
				lTs += sum;
			}

			// Regular TETs
			for (j = 1; j <= layer; j++) {
			for (i = 0, iMax = layer-j; i <= iMax; i++) {
				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+2-k;

				connectOut[row*8+0] = lBs + i + sum;
				connectOut[row*8+1] = connectOut[row*8+0] + 1;
				connectOut[row*8+2] = connectOut[row*8+0] + layer+2-j;

				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+1-k;

				connectOut[row*8+3] = lTs + i + sum;

				row++;
			}}

			// PYRs
			for (j = 1, jMax = layer-1; j <= jMax; j++) {
			for (i = 1, iMax = layer-j; i <= iMax; i++) {
				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+2-k;

				connectOut[row*8+4] = lBs + i + sum;
				connectOut[row*8+0] = connectOut[row*8+4] + layer+1-j;
				connectOut[row*8+1] = connectOut[row*8+0] + 1;

				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+1-k;

				connectOut[row*8+3] = lTs + i-1 + sum;
				connectOut[row*8+2] = connectOut[row*8+3] + 1;

				row++;

				connectOut[row*8+0] = connectOut[(row-1)*8+0];
				connectOut[row*8+1] = connectOut[(row-1)*8+1];
				connectOut[row*8+2] = connectOut[(row-1)*8+2];
				connectOut[row*8+3] = connectOut[(row-1)*8+3];
				connectOut[row*8+4] = connectOut[row*8+3] + layer+1-j;

				row++;
			}}

			// Inverted TETs
			for (j = 2, jMax = layer-1; j <= jMax; j++) {
			for (i = 1, iMax = layer-j; i <= iMax; i++) {
				sum = 0;
				for (k = 1, kMax = j-1; k <= kMax; k++)
					sum += layer+2-k;

				connectOut[row*8+0] = lBs + i +sum;

				sum = 0;
				for (k = 1, kMax = j-2; k <= kMax; k++)
					sum += layer+1-k;

				connectOut[row*8+1] = lTs + i + sum;
				connectOut[row*8+2] = connectOut[row*8+1] + layer+1-j;
				connectOut[row*8+3] = connectOut[row*8+2] + 1;

				row++;
			}}
		}

		row = 0;
		for (i = P; i ; i-- ) {
			// Regular TETs
			sum = 0;
			for (j = 1; j <= i; j++)
				sum += j;
			for (jMax = sum; jMax--; ) {
				typesOut[row] = 10;
				row++;
			}

			// PYRs
			for (jMax = i*(i-1); jMax--; ) {
				typesOut[row] = 14;
				row++;
			}

			// Inverted TETs
			if (i > 2) {
				sum = 0;
				for (j = 1; j <= (i-2); j++)
					sum += j;

				for (jMax = sum; jMax--; ) {
					typesOut[row] = 10;
					row++;
				}
			}
		}

		N = P+1;
		for (i = 0; i < Ne; i++) {
			switch (i) {
				case 0:  j0 = 0;             jS = 1;             jC = 1; break;
				case 1:  j0 = 0;             jS = N;             jC = 1; break;
				case 2:  j0 = 0;             jS = N*(N+1)/2.0;   jC = N; break;
				case 3:  j0 = P;             jS = P;             jC = 1; break;
				case 4:  j0 = P;             jS = N*(N+1)/2.0-1; jC = N; break;
				case 5:  j0 = N*(N+1)/2.0-1; jS = P*N/2.0;       jC = P; break;
				default:
					printf("Error: Unsupported.\n"), EXIT_MSG; break;
			}
			for (j = j0, count = 0; count < N; j += jS, count++) {
				connectEOut[i*N+count] = j;

				if (j != j0 && jS > 1) {
					jS -= jC;
					if (jC > 1)
						jC--;
				}
			}
		}

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connectE = connectEOut;
		*Nn       = NnOut;
		*NE       = NEOut;
	} else if (typeIn == WEDGE) {
		d = 3;

		NEOut = 0;
		for (i = 0; i < P; i++)
			NEOut += 2*i+1;
		NEOut *= P;

		double *rst_TRI, *rst_LINE;
		unsigned int Nn_TRI, Nn_LINE, NE_TRI, NE_LINE,
		             *connect_TRI, *connect_LINE, *dummy_types, *dummy_connectE;

		plotting_element_info(&rst_TRI,&connect_TRI,&dummy_types,&dummy_connectE,&Nn_TRI,&NE_TRI,P,TRI); // free
		free(dummy_types); free(dummy_connectE);
		plotting_element_info(&rst_LINE,&connect_LINE,&dummy_types,&dummy_connectE,&Nn_LINE,&NE_LINE,P,LINE); // free
		free(dummy_types); free(dummy_connectE);

		NnOut = Nn_TRI*Nn_LINE;

		rstOut     = malloc(NnOut*d * sizeof *rstOut);     // keep (requires external free)
		connectOut = calloc(NEOut*8 , sizeof *connectOut); // keep (requires external free)
		typesOut   = malloc(NEOut   * sizeof *typesOut);   // keep (requires external free)

		row = 0;
		for (j = 0; j < Nn_LINE; j++) {
		for (i = 0; i < Nn_TRI; i++) {
			for (k = 0; k < 2; k++)
				rstOut[k*NnOut+row] = rst_TRI[k*Nn_TRI+i];

			rstOut[2*NnOut+row] = rst_LINE[j];

			row++;
		}}

		row = 0;
		for (j = 0, jMax = P; j < jMax; j++) {
		for (i = 0, iMax = NEOut/P; i < iMax; i++) {
			for (k = 0; k < 3; k++)
				connectOut[row*8+k] = connect_TRI[i*8+k] + Nn_TRI*j;
			for (k = 0; k < 3; k++)
				connectOut[row*8+3+k] = connect_TRI[i*8+k] + Nn_TRI*(j+1);

			row++;
		}}

		for (i = 0; i < NEOut; i++)
			typesOut[i] = 13;

		N = P+1;
		for (i = 0; i < Ne; i++) {
			switch (i) {
				case 0:  j0 = P;               jS = P;  jC = 1; break;
				case 1:  j0 = 0;               jS = N;  jC = 1; break;
				case 2:  j0 = 0;               jS = 1;  jC = 1; break;
				case 3:  j0 = N*(N+1)/2.0*P+P; jS = P;  jC = 1; break;
				case 4:  j0 = N*(N+1)/2.0*P;   jS = N;  jC = 1; break;
				case 5:  j0 = N*(N+1)/2.0*P;   jS = 1;  jC = 1; break;
				case 6:  j0 = 0;               jS = N*(N+1)/2.0; jC = 0; break;
				case 7:  j0 = 1;               jS = N*(N+1)/2.0; jC = 0; break;
				case 8:  j0 = 2;               jS = N*(N+1)/2.0; jC = 0; break;
				default:
					printf("Error: Unsupported.\n"), EXIT_MSG; break;
			}
			for (j = j0, count = 0; count < N; j += jS, count++) {
				connectEOut[i*N+count] = j;

				if (j != j0 && jS > 1) {
					jS -= jC;
					if (jC > 1)
						jC--;
				}
			}
		}

		free(rst_TRI);
		free(connect_TRI);

		free(rst_LINE);
		free(connect_LINE);

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connectE = connectEOut;
		*Nn       = NnOut;
		*NE       = NEOut;
	} else if (typeIn == PYR) {
		d = 3;
		Nc = 5;

		NnOut = 0;
		for (i = 1, iMax = P+1; i <= iMax; i++)
			NnOut += pow(i,2);

		NEOut = 0;
		for (i = P; i ; i--) {
			NEOut += pow(i,2);
			NEOut += 2*i*(i-1);
			NEOut += pow(i-1,2);
		}

		double *rst_QUAD;
		unsigned int Nn_QUAD, NE_QUAD,
		             *connect_QUAD, *dummy_types, *dummy_connectE;

		rstOut     = malloc(NnOut*d * sizeof *rstOut);     // keep (requires external free)
		connectOut = calloc(NEOut*8 , sizeof *connectOut); // keep (requires external free)
		typesOut   = malloc(NEOut   * sizeof *typesOut);   // keep (requires external free)

		row = 0;
		for (i = P; i ; i--) {
			di = i;
			plotting_element_info(&rst_QUAD,&connect_QUAD,&dummy_types,&dummy_connectE,&Nn_QUAD,&NE_QUAD,i,QUAD); // free
			free(dummy_types); free(dummy_connectE);

			for (j = 0; j < Nn_QUAD; j++) {
				for (k = 0; k < 2; k++)
					rstOut[k*NnOut+row] = rst_QUAD[k*Nn_QUAD+j]*di/P;

				rstOut[2*NnOut+row] = -1.0/5.0*sqrt(2.0) + sqrt(2.0)/2.0*(1.0+(-1.0*di/P + 1.0*(P-di)/P));
				row++;
			}

			free(rst_QUAD);
			free(connect_QUAD);
		}
		for (k = 0; k < 2; k++)
			rstOut[k*NnOut+row] = 0.0;
		rstOut[2*NnOut+row] = 4.0/5.0*sqrt(2.0);

		row = 0;
		for (layer = P; layer; layer--) {
			lBs = 0;
			for (i = P; i > layer; i--)
				lBs += pow(i+1,2);

			lTs = 0;
			for (i = P; i > (layer-1); i--)
				lTs += pow(i+1,2);

			// Regular PYRs
			for (j = 0; j < layer; j++) {
			for (i = 0; i < layer; i++) {
				connectOut[row*8+0] = lBs + i + j*(layer+1);
				connectOut[row*8+1] = connectOut[row*8+0] + 1;
				connectOut[row*8+2] = connectOut[row*8+1] + layer+1;
				connectOut[row*8+3] = connectOut[row*8+2] - 1;
				connectOut[row*8+4] = lTs + i + j*layer;
				row++;
			}}

			// TETs
			for (j = 0; j < layer; j++) {
			for (i = 1; i < layer; i++) {
				connectOut[row*8+0] = lBs + i + j*(layer+1);
				connectOut[row*8+1] = connectOut[row*8+0] + layer+1;
				connectOut[row*8+2] = lTs + i-1 + j*layer;
				connectOut[row*8+3] = connectOut[row*8+2] + 1;
				row++;
			}}

			for (i = 0; i < layer; i++) {
			for (j = 1; j < layer; j++) {
				connectOut[row*8+0] = lBs + i + j*(layer+1);
				connectOut[row*8+1] = connectOut[row*8+0] + 1;
				connectOut[row*8+2] = lTs + i + (j-1)*layer;
				connectOut[row*8+3] = connectOut[row*8+2] + layer;
				row++;
			}}

			// Inverted PYRs
			for (j = 1; j < layer; j++) {
			for (i = 1; i < layer; i++) {
				connectOut[row*8+0] = lTs + i-1 + (j-1)*layer;
				connectOut[row*8+1] = connectOut[row*8+0] + 1;
				connectOut[row*8+2] = connectOut[row*8+1] + layer;
				connectOut[row*8+3] = connectOut[row*8+2] - 1;
				connectOut[row*8+4] = lBs + i + j*(layer+1);
				row++;
			}}
//array_print_ui(NEOut,8,connectOut,'R');
		}

		row = 0;
		for (i = P; i ; i--) {
			// Regular PYRs
			for (j = 0, jMax = pow(i,2); j < jMax; j++) {
				typesOut[row] = 14;
				row++;
			}

			// TETs
			for (j = 0, jMax = 2*i*(i-1); j < jMax; j++) {
				typesOut[row] = 10;
				row++;
			}

			// Inverted PYRs
			for (j = 0, jMax = pow(i-1,2); j < jMax; j++) {
				typesOut[row] = 14;
				row++;
			}
		}

		N = P+1;
		for (i = 0; i < Ne; i++) {
			switch (i) {
				case 0:  j0 = 0;     jS = 1;     jC = 0; break;
				case 1:  j0 = P*N;   jS = 1;     jC = 0; break;
				case 2:  j0 = 0;     jS = N;     jC = 0; break;
				case 3:  j0 = P;     jS = N;     jC = 0; break;
				case 4:  j0 = 0;     jS = N*N;   jC = 2*N-1; break;
				case 5:  j0 = P;     jS = N*N-1; jC = 2*N-1; break;
				case 6:  j0 = P*N;   jS = P*P+1; jC = 2*P-1; break;
				case 7:  j0 = N*N-1; jS = P*P;   jC = 2*P-1; break;
				default:
					printf("Error: Unsupported.\n"), EXIT_MSG; break;
			}
			for (j = j0, count = 0; count < N; j += jS, count++) {
				connectEOut[i*N+count] = j;

				if (j != j0 && jS > 1) {
					jS -= jC;
					if (jC > 1)
						jC -= 2;
				}
			}
		}

		*rst      = rstOut;
		*connect  = connectOut;
		*types    = typesOut;
		*connectE = connectEOut;
		*Nn       = NnOut;
		*NE       = NEOut;
	}

//array_print_d(NnOut,d,*rst,'C');
//array_print_ui(NEOut,8,*connect,'R');
//array_print_ui(NEOut,1,*types,'R');

}
