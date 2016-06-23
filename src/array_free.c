// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>

#include "database.h"

/*
 *	Purpose:
 *		Free dynamically allocated arrays with more than one level of pointer abstraction.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 *
 */

void array_free2_c(unsigned int iMax, char **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_ui(unsigned int iMax, unsigned int **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_i(unsigned int iMax, int **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_l(unsigned int iMax, long **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_ll(unsigned int iMax, long long **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_f(unsigned int iMax, float **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_d(unsigned int iMax, double **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}

void array_free2_ld(unsigned int iMax, long double **A)
{
	unsigned int i;

	for (i = 0; i < iMax; i++)
		if (A[i] != NULL)
			free(A[i]);
	free(A);
}



void array_free3_c(unsigned int iMax, unsigned int jMax, char ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_ui(unsigned int iMax, unsigned int jMax, unsigned int ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_i(unsigned int iMax, unsigned int jMax, int ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_l(unsigned int iMax, unsigned int jMax, long ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_ll(unsigned int iMax, unsigned int jMax, long long ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_f(unsigned int iMax, unsigned int jMax, float ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_d(unsigned int iMax, unsigned int jMax, double ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free3_ld(unsigned int iMax, unsigned int jMax, long double ***A)
{
	unsigned int i, j;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++)
				if (A[i][j] != NULL)
					free(A[i][j]);
			free(A[i]);
		}
	}
	free(A);
}

void array_free4_d(unsigned int iMax, unsigned int jMax, unsigned int kMax, double ****A)
{
	unsigned int i, j, k;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++) {
				if (A[i][j] != NULL) {
					for (k = 0; k < kMax; k++)
						if (A[i][j][k] != NULL)
							free(A[i][j][k]);
					free(A[i][j]);
				}
			}
			free(A[i]);
		}
	}
	free(A);
}

void array_free5_d(unsigned int iMax, unsigned int jMax, unsigned int kMax, unsigned int lMax, double *****A)
{
	unsigned int i, j, k, l;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++) {
				if (A[i][j] != NULL) {
					for (k = 0; k < kMax; k++) {
						if (A[i][j][k] != NULL) {
							for (l = 0; l < lMax; l++)
								if (A[i][j][k][l] != NULL)
									free(A[i][j][k][l]);
							free(A[i][j][k]);
						}
					}
					free(A[i][j]);
				}
			}
			free(A[i]);
		}
	}
	free(A);
}

void array_free1_CSR_d(struct S_OpCSR *A)
{
	if (A != NULL) {
		free(A->rowIndex);
		free(A->columns);
		free(A->values);
		free(A);
	}
}

void array_free5_CSR_d(unsigned int iMax, unsigned int jMax, unsigned int kMax, unsigned int lMax, struct S_OpCSR *****A)
{
	unsigned int i, j, k, l;

	for (i = 0; i < iMax; i++) {
		if (A[i] != NULL) {
			for (j = 0; j < jMax; j++) {
				if (A[i][j] != NULL) {
					for (k = 0; k < kMax; k++) {
						if (A[i][j][k] != NULL) {
							for (l = 0; l < lMax; l++)
								if (A[i][j][k][l] != NULL) {
									free(A[i][j][k][l]->rowIndex);
									free(A[i][j][k][l]->columns);
									free(A[i][j][k][l]->values);
									free(A[i][j][k][l]);
								}
							free(A[i][j][k]);
						}
					}
					free(A[i][j]);
				}
			}
			free(A[i]);
		}
	}
	free(A);
}
