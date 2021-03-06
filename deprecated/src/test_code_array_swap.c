// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_code_array_swap.h"

#include <stdlib.h>
#include <stdio.h>

/*
 *	Purpose:
 *		Provide functions for testing array_swap.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void array_swap1_ui(register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn,
                    const unsigned int stepIn)
{
	register unsigned int i, N, step;
	register unsigned int tmp;

	for (i = 0, N = NIn, step = stepIn; N-- ; i += step) {
		tmp     = arr1[i];
		arr1[i] = arr2[i];
		arr2[i] = tmp;
	}
}

void array_swap2_ui(register unsigned int *arr1, register unsigned int *arr2, const unsigned int NIn,
                    const unsigned int stepIn)
{
	register unsigned int N, step;
	register unsigned int tmp;

	for (N = NIn, step = stepIn; N-- ; ) {
		tmp   = *arr1;
		*arr1 = *arr2;
		*arr2 = tmp;

		arr1 += step;
		arr2 += step;
	}
}

void array_swap1_d(register double *arr1, register double *arr2, const unsigned int NIn, const unsigned int stepIn)
{
	register unsigned int i, N, step;
	register double tmp;

	for (i = 0, N = NIn, step = stepIn; N-- ; i += step) {
		tmp     = arr1[i];
		arr1[i] = arr2[i];
		arr2[i] = tmp;
	}
}

void array_swap2_d(register double *arr1, register double *arr2, const unsigned int NIn, const unsigned int stepIn)
{
	register unsigned int N, step;
	register double tmp;

	for (N = NIn, step = stepIn; N-- ; ) {
		tmp   = *arr1;
		*arr1 = *arr2;
		*arr2 = tmp;

		arr1 += step;
		arr2 += step;
	}
}
