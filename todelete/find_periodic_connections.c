// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>

#include "functions.h"

#include "petscsys.h"

/*
 *	Purpose:
 *		Given a list of periodic vertex correspondence, find all possible matches between the vertices.
 *
 *	Comments:
 *		Each row of the list is sorted in ascending order. Hence, reversed entries are redundant and not included.
 *
 *	Notation:
 *		PVe : List of periodic vertices (pve x 2 array)
 *		pve : Index of last row in PVe
 *
 *	References:
 *
 */

void find_periodic_connections(unsigned int *PVe, unsigned int *NPVePointer, const unsigned int VeMax)
{
	unsigned int i, j, k,
	             pve, *IndicesDummy, Modified, *PVe1D, NUnique, *PVeUnique, *PVeUniqueOver, *PVeMatches, *IndPVeMatches,
	             IndRow, LenF, match, row, col, n2, *PVePotential, IndPVe;

	pve = *NPVePointer;

// array_print_i(pve,2,PVe,'R');

	// Sort Rows
	for (i = 0; i < pve; i++)
		PetscSortInt(2,(int *)&PVe[i*2+0]);

	// Sort Columns
	IndicesDummy = malloc(pve*2 * sizeof *IndicesDummy); // free
	array_sort_ui(pve,2,PVe,IndicesDummy,'R','T');
	free(IndicesDummy);

// array_print_i(pve,2,PVe,'R');

	Modified = 1;
	while (Modified) {
		Modified = 0;

		PVe1D = malloc(pve*2 * sizeof *PVe1D); // free
		for (i = k = 0; i < pve; i++) {
		for (j = 0; j < 2; j++) {
			PVe1D[k] = PVe[i*2+j];
			k++;
		}}

		PetscSortInt(k,(int *)PVe1D);
// array_print_i(1,k,PVe1D,'R');

		PVeUniqueOver = malloc(pve*2 * sizeof *PVeUniqueOver); // free
		PVeUniqueOver[0] = PVe1D[0];
		for (i = 1, NUnique = 1; i < pve*2; i++) {
			if (PVe1D[i] != PVeUniqueOver[NUnique-1]) {
				PVeUniqueOver[NUnique] = PVe1D[i];
				NUnique++;
			}
		}
		free(PVe1D);

		PVeUnique = malloc(NUnique * sizeof *PVeUnique); // free
		for (i = 0; i < NUnique; i++)
			PVeUnique[i] = PVeUniqueOver[i];
		free(PVeUniqueOver);

// array_print_i(1,NUnique,PVeUnique,'R');

		PVeMatches    = malloc(NUnique*8 * sizeof *PVeMatches); // free
		IndPVeMatches = malloc(NUnique   * sizeof *IndPVeMatches); // free
		for (i = 0; i < NUnique*8; i++) PVeMatches[i]    = VeMax;
		for (i = 0; i < NUnique; i++)   IndPVeMatches[i] = 0;

		for (k = 0; k < 2; k++) {
		for (i = 0; i < pve; i++) {
			array_find_indexo_ui(NUnique,PVeUnique,PVe[i*2+k],&IndRow,&LenF);

			for (j = 0, match = 0; j < IndPVeMatches[IndRow]; j++) {
				if (PVe[i*2+1-k] == PVeMatches[IndRow*8+j])
					match = 1;
			}

			if (!match) {
				PVeMatches[IndRow*8+IndPVeMatches[IndRow]] = PVe[i*2+1-k];

				IndPVeMatches[IndRow]++;
				if (IndPVeMatches[IndRow] > 7)
					printf("Error: Too many periodic connections\n"), exit(1);
			}
		}}

		// sort rows
		for (i = 0; i < NUnique; i++)
			PetscSortInt(IndPVeMatches[i]+1,(int *)&PVeMatches[i*8+0]);

// array_print_i(NUnique,8,PVeMatches,'R');
// array_print_i(1,NUnique,IndPVeMatches,'R');

		PVePotential = malloc(2 * sizeof *PVePotential); // free

		for (row = 0; row < NUnique; row++) {
		for (col = 0; col < IndPVeMatches[row]; col++) {
			array_find_indexo_ui(NUnique,PVeUnique,PVeMatches[row*8+col],&IndRow,&LenF);
			for (i = 0; i < IndPVeMatches[IndRow]; i++) {
				n2 = PVeMatches[IndRow*8+i];
				if (PVeUnique[row] != n2) {
					PVePotential[0] = PVeUnique[row];
					PVePotential[1] = n2;
					PetscSortInt(2,(int *)PVePotential);

					for (IndPVe = 0, match = 0; IndPVe < pve; IndPVe++) {
						if (PVe[IndPVe*2+0] == PVePotential[0] &&
							PVe[IndPVe*2+1] == PVePotential[1]) {
								match = 1;
								break;
						}
					}

					if (!match) {
						for (j = 0; j < 2; j++)
							PVe[pve*2+j] = PVePotential[j];

						pve++;
						Modified = 1;
					}
				}
			}
		}}
		free(PVeUnique);
		free(PVePotential);

		IndicesDummy = malloc(pve*2 * sizeof *IndicesDummy); // free
		array_sort_ui(pve,2,PVe,IndicesDummy,'R','T');
		free(IndicesDummy);

// array_print_i(pve,2,PVe,'R');

		free(PVeMatches);
		free(IndPVeMatches);
	}

	*NPVePointer = pve;
}