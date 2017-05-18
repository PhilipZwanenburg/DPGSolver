// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#include "test_integration_Advection.h"

#include <stdlib.h>
#include <stdio.h>

#include "Parameters.h"

#include "test_code_integration_conv_order.h"
#include "test_code_integration_linearization.h"
#include "test_support.h"

#include "array_free.h"

/*
 *	Purpose:
 *		Test various aspects of the Advection solver implementation:
 *			- Linearization;
 *			- Optimal convergence orders.
 *
 *	Comments:
 *
 *	Notation:
 *
 *	References:
 */

void test_integration_Advection(int nargc, char **argv)
{
	bool const RunTests_linearization = 0,
	           RunTests_conv_order    = 1;

	char **argvNew, *PrintName;

	argvNew    = malloc(2          * sizeof *argvNew);  // free
	argvNew[0] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	argvNew[1] = malloc(STRLEN_MAX * sizeof **argvNew); // free
	PrintName  = malloc(STRLEN_MAX * sizeof *PrintName); // free

	// silence
	strcpy(argvNew[0],argv[0]);

	// **************************************************************************************************** //
	// Linearization Testing
	// **************************************************************************************************** //
	if (RunTests_linearization) {
		struct S_linearization *data_l = calloc(1 , sizeof *data_l); // free

		data_l->nargc     = nargc;
		data_l->argvNew   = argvNew;
		data_l->PrintName = PrintName;

		// 2D (Mixed TRI/QUAD mesh)
		test_linearization(data_l,"Advection_StraightTRI");
		test_linearization(data_l,"Advection_CurvedTRI");
		test_linearization(data_l,"Advection_ToBeCurvedTRI");
		test_linearization(data_l,"Advection_MIXED2D");

		test_print_warning("Advection 3D testing needs to be implemented");

		free(data_l);
	} else {
		test_print_warning("Advection linearization testing currently disabled");
	}

	// **************************************************************************************************** //
	// Convergence Order Testing
	// **************************************************************************************************** //
	if (RunTests_conv_order) {
		struct S_convorder *data_c = calloc(1 , sizeof *data_c); // free

		data_c->nargc     = nargc;
		data_c->argvNew   = argvNew;
		data_c->PrintName = PrintName;

		test_conv_order(data_c,"Advection_n-Cube_Default_TRI");
		test_conv_order(data_c,"Advection_n-Cube_Peterson_TRI");

		free(data_c);
	} else {
		test_print_warning("Advection convergence order testing currently disabled");
	}

	array_free2_c(2,argvNew);
	free(PrintName);
}