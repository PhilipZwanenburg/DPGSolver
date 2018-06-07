/* {{{
This file is part of DPGSolver.

DPGSolver is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or any later version.

DPGSolver is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
Public License for more details.

You should have received a copy of the GNU General Public License along with DPGSolver.  If not, see
<http://www.gnu.org/licenses/>.
}}} */
/** \file
 */


#include "optimize.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>

#include "face.h"
#include "volume.h"

#include "test_integration.h"

#include "definitions_adaptation.h"
#include "definitions_core.h"
#include "definitions_tol.h"

#include "macros.h"

#include "matrix.h"
#include "matrix_constructors.h"

#include "multiarray.h"
#include "complex_multiarray.h"
#include "multiarray_constructors.h"
#include "complex_multiarray_constructors.h"
#include "multiarray_math.h"

#include "intrusive.h"
#include "definitions_intrusive.h"

#include "geometry.h"
#include "geometry_parametric.h"
#include "solve_implicit.h"

#include "simulation.h"

#include "volume_solver.h"
#include "face_solver.h"
#include "test_complex_face_solver.h"
#include "test_complex_volume_solver.h"
#include "visualization.h"
#include "definitions_visualization.h"

#include "optimization_case.h"

#include "adjoint.h"
#include "sensitivities.h"
#include "optimization_minimizers.h"

#include "objective_functions.h"


// TODO: Read from the optimization.data file
#define CONST_L2_GRAD_EXIT 1E-10
#define CONST_OBJECTIVE_FUNC_EXIT 1E-10
#define MAX_NUM_DESIGN_ITERS 250

// Static function declarations ************************************************************************************* //


static struct Optimization_Case* setup_optimization(struct Simulation* sim);

static void copy_data_r_to_c_sim(struct Simulation* sim, struct Simulation* sim_c);

static struct Multiarray_d* compute_gradient(struct Optimization_Case* optimization_case);

static struct Multiarray_d* test_compute_brute_force_gradient(struct Optimization_Case* optimization_case);

static void optimize_line_search_method(struct Optimization_Case* optimization_case);


// NLPQLP Specific Static Functions
static void optimize_NLPQLP(struct Optimization_Case *optimization_case);
static void preprocess_optimize_NLPQLP(struct Optimization_Case *optimization_case);
static void postprocess_optimize_NLPQLP(struct Optimization_Case *optimization_case);
static void compute_function_values_NLPQLP(struct Optimization_Case *optimization_case);
static void compute_gradient_values_NLPQLP(struct Optimization_Case *optimization_case);
static void write_NLPQLP_input_file(struct Optimization_Case *optimization_case);
//static void process_NLPQLP_output_file(struct Optimization_Case *optimization_case);

// Interface functions ********************************************************************************************** //



void optimize(struct Simulation* sim){

	/*
	Optimize the the geometry to minimize a specified functional.

	Arguments:
		sim = The simulation data structure

	Return:
		- 
	*/

	// Preprocessing
	struct Optimization_Case *optimization_case = setup_optimization(sim);
	copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);

	// Optimization routine
	//optimize_line_search_method(optimization_case);
	optimize_NLPQLP(optimization_case);


	// Post Processing
	output_visualization(sim,VIS_GEOM_EDGES);
	output_visualization(sim,VIS_GEOM_VOLUMES);
	output_visualization(sim,VIS_NORMALS);
	output_visualization(sim,VIS_SOLUTION);


	// Clear allocated structures:
	destructor_Optimization_Case(optimization_case);

	
	// Unused functions:
	UNUSED(test_compute_brute_force_gradient);
	UNUSED(optimize_line_search_method);
	UNUSED(optimize_NLPQLP);

}


// Static functions ************************************************************************************************* //


static void optimize_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	Interface with NLPQLP and use it to perform the shape optimization.

	Arguments:
		optimization_case = The data structure with the optimization data
	*/

	struct Simulation *sim = optimization_case->sim;

	preprocess_optimize_NLPQLP(optimization_case);

	// Optimization routine
	int design_iteration = 0;
	
	while (true){

		printf("Design Iteration: %d \n", design_iteration);



	}

	write_NLPQLP_input_file(optimization_case);
	exit(0);


	postprocess_optimize_NLPQLP(optimization_case);


}

static void preprocess_optimize_NLPQLP(struct Optimization_Case *optimization_case){
	
	/*
	Initialize the NLPQLP data and allocate the data structures needed for the NLPQLP
	optimization. Place all data structures in optimization_case.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	// ========================================
	//       Optimization Parameters
	// ========================================

	// Setup the optimization parameters
	int NP = 1;  // Number of processors

	// Number of design variable dofs
	int N = optimization_case->num_design_pts_dofs;
	int NMAX = N + 1;
	
	// Constrained Optimization Problem with
	// constraint that consecutive design variables must be within delta_y 
	// limit
	
	// NOTE: This only handles the case where all design points move only in the 
	//		y direction (1 dof)

	// Consider an unconstrained optimization problem first

	//int CONST_M = CONST_N-2;
	int M = 0;
	int ME = 0;


	// Independent Optimization Parameters:
	optimization_case->NLPQLP_data.NP 		= NP;
	optimization_case->NLPQLP_data.N 		= N;
	optimization_case->NLPQLP_data.NMAX 	= NMAX;
	optimization_case->NLPQLP_data.M 		= M;
	optimization_case->NLPQLP_data.ME 		= ME;



	// Dependent Optimization Parameters:
	optimization_case->NLPQLP_data.IFAIL 	= 0;  // Initial value is 0
	optimization_case->NLPQLP_data.MODE 	= 0;
	optimization_case->NLPQLP_data.MNN2 	= M + N + N + 2;



	// Remaining Optimization Parameters 
	// Set these up using default values (can be overwritten)
	optimization_case->NLPQLP_data.IOUT 	= 6;
	optimization_case->NLPQLP_data.MAXIT 	= 100;
	optimization_case->NLPQLP_data.MAXFUN 	= 10;
	optimization_case->NLPQLP_data.MAXNM 	= 10;
	optimization_case->NLPQLP_data.LQL 		= 1;
	optimization_case->NLPQLP_data.IPRINT 	= 0;

	optimization_case->NLPQLP_data.ACC 		= 1E-8;
	optimization_case->NLPQLP_data.ACCQP 	= 1E-14;
	optimization_case->NLPQLP_data.STPMIN 	= 1E-10;
	optimization_case->NLPQLP_data.RHO 		= 0.0;



	// Memory Objects (initialize with zeros)
	optimization_case->NLPQLP_data.X 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.X, 0.0);	

	optimization_case->NLPQLP_data.XL 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.XL, 0.0);	

	optimization_case->NLPQLP_data.XU 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.XU, 0.0);	



	optimization_case->NLPQLP_data.F = 0;

	optimization_case->NLPQLP_data.G 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, M});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.G, 0.0);	

	optimization_case->NLPQLP_data.dF 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.dF, 0.0);		

	optimization_case->NLPQLP_data.dG 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){M, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.dG, 0.0);	



	optimization_case->NLPQLP_data.C 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){NMAX, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.C, 0.0);	

	optimization_case->NLPQLP_data.U 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, optimization_case->NLPQLP_data.MNN2});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.U, 0.0);	

	optimization_case->NLPQLP_data.D 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, NMAX});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.D, 0.0);	



	optimization_case->NLPQLP_data.LWA = 3*(optimization_case->NLPQLP_data.N+1)*(optimization_case->NLPQLP_data.N+1)/2 + 
			33*(optimization_case->NLPQLP_data.N+1) + 9*optimization_case->NLPQLP_data.M + 150;
	optimization_case->NLPQLP_data.LWA = optimization_case->NLPQLP_data.LWA*10;	
	
	optimization_case->NLPQLP_data.WA 	= constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, optimization_case->NLPQLP_data.LWA});
	set_to_value_Multiarray_d(optimization_case->NLPQLP_data.WA, 0.0);	



	optimization_case->NLPQLP_data.LKWA = N + 25;
	optimization_case->NLPQLP_data.LKWA = optimization_case->NLPQLP_data.LKWA*10;	
	
	optimization_case->NLPQLP_data.KWA 	= constructor_empty_Multiarray_i('C',2,(ptrdiff_t[]){1, optimization_case->NLPQLP_data.LKWA});
	set_to_value_Multiarray_i(optimization_case->NLPQLP_data.KWA, 0.0);	



	optimization_case->NLPQLP_data.LACTIV = 2*optimization_case->NLPQLP_data.M + 10;
	optimization_case->NLPQLP_data.ACTIVE   = constructor_empty_Multiarray_i('C',2,(ptrdiff_t[]){1, optimization_case->NLPQLP_data.LACTIV});
	set_to_value_Multiarray_i(optimization_case->NLPQLP_data.ACTIVE, 0.0);	


	// ========================================
	//         Load Optimization Data
	// ========================================

	// Load the design point data
	struct Multiarray_d* ctrl_pts_and_w = optimization_case->geo_data.control_points_and_weights;
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	struct Multiarray_d* ctrl_pts_lims = optimization_case->geo_data.control_points_optimization_lims;
	int n_pts = (int)ctrl_pts_opt->extents[0];

	int ctrl_pt_index;
	int vec_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			ctrl_pt_index = get_col_Multiarray_i(0, ctrl_pts_opt)[i];

			optimization_case->NLPQLP_data.X->data[vec_index] = get_col_Multiarray_d(j-1, ctrl_pts_and_w)[ctrl_pt_index];
			optimization_case->NLPQLP_data.XL->data[vec_index] = get_col_Multiarray_d(0, ctrl_pts_lims)[vec_index];
			optimization_case->NLPQLP_data.XU->data[vec_index] = get_col_Multiarray_d(1, ctrl_pts_lims)[vec_index];

			vec_index++;
		}
	}

	// Set the initial Hessian approximation to be a large multiple of the 
	// identity matrix to keep the step length small
	for (int i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		get_col_Multiarray_d(i, optimization_case->NLPQLP_data.C)[i] = 500.;
	}
	optimization_case->NLPQLP_data.MODE = 1;


	// Set the initial function values and gradient values
	compute_function_values_NLPQLP(optimization_case);
	compute_gradient_values_NLPQLP(optimization_case);

}


static void postprocess_optimize_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	Postprocess the NLPQLP data and deallocate the data structures created for the NLPQLP
	optimization.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	destructor_Multiarray_d(optimization_case->NLPQLP_data.X);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.XL);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.XU);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.G);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.U);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.D);
	destructor_Multiarray_d(optimization_case->NLPQLP_data.WA);

	destructor_Multiarray_i(optimization_case->NLPQLP_data.KWA);
	destructor_Multiarray_i(optimization_case->NLPQLP_data.ACTIVE);

}

static void compute_function_values_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	Compute the objective function value and constraint function values. Load the data
	into optimization_case NLPQLP data structure.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	// Objective Function:
	optimization_case->NLPQLP_data.F = optimization_case->objective_function(optimization_case->sim);

	// Constraint Functions:

	// TODO: Implement constraint functions to be added to the optimization

}


static void compute_gradient_values_NLPQLP(struct Optimization_Case *optimization_case){

	/*
	Compute the gradient of the objective function and constraint functions. Load the data
	into optimization_case NLPQLP data structure.

	Arguments:
		optimization_case = The data structure with the optimization data

	Return:
		-
	*/

	// =================================
	//    Objective Function Gradient
	// =================================

	setup_adjoint(optimization_case);
	solve_adjoint(optimization_case);


	// Compute dI_dXp and dR_dXp
	compute_sensitivities(optimization_case);

	// Use the adjoint and sensitivities to find the gradient of the objective function
	struct Multiarray_d* grad_I = compute_gradient(optimization_case);  // free

	// Load the gradient data into the NLPQLP data structure
	for (int i = 0; i < optimization_case->num_design_pts_dofs; i++)
		optimization_case->NLPQLP_data.dF->data[i] = grad_I->data[i];

	destructor_Multiarray_d(grad_I);

	// =================================
	//    Constraint Function Gradient
	// =================================

	// TODO: Implement constraint functions for the optimization

}


static void write_NLPQLP_input_file(struct Optimization_Case *optimization_case){

	/*
	Write the file that will be read by the NLPQLP program.

	NOTE: Uses absolute paths to the input and output files here 
		in the optimization directory
	TODO: Read the path to the NLPQLP directory from an optimization.data file

	Arguments:
		optimization_case = The data structure with the optimization information

	Return:
		-
	*/

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../../NLPQLP_Optimizer/");
	strcat(output_name,"INPUT.txt");

	FILE *fp;

	char lql_string[100];
	int i, j;

	if ((fp = fopen(output_name,"w")) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Optimizer Properties
	fprintf(fp, "NP \t N \t NMAX \t M \t ME \t IFAIL \t MODE \n");
	fprintf(fp, "%d \t %d \t %d \t %d \t %d \t %d \t %d \n", 
		optimization_case->NLPQLP_data.NP, optimization_case->NLPQLP_data.N, 
		optimization_case->NLPQLP_data.NMAX, optimization_case->NLPQLP_data.M,
		optimization_case->NLPQLP_data.ME, optimization_case->NLPQLP_data.IFAIL, 
		optimization_case->NLPQLP_data.MODE);
	fprintf(fp, "\n");

	fprintf(fp, "LWA \t LKWA \t LACTIV \t IOUT \t ACC \t ACCQP \n");
	fprintf(fp, "%d \t %d \t %d \t %d \t %e \t %e \n", 
		optimization_case->NLPQLP_data.LWA, optimization_case->NLPQLP_data.LKWA, 
		optimization_case->NLPQLP_data.LACTIV, optimization_case->NLPQLP_data.IOUT, 
		optimization_case->NLPQLP_data.ACC, optimization_case->NLPQLP_data.ACCQP);
	fprintf(fp, "\n");

	if(optimization_case->NLPQLP_data.LQL){
		strcpy(lql_string, "T");
	} else{
		strcpy(lql_string, "F");
	}

	fprintf(fp, "STPMIN \t MAXIT \t MAXFUN \t MAXNM \t RHO \t LQL \t IPRINT \n");
	fprintf(fp, "%e \t %d \t %d \t %d \t %e \t %s \t %d \n", 
		optimization_case->NLPQLP_data.STPMIN, optimization_case->NLPQLP_data.MAXIT, 
		optimization_case->NLPQLP_data.MAXFUN, optimization_case->NLPQLP_data.MAXNM,
		optimization_case->NLPQLP_data.RHO, lql_string, optimization_case->NLPQLP_data.IPRINT);
	fprintf(fp, "\n");

	// Design Variable Values
	fprintf(fp, "X \t XL \t XU \n");
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		fprintf(fp, "%.14e \t %.14e \t %.14e \n", 
			optimization_case->NLPQLP_data.X->data[i], 
			optimization_case->NLPQLP_data.XL->data[i], 
			optimization_case->NLPQLP_data.XU->data[i]);
	}
	fprintf(fp, "\n");

	// Objective Function Evaluation
	fprintf(fp, "F \n");
	fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.F);
	fprintf(fp, "\n");

	// Constraint Function Evaluations
	fprintf(fp, "G \n");
	for (i = 0; i < optimization_case->NLPQLP_data.M; i++)
		fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.G->data[i]);
	fprintf(fp, "\n");

	// Gradient Objective Function
	fprintf(fp, "dF \n");
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++)
		fprintf(fp, "%.14e ", optimization_case->NLPQLP_data.dF->data[i]);
	fprintf(fp, "\n \n");

	// Gradient Constraint Functions (Column Major Ordering)
	fprintf(fp, "dG \n");
	for (i = 0; i < optimization_case->NLPQLP_data.M; i++){
		for (j = 0; j < optimization_case->NLPQLP_data.NMAX; j++){

			fprintf(fp, "%.14e ", get_col_Multiarray_d(j, optimization_case->NLPQLP_data.dG)[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	// Hessian Approximation (C)
	fprintf(fp, "C \n");
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		for (j = 0; j < optimization_case->NLPQLP_data.NMAX; j++){
			fprintf(fp, "%.14e ", get_col_Multiarray_d(j, optimization_case->NLPQLP_data.C)[i]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n");

	// Multipliers (U)
	fprintf(fp, "U \n");
	for (i = 0; i < optimization_case->NLPQLP_data.MNN2; i++){
		fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.U->data[i]);
	}
	fprintf(fp, "\n");

	// D Structure
	fprintf(fp, "D \n");
	for (i = 0; i < optimization_case->NLPQLP_data.NMAX; i++){
		fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.D->data[i]);
	}
	fprintf(fp, "\n");

	// WA Structure
	fprintf(fp, "WA \n");
	for (i = 0; i < optimization_case->NLPQLP_data.LWA; i++){
		fprintf(fp, "%.14e \n", optimization_case->NLPQLP_data.WA->data[i]);
	}
	fprintf(fp, "\n");

	// KWA Structure
	fprintf(fp, "KWA \n");
	for (i = 0; i < optimization_case->NLPQLP_data.LKWA; i++){
		fprintf(fp, "%d \n", optimization_case->NLPQLP_data.KWA->data[i]);
	}
	fprintf(fp, "\n");

	// Active
	fprintf(fp, "ACTIVE \n");
	for (i = 0; i < optimization_case->NLPQLP_data.LACTIV; i++){
		if (optimization_case->NLPQLP_data.ACTIVE->data[i])
			fprintf(fp, "T\n");
		else
			fprintf(fp, "F\n");
	}
	fprintf(fp, "\n");

	fclose(fp);

}


static void optimize_line_search_method(struct Optimization_Case* optimization_case){

	/*
	Use a line search method (gradient descent or BFGS) to find the optimal shape

	Arguments:
		optimization_case = The data structure holding all optimization data

	Return:
		-
	*/

	struct Simulation *sim = optimization_case->sim;

	preprocessor_minimizer(optimization_case);

	int design_iteration = 0;
	double L2_grad, objective_func_value;

	objective_func_value = optimization_case->objective_function(sim);
	printf(" INITIAL -> obj_func : %e  , Cl: %f \n", objective_func_value, compute_Cl(sim));

	/*
	output_visualization(sim,VIS_GEOM_EDGES);
	output_visualization(sim,VIS_GEOM_VOLUMES);
	output_visualization(sim,VIS_NORMALS);
	output_visualization(sim,VIS_SOLUTION);
	*/

	// ================================
	//      Optimization Routine
	// ================================

	while(true){

		printf("------------------------------------\n");
		printf("Design Iteration : %d \n", design_iteration);fflush(stdout);

		// Setup and solve for the adjoint (Chi)
		setup_adjoint(optimization_case);
		solve_adjoint(optimization_case);


		// Compute dI_dXp and dR_dXp
		compute_sensitivities(optimization_case);
		printf("Computed Sensitivities \n");fflush(stdout);


		// Use the adjoint and sensitivities to find the gradient of the objective function
		optimization_case->grad_I = compute_gradient(optimization_case);  // free


		// Find the search direction, step length. Then modify the geometry and solve the flow
		//gradient_descent(optimization_case, design_iteration);
		BFGS_minimizer(optimization_case, design_iteration);


		// Copy the new data from the real to the complex structures (for new complex)
		copy_data_r_to_c_sim(optimization_case->sim, optimization_case->sim_c);


		// Compute optimization values to keep track and for the exit condition
		L2_grad = norm_Multiarray_d(optimization_case->grad_I, "L2");
		objective_func_value = optimization_case->objective_function(sim);
		printf(" L2_grad : %e   obj_func : %e  ", L2_grad, objective_func_value);
		printf(" Cl : %f ", compute_Cl(sim));  // temporary monitoring
		printf("\n");


		// Destruct allocated data structures:
		destruct_sensitivity_structures(optimization_case);
		destructor_Multiarray_d(optimization_case->grad_I);


		// Exit condition
		if (L2_grad < CONST_L2_GRAD_EXIT || 
			objective_func_value < CONST_OBJECTIVE_FUNC_EXIT ||
			design_iteration >= MAX_NUM_DESIGN_ITERS)
			break;

		design_iteration++;

	}

	// ================================
	//         Postprocessing
	// ================================

	postprocessor_minimizer(optimization_case);
}


static void copy_data_r_to_c_sim(struct Simulation* sim, struct Simulation* sim_c){

	/*
	Copy the data from the real sim object to the complex sim object. This method will 
	transfer all the face solver and volume solver data into the complex sim object.

	NOTE: Since the mesh and same control file is read (no adaptation was done) we can 
		assume that the ordering of the volumes and faces in the sim and sim_c structure
		is the same.
	
	Arguments:
		sim = The sim data structure with real Volumes and Faces
		sim_c = The sim data structure with Volumes and Faces that hold complex data.

	Return:
		-
	*/

	// Copy the Volume_Solver data
	struct Intrusive_Link* curr   = sim->volumes->first;
	struct Intrusive_Link* curr_c = sim_c->volumes->first;

	while(true){

		struct Solver_Volume* s_vol 		= (struct Solver_Volume*) curr;
		struct Solver_Volume_c* s_vol_c 	= (struct Solver_Volume_c*) curr_c;
		
		copy_members_r_to_c_Solver_Volume((struct Solver_Volume_c*const)s_vol_c, 
			(const struct Solver_Volume*const)s_vol, sim);

		curr 	= curr->next;
		curr_c 	= curr_c->next;

		if(!curr || !curr_c)
			break;
	}

	// Copy the Face_Solver data
	curr   = sim->faces->first;
	curr_c = sim_c->faces->first;

	while(true){

		struct Solver_Face* s_face 		= (struct Solver_Face*) curr;
		struct Solver_Face_c* s_face_c 	= (struct Solver_Face_c*) curr_c;
		
		// TODO:
		// Needed in order to do the r_to_c conversion. Find a way to fix this issue
		s_face->nf_fc = (const struct const_Multiarray_d*)constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1,1});
		s_face_c->nf_fc = (const struct const_Multiarray_c*)constructor_empty_Multiarray_c('C',2,(ptrdiff_t[]){1,1});

		copy_members_r_to_c_Solver_Face((struct Solver_Face_c*const)s_face_c, 
			(const struct Solver_Face*const)s_face, sim);

		curr 	= curr->next;
		curr_c 	= curr_c->next;

		if(!curr || !curr_c)
			break;
	}

}


static struct Optimization_Case* setup_optimization(struct Simulation* sim){

	/*
	Setup the optimization case. In this method, first the constructor for the 
	Optimization_Case data structure will be called. Following this, a complex 
	sim object will be created. This object will be used for obtaining all linearizations
	using the complex step. The data from the real sim object will simply be copied into
	the complex counterpart between each design cycle.
	
	Arguments:
		sim = The simulation data structure
	
	Return:
		The Optimization_Case data structure holding the data needed for the optimization.

	*/

	struct Optimization_Case *optimization_case = constructor_Optimization_Case(sim);


	// Create the complex simulation object

	struct Integration_Test_Info* int_test_info = constructor_Integration_Test_Info(sim->ctrl_name);
	
	const int* p_ref  = int_test_info->p_ref,
	         * ml_ref = int_test_info->ml;

	struct Simulation* sim_c = NULL;
	const char type_rc = 'c';

	bool ignore_static = false;
	int ml_max = ml_ref[1];

	UNUSED(ml_max);

	int ml_prev = ml_ref[0]-1,
    	p_prev  = p_ref[0]-1;

	int ml = ml_ref[0];
	int p = p_ref[0];

	const int adapt_type = int_test_info->adapt_type;

	structor_simulation(&sim_c, 'c', adapt_type, p, ml, p_prev, ml_prev, sim->ctrl_name, 
		type_rc, ignore_static);
	

	// Store pointers to the sim objects
	optimization_case->sim = sim;
	optimization_case->sim_c = sim_c;

	return optimization_case;

}


static struct Multiarray_d* compute_gradient(struct Optimization_Case* optimization_case){

	/*
	Compute the gradient of the objective function with respect to the design 
	variables. This is equal to 

		grad_I = dI_dXp - Chi_T * dR_dXp

	where Chi_T is the transpose of the adjoint.


	Arguments:
		optimization_case = The data structure holding the optimization information

	Return:
		A multiarray with the gradient of the objective function with respect to the 
		design variables. The multiarray is of dimension [1 x num_design_dofs] and 
		is in column major form.
	*/

	int num_design_dofs = optimization_case->num_design_pts_dofs;

	struct Multiarray_d *Chi_T, *grad_I, *Chi_T_dot_dR_dXp;
	grad_I = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, num_design_dofs});  // returned
	
	Chi_T = constructor_copy_Multiarray_d(optimization_case->Chi);  // destructed
	transpose_Multiarray_d(Chi_T, false);
	struct Matrix_d* Mat_Chi_T = constructor_move_Matrix_d_d('C',
				Chi_T->extents[0], Chi_T->extents[1], false, Chi_T->data); // destructed

	Chi_T_dot_dR_dXp = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, num_design_dofs});  // destructed

	mm_NN1C_Multiarray_d((const struct const_Matrix_d*const)Mat_Chi_T,
		(const struct const_Multiarray_d*const)optimization_case->dR_dXp ,Chi_T_dot_dR_dXp);

	for (int i = 0; i < num_design_dofs; i++){
		grad_I->data[i] = optimization_case->dI_dXp->data[i] - Chi_T_dot_dR_dXp->data[i];
	}

	// Destruct allocated structures
	destructor_Matrix_d(Mat_Chi_T);
	destructor_Multiarray_d(Chi_T_dot_dR_dXp);
	destructor_Multiarray_d(Chi_T);

	return grad_I;

}

static struct Multiarray_d* test_compute_brute_force_gradient(struct Optimization_Case* optimization_case){

	/*
	Compute the gradient using a brute force approach (perturb each design point and 
	run the flow solver). Use this gradient for testing purposes.

	Arguments:
		optimization_case = The data structure holding the optimization information

	Return:
		Multiarray of dimnesion [1 x num_design_dofs] corresponding to the gradient of
		the objective function with respect to the design variables.
	*/

	struct Simulation *sim = optimization_case->sim;

	int num_design_dofs = (int)optimization_case->dI_dXp->extents[1];
	struct Multiarray_d *grad_I = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){1, num_design_dofs});  // returned
	
	struct Multiarray_i* ctrl_pts_opt = optimization_case->geo_data.control_points_optimization;
	int *ctrl_pt_indeces = get_col_Multiarray_i(0, ctrl_pts_opt);
	int n_pts = (int)ctrl_pts_opt->extents[0];

	struct Multiarray_d* ctrl_pts_and_weights = optimization_case->geo_data.control_points_and_weights;
	int control_pt_index;

	double I_0 = optimization_case->objective_function(sim);
	int col_index = 0;

	for (int i = 0; i < n_pts; i++){
		// Loop over the design points

		for (int j = 1; j <= 2; j++){
			// Loop over the degrees of freedom for the design point
			// j = 1 (x degree of freeedom) and j = 2 (y degree of freedom)

			if (!get_col_Multiarray_i(j, ctrl_pts_opt)[i])
				continue;

			control_pt_index = ctrl_pt_indeces[i];
			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] += FINITE_DIFF_STEP;


			update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)ctrl_pts_and_weights);
			set_up_solver_geometry(sim);
			solve_implicit(sim);
			grad_I->data[col_index] = (optimization_case->objective_function(sim) - I_0)/FINITE_DIFF_STEP;

			get_col_Multiarray_d(j-1, ctrl_pts_and_weights)[control_pt_index] -= FINITE_DIFF_STEP;

			col_index++;
			printf("col_index_brute_force : %d \n", col_index); fflush(stdout);
		}
	}

	// Update the geometry one last time with the original control point configurations
	update_geo_data_NURBS_parametric((const struct const_Multiarray_d*)ctrl_pts_and_weights);
	set_up_solver_geometry(sim);
	solve_implicit(sim);

	return grad_I;

}





