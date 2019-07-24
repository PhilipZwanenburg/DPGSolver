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

#include "output_progress.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include <math.h>

#include "macros.h"

#include "simulation.h"
#include "math_functions.h"
#include "definitions_alloc.h"
#include "file_processing.h"

#include "intrusive.h"
#include "volume.h"
#include "volume_solver.h"
#include "face_solver.h"
#include "element_solver.h"
#include "element_solution.h"
#include "operator.h"
#include "compute_face_rlhs.h"
#include "solution_euler.h"
#include "boundary.h"
#include "element_operators.h"
#include "definitions_nodes.h"
#include "nodes.h"

#include "matrix.h"
#include "matrix_constructors.h"
#include "matrix_math.h"
#include "multiarray.h"
#include "multiarray_constructors.h"
#include "multiarray_math.h"
#include "vector.h"

#include "functionals.h"
#include "optimization_case.h"
#include "gradient.h"

// Static function declarations ************************************************************************************* //


/** \brief 	Read the Optimization Value Monitoring section of the optimization.data file and
 * 	determine which parameters should be printed or monitored. Return the string holding the
 *	keywords for which parameters to print.
 */
static char* read_optimization_data(
	);


// Interface functions ********************************************************************************************** //

FILE* constructor_optimization_progress_file(struct Optimization_Case *optimization_case){

	struct Simulation *sim = optimization_case->sim;

	// Create the optimization convergence output file
	char f_name[4*STRLEN_MAX] = { 0, };
	sprintf(f_name,"%s%c%s%c%sML%d_P%d_", sim->pde_name,'/',sim->pde_spec,'/',optimization_case->optimizer_output_files_prefix, 
	sim->ml[0], sim->p_ref[0]);
	strcat(f_name, "Optimization_Convergence.txt");

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../output/");
	strcat(output_name,"optimization/");
	strcat(output_name,f_name);

	FILE* fp;
	if ((fp = fopen_create_dir(output_name)) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);


	// Print the header for the data
	// - Standard data headers
	fprintf(fp, "Design_Iteration CPU_Time(s) L2_Norm_Gradient_Cost_Function Cost_Function ");

	// - Additional headers (based on the test case)
	char *progress_additional_parameters = read_optimization_data(); 
	if (strstr(progress_additional_parameters, "FUNCTIONAL_CL"))
		fprintf(fp, "cl ");
	if (strstr(progress_additional_parameters, "FUNCTIONAL_CM_LE"))
		fprintf(fp, "cm ");
	if (strstr(progress_additional_parameters, "FUNCTIONAL_TARGET_CL"))
		fprintf(fp, "cl_target ");
	if (strstr(progress_additional_parameters, "FUNCTIONAL_MESH_VOLUME"))
		fprintf(fp, "mesh_volume ");
	if (strstr(progress_additional_parameters, "FUNCTIONAL_FRACTIONAL_CHANGE_MESH_VOLUME"))
		fprintf(fp, "volume_fractional_change ");
	if (strstr(progress_additional_parameters, "FUNCTIONAL_INVERSE_PRESSURE_DESIGN"))
		fprintf(fp, "functional_inverse_pressure_design ");
	fprintf(fp, "\n");

	return fp;
}


void destructor_optimization_progress_file(FILE *fp){

	fclose(fp);

}


void progress_file_add_information(FILE *fp, struct Optimization_Case *optimization_case, double L2_grad, 
	double objective_func_value, double time_elapsed, int design_iteration, bool output_to_stdout){


	// \todo
	// In the optimization.data file, have a section where the user can list what data 
	// they want to monitor in the convergence (full names) such as lift coefficient, 
	// pitching moment ... Accordingly, compute those values and output them to the 
	// convergence file.

	// Compute the standard output information and print it to the output file and standard out
	fprintf(fp, "%d %.14e %.14e %.14e ", design_iteration, time_elapsed, L2_grad, objective_func_value);

	char std_out_line[2*STRLEN_MAX];
	int std_out_line_cursor_index = 0;
	if (output_to_stdout){
		
		// Create the information line that will go to stdout
		std_out_line_cursor_index += sprintf(std_out_line + std_out_line_cursor_index, 
			" Design_Iteration: %d   Time : %f   L2_grad : %e   obj_func : %e   ", 
			design_iteration, time_elapsed, L2_grad, objective_func_value);
	}
	
	// Print additional test specific information
	// Get the line of keywords for additional parameters to print and, based on those
	// keywords, print the values.
	char *progress_additional_parameters = read_optimization_data(); 
	if (strstr(progress_additional_parameters, "FUNCTIONAL_CL")){

		double cl = functional_cl(optimization_case->sim);
		
		fprintf(fp, "%.14e ", cl);
		if (output_to_stdout)
			std_out_line_cursor_index += sprintf(std_out_line+std_out_line_cursor_index, "cl: %f ", cl);
	}
	if (strstr(progress_additional_parameters, "FUNCTIONAL_CM_LE")){

		double cm = functional_cm_le(optimization_case->sim);
		
		fprintf(fp, "%.14e ", cm);
		if (output_to_stdout)
			std_out_line_cursor_index += sprintf(std_out_line+std_out_line_cursor_index, "cm: %f ", cm);
	}
	if (strstr(progress_additional_parameters, "FUNCTIONAL_TARGET_CL")){

		double target_cl = functional_target_cl(optimization_case->sim);
		
		fprintf(fp, "%.14e ", target_cl);
		if (output_to_stdout)
			std_out_line_cursor_index += sprintf(std_out_line+std_out_line_cursor_index, "target_cl: %f ", target_cl);
	}
	if (strstr(progress_additional_parameters, "FUNCTIONAL_MESH_VOLUME")){

		double mesh_vol = functional_mesh_volume(optimization_case->sim);
		
		fprintf(fp, "%.14e ", mesh_vol);
		if (output_to_stdout)
			std_out_line_cursor_index += sprintf(std_out_line+std_out_line_cursor_index, "mesh_vol: %f ", mesh_vol);
	}
	if (strstr(progress_additional_parameters, "FUNCTIONAL_FRACTIONAL_CHANGE_MESH_VOLUME")){

		double v_change = functional_mesh_volume_fractional_change(optimization_case->sim);
		
		fprintf(fp, "%.14e ", v_change);
		if (output_to_stdout)
			std_out_line_cursor_index += sprintf(std_out_line+std_out_line_cursor_index, "vfc: %e ", v_change);
	}	
	if (strstr(progress_additional_parameters, "FUNCTIONAL_INVERSE_PRESSURE_DESIGN")){

		double inv_p = functional_inverse_pressure_design(optimization_case->sim);
		
		fprintf(fp, "%.14e ", inv_p);
		if (output_to_stdout)
			std_out_line_cursor_index += sprintf(std_out_line+std_out_line_cursor_index, "inv_p: %f ", inv_p);
	}

	// Finalize print to the output file and to standard out
	fprintf(fp, "\n");

	if (output_to_stdout){

		// Create the divider (using asterixes) for making the data more visible
		int line_len = (int)strlen(std_out_line);
		char std_out_divider_line[2*STRLEN_MAX];

		for (int i = 0; i < line_len; i++)
			sprintf(std_out_divider_line+i, "*");

		// Print to standard out
		printf("\n");
		printf("%s\n", std_out_divider_line);
		printf("%s\n", std_out_line);
		printf("%s\n", std_out_divider_line);
		printf("\n");
		fflush(stdout);
	}	
}


void output_NURBS_patch_information(struct Optimization_Case* optimization_case){

	struct Simulation *sim = optimization_case->sim;

	char f_name[4*STRLEN_MAX] = { 0, };
	sprintf(f_name,"%s%c%s%c%sML%d_P%d_", sim->pde_name,'/',sim->pde_spec,'/',optimization_case->optimizer_output_files_prefix, 
	sim->ml[0], sim->p_ref[0]);
	strcat(f_name, "Optimized_NURBS_Patch.txt");

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../output/");
	strcat(output_name,"optimization/");
	strcat(output_name,f_name);

	FILE* fp;
	if ((fp = fopen_create_dir(output_name)) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Print the patch information
	fprintf(fp, "P(xi_order) %d\n", optimization_case->geo_data.P);
	fprintf(fp, "Q(eta_order) %d\n", optimization_case->geo_data.Q);
	fprintf(fp, "\n");

	fprintf(fp, "knots_xi %d\n", (int)optimization_case->geo_data.knots_xi->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.knots_xi->extents[0]; i++){
		fprintf(fp, "%.14e\n", optimization_case->geo_data.knots_xi->data[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "knots_eta %d\n", (int)optimization_case->geo_data.knots_eta->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.knots_eta->extents[0]; i++){
		fprintf(fp, "%.14e\n", optimization_case->geo_data.knots_eta->data[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Control_Point_Data %d \n", (int)optimization_case->geo_data.control_points->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.control_points->extents[0]; i++){
		fprintf(fp, "%.14e %.14e %.14e\n", 
			get_col_Multiarray_d(0, optimization_case->geo_data.control_points)[i],
			get_col_Multiarray_d(1, optimization_case->geo_data.control_points)[i],
			get_col_Multiarray_d(0, optimization_case->geo_data.control_weights)[i]
			);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Control_Point_Connectivity %d %d\n", 
		(int)optimization_case->geo_data.control_pt_wt_connectivity->extents[0],
		(int)optimization_case->geo_data.control_pt_wt_connectivity->extents[1]);
	for (int i = 0; i < (int)optimization_case->geo_data.control_pt_wt_connectivity->extents[0]; i++){
		for (int j = 0; j < (int)optimization_case->geo_data.control_pt_wt_connectivity->extents[1]; j++){
			fprintf(fp, "%d ", get_row_Multiarray_i(i, optimization_case->geo_data.control_pt_wt_connectivity)[j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);


	// Temporary
	// Write the output to Optimized_NURBS_Patch file in the NURBS_Airfoil directory for easy
	// progress monitoring using the Python script.


	sprintf(f_name,"%s%c%s%c%s", sim->pde_name,'/',sim->pde_spec,'/', "Optimized_NURBS_Patch.txt");

	strcpy(output_name,"../output/");
	strcat(output_name,"optimization/");
	strcat(output_name,f_name);

	if ((fp = fopen_create_dir(output_name)) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	// Print the patch information
	fprintf(fp, "P(xi_order) %d\n", optimization_case->geo_data.P);
	fprintf(fp, "Q(eta_order) %d\n", optimization_case->geo_data.Q);
	fprintf(fp, "\n");

	fprintf(fp, "knots_xi %d\n", (int)optimization_case->geo_data.knots_xi->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.knots_xi->extents[0]; i++){
		fprintf(fp, "%.14e\n", optimization_case->geo_data.knots_xi->data[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "knots_eta %d\n", (int)optimization_case->geo_data.knots_eta->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.knots_eta->extents[0]; i++){
		fprintf(fp, "%.14e\n", optimization_case->geo_data.knots_eta->data[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Control_Point_Data %d \n", (int)optimization_case->geo_data.control_points->extents[0]);
	for (int i = 0; i < (int)optimization_case->geo_data.control_points->extents[0]; i++){
		fprintf(fp, "%.14e %.14e %.14e\n", 
			get_col_Multiarray_d(0, optimization_case->geo_data.control_points)[i],
			get_col_Multiarray_d(1, optimization_case->geo_data.control_points)[i],
			get_col_Multiarray_d(0, optimization_case->geo_data.control_weights)[i]
			);
	}
	fprintf(fp, "\n");

	fprintf(fp, "Control_Point_Connectivity %d %d\n", 
		(int)optimization_case->geo_data.control_pt_wt_connectivity->extents[0],
		(int)optimization_case->geo_data.control_pt_wt_connectivity->extents[1]);
	for (int i = 0; i < (int)optimization_case->geo_data.control_pt_wt_connectivity->extents[0]; i++){
		for (int j = 0; j < (int)optimization_case->geo_data.control_pt_wt_connectivity->extents[1]; j++){
			fprintf(fp, "%d ", get_row_Multiarray_i(i, optimization_case->geo_data.control_pt_wt_connectivity)[j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);
}


void output_gradient(struct Optimization_Case* optimization_case, 
	double *gradient){

	struct Simulation *sim = optimization_case->sim;
	int n_dof = optimization_case->num_design_pts_dofs; // length of gradient

	char f_name[4*STRLEN_MAX] = { 0, };
	//sprintf(f_name,"%s%c%s%c%s%s%s%d%s%d", sim->pde_name,'/',sim->pde_spec,'/',optimization_case->optimizer_output_files_prefix,
	//'/', "ML", sim->ml[0], "_P", sim->p_ref[0]);
	sprintf(f_name,"%s%c%s%c%sML%d_P%d_", sim->pde_name,'/',sim->pde_spec,'/',optimization_case->optimizer_output_files_prefix, 
	sim->ml[0], sim->p_ref[0]);
	strcat(f_name, "Objective_Gradient.txt");

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../output/");
	strcat(output_name,"optimization/");
	strcat(output_name,f_name);

	FILE* fp;
	if ((fp = fopen_create_dir(output_name)) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);

	for (int i = 0; i < n_dof; i++){
		fprintf(fp, "%f\n", gradient[i]);
	}

	fclose(fp);
}


void output_pressure_distribution(struct Optimization_Case* optimization_case){

	assert(DIM == 2);

	printf("print to file\n");

	struct Simulation *sim = optimization_case->sim;

	char f_name[4*STRLEN_MAX] = { 0, };
	sprintf(f_name,"%s%c%s%c%sML%d_P%d_", sim->pde_name,'/',sim->pde_spec,'/',optimization_case->optimizer_output_files_prefix, 
	sim->ml[0], sim->p_ref[0]);
	strcat(f_name, "Pressure_Distribution.txt");

	char output_name[STRLEN_MAX] = { 0, };
	strcpy(output_name,"../output/");
	strcat(output_name,"optimization/");
	strcat(output_name,f_name);

	FILE* fp;
	if ((fp = fopen_create_dir(output_name)) == NULL)
		printf("Error: File %s did not open.\n", output_name), exit(1);


	// First determine the total number of face cubature nodes and print the number
	int total_num_fc_nodes = 0;
	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		const struct Face*const face = (struct Face*) curr;
		const struct Solver_Face*const s_face = (struct Solver_Face*) curr;
		if (!is_face_wall_boundary(face))
			continue;

		if (face->neigh_info[0].ind_lf == 2){
			// This is a eta = -1 face
			total_num_fc_nodes += s_face->p_ref + 1;
		}
	}

	fprintf(fp, "%d\n", total_num_fc_nodes);


	for (struct Intrusive_Link* curr = sim->faces->first; curr; curr = curr->next) {
		// Loop over all the faces that are on the wall boundary

		const struct Face*const face = (struct Face*) curr;
		const struct Solver_Face*const s_face = (struct Solver_Face*) curr;
		if (!is_face_wall_boundary(face))
			continue;


		const struct Volume*const vol = (struct Volume*) face->neigh_info[0].volume;

		const struct const_Element*const el_vol      = vol->element;
		const struct const_Multiarray_d*const xyz_ve = vol->xyz_ve;

		const struct const_Vector_i*const f_ve_lf = el_vol->f_ve->data[face->neigh_info[0].ind_lf];

		const double*const xyz_ve_a[] = { get_row_const_Multiarray_d(f_ve_lf->data[0],xyz_ve),
		                                  get_row_const_Multiarray_d(f_ve_lf->data[1],xyz_ve), };


		// Get the xi limits for this face
		double 	xi_val1  = xyz_ve_a[0][0],
				xi_val2  = xyz_ve_a[1][0],
				eta_val1 = xyz_ve_a[0][1],
				eta_val2 = xyz_ve_a[1][1];

		double xi_range[2], eta_range[2];
		xi_range[0]  = xi_val1  > xi_val2  ? xi_val2  : xi_val1;
		xi_range[1]  = xi_val1  > xi_val2  ? xi_val1  : xi_val2;
		eta_range[0] = eta_val1 > eta_val2 ? eta_val2 : eta_val1;
		eta_range[1] = eta_val1 > eta_val2 ? eta_val1 : eta_val2;


		// Obtain the face cubature (fc) nodes on the parametric knot domain
		const struct const_Nodes *nodes_c = constructor_const_Nodes_tp(face->element->d, s_face->p_ref, NODES_GL);
		const struct const_Multiarray_d *rst_i = constructor_move_const_Multiarray_d_Matrix_d(nodes_c->rst);  // free
		const ptrdiff_t n_fc = rst_i->extents[0];


		// \todo: Map these values using some operator and be able to generalize to any type of
		//	element (not just QUADS)

		// Map the node point locations to knot domain:
		struct Multiarray_d *rst_knots_i = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){n_fc, 2}); // destructed
	
		// Flag is set to true if this is an eta = -1 face
		bool eta_min_1_face = false; 

		switch (face->neigh_info[0].ind_lf){
			case 0:
				// Face 0 : xi is fixed (left face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = xi_range[0];
					get_col_Multiarray_d(1, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(eta_range[1] - eta_range[0])*0.5 
																+ 0.5*(eta_range[1] + eta_range[0]);
				}
				break;

			case 1:
				// Face 1 : xi is fixed (right face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = xi_range[1];
					get_col_Multiarray_d(1, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(eta_range[1] - eta_range[0])*0.5 
																+ 0.5*(eta_range[1] + eta_range[0]);
				}
				break;

			case 2:
				// Face 2 : eta is fixed (bottom face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(xi_range[1] - xi_range[0])*0.5 
																+ 0.5*(xi_range[1] + xi_range[0]);
					get_col_Multiarray_d(1, rst_knots_i)[i] = eta_range[0];
				}
				eta_min_1_face = true;
				break;

			case 3:
				// Face 3 : eta is fixed (top face)
				for (ptrdiff_t i = 0; i < n_fc; i++){
					get_col_Multiarray_d(0, rst_knots_i)[i] = get_col_const_Multiarray_d(0, rst_i)[i]*(xi_range[1] - xi_range[0])*0.5 
																+ 0.5*(xi_range[1] + xi_range[0]);
					get_col_Multiarray_d(1, rst_knots_i)[i] = eta_range[1];
				}
				break;

			default:
				EXIT_ERROR("Unsupported face index"); break;
		}

		// Get the solution (s) at the face cubature (fc) nodes and convert it to primitive variables (pv)
		const struct const_Multiarray_d* pv_fc = constructor_s_fc_interp_d(0, s_face);
		convert_variables((struct Multiarray_d*)pv_fc,'c','p');

		// Pressure (p) at face cubature (fc) nodes as an array
		const double *p_fc_i = get_col_const_Multiarray_d(pv_fc->extents[1]-1,pv_fc);

		if (eta_min_1_face){
			for (int i = 0; i < n_fc; i++){
				// print the xi point for the face cubature node and the pressure at the node
				fprintf(fp, "%.14e %.14e\n", get_col_Multiarray_d(0, rst_knots_i)[i], p_fc_i[i]);
			}
		}


		// Free allocated memory
		destructor_const_Multiarray_d(rst_i);
		destructor_Multiarray_d(rst_knots_i);
		destructor_const_Nodes(nodes_c);
		destructor_const_Multiarray_d(pv_fc);
	}

	fclose(fp);
}


// Static functions ************************************************************************************************* //


static char* read_optimization_data(){

	static bool need_input = true;
	static char progress_monitoring_parameter_keywords[STRLEN_MAX];

	if (need_input){
		need_input = false;

		FILE* input_file = fopen_input('o',NULL,NULL); // closed
		char line[STRLEN_MAX];

		// Read in the information from the file
		while (fgets(line,sizeof(line),input_file)) {

			// Any line with a comment flag should be skipped
			if (strstr(line, "//"))
				continue;

			if (strstr(line, "output_progress_additional_parameters")) read_skip_c(line, progress_monitoring_parameter_keywords);		
		}

		fclose(input_file);
	}

	return progress_monitoring_parameter_keywords;
}


