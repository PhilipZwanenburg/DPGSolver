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

#include "compute_error.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "mpi.h"

#include "macros.h"
#include "definitions_intrusive.h"

#include "multiarray.h"
#include "vector.h"

#include "computational_elements.h"
#include "element_error.h"
#include "volume_solver.h"

#include "const_cast.h"
#include "file_processing.h"
#include "intrusive.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Destructor for a \ref Error_CE container.
static void destructor_Error_CE
	(struct Error_CE* error_ce ///< Standard.
	);

/** \brief Compute the volume of the input \ref Solver_Volume.
 *  \return See brief. */
static double compute_volume
	(const struct Solver_Volume* s_vol ///< \ref Solver_Volume.
	);

/** \brief Constructor for a \ref const_Vector_d\* holding the cubature weights multiplied by the Jacobian determinants.
 *  \return See brief. */
static const struct const_Vector_d* constructor_w_detJ
	(const struct Solver_Volume* s_vol ///< \ref Solver_Volume.
	);

/** \brief Compute the number of 'd'egrees 'o'f 'f'reedom in the volume computational elements.
 *  \return See brief. */
static ptrdiff_t compute_dof_volumes
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Output the errors to the 's'erial/'p'arallel file.
static void output_errors_sp
	(const char sp_type,              ///< 's'erial/'p'arallel specifier.
	 const struct Error_CE* error_ce, ///< \ref Error_CE.
	 const struct Simulation* sim     ///< \ref Simulation.
	);

/// \brief Output the combined errors from all processors.
static void output_errors_global
	(const struct Error_CE* error_ce, ///< \ref Error_CE.
	 const struct Simulation* sim     ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

void output_error (const struct Simulation* sim)
{
	assert(sim->volumes->name == IL_SOLVER_VOLUME);
	assert(sim->faces->name   == IL_SOLVER_FACE);

	constructor_derived_Elements((struct Simulation*)sim,IL_SOLUTION_ELEMENT);
	constructor_derived_Elements((struct Simulation*)sim,IL_ELEMENT_ERROR);

	struct Error_CE* error_ce = sim->test_case->constructor_Error_CE(sim);
	const_cast_ptrdiff(&error_ce->dof,compute_dof(sim));

	output_errors_sp('s',error_ce,sim);
	output_errors_global(error_ce,sim);

	destructor_Error_CE(error_ce);

	destructor_derived_Elements((struct Simulation*)sim,IL_SOLUTION_ELEMENT);
	destructor_derived_Elements((struct Simulation*)sim,IL_ELEMENT);
}

double compute_domain_volume (const struct Simulation* sim)
{
	double domain_volume = 0.0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		domain_volume += compute_volume((struct Solver_Volume*) curr);
	return domain_volume;
}

void increment_vol_errors_l2_2
	(struct Vector_d* errors_l2_2, const struct const_Multiarray_d* err_v, const struct Solver_Volume* s_vol)
{
	assert(errors_l2_2->ext_0 == err_v->extents[1]);

	const struct const_Vector_d* w_detJ = constructor_w_detJ(s_vol); // destructed
	const ptrdiff_t ext_0 = w_detJ->ext_0;

	const ptrdiff_t n_out = errors_l2_2->ext_0;
	for (int i = 0; i < n_out; ++i) {
		const double* err_data = get_col_const_Multiarray_d(i,err_v);
		for (int j = 0; j < ext_0; ++j)
			errors_l2_2->data[i] += w_detJ->data[j]*err_data[j]*err_data[j];
	}
	destructor_const_Vector_d(w_detJ);
}

ptrdiff_t compute_dof (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	switch (sim->method) {
	case METHOD_DG:
		dof += compute_dof_volumes(sim);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}
	return dof;
}

const char* compute_error_file_name (const struct Simulation* sim)
{
	static const char* name_part = "../output/errors/";

	static char output_name[STRLEN_MAX];
	sprintf(output_name,"%s%s%c%s%c%s%s",
	        name_part,sim->pde_name,'/',sim->pde_spec,'/',"l2_errors__",extract_name(sim->ctrl_name_full,true));
	return output_name;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static void destructor_Error_CE (struct Error_CE* error_ce)
{
	destructor_const_Vector_d(error_ce->sol_L2);
	destructor_const_Vector_i(error_ce->expected_order);

	free(error_ce);
}

static double compute_volume (const struct Solver_Volume* s_vol)
{
	const struct const_Vector_d* w_detJ = constructor_w_detJ(s_vol); // destructed
	const ptrdiff_t ext_0 = w_detJ->ext_0;

	double volume = 0.0;
	for (int i = 0; i < ext_0; ++i)
		volume += w_detJ->data[i];

	destructor_const_Vector_d(w_detJ);

	return volume;
}

static const struct const_Vector_d* constructor_w_detJ (const struct Solver_Volume* s_vol)
{
	struct Volume* b_vol = (struct Volume*)s_vol;
	struct Error_Element* e = (struct Error_Element*) b_vol->element;

	const int curved = b_vol->curved,
	          p      = s_vol->p_ref;
	const struct const_Vector_d* w_vc = get_const_Multiarray_Vector_d(e->w_vc[curved],(ptrdiff_t[]){0,0,p,p});
	const struct const_Vector_d jacobian_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
	assert(w_vc->ext_0 == jacobian_det_vc.ext_0);

	const ptrdiff_t ext_0 = w_vc->ext_0;

	struct Vector_d* w_detJ = constructor_empty_Vector_d(ext_0); // returned

	for (int i = 0; i < ext_0; ++i)
		w_detJ->data[i] = w_vc->data[i]*jacobian_det_vc.data[i];

	return (const struct const_Vector_d*) w_detJ;
}

static ptrdiff_t compute_dof_volumes (const struct Simulation* sim)
{
	ptrdiff_t dof = 0;
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next)
		dof += ((struct Solver_Volume*)curr)->sol_coef->extents[0];
	return dof;
}

static void output_errors_sp (const char sp_type, const struct Error_CE* error_ce, const struct Simulation* sim)
{
	const char* output_name = compute_error_file_name(sim);
	FILE* sp_file = fopen_sp_output_file(sp_type,output_name,"txt",sim->mpi_rank); // closed

	const ptrdiff_t n_out = error_ce->sol_L2->ext_0;
	if (sp_type == 'p') {
		fprintf(sp_file,"n_out: %td\n",n_out);
		fprintf(sp_file,"expected_orders: ");

		for (int i = 0; i < n_out; ++i)
			fprintf(sp_file,"%2d ",error_ce->expected_order->data[i]);
		fprintf(sp_file,"\n\n");

		fprintf(sp_file,"%-14s%s\n","h",error_ce->header_spec);
		fprintf(sp_file,"%-14.4e",1.0/pow(error_ce->dof,1.0/(sim->d)));
	} else if (sp_type == 's') {
		fprintf(sp_file,"%-10s%-14s%s\n","dof","vol",error_ce->header_spec);
		fprintf(sp_file,"%-10td%-14.4e",error_ce->dof,error_ce->domain_volume);
	}

	for (int i = 0; i < n_out; ++i)
		fprintf(sp_file,"%-14.4e",error_ce->sol_L2->data[i]);

	fclose(sp_file);
}

static void output_errors_global (const struct Error_CE* error_ce, const struct Simulation* sim)
{
	if (!(sim->mpi_rank == 0))
		return;
	MPI_Barrier(MPI_COMM_WORLD);

	const char* file_name = compute_error_file_name(sim);

	const ptrdiff_t n_out = error_ce->sol_L2->ext_0;

	ptrdiff_t dof_g = 0;
	double domain_volume_g = 0.0;
	struct Vector_d* sol_L2_g = constructor_empty_Vector_d(n_out); // moved
	set_to_value_Vector_d(sol_L2_g,0.0);

	for (int mpi_rank = 0; mpi_rank < sim->mpi_size; ++mpi_rank) {
		FILE* s_file = fopen_sp_input_file('s',file_name,"txt",mpi_rank);

		skip_lines(s_file,1);
		char line[STRLEN_MAX];
		char* line_ptr[1] = {line};
		fgets(line,sizeof(line),s_file);

		int dof = 0;
		read_line_values_i(line_ptr,1,&dof,false);
		dof_g += dof;

		double domain_volume = 0.0;
		read_skip_d_1(line,1,&domain_volume,1);
		domain_volume_g += domain_volume;

		double sol_L2[n_out];
		read_skip_d_1(line,2,sol_L2,n_out);
		for (int i = 0; i < sol_L2_g->ext_0; ++i)
			sol_L2_g->data[i] += sol_L2[i]*sol_L2[i]*domain_volume;

		fclose(s_file);
	}

	for (int i = 0; i < sol_L2_g->ext_0; ++i)
		sol_L2_g->data[i] = sqrt(sol_L2_g->data[i]/domain_volume_g);

	struct Error_CE* error_ce_g = malloc(sizeof *error_ce_g); // destructed

	const_cast_d(&error_ce_g->domain_volume,domain_volume_g);
	error_ce_g->sol_L2 = (const struct const_Vector_d*) sol_L2_g; // keep
	error_ce_g->expected_order = constructor_copy_const_Vector_i_i(n_out,error_ce->expected_order->data); // keep
	const_cast_c1(&error_ce_g->header_spec,error_ce->header_spec);
	const_cast_ptrdiff(&error_ce_g->dof,dof_g);

	output_errors_sp('p',error_ce_g,sim);

	destructor_Error_CE(error_ce_g);
}
