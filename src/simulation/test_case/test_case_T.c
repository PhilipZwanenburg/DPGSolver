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

#include <float.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_core.h"
#include "definitions_dpg.h"
#include "definitions_test_case.h"


#include "def_templates_test_case.h"

#include "def_templates_flux.h"
#include "def_templates_restart.h"
#include "def_templates_solution_advection.h"
#include "def_templates_solution_diffusion.h"
#include "def_templates_solution_euler.h"
#include "def_templates_solution_navier_stokes.h"
#include "def_templates_solution_burgers_inviscid.h"

// Static function declarations ************************************************************************************* //

/// \brief Set associations between `char*` and `int` variables.
static void set_string_associations
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Set pde related parameters.
static void set_pde_related
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Set \ref Simulation::method related parameters.
static void set_method_related
	(struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

/// Set nonlinear solver related parameters.
static void set_nonlinear_solver_related
	(struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

/// \brief Set the function pointer members of \ref Test_Case_T.
static void set_function_pointers
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Read members of \ref Test_Case_T from input file.
static void read_test_case_parameters
	(struct Test_Case_T* test_case,    ///< \ref Test_Case_T.
	 const struct Simulation*const sim ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

struct Test_Case_T* constructor_Test_Case_T (const struct Simulation* sim)
{
	struct Test_Case_T* test_case = calloc(1,sizeof *test_case); // returned

	set_string_associations(test_case,sim);
	set_pde_related(test_case,sim);
	set_method_related(test_case,sim);
	set_nonlinear_solver_related(test_case,sim);

	read_test_case_parameters(test_case,sim);
	get_set_ind_num_flux(test_case->ind_num_flux);
	set_function_pointers(test_case,sim);

	test_case->solver_method_curr = 0;

	return test_case;
}

void destructor_Test_Case_T (const struct Test_Case_T* test_case)
{
	free((void*)test_case);
}

void increment_pointers_T (const int n_ptr, const Type**const ptrs)
{
	for (int i = 0; i < n_ptr; ++i)
		++ptrs[i];
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Correct invalid \ref Test_Case_T parameters if present.
static void correct_invalid_test_case_parameters
	(struct Test_Case_T* test_case, ///< \ref Test_Case_T.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Return a statically allocated array of \ref Flux_Input_T::compute_member flags depending on the pde and
 *         solver method.
 *  \return See brief. */
static const bool* get_compute_member_Flux_Input
	(const char type_ei,                 ///< 'e'xplicit/'i'mplicit type. Options: 'e', 'i'.
	 struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

/** \brief Return a statically allocated array of \ref Boundary_Value_Input_T::compute_member flags depending on the pde
 *         and solver method.
 *  \return See brief. */
static const bool* get_compute_member_Boundary_Value_Input
	(const char type_ei,                 ///< 'e'xplicit/'i'mplicit type. Options: 'e', 'i'.
	 struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	);

/// \brief Set the function pointers relating to the solution used at the start of the simulation.
static void set_function_pointers_start
	(struct Test_Case_T*const test_case ///< \ref Test_Case_T.
	);

static void set_string_associations (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	// pde_index
	if (strstr(sim->pde_name,"advection"))
		const_cast_i(&test_case->pde_index,PDE_ADVECTION);
	else if (strstr(sim->pde_name,"diffusion"))
		const_cast_i(&test_case->pde_index,PDE_DIFFUSION);
	else if (strstr(sim->pde_name,"euler"))
		const_cast_i(&test_case->pde_index,PDE_EULER);
	else if (strstr(sim->pde_name,"navier_stokes"))
		const_cast_i(&test_case->pde_index,PDE_NAVIER_STOKES);
	else if (strstr(sim->pde_name,"burgers_inviscid"))
		const_cast_i(&test_case->pde_index,PDE_BURGERS_INVISCID);
	else
		EXIT_ERROR("Unsupported: %s\n",sim->pde_name);
	get_set_pde_index(&test_case->pde_index);
}

static void set_pde_related (struct Test_Case_T* test_case, const struct Simulation* sim)
{
	switch (test_case->pde_index) {
	case PDE_ADVECTION:
		const_cast_b(&test_case->is_linear,true);
		const_cast_i(&test_case->n_var,1);
		const_cast_i(&test_case->n_eq,1);
		const_cast_b(&test_case->has_1st_order,true);
		const_cast_b(&test_case->has_2nd_order,false);
		break;
	case PDE_DIFFUSION:
		const_cast_b(&test_case->is_linear,true);
		const_cast_i(&test_case->n_var,1);
		const_cast_i(&test_case->n_eq,1);
		const_cast_b(&test_case->has_1st_order,false);
		const_cast_b(&test_case->has_2nd_order,true);
		break;
	case PDE_EULER:
		const_cast_b(&test_case->is_linear,false);
		const_cast_i(&test_case->n_var,DIM+2);
		const_cast_i(&test_case->n_eq,DIM+2);
		const_cast_b(&test_case->has_1st_order,true);
		const_cast_b(&test_case->has_2nd_order,false);
		break;
	case PDE_NAVIER_STOKES:
		const_cast_b(&test_case->is_linear,false);
		const_cast_i(&test_case->n_var,DIM+2);
		const_cast_i(&test_case->n_eq,DIM+2);
		const_cast_b(&test_case->has_1st_order,true);
		const_cast_b(&test_case->has_2nd_order,true);
		break;
	case PDE_BURGERS_INVISCID:
		const_cast_b(&test_case->is_linear,false);
		const_cast_i(&test_case->n_var,DIM);
		const_cast_i(&test_case->n_eq,DIM);
		const_cast_b(&test_case->has_1st_order,true);
		const_cast_b(&test_case->has_2nd_order,false);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}
	get_set_n_var_eq((int[]){test_case->n_var,test_case->n_eq});
	get_set_has_1st_2nd_order((bool[]){test_case->has_1st_order,test_case->has_2nd_order});

	const bool* flux_comp_mem_e = get_compute_member_Flux_Input('e',test_case,sim),
	          * flux_comp_mem_i = get_compute_member_Flux_Input('i',test_case,sim);
	for (int i = 0; i < MAX_FLUX_OUT; ++i) {
		const_cast_b(&test_case->flux_comp_mem_e[i],flux_comp_mem_e[i]);
		const_cast_b(&test_case->flux_comp_mem_i[i],flux_comp_mem_i[i]);
	}

	const bool* boundary_value_comp_mem_e = get_compute_member_Boundary_Value_Input('e',test_case,sim),
	          * boundary_value_comp_mem_i = get_compute_member_Boundary_Value_Input('i',test_case,sim);
	for (int i = 0; i < MAX_BV_OUT; ++i) {
		const_cast_b(&test_case->boundary_value_comp_mem_e[i],boundary_value_comp_mem_e[i]);
		const_cast_b(&test_case->boundary_value_comp_mem_i[i],boundary_value_comp_mem_i[i]);
	}
}

static void set_method_related (struct Test_Case_T*const test_case, const struct Simulation*const sim)
{
	const_cast_b(&test_case->required_unknowns[0],true);
	for (int i = 1; i < MAX_N_UNKNOWNS; ++i)
		const_cast_b(&test_case->required_unknowns[i],false);

	switch (sim->method) {
	case METHOD_DG:
		if (test_case->has_2nd_order)
			const_cast_b(&test_case->required_unknowns[2],true);
		break;
	case METHOD_DPG: // fallthrough
	case METHOD_OPG:
		const_cast_b(&test_case->required_unknowns[1],true);
		if (test_case->has_2nd_order) {
			for (int i = 2; i < MAX_N_UNKNOWNS; ++i)
				const_cast_b(&test_case->required_unknowns[i],true);
		}
		break;
	case METHOD_OPGC0:
		if (get_set_has_1st_2nd_order(NULL)[1])
			EXIT_ADD_SUPPORT;
		break;
	case METHOD_L2_PROJ:
		break; // Do nothing.
	default:
		EXIT_ERROR("Unsupported: %d.",sim->method);
		break;
	}
}

static void set_function_pointers (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	switch (test_case->pde_index) {
		case PDE_ADVECTION:     set_function_pointers_solution_advection_T(test_case,sim);     break;
		case PDE_DIFFUSION:     set_function_pointers_solution_diffusion_T(test_case,sim);     break;
		case PDE_EULER:         set_function_pointers_solution_euler_T(test_case,sim);         break;
		case PDE_NAVIER_STOKES: set_function_pointers_solution_navier_stokes_T(test_case,sim); break;
		case PDE_BURGERS_INVISCID: set_function_pointers_solution_burgers_inviscid_T(test_case,sim); break;
		default: EXIT_ERROR("Unsupported: %d\n",test_case->pde_index); break;
	}

	set_function_pointers_start(test_case);
	if (!using_restart())
		test_case->constructor_Error_CE_restart_test = test_case->constructor_Error_CE;
}

static void read_test_case_parameters (struct Test_Case_T* test_case, const struct Simulation*const sim)
{
	const int count_to_find = 1;

	int count_found = 0,
	    count_tmp = 0;
	char line[STRLEN_MAX];
	FILE* input_file = NULL;

	input_file = fopen_input('t',NULL,NULL); // closed
	while (fgets(line,sizeof(line),input_file)) {
		read_skip_convert_const_i(line,"solver_proc",  &test_case->solver_proc,  &count_found);
		read_skip_convert_const_i(line,"solver_type_e",&test_case->solver_type_e,NULL);
		read_skip_convert_const_i(line,"solver_type_i",&test_case->solver_type_i,NULL);
		read_skip_convert_const_i(line,"lhs_terms",    &test_case->lhs_terms,    NULL);
		read_skip_string_count_const_d("cfl_initial",&count_tmp,line,&test_case->cfl_initial);
		read_skip_string_count_const_d("cfl_max",&count_tmp,line,&test_case->cfl_max);

		read_skip_convert_const_i(line,"geom_parametrization",&test_case->geom_parametrization,NULL);

		read_skip_convert_const_i(line,"num_flux_1st",&test_case->ind_num_flux[0],NULL);
		read_skip_convert_const_i(line,"num_flux_2nd",&test_case->ind_num_flux[1],NULL);
		read_skip_convert_const_i(line,"test_norm",&test_case->ind_test_norm,NULL);
		read_skip_convert_const_i(line,"conservation",&test_case->ind_conservation,NULL);

		if (strstr(line,"time_final")) read_skip_const_d(line,&test_case->time_final,1,false);
		if (strstr(line,"time_step"))  read_skip_const_d(line,&test_case->dt,1,false);

		if (strstr(line,"use_schur_complement")) read_skip_const_b(line,&test_case->use_schur_complement);

		if (strstr(line,"display_progress")) read_skip_const_b(line,&test_case->display_progress);
		if (strstr(line,"has_functional"))   read_skip_const_b(line,&test_case->has_functional);

		read_skip_string_count_const_d("exit_tol_e",  &count_tmp,line,&test_case->exit_tol_e);
		read_skip_string_count_const_d("exit_ratio_e",&count_tmp,line,&test_case->exit_ratio_e);
		read_skip_string_count_const_d("exit_tol_i",  &count_tmp,line,&test_case->exit_tol_i);
		read_skip_string_count_const_d("exit_ratio_i",&count_tmp,line,&test_case->exit_ratio_i);
	}
	fclose(input_file);
	const_cast_b(&test_case->copy_initial_rhs,false);

	correct_invalid_test_case_parameters(test_case,sim);

	if (count_found != count_to_find)
		EXIT_ERROR("Did not find the required number of variables");
}

// Level 1 ********************************************************************************************************** //

static void correct_invalid_test_case_parameters (struct Test_Case_T* test_case, const struct Simulation* sim)
{
	switch (sim->method) {
	case METHOD_DG:  // fallthrough
	case METHOD_OPG: // fallthrough
	case METHOD_OPGC0:
		const_cast_b(&test_case->use_schur_complement,false);
		break;
	case METHOD_DPG:
	case METHOD_L2_PROJ:
		break; // Do nothing.
	default:
		EXIT_ERROR("Unsupported: %d\n",sim->method);
		break;
	}

	switch (test_case->lhs_terms) {
	case 0:
		const_cast_i(&test_case->lhs_terms,LHS_FULL_NEWTON);
		break;
	case LHS_FULL_NEWTON:
		break; // Do nothing.
	case LHS_CFL_RAMPING:
		switch (test_case->pde_index) {
		case PDE_ADVECTION: // fallthrough
		case PDE_DIFFUSION:
			EXIT_ERROR("Unsupported: Use full newton for linear problems.\n");
			break;
		case PDE_EULER:         // fallthrough
		case PDE_NAVIER_STOKES:
			break; // do nothing.
		default:
			EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
			break;
		}
		assert(test_case->cfl_initial != 0.0);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->lhs_terms);
		break;
	}
}

static const bool* get_compute_member_Flux_Input
	(const char type_ei, struct Test_Case_T*const test_case, const struct Simulation*const sim)
{
	static const bool cm_100000[] = {1,0,0,0,0,0,},
	                  cm_110000[] = {1,1,0,0,0,0,},
	                  cm_101000[] = {1,0,1,0,0,0,},
	                  cm_111000[] = {1,1,1,0,0,0,},
	                  cm_110100[] = {1,1,0,1,0,0,};

	assert(type_ei == 'e' || type_ei == 'i');
	switch (test_case->pde_index) {
	case PDE_ADVECTION:
		if (type_ei == 'e')
			return cm_100000;
		else if (type_ei == 'i')
			return cm_110000;
		break;
	case PDE_DIFFUSION:
		if (type_ei == 'e')
			return cm_100000;
		else if (type_ei == 'i')
			return cm_101000;
		break;
	case PDE_EULER:
		switch (sim->method) {
		case METHOD_DG:
			if (type_ei == 'e')
				return cm_100000;
			else if (type_ei == 'i')
				return cm_110000;
			break;
		case METHOD_DPG:
			if (type_ei == 'e')
				return cm_110000;
			else if (type_ei == 'i')
				return cm_110100;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",sim->method);
			break;
		}
		break;
	case PDE_NAVIER_STOKES:
		switch (sim->method) {
		case METHOD_DG:
			if (type_ei == 'e')
				return cm_100000;
			else if (type_ei == 'i')
				return cm_111000;
			break;
		case METHOD_DPG:
			EXIT_ADD_SUPPORT;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",sim->method);
			break;
		}
		break;
	case PDE_BURGERS_INVISCID:
		switch (sim->method) {
		case METHOD_DG:
			if (type_ei == 'e')
				return cm_100000;
			else if (type_ei == 'i')
				return cm_110000;
			break;
		default:
			EXIT_ERROR("Unsupported: %d\n",sim->method);
			break;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}
	EXIT_ERROR("Should not have reached this point.\n");
}

static const bool* get_compute_member_Boundary_Value_Input
	(const char type_ei, struct Test_Case_T*const test_case, const struct Simulation*const sim)
{
	UNUSED(sim);
	static const bool cm_100000[] = {1,0,0,0,0,0,},
	                  cm_110000[] = {1,1,0,0,0,0,},
	                  cm_101000[] = {1,0,1,0,0,0,},
	                  cm_111100[] = {1,1,1,1,0,0,};

	assert(type_ei == 'e' || type_ei == 'i');
	switch (test_case->pde_index) {
	case PDE_ADVECTION:
		if (type_ei == 'e')
			return cm_100000;
		else if (type_ei == 'i')
			return cm_110000;
		break;
	case PDE_DIFFUSION:
		if (type_ei == 'e')
			return cm_101000;
		else if (type_ei == 'i')
			return cm_111100;
		break;
	case PDE_EULER:
		if (type_ei == 'e')
			return cm_100000;
		else if (type_ei == 'i')
			return cm_110000;
		break;
	case PDE_NAVIER_STOKES:
		if (type_ei == 'e')
			return cm_101000;
		else if (type_ei == 'i')
			return cm_111100; // May potentially require cm_111110 in future.
		break;
	case PDE_BURGERS_INVISCID:
		if (type_ei == 'e')
			return cm_100000;
		else if (type_ei == 'i')
			return cm_110000;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->pde_index);
		break;
	}
	EXIT_ERROR("Should not have reached this point.\n");
}

static void set_nonlinear_solver_related
	(struct Test_Case_T*const test_case, ///< \ref Test_Case_T.
	 const struct Simulation*const sim   ///< \ref Simulation.
	)
{
	UNUSED(sim);
	test_case->cfl_max = 1e30;
}

static void set_function_pointers_start (struct Test_Case_T*const test_case)
{
	if (!using_restart()) {
		test_case->set_sol_start         = test_case->set_sol;
		test_case->constructor_sol_start = test_case->constructor_sol;
	} else {
		test_case->set_sol_start         = set_sol_restart_T;
		test_case->constructor_sol_start = constructor_const_sol_restart_T;
	}

	if (get_set_method(NULL) == METHOD_L2_PROJ)
		const_cast_b(&test_case->copy_initial_rhs,false);
}

#include "undef_templates_test_case.h"

#include "undef_templates_flux.h"
#include "undef_templates_restart.h"
#include "undef_templates_solution_advection.h"
#include "undef_templates_solution_diffusion.h"
#include "undef_templates_solution_euler.h"
#include "undef_templates_solution_navier_stokes.h"
#include "undef_templates_solution_burgers_inviscid.h"
