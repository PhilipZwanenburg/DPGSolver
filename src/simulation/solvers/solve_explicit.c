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

#include "solve_explicit.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>

#include "macros.h"
#include "definitions_intrusive.h"
#include "definitions_test_case.h"

#include "computational_elements.h"
#include "volume_solver.h"
#include "volume_solver_dg.h"

#include "multiarray.h"

#include "intrusive.h"
#include "simulation.h"
#include "solve.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Function pointer to a time-stepping function.
 *  \return The absolute value of the maximum rhs at the current time.
 *  \param dt  The time step.
 *  \param sim \ref Simulation.
 */
typedef double (*time_step_fptr)
	(const double dt,
	 const struct Simulation* sim
	);

/** \brief Set the function pointer to the appropriate time-stepping function.
 *  \return See brief. */
time_step_fptr set_time_step
	(const struct Simulation* sim ///< \ref Simulation
	);

/// \brief Display the solver progress.
void display_progress
	(const struct Test_Case* test_case, ///< \ref Test_Case.
	 const int t_step,                  ///< The current time step.
	 const double max_rhs,              ///< The current maximum value of the rhs term.
	 const double max_rhs0              ///< The initial maximum value of the rhs term.
	);

/** \brief Check the exit conditions.
 *  \return `true` if exit conditions are satisfied, `false` otherwise. */
bool check_exit
	(const struct Test_Case* test_case, ///< \ref Test_Case.
	 const double max_rhs,              ///< The current maximum value of the rhs term.
	 const double max_rhs0              ///< The initial maximum value of the rhs term.
	);

// Interface functions ********************************************************************************************** //

void solve_explicit (struct Simulation* sim)
{
	assert(sim->method == METHOD_DG); // Can be made flexible in future.

	sim->test_case->solver_method_curr = 'e';
	constructor_derived_Elements(sim,IL_ELEMENT_SOLVER_DG);       // destructed
	constructor_derived_computational_elements(sim,IL_SOLVER_DG); // destructed

	time_step_fptr time_step = set_time_step(sim);

	struct Test_Case* test_case = sim->test_case;
	const double time_final = test_case->time_final;
	double dt = test_case->dt;
	assert(time_final > 0.0);

	double max_rhs0 = 0.0;

	test_case->time = 0.0;
	for (int t_step = 0; ; ++t_step) {
		if (test_case->time + dt > time_final)
			dt = time_final - test_case->time;
		test_case->time += dt;

		const double max_rhs = time_step(dt,sim);
		if (t_step == 0)
			max_rhs0 = max_rhs;

		display_progress(test_case,t_step,max_rhs,max_rhs0);
		if (check_exit(test_case,max_rhs,max_rhs0))
			break;
	}

	destructor_derived_computational_elements(sim,IL_SOLVER);
	destructor_derived_Elements(sim,IL_ELEMENT);
	sim->test_case->solver_method_curr = 0;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Perform time-stepping using the forward Euler method.
 *  \return See \ref time_step_fptr. */
static double time_step_euler
	(const double dt,             ///< The time step.
	 const struct Simulation* sim ///< Defined for \ref time_step_fptr.
	);

/** \brief Perform time-stepping using the strong stability preserving Runge-Kutta 3-stage, 3rd order method.
 *  \return See \ref time_step_fptr.
 *
 *  Reference: eq. (4.2) \cite Gottlieb2001.
 */
static double time_step_ssp_rk_33
	(const double dt,             ///< The time step.
	 const struct Simulation* sim ///< Defined for \ref time_step_fptr.
	);

/** \brief Perform time-stepping using the low-storage Runge-Kutta 5-stage, 4rd order method.
 *  \return See \ref time_step_fptr. */
static double time_step_ls_rk_54
	(const double dt,             ///< The time step.
	 const struct Simulation* sim ///< Defined for \ref time_step_fptr.
	);

time_step_fptr set_time_step (const struct Simulation* sim)
{
	const struct Test_Case* test_case = sim->test_case;
	assert((test_case->solver_proc == SOLVER_E) || (test_case->solver_proc == SOLVER_EI));

	switch (test_case->solver_type_e) {
	case SOLVER_E_EULER:
		return time_step_euler;
		break;
	case SOLVER_E_SSP_RK_33:
		return time_step_ssp_rk_33;
		break;
	case SOLVER_E_LS_RK_54:
		return time_step_ls_rk_54;
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_type_e);
		break;
	}
}

void display_progress (const struct Test_Case* test_case, const int t_step, const double max_rhs, const double max_rhs0)
{
	if (!test_case->display_progress)
		return;

	switch (test_case->solver_proc) {
	case SOLVER_E:
		printf("Complete: % 7.2f%%, tstep: %8d, maxRHS: % .3e\n",
		       100*(test_case->time)/(test_case->time_final),t_step,max_rhs);
		break;
	case SOLVER_EI:
		printf("Exit Conditions (tol, ratio): % .3e (% .3e), % .3e (% .3e), tstep: %8d\n",
		       max_rhs,test_case->exit_tol_e,max_rhs/max_rhs0,test_case->exit_ratio_e,t_step);
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}
}

bool check_exit (const struct Test_Case* test_case, const double max_rhs, const double max_rhs0)
{
	bool exit_now = false;
	switch (test_case->solver_proc) {
	case SOLVER_E:
		if (test_case->time == test_case->time_final)
			exit_now = true;
		break;
	case SOLVER_EI:
		if (max_rhs < test_case->exit_tol_e) {
			printf("Complete: max_rhs is below the exit tolerance.\n");
			exit_now = true;
		}

		if (max_rhs/max_rhs0 < test_case->exit_ratio_e) {
			printf("Complete: max_rhs dropped by % .2e orders.\n",log10(max_rhs/max_rhs0));
			exit_now = true;
		}
		break;
	default:
		EXIT_ERROR("Unsupported: %d\n",test_case->solver_proc);
		break;
	}
	return exit_now;
}

// Level 1 ********************************************************************************************************** //

static double time_step_euler (const double dt, const struct Simulation* sim)
{
	assert(sim->method == METHOD_DG);
UNUSED(dt);
UNUSED(sim);
EXIT_ADD_SUPPORT;
}

static double time_step_ssp_rk_33 (const double dt, const struct Simulation* sim)
{
	assert(sim->method == METHOD_DG);
	assert(sim->volumes->name == IL_VOLUME_SOLVER_DG);
	assert(sim->faces->name   == IL_FACE_SOLVER_DG);

	double max_rhs = 0.0;
	for (int rk = 0; rk < 3; rk++) {
		max_rhs = compute_rhs(sim);
		for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
			struct Solver_Volume*    s_vol    = (struct Solver_Volume*) curr;
			struct DG_Solver_Volume* s_vol_dg = (struct DG_Solver_Volume*) curr;

			struct Multiarray_d* sol_coef   = s_vol->sol_coef,
			                   * sol_coef_p = s_vol_dg->sol_coef_p,
			                   * rhs        = s_vol_dg->rhs;

			double* data_s   = sol_coef->data,
			      * data_sp  = sol_coef_p->data,
			      * data_rhs = rhs->data;

			const ptrdiff_t i_max = compute_size(sol_coef->order,sol_coef->extents);
			assert(i_max == compute_size(sol_coef_p->order,sol_coef_p->extents));
			assert(i_max == compute_size(rhs->order,rhs->extents));

			switch (rk) {
			case 0:
				for (ptrdiff_t i = 0; i < i_max; ++i) {
					data_sp[i]  = data_s[i];
					data_s[i]  += dt*data_rhs[i];
				}
				break;
			case 1:
				for (ptrdiff_t i = 0; i < i_max; ++i)
					data_s[i] = (1.0/4.0)*(3.0*data_sp[i] + data_s[i] + dt*data_rhs[i]);
				break;
			case 2:
				for (ptrdiff_t i = 0; i < i_max; ++i)
					data_s[i] = (1.0/3.0)*(data_sp[i] + 2.0*data_s[i] + 2.0*dt*data_rhs[i]);
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",rk);
				break;
			}
//			enforce_positivity_highorder(sol_coef);
		}
	}

	return max_rhs;
EXIT_ADD_SUPPORT;
}

static double time_step_ls_rk_54 (const double dt, const struct Simulation* sim)
{
// use 26 decimal precision from carpenter for coefficients.
	assert(sim->method == METHOD_DG);
UNUSED(dt);
UNUSED(sim);
EXIT_ADD_SUPPORT;
}