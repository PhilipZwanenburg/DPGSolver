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

#include "solve_dg.h"

#include <assert.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "gsl/gsl_math.h"
#include "petscvec.h"

#include "macros.h"
#include "definitions_test_case.h"
#include "definitions_tol.h"

#include "face.h"
#include "face_solver.h"
#include "element.h"
#include "element_solver_dg.h"
#include "volume.h"
#include "volume_solver_dg.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "compute_face_rlhs.h"
#include "compute_grad_coef_dg.h"
#include "compute_volume_rlhs_dg.h"
#include "compute_face_rlhs_dg.h"
#include "compute_rlhs.h"
#include "compute_source_rlhs_dg.h"
#include "const_cast.h"
#include "intrusive.h"
#include "math_functions.h"
#include "multiarray_operator.h"
#include "operator.h"
#include "simulation.h"
#include "solution_euler.h"
#include "solution_navier_stokes.h"
#include "solve.h"
#include "solve_implicit.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

/// \brief Call the common functions when computing the rlhs values for explicit and implicit dg schemes.
static void compute_rlhs_common_dg
	(const struct Simulation*const sim,       ///< \ref Simulation.
	 struct Solver_Storage_Implicit*const ssi ///< \ref Solver_Storage_Implicit.
	);

/// \brief Scale the rhs terms by the the inverse mass matrices (for explicit schemes).
static void scale_rhs_by_m_inv
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Fill \ref Solver_Storage_Implicit::b with the negated rhs values.
 *
 *  See comments in \ref solve_implicit.h for why the negated values are used here.
 */
static void fill_petsc_Vec_b_dg
	(const struct Simulation*const sim,       ///< Defined for \ref compute_rlhs_dg.
	 struct Solver_Storage_Implicit*const ssi ///< Defined for \ref compute_rlhs_dg.
	);

/** \brief Compute and add diagonal block terms to the LHS matrix in case \ref LHS_CFL_RAMPING was selected.
 *
 *  The addition of these terms corresponds to the use of the backwards (implicit) Euler method with bounded time step,
 *  as opposed to the full Newton method which is the limit of the backwards Euler as the time step goes to infinity.
 *  The 'C'ourant'F'riedrichs'L'ewy number is used to determine the time step and is increased as the residual
 *  decreases.
 *
 *  The explanation for what term must be added to the `A` matrix to achieve the desired effect follows.
 *
 *  The equation to be solved for the backwards Euler method for each \ref Volume is:
 *  \f[
 *  	M \frac{s_{coef}-s_{coef,curr}}{\Delta t} - \text{rhs}(s_{coef}) = 0,
 *  \f]
 *  where \f$ M \f$ is \ref DG_Solver_Volume_T::m.
 *
 *  Using Newton's method to solve the above equation taking \f$ s_{coef,curr} \f$ as the initial guess yields
 *  \f[
 *  	0 = M \frac{s_{coef,curr}-s_{coef,curr}}{\Delta t} - \text{rhs}(s_{coef}) = 0,
 *  \f]
 *  \f{eqnarray*}{
 *  	0 &=& M \frac{s_{coef,curr}-s_{coef,curr}}{\Delta t} - \text{rhs}(s_{coef,curr})
 *  	    + \left( M \frac{1}{\Delta t} - \text{lhs}(s_{coef,curr}) \right)\ \Delta(s_{coef}) \\
 *  	  &=& 0 - \text{rhs}(s_{coef,curr})
 *  	    + \left( M \frac{1}{\Delta t} - \text{lhs}(s_{coef,curr}) \right)\ \Delta(s_{coef}).
 *  \f}
 *
 *  Rearranging terms such that the system corresponds to that described in \ref solve_implicit.h :
 *  \f[
 *  	\left( - M \frac{1}{\Delta t} + \text{lhs}(s_{coef,curr}) \right)\ \Delta(s_{coef}) = - \text{rhs}(s_{coef,curr}).
 *  \f]
 *
 *  Thus the mass matrix contributions should be scaled by \f$ -\frac{1}{\Delta t} \f$.
 */
static void compute_CFL_ramping
	(struct Solver_Storage_Implicit*const ssi, ///< \ref Solver_Storage_Implicit.
	 const struct Simulation*const sim         ///< \ref Simulation.
	);

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solve_dg_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solve_dg_T.c"
#include "undef_templates_type.h"

double compute_rhs_dg (const struct Simulation* sim)
{
	compute_rlhs_common_dg(sim,NULL);
	scale_rhs_by_m_inv(sim);

	return compute_max_rhs_dg_like(sim);
}

double compute_rlhs_dg (const struct Simulation* sim, struct Solver_Storage_Implicit* ssi)
{
	compute_rlhs_common_dg(sim,ssi);
	compute_CFL_ramping(ssi,sim);
	fill_petsc_Vec_b_dg(sim,ssi);

	return compute_max_rhs_dg_like(sim);
}

void set_petsc_Mat_row_col_dg
	(struct Solver_Storage_Implicit*const ssi, const struct Solver_Volume* v_l, const int eq,
	 const struct Solver_Volume* v_r, const int vr)
{
	ssi->row = (int)(v_l->ind_dof+v_l->sol_coef->extents[0]*eq);
	ssi->col = (int)(v_r->ind_dof+v_r->sol_coef->extents[0]*vr);
}

void compute_flux_imbalances_dg (const struct Simulation*const sim)
{
	compute_flux_imbalances_faces_dg(sim);
	compute_flux_imbalances_source_dg(sim);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Version of \ref scale_rhs_by_m_inv used for non-collocated runs.
static void scale_rhs_by_m_inv_std
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/// \brief Version of \ref scale_rhs_by_m_inv used for collocated runs.
static void scale_rhs_by_m_inv_col
	(const struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Return the time step relating to the CFL ramping.
 *  \return See brief. */
static double compute_dt_cfl_constrained
	(const double max_rhs,                   ///< The maximum of the RHS terms (negated steady residual).
	 const struct Solver_Volume*const s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim       ///< \ref Simulation.
	);

static void compute_rlhs_common_dg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	initialize_zero_memory_volumes(sim->volumes);
	compute_grad_coef_dg(sim,sim->volumes,sim->faces);
	compute_volume_rlhs_dg(sim,ssi,sim->volumes);
	compute_face_rlhs_dg(sim,ssi,sim->faces);
	compute_source_rhs_dg_like(sim);
}

static void scale_rhs_by_m_inv (const struct Simulation*const sim)
{
	struct Test_Case* test_case = (struct Test_Case*)sim->test_case_rc->tc;
	assert((test_case->solver_proc == SOLVER_E) ||
	       (test_case->solver_proc == SOLVER_EI));

	if (!sim->collocated)
		scale_rhs_by_m_inv_std(sim);
	else
		scale_rhs_by_m_inv_col(sim);
}

static void fill_petsc_Vec_b_dg (const struct Simulation*const sim, struct Solver_Storage_Implicit*const ssi)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		const int ind_dof        = (int)((struct Solver_Volume*)curr)->ind_dof;
		struct Multiarray_d* rhs = ((struct Solver_Volume*)curr)->rhs;

		const int ni = (int)compute_size(rhs->order,rhs->extents);

		PetscInt    ix[ni];
		PetscScalar y[ni];

		for (int i = 0; i < ni; ++i) {
			ix[i] = ind_dof+i;
			y[i]  = -(rhs->data[i]);
		}

		VecSetValues(ssi->b,ni,ix,y,INSERT_VALUES);
	}
}

static void compute_CFL_ramping (struct Solver_Storage_Implicit*const ssi, const struct Simulation*const sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;

	if (test_case->lhs_terms != LHS_CFL_RAMPING)
		return;

	const double max_rhs = compute_max_rhs_dg_like(sim);
	const int n_eq = test_case->n_eq;

	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume* s_vol = (struct Solver_Volume*) curr;
		const double dt = compute_dt_cfl_constrained(max_rhs,s_vol,sim);

		struct DG_Solver_Volume* dg_s_vol = (struct DG_Solver_Volume*) curr;

		/// \note See the comments for this function for the justification of the scaling.
		const struct const_Matrix_d*const m_dt = constructor_copy_scale_const_Matrix_d(dg_s_vol->m,-1.0/dt); // dest.

		for (int eq = 0; eq < n_eq; ++eq) {
			set_petsc_Mat_row_col_dg(ssi,s_vol,eq,s_vol,eq);
			add_to_petsc_Mat(ssi,(struct const_Matrix_d*)m_dt);
		}
		destructor_const_Matrix_d(m_dt);
	}
}

// Level 1 ********************************************************************************************************** //

/** \brief Return the minimum measure of the length of the input volume.
 *  \return See brief. */
static double compute_min_length_measure
	(const struct Solver_Volume*const s_vol, ///< \ref Solver_Volume_T.
	 const struct Simulation*const sim       ///< \ref Simulation.
	);

/** \brief Return the ratio of the initial max_rhs divided by the current max_rhs.
 *  \return See brief. */
static double compute_max_rhs_ratio
	(const double max_rhs ///< The current maximum rhs value.
	);

static void scale_rhs_by_m_inv_std (const struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume*const s_vol       = (struct Solver_Volume*) curr;
		struct DG_Solver_Volume*const dg_s_vol = (struct DG_Solver_Volume*) curr;
		mm_NN1C_overwrite_Multiarray_d(dg_s_vol->m_inv,&s_vol->rhs);
	}
}

static void scale_rhs_by_m_inv_col (const struct Simulation*const sim)
{
	for (struct Intrusive_Link* curr = sim->volumes->first; curr; curr = curr->next) {
		struct Solver_Volume*const s_vol = (struct Solver_Volume*) curr;

		const struct const_Vector_d jac_det_vc = interpret_const_Multiarray_as_Vector_d(s_vol->jacobian_det_vc);
		scale_Multiarray_by_Vector_d('L',1.0,s_vol->rhs,&jac_det_vc,true);
	}
}

static double compute_dt_cfl_constrained
	(const double max_rhs, const struct Solver_Volume*const s_vol, const struct Simulation*const sim)
{
	const struct Test_Case*const test_case = (struct Test_Case*) sim->test_case_rc->tc;

	struct Multiarray_d* s_coef_b = constructor_s_coef_bezier(s_vol,sim); // destructed

	const ptrdiff_t ext_0 = s_coef_b->extents[0];

	// Max wave speed.
	struct Multiarray_d* V_p_c = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){ext_0,1}); // destructed
	compute_max_wavespeed(V_p_c,(struct const_Multiarray_d*)s_coef_b,'c');

	const double max_wave_speed = maximum_dd(V_p_c->data,ext_0);
	destructor_Multiarray_d(V_p_c);

	// Max viscosity.
	double max_viscosity = EPS;
	if (test_case->has_2nd_order) {
		struct Multiarray_d* mu = constructor_empty_Multiarray_d('C',2,(ptrdiff_t[]){ext_0,1}); // destructed
		compute_viscosity(mu,(struct const_Multiarray_d*)s_coef_b,'c');

		max_viscosity = maximum_dd(mu->data,ext_0);
		destructor_Multiarray_d(mu);
	}
	destructor_Multiarray_d(s_coef_b);

	// Min length measure.
	double dx = compute_min_length_measure(s_vol,sim);

	const double max_rhs_ratio = compute_max_rhs_ratio(max_rhs);
	double cfl = test_case->cfl_initial * ( max_rhs_ratio < 1.0 ? 1.0 : max_rhs_ratio );
	if (cfl > test_case->cfl_max)
		cfl = test_case->cfl_max;
//printf("cfl: %f %f %f\n",cfl,dx,max_wave_speed);

	return cfl*GSL_MIN(dx/max_wave_speed,( !test_case->has_2nd_order ? DBL_MAX : dx*dx/max_viscosity ));
}

// Level 2 ********************************************************************************************************** //

static double compute_min_length_measure (const struct Solver_Volume*const s_vol, const struct Simulation*const sim)
{
	/** In the 1D case, the length of the volume is being returned. In higher dimensional cases, the minimum of the
	 *  lengths of each of the edges is returned.
	 *
	 *  Note that the ideal length measure is the distance from the outgoing characteristic wave to the edge of the
	 *  neighbouring volume (relating to the physics of the CFL constraint). The approximation used here is chosen for
	 *  efficiency/convenience as it can be computed locally to each volume but that this would require that a smaller
	 *  CFL number than the maximum allowable for linear stability be used.
	 */
	double min_length_measure = DBL_MAX;

	// sim may be used to store a parameter establishing which type of operator to use for the computation.
	UNUSED(sim);
	const char op_format = 'd';

	const struct Volume*const vol    = (struct Volume*) s_vol;
	const struct const_Element* el   = (struct const_Element*) vol->element;
	const struct Solver_Element* s_e = (struct Solver_Element*) el;

	const bool curved = vol->curved;
	const int p = s_vol->p_ref;
	const struct const_Multiarray_d*const g_coef = s_vol->geom_coef;

	struct Multiarray_d*const v_vals = constructor_empty_Multiarray_d('R',2,(ptrdiff_t[]){2,DIM}); // destructed
	if (el->d == 1) {
		const struct Operator*const cv0_vg_vv =
			( curved ? get_Multiarray_Operator(s_e->cv0_vg_vv[curved],(ptrdiff_t[]){0,0,1,p})
			         : get_Multiarray_Operator(s_e->cv0_vg_vv[curved],(ptrdiff_t[]){0,0,1,1}) );

		mm_NN1_Operator_Multiarray_d(cv0_vg_vv,g_coef,v_vals,op_format,2,NULL,NULL);

		assert(v_vals->extents[0] == 2);
		assert(v_vals->extents[1] == 1);
		min_length_measure = fabs(v_vals->data[0]-v_vals->data[1]);
	} else {
		const int n_e = el->n_e;
		for (int le = 0; le < n_e; ++le) {
			const struct Operator*const cv0_vg_ev =
				( curved ? get_Multiarray_Operator(s_e->cv0_vg_ev[curved],(ptrdiff_t[]){le,0,0,1,p})
				         : get_Multiarray_Operator(s_e->cv0_vg_ev[curved],(ptrdiff_t[]){le,0,0,1,1}) );

			mm_NN1_Operator_Multiarray_d(cv0_vg_ev,g_coef,v_vals,op_format,2,NULL,NULL);

			assert(v_vals->extents[0] == 2);
			double*const data[2] = { get_row_Multiarray_d(0,v_vals), get_row_Multiarray_d(1,v_vals), };
			for (int d = 0; d < DIM; ++d)
				data[0][d] -= data[1][d];

			const double edge_len = norm_d(DIM,data[0],"L2");
			if (edge_len < min_length_measure)
				min_length_measure = edge_len;
		}
	}
	destructor_Multiarray_d(v_vals);
	return min_length_measure;
}

static double compute_max_rhs_ratio (const double max_rhs)
{
	static double max_rhs0 = 0.0;
	if (max_rhs0 == 0.0)
		max_rhs0 = max_rhs;
	return max_rhs0/max_rhs;
}
