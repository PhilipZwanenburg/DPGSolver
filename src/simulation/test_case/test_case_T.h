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
 *  \brief Provides templated container(s) and functions relating to the test cases.
 */

#include <stdbool.h>
#include "definitions_test_case.h"

#include "def_templates_boundary.h"
#include "def_templates_flux.h"
#include "def_templates_geometry.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_solution.h"
#include "def_templates_test_case.h"

struct Simulation;

/** \brief Container for test case specific information.
 *
 *  This container is used to hold test case specific variables (such as those related to the pde under consideration)
 *  as well as function pointers to various functions such that the control flow is not broken during run-time.
 *
 *  As this container is a general template for many supported test cases, it is frequently the case that some of the
 *  function pointers are not necessary. In this case, they are set to point to functions which simply return
 *  immediately.
 *
 *  Regarding the solution related function pointers, different pointers are provided for the starting solution and the
 *  run-time solution to avoid inefficient use of the restarted solution during run-time if applicable. For example,
 *  when using the Riemann invariant boundary condition, the solution is constructed while solving for the boundary
 *  numerical flux and the free-stream solution should be used despite this solution not being used for the
 *  initialization or error computation. If no restart file is being used, the function pointers are set to be equal.
 *
 * \todo Remove: n_var, n_eq, ind_num_flux from Test_Case_T in favor of using the static function for retrieval.
 *       This should eliminate a large number of instances where Simulation/Test_Case_T are passed as function
 *       arguments.
 */
struct Test_Case_T {
	const int pde_index;  ///< Index corresponding to \ref Simulation::pde_name.
	const bool is_linear; ///< Flag for whether the pde under consideration is linear.

	const bool has_1st_order, ///< Flag for whether the pde under consideration has 1st order terms.
	           has_2nd_order; ///< Flag for whether the pde under consideration has 2nd order terms.

	/** The required unknowns for the solver (i.e. which must be computed **before** the solving stage) for the
	 *  \ref Simulation::method.
	 *
	 *  The indices correspond to:
	 *  0: \ref Solver_Volume_T::sol_coef.
	 *  1: \ref Solver_Face_T::nf_coef.
	 *  2: \ref Solver_Volume_T::grad_coef.
	 *  3: \ref Solver_Face_T::s_coef.
	 */
	const bool required_unknowns[MAX_N_UNKNOWNS];

	const int n_var, ///< Number of variables in the PDE under consideration.
	          n_eq;  ///< Number of equations in the PDE under consideration.

	set_sol_fptr_T set_sol;  ///< Function pointer to the function used to set data relating to the solution.
	set_sol_fptr_T set_grad; ///< Function pointer to the function used to set data relating to the solution gradients.

	/** Pointer to function constructing the physical xyz coordinates from either parametric space or from straight to
	 *  curved elements. */
	constructor_xyz_fptr_T constructor_xyz;

	constructor_sol_fptr_T constructor_sol;  ///< Pointer to the function used to construct the solution.
	constructor_sol_fptr_T constructor_grad; ///< Pointer to the function used to construct the solution gradient.

	set_sol_fptr_T         set_sol_start;         ///< Pointer to function used to set the starting solution data.
	constructor_sol_fptr_T constructor_sol_start; ///< Pointer to function used to construct the starting solution data.

	// Geometry related parameters
	/// The surface parametrization type to use for geometry blending. Options: see \ref definitions_geometry.h.
	const int geom_parametrization;

	// Solver related parameters
	char solver_method_curr; ///< The current solver method. Options: 'e'xplicit, 'i'mplicit.

	/// The type of solver procedure to be used for the simulation. Options: See definitions_test_case.h.
	const int solver_proc;

	// Parameters for explicit simulations.
	const int solver_type_e; ///< The explicit solver type. Options: See definitions_test_case.h.

	double time;             ///< The current time.
	const double time_final; ///< The final time.
	const double dt;         ///< The time increment at each stage of the explicit solve.

	// Parameters for implicit simulations.
	/** Flag for whether the Schur complement should be used for the global system solve. This option is available
	 *  whenever it is possible for certain degrees of freedom to be statically condensed out of the global
	 *  system. */
	const bool use_schur_complement;

	// Parameters for explicit/implicit simulations.
	const int solver_type_i; ///< The implicit solver type. Options: See definitions_test_case.h.

	/// Parameter relating to the terms to be included in the LHS matrix. Options: See definitions_test_case.h.
	const int lhs_terms;

	const double cfl_initial; ///< Initial CFL number in case \ref LHS_CFL_RAMPING is selected.

	double cfl_max;           ///< Maximum CFL number in case \ref LHS_CFL_RAMPING is selected.

	/// Integer indices of the 1st/2nd order numerical fluxes. See \ref definitions_test_case.h
	const int ind_num_flux[2];

	const int ind_test_norm; ///< Integer index of the test_norm. See \ref definitions_dpg.h.

	/// Integer index of the conservation enforcement procedure. See \ref definitions_dpg.h.
	const int ind_conservation;

	const double exit_tol_e,   ///< The exit tolerance for the residual during the explicit solver stage.
	             exit_ratio_e, ///< The exit ratio for the residual during the explicit solver stage.
	             exit_tol_i,   ///< The exit tolerance for the residual during the implicit solver stage.
	             exit_ratio_i; ///< The exit ratio for the residual during the implicit solver stage.

	const bool flux_comp_mem_e[MAX_FLUX_OUT],         ///< \ref Flux_Input_T::compute_member (explicit).
	           flux_comp_mem_i[MAX_FLUX_OUT],         ///< \ref Flux_Input_T::compute_member (implicit).
	           boundary_value_comp_mem_e[MAX_BV_OUT], ///< \ref Boundary_Value_Input_T::compute_member (explicit).
	           boundary_value_comp_mem_i[MAX_BV_OUT]; ///< \ref Boundary_Value_Input_T::compute_member (implicit).

	/// Pointer to the function used to call the combination of 1st and 2nd order flux functions.
	compute_Flux_fptr_T compute_Flux;

	/// Pointers to the functions used to compute the 1st order (inviscid) and 2nd order (viscous) fluxes.
	compute_Flux_fptr_T compute_Flux_iv[2];

	/// Function pointer to the function used to call the combination of 1st and 2nd order numerical flux functions.
	compute_Numerical_Flux_fptr_T compute_Numerical_Flux;

/// \todo Merge the numerical flux computation functions (without and with Jacobian) and remove one of these function pointers.
	/// Function pointers to the functions used to compute the 1st/2nd order numerical fluxes for the explicit solver.
	compute_Numerical_Flux_fptr_T compute_Numerical_Flux_e[2];

	/** Function pointers to the functions used to compute the 1st/2nd order numerical fluxes (and optionally
	 *  Jacobians) for the implicit solver. */
	compute_Numerical_Flux_fptr_T compute_Numerical_Flux_i[2];

	/** Constructor for solution and gradient from the left volume at face cubature nodes as seen from the left
	 *  volume. */
	constructor_Boundary_Value_Input_face_fptr_T constructor_Boundary_Value_Input_face_fcl;

	/** \note constructor_*_r_fcl are provided as part of the Solver_Face_* container as they are dependent upon
	 *        whether the face is on a boundary or not. */

	compute_source_rhs_fptr_T compute_source_rhs; ///< Function pointer to the rhs source computation function.

	/// Function pointer to the source flux imbalance computation function.
	compute_source_rhs_fptr_T add_to_flux_imbalance_source;

	constructor_Error_CE_fptr constructor_Error_CE;         ///< Pointer to the function computing the error.

	/// Pointer to the function computing the error when testing restart functionality.
	constructor_Error_CE_fptr constructor_Error_CE_restart_test;

	/// Pointer to the function computing the functional errors for the computational elements.
	constructor_Error_CE_fptr constructor_Error_CE_functionals;

	// Miscellaneous parameters
	const bool display_progress; ///< Flag for whether the solver progress should be displayed (in stdout).
	const bool has_functional;   ///< Flag for whether a functional error should be computed for the test case.
	const bool has_analytical;   ///< Flag for whether the test case has an analytical solution.

	/** Flag for whether the rhs term computed from the initial state should be copied to the \ref Solver_Volume_T\*s.
	 *  This is used to check free-stream preservation, for example. */
	const bool copy_initial_rhs;
};

/** \brief Constructor for a \ref Test_Case_T.
 *  \return See brief. */
struct Test_Case_T* constructor_Test_Case_T
	(const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Test_Case_T.
void destructor_Test_Case_T
	(const struct Test_Case_T* test_case ///< \ref Test_Case_T.
	);

/// \brief Increment input number of pointers by one.
void increment_pointers_T
	(const int n_ptr,       ///< The number of pointers.
	 const Type**const ptrs ///< Pointer to the array of pointers.
	);

#include "undef_templates_boundary.h"
#include "undef_templates_flux.h"
#include "undef_templates_geometry.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_solution.h"
#include "undef_templates_test_case.h"
