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
 *  \brief Provides the interface to functions used for solution specification (initialization).
 */

#include <stdbool.h>

#include "def_templates_solution.h"
#include "def_templates_flux.h"
#include "def_templates_multiarray.h"
#include "def_templates_face_solver.h"
#include "def_templates_volume_solver.h"

struct const_Multiarray_T;
struct Flux_Input_T;
struct Simulation;
struct Solver_Volume_T;
struct Solver_Face_T;
struct Solution_Container_T;

/** \brief Function pointer to a function which constructs the solution given the input xyz coordinates.
 *  \return A \ref const_Multiarray_T\* container holding the solution.
 *
 *  \param xyz Input xyz coordinates.
 *  \param sim \ref Simulation.
 */
typedef const struct const_Multiarray_T* (*constructor_sol_fptr_T)
	(const struct const_Multiarray_T* xyz,
	 const struct Simulation* sim
	);

/** \brief `mutable` version of \ref constructor_sol_fptr_T.
 *  \return See brief.
 *
 *  \param xyz See brief.
 *  \param sim See brief.
 */
typedef struct Multiarray_T* (*mutable_constructor_sol_fptr_T)
	(const struct const_Multiarray_T* xyz,
	 const struct Simulation* sim
	);

/** \brief Function pointer to volume solution setting function.
 *  \param sim      \ref Simulation.
 *  \param sol_cont \ref Solution_Container_T.
 */
typedef void (*set_sol_fptr_T)
	(const struct Simulation* sim,
	 struct Solution_Container_T sol_cont
	);

/** \brief Function pointer to the function setting the source contribution of the rhs term.
 *  \param sim   \ref Simulation.
 *  \param s_vol \ref Solver_Volume_T.
 *  \param rhs   Memory in which to add the rhs contribution.
 */
typedef void (*compute_source_rhs_fptr_T)
	(const struct Simulation* sim,
	 const struct Solver_Volume_T* s_vol,
	 struct Multiarray_T* rhs
	);

/// Container for members relating to the solution computation.
struct Solution_Container_T {
	char ce_type,   ///< The type of computational element associated with the solution data being set.
	     cv_type,   ///< The format in which to return the solution. Options: 'c'oefficients, 'v'alues.
	     node_kind; ///< The kind of nodes to be used. Options: 's'olution, 'c'ubature.

	bool using_restart; ///< Flag for whether the solution is being computed from a restart file.

	struct Solver_Volume_T* volume; ///< \ref Solver_Volume_T.
	struct Solver_Face_T* face;     ///< \ref Solver_Face_T.

	struct Multiarray_T* sol; ///< The container for the computed solution.
};

// Interface functions ********************************************************************************************** //

/** \brief Version of \ref constructor_sol_fptr_T to be used what a call to this function should be invalid.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_const_sol_invalid_T
	(const struct const_Multiarray_T* xyz, ///< Defined for \ref constructor_sol_fptr_T.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr_T.
	);

/** \brief Set up the initial solution for the simulation. Sets:
 *	- \ref Solver_Volume_T::sol_coef;
 *	- \ref Solver_Volume_T::grad_coef (if applicable);
 *	- \ref Solver_Face_T::nf_coef     (if applicable);
 */
void set_initial_solution_T
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Function pointer to be used for \ref Test_Case_T::set_sol or \ref Test_Case_T::set_grad when this solution is
 *         not required. */
void set_sg_do_nothing_T
	(const struct Simulation* sim,      ///< Defined for \ref set_sol_fptr_T.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr_T.
	);

/** \brief Function pointer to be used for \ref Test_Case_T::set_grad setting the multiarray of the appropriate size to
 *         zero. */
void set_sg_zero_T
	(const struct Simulation*const sim,   ///< Defined for \ref set_sol_fptr_T.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr_T.
	);

/** \brief Contructor for a \ref const_Multiarray_T\* holding the xyz coordinates associated with the
 *         \ref Solution_Container_T.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_xyz_sol_T
	(const struct Simulation* sim, ///< \ref Simulation.
	 const struct Solution_Container_T* sol_cont ///< \ref Solution_Container_T.
	);

/** \brief Contructor for a \ref Multiarray_T\* holding the solution at volume nodes of input kind.
 *  \return See brief. */
struct Multiarray_T* constructor_sol_v_T
	(const struct Simulation* sim, ///< \ref Simulation.
	 struct Solver_Volume_T* s_vol, ///< \ref Solver_Volume_T.
	 const char node_kind           ///< The kind of node. Options: 's'olution, 'c'ubature.
	);

/// \brief Function pointer to be used for \ref Test_Case_T::compute_source_rhs when there is no source term.
void compute_source_rhs_do_nothing_T
	(const struct Simulation* sim,      ///< See brief.
	 const struct Solver_Volume_T* s_vol, ///< See brief.
	 struct Multiarray_T* rhs           ///< See brief.
	);

/// \brief Function pointer to be used for \ref Test_Case_T::add_to_flux_imbalance_source when there is no source term.
void add_to_flux_imbalance_source_do_nothing_T
	(const struct Simulation* sim,      ///< See brief.
	 const struct Solver_Volume_T* s_vol, ///< See brief.
	 struct Multiarray_T* rhs           ///< See brief.
	);

/// \brief Update \ref Solution_Container_T::sol based on the input solution values.
void update_Solution_Container_sol_T
	(struct Solution_Container_T*const sol_cont, ///< Defined for \ref set_sol_fptr_T.
	 struct Multiarray_T*const sol,              ///< The solution values.
	 const struct Simulation*const sim           ///< \ref Simulation.
	);

/// \brief Update \ref Solution_Container_T::sol based on the input solution gradient values.
void update_Solution_Container_grad_T
	(struct Solution_Container_T*const sol_cont, ///< Defined for \ref set_sol_fptr_T.
	 struct Multiarray_T*const grad,             ///< The gradient values.
	 const struct Simulation*const sim           ///< \ref Simulation.
	);

/** \brief Constructor for the xyz coordinates evaluated at the volume cubature nodes using interpolation.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_xyz_vc_interp_T
	(const struct Solver_Volume_T* s_vol, ///< The current volume.
	 const struct Simulation* sim       ///< \ref Simulation.
	);

/// \brief Constructor for \ref Solver_Face_T::nf_coef using the solution of the dominant neighbouring volume.
void constructor_Solver_Face__nf_coef_T
	(struct Solver_Face_T*const s_face, ///< \ref Solver_Face_T.
	 struct Flux_Input_T*const flux_i,  ///< \ref Flux_Input_T.
	 const struct Simulation*const sim, ///< \ref Simulation.
	 const char method                  /**< Method to use to get the coefficients. Options:
	                                     *   - Compute from 'i'nitial solution;
	                                     *   - Compute from 'c'urrent solution.
	                                     */
	);

/// \brief Function to be used for \ref Test_Case_T::set_grad for the test cases with zero gradients.
void set_grad_zero_T
	(const struct Simulation* sim,        ///< Defined for \ref set_sol_fptr_T.
	 struct Solution_Container_T sol_cont ///< Defined for \ref set_sol_fptr_T.
	);

/** \brief Function to be used for \ref Test_Case_T::constructor_grad for test cases with zero gradients.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_const_grad_zero_T
	(const struct const_Multiarray_T* xyz, ///< Defined for \ref constructor_sol_fptr_T.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr_T.
	);

/// \brief Set the \ref Solver_Volume_T::rhs terms to zero.
void set_to_zero_residual_T
	(const struct Simulation*const sim ///< Standard.
	);

/// \brief Set up the initial \ref Solver_Volume_T::test_s_coef and \todo ref Solver_Volume_T::test_g_coef.
void set_initial_v_test_sg_coef_T
	(struct Simulation*const sim ///< \ref Simulation.
	);

/** \brief Contructor for a \ref const_Multiarray_T\* holding the xyz coordinates at volume nodes of input kind.
 *  \return See brief. */
const struct const_Multiarray_T* constructor_xyz_v
	(const struct Simulation*const sim,        ///< Standard.
	 const struct Solver_Volume_T*const s_vol, ///< Standard.
	 const char node_kind,                     ///< The kind of node. Options: 's'olution, 'c'ubature.
	 const bool using_restart                  ///< Flag for whether a restart solution is being used.
		);

#include "undef_templates_solution.h"
#include "undef_templates_flux.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_face_solver.h"
#include "undef_templates_volume_solver.h"
