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

#ifndef DPG__solution_h__INCLUDED
#define DPG__solution_h__INCLUDED
/** \file
 *  \brief Provides the interface to functions used for solution specification (initialization).
 */

struct const_Multiarray_d;
struct Simulation;
struct Solver_Volume;
struct Solver_Face;
struct Solution_Container;

/** \brief Function pointer to a function which constructs the solution given the input xyz coordinates.
 *  \return A \ref const_Multiarray_d\* container holding the solution.
 *
 *  \param xyz \ref Input xyz coordinates.
 *  \param sim \ref Simulation.
 */
typedef const struct const_Multiarray_d* (*constructor_sol_fptr)
	(const struct const_Multiarray_d* xyz,
	 const struct Simulation* sim
	);

/** \brief Function pointer to volume solution setting function.
 *  \param sim      \ref Simulation.
 *  \param sol_cont \ref Solution_Container.
 */
typedef void (*set_sol_fptr)
	(const struct Simulation* sim,
	 struct Solution_Container sol_cont
	);

/** \brief Function pointer to the function setting the source contribution of the rhs term.
 *  \param sim    \ref Simulation.
 *  \param volume \ref Solver_Volume.
 */
typedef void (*compute_source_fptr)
	(const struct Simulation* sim,
	 struct Solver_Volume* volume
	);

/// Container for members relating to the solution computation.
struct Solution_Container {
	const char ce_type,   ///< The type of computational element associated with the solution data being set.
	           cv_type,   ///< The format in which to return the solution. Options: 'c'oefficients, 'v'alues.
	           node_kind; ///< The kind of nodes to be used. Options: 's'olution, 'c'ubature.

	struct Solver_Volume* volume; ///< \ref Solver_Volume.
	struct Solver_Face* face;     ///< \ref Solver_Face.

	struct Multiarray_d* sol; ///< The container for the computed solution.
};

// Interface functions ********************************************************************************************** //

/** \brief Version of \ref constructor_sol_fptr to be used what a call to this function should be invalid.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_const_sol_invalid
	(const struct const_Multiarray_d* xyz, ///< Defined for \ref constructor_sol_fptr.
	 const struct Simulation* sim          ///< Defined for \ref constructor_sol_fptr.
	);

/** \brief Set up the initial solution for the simulation. Sets:
 *	- \ref Solver_Volume::sol_coef;
 *	- \ref Solver_Volume::grad_coef (if applicable);
 *	- \ref Solver_Face::sol_coef    (if applicable);
 *	- \ref Solver_Face::grad_coef   (if applicable).
 */
void set_initial_solution
	(struct Simulation* sim ///< \ref Simulation.
	);

/** \brief Function pointer to be used for \ref Test_Case::set_sol or \ref Test_Case::set_grad when this solution is not
 *         required. */
void set_sg_do_nothing
	(const struct Simulation* sim,      ///< Defined for \ref set_sol_fptr.
	 struct Solution_Container sol_cont ///< Defined for \ref set_sol_fptr.
	);

/** \brief Contructor for a \ref const_Multiarray_d\* holding the xyz coordinates at volume nodes of input kind.
 *  \return See brief. */
const struct const_Multiarray_d* constructor_xyz_v
	(const struct Simulation* sim, ///< \ref Simulation.
	 struct Solver_Volume* volume, ///< \ref Solver_Volume.
	 const char node_kind          ///< The kind of node. Options: 's'olution, 'c'ubature.
	);

/// \brief Compute the coefficients associated with the values of the volume solution.
void compute_coef_from_val_vs
	(const struct Solver_Volume* s_vol,        ///< \ref Solver_Volume.
	 const struct const_Multiarray_d* sol_val, ///< The solution values.
	 struct Multiarray_d* sol_coef             ///< To hold the solution coefficients.
	);

/** \brief Contructor for a \ref Multiarray_d\* holding the solution at volume nodes of input kind.
 *  \return See brief. */
struct Multiarray_d* constructor_sol_v
	(const struct Simulation* sim, ///< \ref Simulation.
	 struct Solver_Volume* volume, ///< \ref Solver_Volume.
	 const char node_kind          ///< The kind of node. Options: 's'olution, 'c'ubature.
	);

/// \brief Function pointer to be used for \ref Test_Case::compute_source when there is no source term.
void compute_source_do_nothing
	(const struct Simulation* sim, ///< Defined for \ref compute_source_fptr.
	 struct Solver_Volume* volume  ///< Defined for \ref compute_source_fptr.
	);

#endif // DPG__solution_h__INCLUDED