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
 *  \brief Provides the interface for the templated \ref DG_Solver_Face container and associated functions.
 *
 *  These faces are needed by the 'D'iscontinuous 'G'alerkin solver functions.
 */

/// \brief Container for data relating to the DG solver faces.
struct DG_Solver_Face_T {
	struct Solver_Face_T face; ///< The base \ref Solver_Face.

	// Members required for 2nd order PDE terms.

	/// The face contributions to the solution gradient coefficients in each of the neighbouring volumes.
	struct Multiarray_T* grad_coef_f[2];
};

/// \brief Constructor for a derived \ref DG_Solver_Face.
void constructor_derived_DG_Solver_Face_T
	(struct Face* face_ptr,       ///< Pointer to the face.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a derived \ref DG_Solver_Face.
void destructor_derived_DG_Solver_Face_T
	(struct Face* face_ptr ///< Pointer to the face.
	);