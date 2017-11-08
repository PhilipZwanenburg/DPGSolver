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

#ifndef DPG__geometry_h__INCLUDED
#define DPG__geometry_h__INCLUDED
/**	\file
 *	\brief Provides the interface to functions used for geometry processing.
 */

struct Simulation;
struct Intrusive_List;
struct Solver_Volume;
struct Multiarray_d;
struct const_Multiarray_d;

/**	\brief Set up the solver geometry:
 *	- \ref Solver_Volume::metrics_vm;
 *	- \ref Solver_Volume::metrics_vc;
 *	- \ref Solver_Volume::jacobian_det_vc;
 *	- \todo [ref here] Solver_Face::normals_fc;
 *	- \todo [ref here] Solver_Face::jacobian_det_fc;
 *
 *	Requires that:
 *	- \ref Simulation::volumes points to a list of \ref Solver_Volume\*s;
 *	- \ref Simulation::faces   points to a list of \todo [ref here] Solver_Face\*s.
 */
void set_up_solver_geometry
	(struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Compute the face unit normal vectors at the nodes corresponding to the given face metrics.
void compute_unit_normals
	(const int ind_lf,                             ///< Defined for \ref compute_unit_normals_and_det.
	 const struct const_Multiarray_d* normals_ref, ///< Defined for \ref compute_unit_normals_and_det.
	 const struct const_Multiarray_d* metrics_f,   ///< Defined for \ref compute_unit_normals_and_det.
	 struct Multiarray_d* normals_f                ///< Defined for \ref compute_unit_normals_and_det.
	);

#endif // DPG__geometry_h__INCLUDED
