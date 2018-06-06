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

#ifndef DPG__compute_face_rlhs_opg_h__INCLUDED
#define DPG__compute_face_rlhs_opg_h__INCLUDED
/** \file
 *  \brief Provides real functions used for computing the face contributions to the right and left-hand side (rlhs)
 *         terms of the OPG scheme.
 */

#include "def_templates_type_d.h"
#include "def_templates_compute_face_rlhs_opg.h"
#include "def_templates_numerical_flux.h"
#include "def_templates_face_solver_opg.h"
#include "compute_face_rlhs_opg_T.h"
#include "undef_templates_type.h"
#include "undef_templates_compute_face_rlhs_opg.h"
#include "undef_templates_numerical_flux.h"
#include "undef_templates_face_solver_opg.h"

/** \brief Update the values of \ref Solver_Face_T::nf_coef based on the updated \ref Solver_Volume_T::test_s_coef
 *         values. */
void update_coef_nf_f_opg
	(const struct Simulation*const sim ///< Standard.
	);

#endif // DPG__compute_face_rlhs_opg_h__INCLUDED