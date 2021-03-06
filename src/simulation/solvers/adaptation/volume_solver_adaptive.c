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

#include "volume_solver_adaptive.h"

#include "macros.h"
#include "definitions_adaptation.h"

#include "simulation.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

void constructor_derived_Adaptive_Solver_Volume (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) volume_ptr;
	struct Solver_Volume* s_vol            = (struct Solver_Volume*) volume_ptr;

	a_s_vol->p_ref_prev = s_vol->p_ref;
	initialize_Adaptive_Solver_Volume(a_s_vol);
	UNUSED(sim);
}

void destructor_derived_Adaptive_Solver_Volume (struct Volume* volume_ptr)
{
	struct Adaptive_Solver_Volume* a_s_vol = (struct Adaptive_Solver_Volume*) volume_ptr;
	UNUSED(a_s_vol);
}

void initialize_Adaptive_Solver_Volume (struct Adaptive_Solver_Volume*const a_s_vol)
{
	a_s_vol->adapt_type = ADAPT_NONE;
	a_s_vol->ind_h      = -1;
	a_s_vol->updated              = false;
	a_s_vol->updated_geom_nc      = false;
	a_s_vol->updated_geom_nc_face = false;

	a_s_vol->child_0 = NULL;
	a_s_vol->parent  = NULL;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
