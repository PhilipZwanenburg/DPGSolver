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

#include "face_solver.h"

#include "volume.h"

#include "multiarray.h"

#include "boundary_advection.h"
#include "boundary_euler.h"
#include "const_cast.h"
#include "geometry.h"
#include "simulation.h"
#include "test_case.h"

// Static function declarations ************************************************************************************* //

// Templated functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "face_solver_T.c"

// Interface functions ********************************************************************************************** //

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

