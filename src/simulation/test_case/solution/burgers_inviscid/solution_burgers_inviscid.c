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

#include "solution_burgers_inviscid.h"

#include "macros.h"
#include "definitions_solution.h"

#include "boundary.h"
#include "compute_error.h"
#include "compute_error_burgers_inviscid.h"
#include "const_cast.h"
#include "file_processing.h"
#include "flux_burgers_inviscid.h"
#include "numerical_flux_burgers_inviscid.h"
#include "simulation.h"
#include "solution.h"
#include "test_case.h"

#include "trigonometric/solution_trigonometric.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "solution_burgers_inviscid_T.c"
#include "undef_templates_type.h"

#include "def_templates_type_dc.h"
#include "solution_burgers_inviscid_T.c"
#include "undef_templates_type.h"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
