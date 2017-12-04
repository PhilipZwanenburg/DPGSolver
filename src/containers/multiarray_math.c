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

#include "multiarray_math.h"

#include <assert.h>

#include "macros.h"

#include "matrix.h"
#include "multiarray.h"
#include "vector.h"

#include "math_functions.h"

// Static function declarations ************************************************************************************* //

// Interface functions ********************************************************************************************** //

#include "def_templates_type_d.h"
#include "def_templates_matrix_d.h"
#include "def_templates_multiarray_d.h"
#include "def_templates_vector_d.h"
#include "multiarray_math_T.c"
#include "undef_templates_type.h"
#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //
