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
 *  \brief Undefine macro definitions for c-style templated containers/functions relating to \ref Test_Case_T.
 */

///\{ \name Data types
#undef Test_Case_T
///\}

///\{ \name Function names
#undef constructor_Test_Case_T
#undef destructor_Test_Case_T
#undef increment_pointers_T
///\}

#undef set_string_associations
#undef set_pde_related
#undef set_method_related
#undef set_nonlinear_solver_related
#undef set_function_pointers
#undef read_test_case_parameters
#undef correct_invalid_test_case_parameters
#undef get_compute_member_Flux_Input
#undef get_compute_member_Boundary_Value_Input
#undef set_function_pointers_start
