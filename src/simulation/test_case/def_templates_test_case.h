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
 *  \brief Provides the macro definitions used for c-style templating related to the \ref Test_Case_T functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Test_Case_T Test_Case
///\}

///\{ \name Function names
#define constructor_Test_Case_T constructor_Test_Case
#define destructor_Test_Case_T  destructor_Test_Case
#define increment_pointers_T    increment_pointers
///\}

///\{ \name Static names
#define set_string_associations set_string_associations
#define set_pde_related set_pde_related
#define set_method_related set_method_related
#define set_nonlinear_solver_related set_nonlinear_solver_related
#define set_function_pointers set_function_pointers
#define read_test_case_parameters read_test_case_parameters
#define correct_invalid_test_case_parameters correct_invalid_test_case_parameters
#define get_compute_member_Flux_Input get_compute_member_Flux_Input
#define get_compute_member_Boundary_Value_Input get_compute_member_Boundary_Value_Input
#define set_function_pointers_start set_function_pointers_start
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Test_Case_T Test_Case_c
///\}

///\{ \name Function names
#define constructor_Test_Case_T constructor_Test_Case_c
#define destructor_Test_Case_T  destructor_Test_Case_c
#define increment_pointers_T    increment_pointers_c
///\}

///\{ \name Static names
#define set_string_associations set_string_associations_c
#define set_pde_related set_pde_related_c
#define set_method_related set_method_related_c
#define set_nonlinear_solver_related set_nonlinear_solver_related_c
#define set_function_pointers set_function_pointers_c
#define read_test_case_parameters read_test_case_parameters_c
#define correct_invalid_test_case_parameters correct_invalid_test_case_parameters_c
#define get_compute_member_Flux_Input get_compute_member_Flux_Input_c
#define get_compute_member_Boundary_Value_Input get_compute_member_Boundary_Value_Input_c
#define set_function_pointers_start set_function_pointers_start_c
///\}

#endif
