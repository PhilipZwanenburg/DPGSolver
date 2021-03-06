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
 *  \brief Provides the macro definitions used for c-style templating related to the \ref Solver_Volume_T
 *         containers/functions.
 */

#if TYPE_RC == TYPE_REAL

///\{ \name Data types
#define Solver_Volume_T Solver_Volume
///\}

///\{ \name Function names
#define constructor_derived_Solver_Volume_T constructor_derived_Solver_Volume
#define destructor_derived_Solver_Volume_T  destructor_derived_Solver_Volume
#define get_operator__w_vc__s_e_T           get_operator__w_vc__s_e
#define constructor_mass_T                  constructor_mass
#define constructor_inverse_mass_T          constructor_inverse_mass
#define constructor_l2_proj_operator_s_T constructor_l2_proj_operator_s
#define get_operator__cv0_vs_vs_T get_operator__cv0_vs_vs_d
#define get_operator__cv0_vr_vs_T get_operator__cv0_vr_vs_d
#define get_operator__cv0_vg_vs_T get_operator__cv0_vg_vs_d
#define get_operator__cv0_vs_vc_T get_operator__cv0_vs_vc
#define get_operator__cv0_vt_vc_T get_operator__cv0_vt_vc
#define get_operator__cv0_vr_vc_T get_operator__cv0_vr_vc
#define get_operator__tw1_vt_vc_T get_operator__tw1_vt_vc
#define get_operator__cv1_vt_vc_T get_operator__cv1_vt_vc
///\}

///\{ \name Static names
#define set_function_pointers_constructor_xyz_surface set_function_pointers_constructor_xyz_surface
///\}

#elif TYPE_RC == TYPE_COMPLEX

///\{ \name Data types
#define Solver_Volume_T Solver_Volume_c
///\}

///\{ \name Function names
#define constructor_derived_Solver_Volume_T constructor_derived_Solver_Volume_c
#define destructor_derived_Solver_Volume_T  destructor_derived_Solver_Volume_c
#define get_operator__w_vc__s_e_T           get_operator__w_vc__s_e_c
#define constructor_mass_T                  constructor_mass_c
#define constructor_inverse_mass_T          constructor_inverse_mass_c
#define constructor_l2_proj_operator_s_T constructor_l2_proj_operator_s_c
#define get_operator__cv0_vs_vs_T get_operator__cv0_vs_vs_c
#define get_operator__cv0_vr_vs_T get_operator__cv0_vr_vs_c
#define get_operator__cv0_vg_vs_T get_operator__cv0_vg_vs_c
#define get_operator__cv0_vs_vc_T get_operator__cv0_vs_vc_c
#define get_operator__cv0_vt_vc_T get_operator__cv0_vt_vc_c
#define get_operator__cv0_vr_vc_T get_operator__cv0_vr_vc_c
#define get_operator__tw1_vt_vc_T get_operator__tw1_vt_vc_c
#define get_operator__cv1_vt_vc_T get_operator__cv1_vt_vc_c
///\}

///\{ \name Static names
#define set_function_pointers_constructor_xyz_surface set_function_pointers_constructor_xyz_surface_c
///\}

#endif
