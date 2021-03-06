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
 *  \brief Undefine macro definitions for c-style templated containers/functions relating to solution for the
 *         diffusion equation (test case: default_steady).
 */

#undef set_sol_diffusion_default_steady_T
#undef set_grad_diffusion_default_steady_T
#undef constructor_const_sol_diffusion_default_steady_T
#undef constructor_const_grad_diffusion_default_steady_T
#undef compute_source_rhs_diffusion_default_steady_T
#undef add_to_flux_imbalance_source_diffusion_default_steady_T

#undef constructor_sol_diffusion_default_steady
#undef constructor_grad_diffusion_default_steady
#undef constructor_source_diffusion_default_steady
#undef constructor_sol_diffusion_default_steady_1d
#undef constructor_grad_diffusion_default_steady_1d
#undef constructor_source_diffusion_default_steady_1d
#undef constructor_sol_diffusion_default_steady_2d
#undef constructor_grad_diffusion_default_steady_2d
#undef constructor_source_diffusion_default_steady_2d
#undef Sol_Data__dd
#undef get_sol_data
#undef read_data_default_diffusion
