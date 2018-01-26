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
 *  \brief Provides templated functions relating to diffusion numerical fluxes.
 */

struct Numerical_Flux_Input_T;
struct mutable_Numerical_Flux_T;

/** \brief Version of \ref compute_Numerical_Flux_fptr_T computing the numerical fluxes using the BR2 method (Table 1,
 *         \cite Brdar2012). */
void compute_Numerical_Flux_T_diffusion_br2
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	);

/** \brief Version of \ref compute_Numerical_Flux_fptr_T computing the numerical fluxes and Jacobians using the BR2
 *         method.
 *  See comments for \ref compute_Numerical_Flux_T_diffusion_br2. */
void compute_Numerical_Flux_T_diffusion_br2_jacobian
	(const struct Numerical_Flux_Input_T* num_flux_i, ///< See brief.
	 struct mutable_Numerical_Flux_T* num_flux        ///< See brief.
	);