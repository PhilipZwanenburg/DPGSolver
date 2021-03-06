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
 *  \brief Provides several standard math functions (templated).
 */

#include <stdbool.h>
#include <stddef.h>
#include <complex.h>

#include "def_templates_math_functions.h"

#define POW2(a) ((a)*(a)) ///< General pow with exponent 2.

/** \brief Compares input values for approximate equality using the relative infinity norm.
 *  \return `true` if the difference is less than the tolerance. */
bool equal_T
	(const Type x0, ///< Input 0.
	 const Type x1, ///< Input 1.
	 const Real tol ///< The tolerance.
		);

/** \brief Version of \ref equal_T without specified denominator for the relative norm.
 *  \return See brief. */
bool equal_spec_rel_T
	(const Type x0,  ///< See brief.
	 const Type x1,  ///< See brief.
	 const Real tol, ///< See brief.
	 const Type den  ///< The value of the denominator to be used for the relative norm.
		);

/** \brief Computes the norm of the input `Type*` data with the specified norm type.
 *  \return See brief. */
Type norm_T
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const Type*const data,     ///< The data.
	 const char*const norm_type ///< The norm type. Options: "L2", "Inf".
		);

/** \brief `version` of \ref norm_T computing the 'R'eal norm from an input array of input 'T'ype.
 *  \return See brief. */
Real norm_R_from_T
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const Type*const data,     ///< The data.
	 const char*const norm_type ///< The norm type. Options: "L2", "Inf".
	);

/** \brief Computes the relative norm of the difference between the input `Type*` data with the specified norm type.
 *  \return See brief. */
Type norm_diff_T
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const Type*const data_0,   ///< The data for input 0.
	 const Type*const data_1,   ///< The data for input 1.
	 const char*const norm_type ///< The norm type. Options: "Inf".
		);

/** \brief `version` of \ref norm_diff_T taking one 'R'eal input and one 'T'ype input acting on the real components.
 *  \return See brief. */
Real norm_diff_RT
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const Real*const data_0,   ///< The data for input 0.
	 const Type*const data_1,   ///< The data for input 1.
	 const char*const norm_type ///< The norm type. Options: "Inf".
	);

/** \brief Computes the infinity norm (not relative) of the difference between the input `Type*` data.
 *  \return See brief. */
Type norm_diff_inf_no_rel_T
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const Type*const data_0,   ///< The data for input 0.
	 const Type*const data_1    ///< The data for input 1.
		);

/** \brief Compute the maximum absolute value of the two inputs and return it.
 *  \return See brief. */
Type max_abs_T
	(const Type a, ///< Input 0.
	 const Type b  ///< Input 1.
	);

/// \brief Variation on the axpy (BLAS 1) function: z = y*x + z.
void z_yxpz_T
	(const int n,   ///< The number of entries.
	 const Type* x, ///< Input x.
	 const Type* y, ///< Input y.
	 Type* z        ///< Location to store the sum.
	);

/// \brief Version of \ref z_yxpz_T with input types: `Real`, `Type`, `Type`.
void z_yxpz_RTT
	(const int n,   ///< The number of entries.
	 const Real* x, ///< Input x.
	 const Type* y, ///< Input y.
	 Type* z        ///< Location to store the sum.
	);

/** \brief Compute the average of the input array entries.
 *  \return See brief. */
Type average_T
	(const Type*const data,    ///< The array of data.
	 const ptrdiff_t n_entries ///< The number of entries.
	);

/** \brief Compute the minimum of the input array entries.
 *  \return See brief. */
Type minimum_T
	(const Type*const data,    ///< The array of data.
	 const ptrdiff_t n_entries ///< The number of entries.
	);

/** \brief Compute the maximum of the absolute values of the input array entries.
 *  \return See brief. */
Type maximum_abs_T
	(const Type*const data,    ///< The array of data.
	 const ptrdiff_t n_entries ///< The number of entries.
	);

/** \brief Compute the 'R'eal maximum of the 'T'ype values of the input array entries.
 *  \return See brief. */
Real maximum_RT
	(const Type*const data,    ///< The array of data.
	 const ptrdiff_t n_entries ///< The number of entries.
	);

/// \brief Add the constant to the entries of the array.
void add_to_T
	(Type*const data,          ///< The array of data.
	 const Type c_add,         ///< The 'c'onstant to add to the data.
	 const ptrdiff_t n_entries ///< The number of entries.
	);

/** \brief Return the dot product of the two input arrays.
 *  \return See brief. */
Type dot_T
	(const ptrdiff_t n,  ///< The number of entries.
	 const Type*const a, ///< The 1st input.
	 const Type*const b  ///< The 2nd input.
		);

/** \brief `version` of \ref dot_T return the 'R'eal dot product 'R'eal and 'T'ype input arrays.
 *  \return See brief. */
Real dot_R_from_RT
	(const ptrdiff_t n,  ///< The number of entries.
	 const Real*const a, ///< The 1st input.
	 const Type*const b  ///< The 2nd input.
	);

/** \brief Return the positive maximum value of the two inputs.
 *  \return See brief. */
Type max_abs_real_T
	(const Type a, ///< First value in the comparison.
	 const Type b  ///< Second value in the comparison.
		);

/** \brief Return the positive minimum value of the two inputs.
 *  \return See brief. */
Type min_abs_real_T
	(const Type a, ///< First value in the comparison.
	 const Type b  ///< Second value in the comparison.
		);

#include "undef_templates_math_functions.h"
