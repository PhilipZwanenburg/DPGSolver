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

#ifndef DPG__multiarray_print_h__INCLUDED
#define DPG__multiarray_print_h__INCLUDED
/** \file
 *  \brief Provides Multiarray_\* printing functions.
 */

struct Multiarray_d;
struct Multiarray_Vector_i;
struct Multiarray_Matrix_d;
struct const_Multiarray_d;
struct const_Multiarray_Vector_i;
struct const_Multiarray_Matrix_d;

/// \brief Print a \ref Multiarray_Vector_i\* to the terminal.
void print_Multiarray_Vector_i
	(const struct Multiarray_Vector_i*const a ///< Standard.
	);

/// \brief `const` version of \ref print_Multiarray_Vector_i.
void print_const_Multiarray_Vector_i
	(const struct const_Multiarray_Vector_i*const a ///< Defined for \ref print_Multiarray_Vector_i.
	);

/// \brief Print a \ref Multiarray_d\* to the terminal displaying entries below the tolerance as 0.0.
void print_Multiarray_d
	(const struct Multiarray_d*const a, ///< Standard.
	 const double tol                   ///< The tolerance.
	);

/// \brief `const` version of \ref print_Multiarray_d.
void print_const_Multiarray_d
	(const struct const_Multiarray_d*const a, ///< Defined for \ref print_Multiarray_d.
	 const double tol                         ///< Defined for \ref print_Multiarray_d.
	);

/// \brief Print a \ref Multiarray_Matrix_d\* to the terminal.
void print_Multiarray_Matrix_d
	(const struct Multiarray_Matrix_d*const a, ///< Standard.
	 const double tol                          ///< The tolerance.
	);

/// \brief `const` version of \ref print_Multiarray_Matrix_d.
void print_const_Multiarray_Matrix_d
	(const struct const_Multiarray_Matrix_d*const a, ///< Defined for \ref print_Multiarray_Matrix_d.
	 const double tol                                ///< Defined for \ref print_Multiarray_Matrix_d.
	);

#endif // DPG__multiarray_print_h__INCLUDED