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
 *  \brief Provides Vector_\* math functions.
 */

struct Vector_T;
struct Vector_R;

/// \brief Invert each of the entries of the input \ref Vector_T\*.
void invert_Vector_R
	(struct Vector_R* a ///< Input vector.
	);
#if TYPE_RC == TYPE_COMPLEX
/// \brief Invert each of the entries of the input \ref Vector_T\*.
void invert_Vector_T
	(struct Vector_T* a ///< Input vector.
	);
#endif
/// \brief Add to a \ref Vector_T\*.
void add_to_Vector_T_T
	(struct Vector_T* a, ///< To be added to.
	 const Type* b     ///< Data to add.
	);