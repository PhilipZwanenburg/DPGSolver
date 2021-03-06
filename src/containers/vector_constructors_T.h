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
 *  \brief Provides Vector_\* constructors and destructors.
 */

#include <stddef.h>
#include <stdbool.h>
#include <stdio.h>

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

struct Vector_R;
struct Matrix_T;
struct Multiarray_T;
struct const_Vector_R;
struct const_Vector_T;
struct const_Matrix_T;
struct const_Multiarray_T;

// Default constructors ********************************************************************************************* //

/** \brief Constructs a default \ref Vector_T\*.
 *  \return Standard. */
struct Vector_T* constructor_default_Vector_T ();

/** \brief Constructor for a default \ref const_Vector_T\*.
 *  \return Standard. */
const struct const_Vector_T* constructor_default_const_Vector_T ();

/** \brief Constructs a default \ref Vector_T\*\*.
 *  \return Standard. */
struct Vector_T** constructor_default_Vector_T_2
	(const ptrdiff_t n_dest ///< The number of \ref Vector_T\* components.
	);

// Empty constructors *********************************************************************************************** //

/** \brief Constructs an empty \ref Vector_T\*.
 *  \return Standard. */
struct Vector_T* constructor_empty_Vector_T
	(const ptrdiff_t ext_0 ///< Defined in \ref Vector_T.
	);

// Zero constructors ************************************************************************************************ //

/** \brief The same as \ref constructor_empty_Vector_T but allocating using calloc.
 *  \return Standard. */
struct Vector_T* constructor_zero_Vector_T
	(const ptrdiff_t ext_0 ///< Defined for \ref constructor_empty_Vector_T.
	);

// Copy constructors ************************************************************************************************ //

/** \brief Copy constructor for a \ref Vector_T\* from a `Vector_T*`.
 *  \return Standard. */
struct Vector_T* constructor_copy_Vector_T
	(const struct Vector_T*const src ///< The source data.
	);

/** \brief `const` version of \ref constructor_copy_Vector_T.
 *  \return Standard. */
const struct const_Vector_T* constructor_copy_const_Vector_T
	(const struct const_Vector_T*const src ///< See brief.
	);

/** \brief Copy constructor for a \ref Vector_T\* from a `const Type*`.
 *  \return Standard. */
struct Vector_T* constructor_copy_Vector_T_T
	(const ptrdiff_t ext_0,      ///< The value of ext_0.
	 const Type*const data_src ///< The source data.
	);

/** \brief `const` version of \ref constructor_copy_Vector_T_T.
 *  \return Standard. */
const struct const_Vector_T* constructor_copy_const_Vector_T_T
	(const ptrdiff_t ext_0,      ///< See brief.
	 const Type*const data_src ///< See brief.
	);

/** \brief Copy constructor for a \ref Vector_T\* from a \ref Vector_T\*.
 *  \return See brief. */
struct Vector_T* constructor_copy_Vector_T_Vector_R
	(struct Vector_R* src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_copy_Vector_T_Vector_R.
 *  \return See brief. */
const struct const_Vector_T* constructor_copy_const_Vector_T_Vector_R
	(const struct const_Vector_R* src ///< See brief.
	);

// Move constructors ************************************************************************************************ //

/** \brief Move constructor for a \ref Vector_T\* from a `Type*`.
 *  \return Standard. */
struct Vector_T* constructor_move_Vector_T_T
	(const ptrdiff_t ext_0, ///< Standard.
	 const bool owns_data,  ///< Standard.
	 Type*const data      ///< Standard.
	);

/** \brief Move constructor for a \ref const_Vector_T\* from a `const Type*`.
 *  \return Standard. */
struct const_Vector_T* constructor_move_const_Vector_T_T
	(const ptrdiff_t ext_0, ///< The value of ext_0.
	 const bool owns_data,  ///< Standard.
	 const Type*const data   ///< Standard.
	);

/** \brief Move constructor for a \ref const_Vector_T\* from a row of a \ref const_Matrix_T\*.
 *  \return Standard. */
const struct const_Vector_T* constructor_move_const_Vector_Matrix_row_T
	(const int row,                    ///< The row index.
	 const struct const_Matrix_T* src, ///< The source matrix.
	 const int owns_data               ///< Standard.
	);

/** \brief Move constructor for a \ref Vector_T\* from a \ref Matrix_T\*.
 *  \return Standard. */
struct Vector_T* constructor_move_Vector_T_Matrix_T
	(struct Matrix_T* src ///< The source matrix.
	);

/// \brief Move constructor for a `const` \ref const_Vector_T `*const`.
void const_constructor_move_Vector_T
	(const struct const_Vector_T*const* dest, ///< Destination.
	 struct Vector_T* src                     ///< Source.
	);

/// \brief Move constructor for a `const` \ref const_Vector_T `*const` from a `const` \ref const_Vector_T `*const`.
void const_constructor_move_const_Vector_T
	(const struct const_Vector_T*const* dest, ///< Destination.
	 const struct const_Vector_T* src         ///< Source.
	);

// Set constructors ************************************************************************************************* //
#ifdef TYPE_RC
/** \brief Constructor for a \ref Vector_T\* from a sub range of a \ref Multiarray_T\*.
 *  \return Standard. */
struct Vector_T* constructor_set_Vector_T_Multiarray_T
	(struct Multiarray_T* src,         ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/** \brief `const` version of \ref constructor_set_Vector_T_Multiarray_T.
 *  \return Standard. */
const struct const_Vector_T* constructor_set_const_Vector_T_Multiarray_T
	(const struct const_Multiarray_T* src, ///< Defined for \ref constructor_set_Vector_T_Multiarray_T.
	 const ptrdiff_t*const sub_indices     ///< Defined for \ref constructor_set_Vector_T_Multiarray_T.
	);

// Special constructors ********************************************************************************************* //

/** \brief Constructor for a \ref Vector_T\* as the inverse of the source vector.
 *  \return See brief. */
struct Vector_T* constructor_inverse_Vector_T
	(const struct const_Vector_T* src ///< The source vector.
	);

/** \brief `const` version of \ref constructor_inverse_Vector_T.
 *  \return See brief. */
const struct const_Vector_T* constructor_inverse_const_Vector_T
	(const struct const_Vector_T* src ///< See brief.
	);

/** \brief Constructor for a \ref const_Vector_T\* which is defined as the input vector appended to itself the input
 *         number of times.
 *  \return See brief.
 */
const struct const_Vector_T* constructor_repeated_const_Vector_T
	(const Type alpha,                    ///< Scaling constant.
	 const struct const_Vector_T*const a, ///< The input vector.
	 const int n_repeated                 ///< The number of times the sub-vector should be repeated.
	);

/** \brief Constructor for a \ref const_Vector_T\* using an element-wise multiplication of each entry of the inputs.
 *  \return See brief. */
const struct const_Vector_T* constructor_dot_mult_const_Vector_T
	(const Type alpha,               ///< Scaling constant.
	 const struct const_Vector_T* a, ///< The 1st input.
	 const struct const_Vector_T* b, ///< The 2nd input.
	 const int n_repeated            ///< The number of times the sub-vector should be repeated.
		);

/** \brief Version of \ref constructor_dot_mult_const_Vector_T taking 'R'eal and 'T'ype inputs.
 *  \return See brief. */
const struct const_Vector_T* constructor_dot_mult_const_Vector_T_RT
	(const Type alpha,               ///< See brief.
	 const struct const_Vector_R* a, ///< See brief.
	 const struct const_Vector_T* b, ///< See brief.
	 const int n_repeated            ///< See brief.
	);

/** \brief Version of \ref constructor_dot_mult_const_Vector_T which inverts the entries of the 2nd input while
 *         performing the dot product.
 *  \return See brief. */
const struct const_Vector_T* constructor_dot_mult_inverse_2nd_const_Vector_T
	(const Type alpha,                    ///< See brief.
	 const struct const_Vector_T*const a, ///< See brief.
	 const struct const_Vector_T*const b, ///< The 2nd input (inverted in dot product).
	 const int n_repeated                 ///< See brief.
		);

/** \brief Constructs a \ref Vector_T\* as the sum of the entries of the two inputs scaled by the input constants.
 *  \return See brief. */
struct Vector_T* constructor_sum_Vectors_Vector_T
	(const Type alpha_0,          ///< The 1st scaling constant.
	 struct Vector_T*const src_0, ///< The 1st source vector.
	 const Type alpha_1,          ///< The 2nd scaling constant.
	 struct Vector_T*const src_1  ///< The 2nd source vector.
	);

/** \brief `const` version of \ref constructor_sum_Vectors_Vector_T.
 *  \return See brief. */
const struct const_Vector_T* constructor_sum_Vectors_const_Vector_T
	(const Type alpha_0,                      ///< See brief.
	 const struct const_Vector_T*const src_0, ///< See brief.
	 const Type alpha_1,                      ///< See brief.
	 const struct const_Vector_T*const src_1  ///< See brief.
	);

/** \brief Constructs a \ref Vector_T\* as the sum of the rows/columns of the input Matrix in the specified direction.
 *  \return Standard. */
struct Vector_T* constructor_sum_Vector_T_const_Matrix_T
	(const char sum_dir,                   ///< The direction in which to sum the entries. Options: 'R'ow, 'C'olumn.
	 const struct const_Matrix_T*const src ///< The source matrix.
	);

/** \brief `const` version of \ref constructor_sum_Vector_T_const_Matrix_T.
 *  \return Standard. */
const struct const_Vector_T* constructor_sum_const_Vector_T_const_Matrix_T
	(const char sum_dir,                   ///< Defined for \ref constructor_sum_Vector_T_const_Matrix_T.
	 const struct const_Matrix_T*const src ///< Defined for \ref constructor_sum_Vector_T_const_Matrix_T.
	);

/** \brief Constructor for a \ref Vector_T\* from a matrix-vector multiplication.
 *  \return Standard. */
struct Vector_T* constructor_mv_Vector_T
	(const char trans_a_i,                ///< Defined for \ref mv_T.
	 const Type alpha,                    ///< Defined for \ref mv_T.
	 const struct const_Matrix_T*const a, ///< Defined for \ref mv_T.
	 const struct const_Vector_T*const b  ///< Defined for \ref mv_T.
	);

/** \brief `const` version of \ref constructor_mv_Vector_T.
 *  \return Standard. */
const struct const_Vector_T* constructor_mv_const_Vector_T
	(const char trans_a_i,                ///< See brief.
	 const Type alpha,                    ///< See brief.
	 const struct const_Matrix_T*const a, ///< See brief.
	 const struct const_Vector_T*const b  ///< See brief.
	);

/** \brief Constructor for a \ref Vector_T\* from the solution of a linear system using [LAPACKE_dsgesv][dsgesv].
 *  \return Standard.
 *
 *  <!-- References: -->
 *  [dsgesv]: https://software.intel.com/en-us/mkl-developer-reference-c-gesv
 */
struct Vector_T* constructor_sgesv_Vector_T
	(struct Matrix_T* A_i, ///< The LHS input matrix A.
	 struct Vector_T* B_i  ///< The RHS input vector B.
	);

/** \brief `const` version of \ref constructor_sgesv_Vector_T.
 *  \return Standard. */
const struct const_Vector_T* constructor_sgesv_const_Vector_T
	(const struct const_Matrix_T* A_i, ///< The LHS input matrix A.
	 const struct const_Vector_T* B_i  ///< The RHS input vector B.
	);

/// \brief Set a \ref Vector_T\* from a sub range of a \ref Matrix_T\*.
void set_Vector_from_Matrix_T
	(struct Vector_T* dest,            ///< The destination.
	 struct Matrix_T* src,             ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Vector_from_Matrix_T.
void set_const_Vector_from_Matrix_T
	(const struct const_Vector_T* dest, ///< Defined for \ref set_Vector_from_Matrix_T.
	 const struct const_Matrix_T* src,  ///< Defined for \ref set_Vector_from_Matrix_T.
	 const ptrdiff_t*const sub_indices  ///< Defined for \ref set_Vector_from_Matrix_T.
	);

/// \brief Set a \ref Vector_T\* from a sub range of a \ref Multiarray_T\*.
void set_Vector_from_Multiarray_T
	(struct Vector_T* dest,            ///< The destination.
	 struct Multiarray_T* src,         ///< The source.
	 const ptrdiff_t*const sub_indices ///< The sub-indices used to specify which part of the source to extract.
	);

/// \brief `const` version of \ref set_Vector_from_Multiarray_T.
void set_const_Vector_from_Multiarray_T
	(const struct const_Vector_T* dest,    ///< Defined for \ref set_Vector_from_Multiarray_T.
	 const struct const_Multiarray_T* src, ///< Defined for \ref set_Vector_from_Multiarray_T.
	 const ptrdiff_t*const sub_indices     ///< Defined for \ref set_Vector_from_Multiarray_T.
	);

#endif
// File constructors ************************************************************************************************ //

#if (defined(TYPE_RC) && TYPE_RC == TYPE_REAL) || !defined(TYPE_RC)
/** \brief Constructor for a \ref Vector_T\* from data in the input file of the given name.
 *  \return Standard. */
struct Vector_T* constructor_file_name_Vector_T
	(const char*const var_name,      ///< The name of the variable to be read in from the file.
	 const char*const file_name_full ///< The name of the file (including the full path).
	);

/** \brief `const` version of \ref constructor_file_name_Vector_T.
 *  \return Standard. */
const struct const_Vector_T* constructor_file_name_const_Vector_T
	(const char*const var_name,      ///< See brief.
	 const char*const file_name_full ///< See brief.
	);

/** \brief Constructor for a \ref Vector_T\* from the current line in the input file.
 *  \return Standard. */
struct Vector_T* constructor_file_Vector_T
	(FILE* data_file,           ///< The pointer to the file from which to read the data.
	 const bool check_container ///< Flag for whether the container type should be checked.
	);

/** \brief `const` version of \ref constructor_file_Vector_T.
 *  \return See brief. */
const struct const_Vector_T* constructor_file_const_Vector_T
	(FILE* data_file,           ///< See brief.
	 const bool check_container ///< See brief.
	);

#endif

// Destructors ****************************************************************************************************** //

/// \brief Destructs a \ref Vector_T\*.
void destructor_Vector_T
	(struct Vector_T* a ///< Standard.
	);

/// \brief `const` version of \ref destructor_Vector_T.
void destructor_const_Vector_T
	(const struct const_Vector_T* a ///< Standard.
	);

/// \brief Destructs a \ref Vector_T\* if it is not `NULL`.
void destructor_conditional_Vector_T
	(struct Vector_T* a ///< Standard.
	);

/// \brief `const` version of \ref destructor_conditional_Vector_T.
void destructor_conditional_const_Vector_T
	(const struct const_Vector_T* a ///< Standard.
	);

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"
