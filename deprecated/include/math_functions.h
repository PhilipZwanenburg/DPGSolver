// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__math_functions_h__INCLUDED
#define DPG__math_functions_h__INCLUDED
/**	\file
 *	\brief Provide several standard math functions.
 *
 *	<!-- References: -->
 *	Press(1992,2nd edition) - Numerical Recipes in C - The Art of Scientific Computing (Ch. 6.1)
 */

#include <stdbool.h>
#include <stddef.h>

/** \brief Evaluates the factorial of the input returning an `double` result.
 *	\return See brief. */
double factorial_d
	(const unsigned int n ///< Input.
	);

/** \brief Evaluates the factorial of the input returning an `long long unsigned int` result.
 *	\return See brief. */
long long unsigned int factorial_ull
	(const unsigned int n ///< Input.
	);

/** \brief Evaluates the gamma function at the input value.
 *	\return See brief. */
double gamma_d
	(const double x ///< Input.
	);

/** \brief Compares input values for approximate equality.
 *	\return `true` if the difference is less than the tolerance. */
bool equal_d
	(const double x0, ///< Input 0.
	 const double x1, ///< Input 1.
	 const double tol ///< The tolerance.
	);

/** \brief Computes the norm of the input `double*` data with the specified norm type.
 *	\return See brief. */
double norm_d
	(const ptrdiff_t n_entries, ///< The number of entries.
	 const double*const data,   ///< The data.
	 const char*const norm_type ///< The norm type. Options: "L2".
	);

#endif // DPG__math_functions_h__INCLUDED
