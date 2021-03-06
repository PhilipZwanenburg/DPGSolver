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

#ifndef DPG__test_base_h__INCLUDED
#define DPG__test_base_h__INCLUDED
/** \file
 *  \brief Provides base functions/structures for general testing.
 */

#include <stdbool.h>
#include <time.h>

#include "definitions_alloc.h"

// Containers ******************************************************************************************************* //

/// Container for test related information.
struct Test_Info {
	char name[STRLEN_MAX]; ///< The test name.

	int n_warn; ///< The number of warnings generated.
};

// Interface functions ********************************************************************************************** //

/// \brief Check whether the test condition is satisfied and abort if not.
void assert_condition
	(const bool cond ///< The condition to check.
	);

/// \brief Call \ref assert_condition, printing the string to the terminal before aborting if relevant.
void assert_condition_message
	(const bool cond,     ///< The condition to check.
	 const char* cond_str ///< A string to be printed to the terminal if the condition is not satisfied.
	);

/// \brief Check whether the test conditions is satisfied and printing the condition string to the terminal if not.
void expect_condition
	(const bool cond,     ///< The condition to check.
	 const char* cond_str ///< A string to be printed to the terminal if the condition is not satisfied.
	);

/// \brief Print a warning and increment \ref Test_Info::n_warn.
void test_print_warning
	(struct Test_Info*const test_info, ///< \ref Test_Info.
	 const char*const warn_name        ///< The warning message string.
	);

/// \brief Output the warning count if warnings were generated.
void output_warning_count
	(struct Test_Info*const test_info ///< \ref Test_Info.
	);

/** \brief Sets the name of the data file in a `static char` (no free necessary).
 *  \return See brief. */
const char* set_data_file_name_unit
	(const char*const file_name_spec /**< The specific name of the data file (`file_name` without extension but with
	                                  *   relative path from '$PROJECT_INPUT_DIR/testing/unit'). */
	);

/** \brief Sets the name of the integration test data file.
 *  \return See brief. */
const char* set_data_file_name_integration
	(const char*const ctrl_name,    ///< The name of the control file.
	 const char*const int_test_type ///< The type of integration test.
	);

#endif // DPG__test_base_h__INCLUDED
