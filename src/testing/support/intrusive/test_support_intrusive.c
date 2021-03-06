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
 *  \todo Likely name change to 'test_support_computational_elements'.
 */

#include "test_support_intrusive.h"

#include "test_support_volume.h"
#include "test_support_face.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "definitions_alloc.h"
#include "definitions_intrusive.h"

#include "volume.h"
#include "face.h"
#include "file_processing.h"

// Static function declarations ************************************************************************************* //

/// Container for intrusive list related information.
struct IL_Info {
	int list_name; ///< \ref Intrusive_List::name.
};

/**	\brief Set the \ref IL_Info.
 *	\return A copy of the constructed \ref IL_Info. */
struct IL_Info set_IL_Info
	(const char*const list_name ///< The `char*` list name.
	);

// Interface functions ********************************************************************************************** //

struct Intrusive_List* constructor_file_name_IL
	(const char*const list_name, const char*const file_name, const struct const_Intrusive_List*const elements,
	 const struct Intrusive_List*const volumes)
{
	FILE* file = fopen_checked(file_name); // closed

	struct IL_Info il_info = set_IL_Info(list_name);

	struct Intrusive_List* intrusive_list = constructor_empty_IL(il_info.list_name,NULL);

	bool found_var = false;
	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),file)) {
		if (strstr(line,list_name)) {
			found_var = true;
			if (il_info.list_name == IL_VOLUME)
				push_back_IL(intrusive_list,(struct Intrusive_Link*) constructor_Volume(file,line,elements));
			else if (il_info.list_name == IL_FACE)
				push_back_IL(intrusive_list,(struct Intrusive_Link*) constructor_Face_file(file,line,elements,volumes));
			else
				EXIT_UNSUPPORTED;
		}
	}

	fclose(file);

	if (!found_var)
		EXIT_ERROR("Did not find a '%s' member in the file: %s.",list_name,file_name);

	return intrusive_list;
}

void constructor_file_name_derived_Faces (struct Intrusive_List* faces, const char*const file_name)
{
	const int list_name = faces->name;

	FILE* file = fopen_checked(file_name); // closed

	struct Intrusive_Link* curr = faces->first;
	char var_name[STRLEN_MIN] = { 0, };
	bool found_last = 0;

	char line[STRLEN_MAX];
	while (fgets(line,sizeof(line),file)) {

		struct Face* face = (struct Face*)curr;
		sprintf(var_name,"%s %d","Face",face->index);

		if (strstr(line,var_name)) {
			if (found_last)
				EXIT_ERROR("Already found the last face.");

			switch (list_name) {
			case IL_FACE_SOLVER:
				constructor_file_Solver_Face(file,line,face);
				break;
			default:
				EXIT_ERROR("Unsupported: %d\n",list_name);
				break;
			}

			if (curr == faces->last)
				found_last = true;

			curr = curr->next;
		}
	}
	if (!found_last)
		EXIT_ERROR("Did not find all faces in the file: %s.",file_name);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

struct IL_Info set_IL_Info (const char*const list_name)
{
	struct IL_Info il_info;

	if (strstr(list_name,"Volume")) {
		il_info.list_name = IL_VOLUME;
	} else if (strstr(list_name,"Face")) {
		il_info.list_name = IL_FACE;
	} else {
		EXIT_UNSUPPORTED;
	}

	return il_info;
}
