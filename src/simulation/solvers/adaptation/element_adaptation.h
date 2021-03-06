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

#ifndef DPG__element_adaptation_h__INCLUDED
#define DPG__element_adaptation_h__INCLUDED
/** \file
 *  \brief Provides the interface for the derived \ref Adaptation_Element container and associated functions.
 */

#include "element.h"

struct Simulation;

/// \brief Container for data relating to the dg solver element.
struct Adaptation_Element {
	const struct const_Element element; ///< Base \ref const_Element.

	const struct Multiarray_Operator* cc0_vs_vs; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cc0_vr_vr; ///< See notation in \ref element_operators.h.
	const struct Multiarray_Operator* cc0_ff_ff; ///< See notation in \ref element_operators.h.

	const struct Multiarray_Operator* vv0_vv_vv; ///< See notation in \ref element_operators.h.

	const struct const_Multiarray_Vector_i* nc_ff; ///< Node correspondence for 'f'ace 'f'lux.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructor for a derived \ref Adaptation_Element.
void constructor_derived_Adaptation_Element
	(struct Element* element,     ///< \ref Adaptation_Element.
	 const struct Simulation* sim ///< \ref Simulation.
	);

/// \brief Destructor for a \ref Adaptation_Element.
void destructor_derived_Adaptation_Element
	(struct Element* element ///< Standard.
	);

#endif // DPG__element_adaptation_h__INCLUDED
