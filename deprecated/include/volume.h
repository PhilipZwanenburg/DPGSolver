// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Volume_h__INCLUDED
#define DPG__Volume_h__INCLUDED
/**	\file
 *	Provides the interface for the base \ref Volume container and associated functions.
 *
 *	A \ref Volume is a `d` dimensional finite element.
 */

#include <stdbool.h>
#include "constants_elements.h"
#include "intrusive.h"
#include "simulation.h"
#include "mesh.h"

// ****************************************************************************************************************** //

/// \brief Container for data relating to the base Volumes.
struct Volume {
	struct Intrusive_Link lnk; ///< The \ref Intrusive_Link.

	const bool boundary, ///< Flag for whether the volume is on a domain boundary.
	           curved;   ///< Flag for whether the volume is curved.

	const struct const_Matrix_d*const xyz_ve; ///< The xyz coordinates of the volume vertices.

	const struct Face*const faces[NFMAX][NSUBFMAX]; ///< Array of pointers to the neighbouring \ref Face containers.

	const struct const_Element*const element; ///< Pointer to the associated \ref const_Element.
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs the base \ref Volume \ref Intrusive_List.
struct Intrusive_List* constructor_Volume_List
	(struct Simulation*const sim, ///< The \ref Simulation.
	 const struct Mesh*const mesh ///< The \ref Mesh.
	);

/// \brief Destructs the base \ref Volume \ref Intrusive_List.
void destructor_Volumes
	(struct Intrusive_List* Volumes ///< Standard.
	);

// Helper functions ************************************************************************************************* //

/** \brief Check if a sufficient number of vertices (2) satisfy the condition on a domain boundary.
 *
 *	When two vertices satisfy the condition, this implies that a volume edge satisfies the condition and the
 *	qualification is passed to the associated volume.
 *
 *	This function can be used to check:
 *		- faces of volumes: `ve_inds` = volume vertex indices, `f_ve` = \ref Volume::element;
 *		- edges of faces:   `ve_inds` = face vertex indices,   `f_ve` = \ref Face::element.
 *
 *	Currently used to check:
 *		- curved (curved_only = `true`);
 *		- boundary (curved only = `false`).
 *
 *	\return See brief. */
bool check_ve_condition
	(const struct const_Multiarray_Vector_i*const f_ve,  ///< Defined in \ref Element.
	 const struct const_Vector_i*const ve_inds,          ///< The vertex indices for the volume.
	 const struct const_Vector_i*const ve_condition,     ///< The vertex condition to check for.
	 const struct const_Multiarray_Vector_i*const ve_bc, ///< Defined in \ref Mesh_Vertices.
	 const bool curved_only                              /**< Flag indicating whether only curved boundary conditions
	                                                      *   should be considered. */
	);

#endif // DPG__Volume_h__INCLUDED
