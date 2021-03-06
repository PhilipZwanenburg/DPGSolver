// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__Face_h__INCLUDED
#define DPG__Face_h__INCLUDED
/**	\file
 *	Provides the interface for the base \ref Face container and associated functions.
 *
 *	A \ref Face is a `d-1` dimensional finite element which is found on the face of a \ref Volume.
 */

#include <stdbool.h>
#include "intrusive.h"
#include "simulation.h"
#include "mesh.h"

/// \brief Container for data relating to the base Faces.
struct Face {
	struct Intrusive_Link lnk; ///< \ref Intrusive_Link.

	const bool boundary, ///< Flag for whether the face is on a domain boundary.
	           curved;   ///< Flag for whether the face is curved.

	const int bc;        ///< The boundary condition associated with the face (if relevant).

/// \todo Add this member to Face_Solver.
//	const char cub_type; /**< Type of cubature to be used for the face. Options: 's'traight, 'c'urved. This should be
//	                      *   curved whenever an adjacent volume is curved. */

	const struct const_Element*const element; ///< Pointer to the associated \ref const_Element.

	/** \brief Container for information relating to the neighbouring \ref Volume on either side of the face.
	 *
	 *	The information for the first index `neigh_info[0]` relates to the Volume whose outward normal vector on the
	 *	current face coincides with that stored as part of the Face; this is generally referred to as the left volume.
	 */
	struct Neighbour_Info {
		int ind_lf,   ///< Local face index in relation to the neighbouring volume.
		    ind_href, /**< Local face h-refinement index.
		               *   Range: [0:NFREFMAX).
		               *   In the case of a conforming mesh (h-refinement disabled), ind_href = 0. */
		    ind_sref, /**< Local sub-face h-refinement index.
		               *   Range: [0:NFSUBFMAX).
		               *
		               *   Example:
		               *
		               *   For a Face of type QUAD, there are nine possible values for ind_sref, related to the allowed
		               *   h-refinements:
		               *   	None       (1 -> 1): ind_href = [0]
		               *   	Isotropic  (1 -> 4): ind_href = [1,4]
		               *   	Horizontal (1 -> 2): ind_href = [5,6]
		               *   	Vertical   (1 -> 2): ind_href = [7,8]
		               *
		               *   	ind_href = 0 : Indsfh = 0 (This is a conforming Face)
		               *   	           1 : Indsfh = 0 (This is the first  sub-face of the h-refined macro Face)
		               *   	           2 : Indsfh = 1 (This is the second sub-face of the h-refined macro Face)
		               *   	           3 : Indsfh = 2 (This is the third  sub-face of the h-refined macro Face)
		               *   	           4 : Indsfh = 3 (This is the fourth sub-face of the h-refined macro Face)
		               *   	           5 : Indsfh = 0
		               *   	           6 : Indsfh = 1
		               *   	           7 : Indsfh = 0
		               *   	           8 : Indsfh = 1
		               */
		    ind_ord;  /**< The local face ordering index.
		               *   Specifies the index for the ordering between Faces as seen from adjacent Volumes. When
		               *   interpolating values from Volumes to Faces, it is implicitly assumed that the Face is in the
		               *   reference configuration of the current Volume. When these interpolated values must be used in
		               *   relation to the neighbouring Face, their orientation is likely to be incorrect and the
		               *   reordering necessary to be seen correctly is required. This variable provides the index to
		               *   the ordering array which transfers the current ordering to that required as if the Face were
		               *   seen from the opposite Volume.
		               */

		struct Volume *volume; ///< Pointer to the neighbouring \ref Volume.
	} neigh_info[2];
};

// Constructor/Destructor functions ********************************************************************************* //

/// \brief Constructs the base \ref Volume \ref Intrusive_List.
struct Intrusive_List* constructor_Face_List
	(struct Simulation*const sim, ///< The \ref Simulation.
	 const struct Mesh*const mesh ///< The \ref Mesh.
	);

/// \brief Destructs the base \ref Face \ref Intrusive_List.
void destructor_Faces
	(struct Intrusive_List* Faces ///< Standard.
	);

/// \brief Cast from \ref Face\* to `const` \ref Face `*const`.
void const_cast_Face
	(const struct Face*const* dest, ///< Destination.
	 const struct Face*const src    ///< Source.
	);

// Helper functions ************************************************************************************************* //

#endif // DPG__Face_h__INCLUDED
