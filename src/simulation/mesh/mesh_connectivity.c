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
/// \file

#include "mesh_connectivity.h"
#include "mesh_readers.h"
#include "intrusive.h"

#include <assert.h>
#include <limits.h>

#include "macros.h"
#include "definitions_bc.h"
#include "definitions_core.h"
#include "definitions_elements.h"
#include "definitions_mesh.h"

#include "multiarray.h"
#include "matrix.h"
#include "vector.h"

#include "element.h"
#include "mesh.h"
#include "mesh_periodic.h"
#include "const_cast.h"

// Static function declarations ************************************************************************************* //

/// \brief Container for locally computed \ref Mesh_Connectivity members.
struct Mesh_Connectivity_l {
	struct Multiarray_Vector_i* v_to_v;     ///< Local version of \ref Mesh_Connectivity::v_to_v.
	struct Multiarray_Vector_i* v_to_lf;    ///< Local version of \ref Mesh_Connectivity::v_to_lf.
	struct Multiarray_Vector_i* v_to_lf_wp; ///< Local version of \ref Mesh_Connectivity::v_to_lf_wp.
};

/** \brief Constructor for \ref Conn_info.
 *  \return Standard. */
static struct Conn_info* constructor_Conn_info
	(const struct Mesh_Data*const mesh_data,          ///< Standard.
	 const struct const_Intrusive_List*const elements ///< Standard.
	);

/// \brief Destructor for \ref Conn_info.
static void destructor_Conn_info
	(struct Conn_info* conn_info ///< Standard.
	);

/** \brief Compute the list of (f)ace (ve)rtices for each face.
 *  \return See brief. */
static void compute_f_ve
	(const struct Mesh_Data*const mesh_data,           ///< Standard.
	 const struct const_Intrusive_List*const elements, ///< Standard.
	 struct Conn_info* conn_info                       ///< The \ref Conn_info.
	);

/** \brief Compute the volume to (volume, local face) correspondence.
 *
 *	This function works by finding the volume and local face indices corresponding to the faces and then establishing
 *	connections by checking for adjacent matching face vertices in the sorted list of face vertices. If no connection
 *	is found, the entries are set to -1.
 */
static void compute_v_to__v_lf
	(const struct Conn_info*const conn_info,      ///< The \ref Conn_info.
	 struct Mesh_Connectivity_l*const mesh_conn_l ///< The \ref Mesh_Connectivity_l.
	);

/** \brief Add the boundary condition information to \ref Mesh_Connectivity::v_to_lf.
 *
 *	This function replaces all invalid entries in the volume to local face container with the number corresponding to
 *	the appropriate boundary condition.
 */
static void add_bc_info
	(const struct Mesh_Data*const mesh_data,       ///< Standard.
	 const struct Conn_info*const conn_info,       ///< The \ref Conn_info.
	 struct Mesh_Connectivity_l*const mesh_conn_l, ///< The \ref Mesh_Connectivity_l.
	 const bool include_periodic                   ///< Flag for whether periodic boundaries should be included.
	);

// Interface functions ********************************************************************************************** //

struct Mesh_Connectivity* constructor_Mesh_Connectivity
	(const struct Mesh_Data*const mesh_data, const struct const_Intrusive_List* elements)
{
	struct Mesh_Connectivity_l mesh_conn_l;
	struct Conn_info* conn_info = constructor_Conn_info(mesh_data,elements); // destructed

	compute_f_ve(mesh_data,elements,conn_info);
	compute_v_to__v_lf(conn_info,&mesh_conn_l);
	add_bc_info(mesh_data,conn_info,&mesh_conn_l,false);
	add_bc_info(mesh_data,conn_info,&mesh_conn_l,true);

	destructor_Multiarray_Vector_i(conn_info->f_ve);
	destructor_Vector_i(conn_info->ind_f_ve);
	destructor_conditional_Multiarray_Vector_i(conn_info->f_ve_per);
	destructor_Conn_info(conn_info);

	struct Mesh_Connectivity* mesh_conn = calloc(1,sizeof *mesh_conn); // returned

	const_constructor_move_Multiarray_Vector_i(&mesh_conn->v_to_v,mesh_conn_l.v_to_v);
	const_constructor_move_Multiarray_Vector_i(&mesh_conn->v_to_lf,mesh_conn_l.v_to_lf);
	const_constructor_move_Multiarray_Vector_i(&mesh_conn->v_to_lf_wp,mesh_conn_l.v_to_lf_wp);

	return mesh_conn;
}

void destructor_Mesh_Connectivity (struct Mesh_Connectivity* mesh_conn)
{
	destructor_const_Multiarray_Vector_i(mesh_conn->v_to_v);
	destructor_const_Multiarray_Vector_i(mesh_conn->v_to_lf);
	destructor_conditional_const_Multiarray_Vector_i(mesh_conn->v_to_lf_wp);

	free(mesh_conn);
}

void set_f_node_nums (struct Vector_i**const f_node_nums, const struct const_Vector_i*const node_nums)
{
	*f_node_nums = constructor_copy_Vector_i_i(node_nums->ext_0,node_nums->data);
	sort_Vector_i(*f_node_nums);
}

bool check_pfe_boundary (const int bc, const bool include_periodic)
{
	const int bc_base = bc % BC_STEP_SC;
	switch (bc_base) {
	case BC_INFLOW: case BC_INFLOW_ALT1: case BC_INFLOW_ALT2: // deprecated
	case BC_OUTFLOW: case BC_OUTFLOW_ALT1: case BC_OUTFLOW_ALT2: // deprecated
	case BC_UPWIND: // Advection
	case BC_UPWIND_ALT1:
	case BC_UPWIND_ALT2:
	case BC_UPWIND_ALT3:
	case BC_UPWIND_ALT4:
	case BC_UPWIND_ALT5:
	case BC_DIRICHLET:        // Diffusion
	case BC_DIRICHLET_ALT1:
	case BC_NEUMANN:
	case BC_NEUMANN_ALT1:
	case BC_RIEMANN:          // Euler
	case BC_SLIPWALL:
	case BC_BACKPRESSURE:
	case BC_TOTAL_TP:
	case BC_SUPERSONIC_IN:
	case BC_SUPERSONIC_OUT:
	case BC_NOSLIP_ADIABATIC:    // Navier-Stokes
	case BC_NOSLIP_DIABATIC:
	case BC_NOSLIP_ALL_ROTATING:
		return true;
		break;
	case PERIODIC_XL: case PERIODIC_XR:
	case PERIODIC_YL: case PERIODIC_YR:
	case PERIODIC_ZL: case PERIODIC_ZR:
	case PERIODIC_XL_REFLECTED_Y: case PERIODIC_XR_REFLECTED_Y:
		return include_periodic;
	default:
		EXIT_ERROR("Unsupported: %d\n",bc_base);
		break;
	}
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/// \brief Container for the list of boundary face information.
struct Boundary_Face_Info {
	ptrdiff_t n_pfe;                ///< The number of physical face elements.
	ptrdiff_t n_bf;                 /**< The number of boundary faces (n_bf = n_pfe - n_pf, where n_pf: number of
	                                 *   periodic faces). */
	struct Boundary_Face** b_faces; ///< The list of \ref Boundary_Face entities.
};

/// \brief Container for boundary face information.
struct Boundary_Face {
	int bc;                     ///< The value of the boundary condition.
	struct Vector_i* node_nums; ///< The node numbers of the face vertices.
};

/** \brief Constructor for \ref Boundary_Face_Info.
 *  \return Standard. */
static struct Boundary_Face_Info* constructor_Boundary_Face_Info
	(ptrdiff_t const n_pfe, ///< \ref Boundary_Face_Info::n_pfe.
	 ptrdiff_t const n_bf   ///< \ref Boundary_Face_Info::n_bf.
	);

/// \brief Destructor for \ref Boundary_Face_Info.
static void destructor_Boundary_Face_Info
	(struct Boundary_Face_Info* bf_info ///< Standard.
	);

/** \brief See return.
 *  \return The sum of the number of faces of all volumes.
 */
static ptrdiff_t compute_sum_n_f
	(const struct const_Intrusive_List*const elements, ///< Standard.
	 const struct const_Vector_i*const volume_types    ///< Defined in \ref Conn_info.
	);

/// \brief Set boundary face info for all entries in the list.
static void set_bf_info
	(struct Boundary_Face_Info* bf_info,     ///< The \ref Boundary_Face_Info.
	 const ptrdiff_t ind_pfe,                ///< Index of the first physical face element in the mesh element list.
	 const struct Mesh_Data*const mesh_data, ///< The \ref Mesh_Data.
	 const bool include_periodic             ///< Flag for whether periodic boundaries should be included.
	);

/** \brief Count the number of boundary faces.
 *  \return See brief. */
static ptrdiff_t count_boundary_faces
	(const ptrdiff_t ind_pfe,                     /**< Index of the first physical face element in the mesh element
	                                               *   list. */
	 const ptrdiff_t n_pfe,                       ///< The number of physical face elements.
	 const struct const_Matrix_i*const elem_tags, ///< \ref Mesh_Data::elem_tags.
	 const bool include_periodic                  ///< Flag for whether periodic faces should be included.
	);

/** \brief Reorder the boundary face entries in the list according to the node numbering.
 *	\warning This is not currently done in place.
 */
static void reorder_b_faces
	(struct Boundary_Face**const b_faces, ///< \ref Boundary_Face_Info::b_faces.
	 struct Vector_i* ordering            ///< The new ordering
	);

/** \brief Comparison function for std::bsearch between \ref Boundary_Face\*\* `a` and `b`.
 *  \return The \ref cmp_Vector_T of the `node_nums` of `a` and `b`.
 *
 *  \note Input Vectors must be have sorted data.
 */
static int cmp_Boundary_Face
	(const void *a, ///< Variable 1.
	 const void *b  ///< Variable 2.
	);

/// \brief Update v_to_lf replacing unused boundary entries with the boundary values.
static void update_v_to_lf_bc
	(struct Multiarray_Vector_i*const v_to_lf, ///< \ref Mesh_Connectivity_l::v_to_lf.
	 const int*const v_to_lf_i                 ///< The `v_to_lf` data with bc information included.
	);

static struct Conn_info* constructor_Conn_info
	(const struct Mesh_Data*const mesh_data, const struct const_Intrusive_List*const elements)
{
	const int d = (int)mesh_data->nodes->ext_1;

	const struct const_Vector_i*const elem_per_dim = mesh_data->elem_per_dim;

	const ptrdiff_t n_v   = mesh_data->elem_per_dim->data[d],
	                ind_v = get_first_volume_index(mesh_data->elem_per_dim,d);

	struct const_Vector_i* volume_types =
		constructor_move_const_Vector_i_i(n_v,false,&mesh_data->elem_types->data[ind_v]); // keep

	struct Vector_i* v_n_lf = constructor_empty_Vector_i(n_v);
	for (ptrdiff_t v = 0; v < n_v; ++v) {
		const struct const_Element*const element = get_element_by_type(elements,volume_types->data[v]);
		v_n_lf->data[v] = element->n_f;
	}

	struct Conn_info* conn_info = calloc(1,sizeof *conn_info); // returned;
	const_cast_i(&conn_info->d,d);
	*(const struct const_Vector_i**)&conn_info->elem_per_dim = elem_per_dim;
	conn_info->volume_types = volume_types;
	conn_info->v_n_lf       = v_n_lf;

	return conn_info;
}

static void destructor_Conn_info (struct Conn_info* conn_info)
{
	destructor_const_Vector_i(conn_info->volume_types);
	destructor_Vector_i(conn_info->v_n_lf);
	free(conn_info);
}

static void compute_f_ve
	(const struct Mesh_Data*const mesh_data, const struct const_Intrusive_List*const elements,
	 struct Conn_info* conn_info)
{
	const int d = conn_info->d;
	const ptrdiff_t ind_v = get_first_volume_index(conn_info->elem_per_dim,d),
	                n_v   = conn_info->elem_per_dim->data[d];

	struct const_Vector_i* volume_types = conn_info->volume_types;

	const ptrdiff_t sum_n_f = compute_sum_n_f(elements,volume_types);
	struct Multiarray_Vector_i* f_ve = constructor_empty_Multiarray_Vector_i(true,1,&sum_n_f); // returned

	const struct const_Vector_i*const*const volume_nums = &mesh_data->node_nums->data[ind_v];
	for (ptrdiff_t v = 0, ind_f = 0; v < n_v; ++v) {
		const struct const_Element*const element = get_element_by_type(elements,volume_types->data[v]);
		for (ptrdiff_t f = 0, f_max = conn_info->v_n_lf->data[v]; f < f_max; ++f) {
			const struct const_Vector_i*const f_ve_f = element->f_ve->data[f];
			const ptrdiff_t n_n = f_ve_f->ext_0;

			struct Vector_i* f_ve_curr = f_ve->data[ind_f];
			resize_Vector_i(f_ve_curr,n_n);
			for (ptrdiff_t n = 0; n < n_n; ++n)
				f_ve_curr->data[n] = volume_nums[v]->data[f_ve_f->data[n]];
			++ind_f;
		}
	}

	conn_info->f_ve     = f_ve; // keep
	conn_info->ind_f_ve = sort_Multiarray_Vector_i(f_ve,true); // keep

	correct_f_ve_for_periodic(mesh_data,conn_info);
}

static void compute_v_to__v_lf (const struct Conn_info*const conn_info, struct Mesh_Connectivity_l*const mesh_conn_l)
{
	struct Vector_i* ind_f_ve_V = conn_info->ind_f_ve;

	const ptrdiff_t d   = conn_info->d,
	                n_v = conn_info->elem_per_dim->data[d],
	                n_f = ind_f_ve_V->ext_0;

	struct Vector_i* v_n_lf = conn_info->v_n_lf;

	// Store global volume and local face indices corresponding to each global face (reordered).
	int* ind_v_i  = malloc((size_t)n_f * sizeof *ind_v_i);  // moved
	int* ind_lf_i = malloc((size_t)n_f * sizeof *ind_lf_i); // moved

	for (int ind_vf = 0, v = 0; v < n_v; ++v) {
		const ptrdiff_t lf_max = v_n_lf->data[v];
		for (int lf = 0; lf < lf_max; ++lf) {
			ind_v_i[ind_vf]  = v;
			ind_lf_i[ind_vf] = lf;
			++ind_vf;
		}
	}

	struct Vector_i* ind_v_V  = constructor_move_Vector_i_i(n_f,true,ind_v_i);  // destructed
	struct Vector_i* ind_lf_V = constructor_move_Vector_i_i(n_f,true,ind_lf_i); // destructed

	reorder_Vector_i(ind_v_V,ind_f_ve_V->data);
	reorder_Vector_i(ind_lf_V,ind_f_ve_V->data);

	// Compute v_to_v, v_to_lf, and (optionally) v_to_lf_wp
	int*const v_to_v_i     = malloc((size_t)n_f * sizeof *v_to_v_i);  // free
	int*const v_to_lf_i    = malloc((size_t)n_f * sizeof *v_to_lf_i); // free
	int*const v_to_lf_wp_i = malloc((size_t)n_f * sizeof *v_to_lf_i); // free
	for (int i = 0; i < n_f; ++i) {
		v_to_v_i[i]     = -1;
		v_to_lf_i[i]    = -1;
		v_to_lf_wp_i[i] = -1;
	}

	const int*const ind_f_ve_i = ind_f_ve_V->data;
	struct Multiarray_Vector_i*const f_ve     = conn_info->f_ve;
	struct Multiarray_Vector_i*const f_ve_per = conn_info->f_ve_per;

	const bool include_periodic = (f_ve_per != NULL);

	for (ptrdiff_t f = 0; f < n_f; ++f) {
		const ptrdiff_t ind_0 = ind_f_ve_i[f];
		if ((f+1) < n_f && check_equal_Vector_i(f_ve->data[f],f_ve->data[f+1])) {
			const ptrdiff_t ind_1 = ind_f_ve_i[f+1];
			v_to_v_i[ind_0]  = ind_v_i[f+1];
			v_to_lf_i[ind_0] = ind_lf_i[f+1];
			v_to_v_i[ind_1]  = ind_v_i[f];
			v_to_lf_i[ind_1] = ind_lf_i[f];

			if (include_periodic && (check_equal_Vector_i(f_ve_per->data[f],f_ve_per->data[f+1]))) {
				v_to_lf_wp_i[ind_0] = v_to_lf_i[ind_0];
				v_to_lf_wp_i[ind_1] = v_to_lf_i[ind_1];
			}
			++f;
		}
	}
	destructor_Vector_i(ind_v_V);
	destructor_Vector_i(ind_lf_V);

	mesh_conn_l->v_to_v  = constructor_copy_Multiarray_Vector_i_i(v_to_v_i,conn_info->v_n_lf->data,1,&n_v);  // keep
	mesh_conn_l->v_to_lf = constructor_copy_Multiarray_Vector_i_i(v_to_lf_i,conn_info->v_n_lf->data,1,&n_v); // keep
	mesh_conn_l->v_to_lf_wp = ( !include_periodic ? NULL :
		constructor_copy_Multiarray_Vector_i_i(v_to_lf_wp_i,conn_info->v_n_lf->data,1,&n_v)); // keep

	free(v_to_v_i);
	free(v_to_lf_i);
	free(v_to_lf_wp_i);
}

static void add_bc_info
	(const struct Mesh_Data*const mesh_data, const struct Conn_info*const conn_info,
	 struct Mesh_Connectivity_l*const mesh_conn_l, const bool include_periodic)
{
	const int d = conn_info->d;
	const ptrdiff_t ind_pfe = get_first_volume_index(conn_info->elem_per_dim,d-1),
	                n_pfe   = conn_info->elem_per_dim->data[d-1],
	                n_bf    = count_boundary_faces(ind_pfe,n_pfe,mesh_data->elem_tags,include_periodic);

	if ((n_bf == 0) || (include_periodic && conn_info->f_ve_per == NULL))
		return;

	// Set the boundary face info and sort it by node_nums.
	struct Boundary_Face_Info* bf_info = constructor_Boundary_Face_Info(n_pfe,n_bf); // destructed
	set_bf_info(bf_info,ind_pfe,mesh_data,include_periodic);

	// Copy the pointers to the node_nums into a Multiarray_Vector_i (for sorting).
	struct Multiarray_Vector_i* bf_ve = constructor_empty_Multiarray_Vector_i(true,1,&n_bf); // destructed
	bf_ve->owns_data = false;
	for (ptrdiff_t i = 0; i < n_bf; ++i) {
		destructor_Vector_i(bf_ve->data[i]);
		bf_ve->data[i] = bf_info->b_faces[i]->node_nums;
	}

	struct Vector_i* ind_bf_ve = sort_Multiarray_Vector_i(bf_ve,true); // destructed
	destructor_Multiarray_Vector_i(bf_ve);

	reorder_b_faces(bf_info->b_faces,ind_bf_ve);
	destructor_Vector_i(ind_bf_ve);

	// Set the unused entries in v_to_lf to the values of the associated boundary conditions.
	struct Multiarray_Vector_i* v_to_lf = NULL;
	struct Vector_i*const* f_ve_V = NULL;

	if (!include_periodic) {
		v_to_lf = mesh_conn_l->v_to_lf;
		f_ve_V  = conn_info->f_ve->data;
	} else {
		v_to_lf = mesh_conn_l->v_to_lf_wp;
		f_ve_V  = conn_info->f_ve_per->data;
	}
	struct Vector_i*const v_to_lf_V = collapse_Multiarray_Vector_i(v_to_lf); // destructed

	int*const v_to_lf_i  = v_to_lf_V->data,
	   *const ind_f_ve_i = conn_info->ind_f_ve->data;

	struct Boundary_Face*const f_curr = calloc(1,sizeof *f_curr); // free
	f_curr->bc = -1;

	const ptrdiff_t n_f = v_to_lf_V->ext_0;
	for (ptrdiff_t n = 0; n < n_f; ++n) {
		const ptrdiff_t f = ind_f_ve_i[n];
		if (v_to_lf_i[f] != -1)
			continue;

		f_curr->node_nums = f_ve_V[n];
		struct Boundary_Face*const*const bf_curr =
			bsearch(&f_curr,bf_info->b_faces,(size_t)n_bf,sizeof(bf_info->b_faces[0]),cmp_Boundary_Face);

		assert(bf_curr != NULL);
		v_to_lf_i[f] = (*bf_curr)->bc;
	}
	free(f_curr);

	update_v_to_lf_bc(v_to_lf,v_to_lf_i);

	destructor_Vector_i(v_to_lf_V);
	destructor_Boundary_Face_Info(bf_info);
}

// Level 1 ********************************************************************************************************** //

/** \brief Constructor for \ref Boundary_Face.
 *  \return Standard. */
static struct Boundary_Face* constructor_Boundary_Face ( );

/** \brief Copy constructor for \ref Boundary_Face_Info::b_faces.
 *  \return Standard. */
static struct Boundary_Face** constructor_b_faces
	(const ptrdiff_t n_bf,            ///< \ref Boundary_Face_Info::n_bf.
	 struct Boundary_Face**const src  ///< The source \ref Boundary_Face_Info::b_faces data.
	);

/// Destructor for \ref Boundary_Face_Info::b_faces.
static void destructor_b_faces
	(const ptrdiff_t n_bf,            ///< \ref Boundary_Face_Info::n_bf.
	 struct Boundary_Face**const src  ///< Standard.
	);

/// \brief Destructor for \ref Boundary_Face.
static void destructor_Boundary_Face
	(struct Boundary_Face* bf ///< Standard.
	);

static ptrdiff_t compute_sum_n_f
	(const struct const_Intrusive_List*const elements, const struct const_Vector_i*const volume_types)
{
	ptrdiff_t sum_n_f = 0;
	for (ptrdiff_t v = 0, v_max = volume_types->ext_0; v < v_max; ++v) {
		const struct const_Element*const element = get_element_by_type(elements,volume_types->data[v]);
		sum_n_f += element->n_f;
	}
	return sum_n_f;
}

static struct Boundary_Face_Info* constructor_Boundary_Face_Info (const ptrdiff_t n_pfe, const ptrdiff_t n_bf)
{
	struct Boundary_Face_Info* bf_info = calloc(1,sizeof *bf_info); // returned

	struct Boundary_Face** b_faces = malloc((size_t)n_bf * sizeof *b_faces); // keep
	for (ptrdiff_t i = 0; i < n_bf; ++i)
		b_faces[i] = constructor_Boundary_Face(); // destructed

	bf_info->n_pfe   = n_pfe;
	bf_info->n_bf    = n_bf;
	bf_info->b_faces = b_faces;

	return bf_info;
}

static void destructor_Boundary_Face_Info (struct Boundary_Face_Info* bf_info)
{
	const ptrdiff_t n_bf = bf_info->n_bf;
	for (ptrdiff_t i = 0; i < n_bf; ++i)
		destructor_Boundary_Face(bf_info->b_faces[i]);
	free(bf_info->b_faces);

	free(bf_info);
}

static void set_bf_info
	(struct Boundary_Face_Info* bf_info, const ptrdiff_t ind_pfe, const struct Mesh_Data*const mesh_data,
	 const bool include_periodic)
{
	const struct const_Matrix_i*const            elem_tags = mesh_data->elem_tags;
	const struct const_Multiarray_Vector_i*const node_nums = mesh_data->node_nums;

	ptrdiff_t count_bf = 0;

	const ptrdiff_t n_max = ind_pfe+bf_info->n_pfe;
	for (ptrdiff_t n = ind_pfe; n < n_max; ++n) {
		const int bc = get_val_const_Matrix_i(n,0,elem_tags);

		if (!check_pfe_boundary(bc,include_periodic))
			continue;

		struct Boundary_Face*const bf = bf_info->b_faces[count_bf];
		bf->bc = bc;

		set_f_node_nums(&bf->node_nums,node_nums->data[n]);

		++count_bf;
	}

	if (count_bf != bf_info->n_bf)
		EXIT_ERROR("Did not find the correct number of boundary face entities: %td %td",count_bf,bf_info->n_bf);
}

static ptrdiff_t count_boundary_faces
	(const ptrdiff_t ind_pfe, const ptrdiff_t n_pfe, const struct const_Matrix_i*const elem_tags,
	 const bool include_periodic)
{
	ptrdiff_t count = 0;

	const ptrdiff_t n_max = ind_pfe+n_pfe;
	for (ptrdiff_t n = ind_pfe; n < n_max; ++n) {
		if (check_pfe_boundary(get_val_const_Matrix_i(n,0,elem_tags),include_periodic))
			++count;
	}
	return count;
}

static void reorder_b_faces (struct Boundary_Face**const b_faces, struct Vector_i* ordering)
{
	const ptrdiff_t n_bf = ordering->ext_0;

	struct Boundary_Face** b_faces_cpy = constructor_b_faces(n_bf,b_faces); // destructed

	const int*const ordering_i = ordering->data;
	for (ptrdiff_t i = 0; i < n_bf; ++i) {
		const int ind_bf = ordering_i[i];
		b_faces_cpy[i]->bc = b_faces[ind_bf]->bc;
		set_to_data_Vector_i(b_faces_cpy[i]->node_nums,b_faces[ind_bf]->node_nums->data);
	}

	for (ptrdiff_t i = 0; i < n_bf; ++i) {
		b_faces[i]->bc = b_faces_cpy[i]->bc;
		set_to_data_Vector_i(b_faces[i]->node_nums,b_faces_cpy[i]->node_nums->data);
	}

	destructor_b_faces(n_bf,b_faces_cpy);
}

static int cmp_Boundary_Face (const void *a, const void *b)
{
	const struct Boundary_Face*const*const ia = (const struct Boundary_Face*const*const) a,
	                          *const*const ib = (const struct Boundary_Face*const*const) b;

	return cmp_Vector_i(&(*ia)->node_nums,&(*ib)->node_nums);
}

static void update_v_to_lf_bc (struct Multiarray_Vector_i*const v_to_lf, const int*const v_to_lf_i)
{
	ptrdiff_t count = 0;
	const ptrdiff_t i_max = compute_size(v_to_lf->order,v_to_lf->extents);
	for (ptrdiff_t i = 0; i < i_max; ++i) {
		struct Vector_i* v_to_lf_curr = v_to_lf->data[i];
		int*const data = v_to_lf_curr->data;
		const ptrdiff_t j_max = v_to_lf_curr->ext_0;
		for (ptrdiff_t j = 0; j < j_max; ++j) {
			if (data[j] == -1)
				data[j] = v_to_lf_i[count];
			count++;
		}
	}
}

// Level 2 ********************************************************************************************************** //

/** \brief Constructor for \ref Boundary_Face where the memory is allocated for the node_nums but not set.
 *  \return Standard. */
static struct Boundary_Face* constructor_empty_Boundary_Face
	(struct Boundary_Face* src ///< The source data.
	);

static struct Boundary_Face* constructor_Boundary_Face ( )
{
	struct Boundary_Face* bf = calloc(1,sizeof *bf); // returned

	return bf;
}

static void destructor_Boundary_Face (struct Boundary_Face* bf)
{
	destructor_Vector_i(bf->node_nums);
	free(bf);
}

static struct Boundary_Face** constructor_b_faces (const ptrdiff_t n_bf, struct Boundary_Face**const src)
{
	struct Boundary_Face** dest = malloc((size_t)n_bf * sizeof *dest); // free
	for (ptrdiff_t i = 0; i < n_bf; ++i)
		dest[i] = constructor_empty_Boundary_Face(src[i]); // destructed

	return dest;
}

static void destructor_b_faces (const ptrdiff_t n_bf, struct Boundary_Face** src)
{
	for (ptrdiff_t i = 0; i < n_bf; ++i)
		destructor_Boundary_Face(src[i]);
	free(src);
}

// Level 3 ********************************************************************************************************** //

static struct Boundary_Face* constructor_empty_Boundary_Face (struct Boundary_Face* src)
{
	struct Boundary_Face* bf = calloc(1,sizeof *bf); // returned

	bf->bc        = src->bc;
	bf->node_nums = constructor_empty_Vector_i(src->node_nums->ext_0);

	return bf;
}
