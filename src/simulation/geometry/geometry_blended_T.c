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
 */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "gsl/gsl_math.h"

#include "macros.h"
#include "definitions_bc.h"
#include "definitions_core.h"
#include "definitions_geometry.h"
#include "definitions_mesh.h"
#include "definitions_tol.h"


#include "def_templates_geometry_blended.h"

#include "def_templates_volume_solver.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_operators.h"
#include "def_templates_test_case.h"

// Static function declarations ************************************************************************************* //

/** \brief Version of \ref constructor_xyz_fptr_T used for the blended curved surface geometry corrections.
 *  \return See brief. */
static const struct const_Multiarray_R* constructor_xyz_blended_ce
	(const char ce_type,                     ///< The surface computational element type.
	 const char n_type,                      ///< Defined for \ref constructor_xyz_fptr_T.
	 const struct const_Multiarray_R* xyz_i, ///< Defined for \ref constructor_xyz_fptr_T.
	 const struct Solver_Volume_T* s_vol,    ///< Defined for \ref constructor_xyz_fptr_T.
	 const struct Simulation* sim            ///< Defined for \ref constructor_xyz_fptr_T.
	);

/** \brief Compute the radius of the curved boundary as that of the first vertex on the boundary.
 *  \return See brief. */
static double compute_radius_from_xyz_ve
	(const int n_component_radius,                      /**< The number of components in xyz_ve to use for the
	                                                     *   computation of the radius. */
	 const struct const_Matrix_R*const xyz,             ///< The xyz coordinates in which to search.
	 const struct const_Matrix_R*const xyz_ve,          ///< The xyz coordinates of the vertices.
	 const struct Blended_Parametric_Data_T*const b_p_d ///< \ref Blended_Parametric_Data_T.
	);

// Interface functions ********************************************************************************************** //

const struct const_Multiarray_R* constructor_xyz_blended_T
	(const char n_type, const struct const_Multiarray_R* xyz_i, const struct Solver_Volume_T* s_vol,
	 const struct Simulation* sim)
{
	UNUSED(s_vol);
	UNUSED(sim);
	assert(DIM >= 2);
	assert(DIM == xyz_i->extents[1]);

	const struct const_Multiarray_R*const xyz_e =
		( DIM == DMAX ? constructor_xyz_blended_ce('e',n_type,xyz_i,s_vol,sim) : xyz_i ); // dest. (if required)

	const struct const_Multiarray_R*const xyz = constructor_xyz_blended_ce('f',n_type,xyz_e,s_vol,sim); // returned
	if (DIM == DMAX)
		destructor_const_Multiarray_R(xyz_e);

	return xyz;
}

struct Boundary_Comp_Elem_Data_T constructor_static_Boundary_Comp_Elem_Data_T
	(const char ce_type, const char n_type, const int p_geom, const struct Solver_Volume_T*const s_vol)
{
	assert(DIM == DMAX || ce_type != 'e');
	assert(ce_type == 'e' || ce_type == 'f');

	struct Boundary_Comp_Elem_Data_T b_ce_d;

	struct Volume* vol = (struct Volume*) s_vol;
	const struct const_Element*const e = vol->element;

	b_ce_d.s_type = e->s_type;
	b_ce_d.xyz_ve = vol->xyz_ve;

	if (ce_type == 'f') {
		b_ce_d.bc_boundaries = vol->bc_faces;
		b_ce_d.b_ve          = e->f_ve;
	} else {
		b_ce_d.bc_boundaries = vol->bc_edges;
		b_ce_d.b_ve          = e->e_ve;
	}

	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;
	switch (n_type) {
	case 'g':
		b_ce_d.vv0_vv_vX = get_Multiarray_Operator(g_e->vv0_vv_vg[1],(ptrdiff_t[]){0,0,p_geom,1});
		if (ce_type == 'f') {
			b_ce_d.vv0_bv_vX = alloc_MO_from_MO(g_e->vv0_fv_vgc,1,(ptrdiff_t[]){0,0,p_geom,1}); // free
		} else {
//			b_ce_d.vv0_bv_vX = set_MO_from_MO(g_e->vv0_ev_vg[1],1,(ptrdiff_t[]){0,0,p_geom,1});
			EXIT_ADD_SUPPORT; // Required for 3D.
		}
		break;
	case 'v':
//		b_ce_d.vv0_vv_vX = get_Multiarray_Operator(g_e->vv0_vv_vv,(ptrdiff_t[]){0,0,p_geom,1});
		EXIT_ADD_SUPPORT; // Required for h-adaptation.
		break;
	case 'c':
		break; // do nothing.
	default:
		EXIT_ERROR("Unsupported: %c",n_type);
		break;
	}

	return b_ce_d;
}

void destructor_static_Boundary_Comp_Elem_Data_T (struct Boundary_Comp_Elem_Data_T*const b_ce_d)
{
	free_MO_from_MO(b_ce_d->vv0_bv_vX);
}

void set_Boundary_Comp_Elem_operators_T
	(struct Boundary_Comp_Elem_Data_T*const b_ce_d, const struct Solver_Volume_T*const s_vol, const char ce_type,
	 const char n_type, const int p, const int b)
{
	assert(ce_type == 'e' || ce_type == 'f');

	struct Volume* vol = (struct Volume*) s_vol;
	const struct Geometry_Element* g_e = &((struct Solver_Element*)vol->element)->g_e;
	switch (n_type) {
	case 'g':
		if (ce_type == 'f') {
			b_ce_d->vv0_vX_bX = get_Multiarray_Operator(g_e->vv0_vgc_fgc,(ptrdiff_t[]){b,0,0,p,p});
			b_ce_d->vv0_bX_vX = get_Multiarray_Operator(g_e->vv0_fgc_vgc,(ptrdiff_t[]){b,0,0,p,p});
			b_ce_d->vv0_vv_bX = get_Multiarray_Operator(g_e->vv0_vv_fgc, (ptrdiff_t[]){b,0,0,p,1});
			b_ce_d->vv0_vv_bv = get_Multiarray_Operator(g_e->vv0_vv_fv,  (ptrdiff_t[]){b,0,0,1,1});
			b_ce_d->vv0_bv_bX = get_Multiarray_Operator(g_e->vv0_fv_fgc, (ptrdiff_t[]){0,0,0,0,p,1});
		} else {
			EXIT_ADD_SUPPORT;
		}
		break;
	case 'v':
		EXIT_ADD_SUPPORT;
		break;
	case 'c':
		assert(ce_type == 'f');
/// \todo Can probably use 'v' operators here when implemented.
		b_ce_d->vv0_vv_bX  = get_Multiarray_Operator(g_e->vv0_vv_fgc,(ptrdiff_t[]){b,0,0,p,1});
		b_ce_d->vv0_vv_fcc = get_Multiarray_Operator(g_e->vv0_vv_fcc,(ptrdiff_t[]){b,0,0,p,1});
		break;
	default:
		EXIT_ERROR("Unsupported: %c",n_type);
		break;
	}
}

const struct const_Matrix_R* constructor_xyz_surf_diff_T
	(const struct Boundary_Comp_Elem_Data_T*const b_ce_d, const struct const_Matrix_R*const xyz_i,
	 const struct Solver_Volume_T*const s_vol, const char n_type, const struct Simulation*const sim)
{
	const struct Volume*const vol = (struct Volume*) s_vol;
	const struct Test_Case_T*const test_case = (struct Test_Case_T*) sim->test_case_rc->tc;
	struct Blended_Parametric_Data_T b_p_d =
		{ .xyz_ve = vol->xyz_ve,
		  .vv0_vv_bX   = b_ce_d->vv0_vv_bX,
		  .vv0_vv_bv   = b_ce_d->vv0_vv_bv,
		  .vv0_vv_fcc  = b_ce_d->vv0_vv_fcc,
		  .vv0_bv_bX   = b_ce_d->vv0_bv_bX,
		  .normal      = NULL,
		  .n_type      = n_type,
		  .domain_type = sim->domain_type,
		  .constructor_xyz = test_case->constructor_xyz,
		};

	assert(s_vol->constructor_xyz_surface != NULL);
	const struct const_Matrix_R*const xyz_b = s_vol->constructor_xyz_surface(&b_p_d);

	struct Matrix_R* xyz_surf_diff = NULL;
	switch (n_type) {
	case 'g': // fallthrough
	case 'v':
		xyz_surf_diff = constructor_mm_Matrix_R('N','N',-1.0,b_ce_d->vv0_vX_bX->op_std,xyz_i,'C'); // returned
		break;
	case 'c':
		xyz_surf_diff = constructor_copy_Matrix_R((struct Matrix_R*)xyz_i); // returned
		break;
	default:
		EXIT_ERROR("Unsupported: %c",n_type);
		break;
	}
	add_in_place_Matrix_R(1.0,xyz_surf_diff,xyz_b);

	return (struct const_Matrix_R*) xyz_surf_diff;
}


/// \todo Consider moving geometry specific functions out of the general blending function.

constructor_xyz_surface_fptr_T set_constructor_xyz_surface_fptr_T (const char*const geom_type, const int geom_prm_type)
{
	constructor_xyz_surface_fptr_T ptr = NULL;

	if (strcmp(geom_type,"n-cylinder") == 0) {
		if      (geom_prm_type == GEOM_PRM_RADIAL_PROJ) ptr = constructor_xyz_surface_cylinder_radial_proj_T;
		else if (geom_prm_type == GEOM_PRM_ARC_LENGTH ) ptr = constructor_xyz_surface_cylinder_arc_length_T;
		else if (geom_prm_type == GEOM_PRM_NORMAL_PROJ) ptr = constructor_xyz_surface_cylinder_normal_proj_T;
		else                                            EXIT_ERROR("Unsupported: %d\n",geom_prm_type);
	} else {
		EXIT_ERROR("Unsupported: %s\n",geom_type);
	}

	return ptr;
}

const struct const_Matrix_R* constructor_xyz_surface_cylinder_radial_proj_T
	(const struct Blended_Parametric_Data_T*const b_p_d)
{
	const struct const_Matrix_R xyz_ve_M = interpret_const_Multiarray_as_Matrix_R(b_p_d->xyz_ve);
	const struct const_Matrix_R*const xyz_ve_bX =
		constructor_mm_const_Matrix_R('N','N',1.0,b_p_d->vv0_vv_bX->op_std,&xyz_ve_M,'R'); // destructed

// make external
	const struct const_Matrix_R* xyz_ve_b = NULL;
	switch (b_p_d->n_type) {
	case 'g': // fallthrough
	case 'v':
		xyz_ve_b = xyz_ve_bX;
		break;
	case 'c': {
		const struct const_Multiarray_R*const xyz_ve_fcc = constructor_mm_NN1_Operator_const_Multiarray_R
			(b_p_d->vv0_vv_fcc,b_p_d->xyz_ve,'C','d',b_p_d->xyz_ve->order,NULL); // destructed

		const struct const_Multiarray_R*const xyz_fcc_Ma = b_p_d->constructor_xyz(0,xyz_ve_fcc,NULL,NULL); // dest.
		destructor_const_Multiarray_R(xyz_ve_fcc);

		const ptrdiff_t ext_0 = xyz_fcc_Ma->extents[0],
		                ext_1 = xyz_fcc_Ma->extents[1];
		xyz_ve_b = constructor_empty_const_Matrix_R(xyz_fcc_Ma->layout,ext_0,ext_1); // destructed
		reinterpret_const_Multiarray_as_Matrix_R(xyz_fcc_Ma,xyz_ve_b,ext_0,ext_1);
		const_cast_b(&xyz_fcc_Ma->owns_data,false);
		destructor_const_Multiarray_R(xyz_fcc_Ma);
/// \todo Check that computation of xyz_surf is not redundant below (i.e. exact boundary values are not already given in
///        xyz_ve_b).

		transpose_Matrix_R((struct Matrix_R*)xyz_ve_b,true);
		break;
	} default:
		EXIT_ERROR("Unsupported: %c.\n",b_p_d->n_type);
		break;
	}

	const double r = compute_radius_from_xyz_ve(2,&xyz_ve_M,xyz_ve_bX,b_p_d);

	const ptrdiff_t n_n = xyz_ve_b->ext_0,
	                dim = xyz_ve_b->ext_1;
	assert(dim >= 2);

	struct Matrix_R* xyz_surf = constructor_empty_Matrix_R('R',n_n,dim); // returned
	for (int n = 0; n < n_n; ++n) {
		const double*const xyz_b = get_row_const_Matrix_R(n,xyz_ve_b);
		const double theta = atan2(xyz_b[1],xyz_b[0]);

		double*const xyz_s = get_row_Matrix_R(n,xyz_surf);
		xyz_s[0] = r*cos(theta);
		xyz_s[1] = r*sin(theta);
		if (dim == DMAX)
			xyz_s[2] = xyz_b[2];
	}
	if (xyz_ve_b != xyz_ve_bX)
		destructor_const_Matrix_R(xyz_ve_b);
	destructor_const_Matrix_R(xyz_ve_bX);

	transpose_Matrix_R(xyz_surf,true);
	return (struct const_Matrix_R*) xyz_surf;
}

const struct const_Matrix_R* constructor_xyz_surface_cylinder_arc_length_T
	(const struct Blended_Parametric_Data_T*const b_p_d)
{
// test for arc length: Should be the same as \theta parametrization for circle.
//	1. Interpolate xyz_ve to boundary (vv0_vv_bv);
//	2. Compute parametric coordinates of vertices
//	3. Interpolate parametric coordinates to curved boundary nodes (vv0_bv_bgc)
//	4. Compute xyz curved
UNUSED(b_p_d); EXIT_ADD_SUPPORT;
}

const struct const_Matrix_R* constructor_xyz_surface_cylinder_normal_proj_T
	(const struct Blended_Parametric_Data_T*const b_p_d)
{
//	1. Interpolate xyz_ve to boundary (vv0_vv_bgc);
//	2. Compute normal vector to surface (pointing outwards)
//	3. Compute xyz curved
UNUSED(b_p_d); EXIT_ADD_SUPPORT;
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

/** \brief Compute the minimum polynomial degree to use as a base for the correction of the next approximation of the
 *         blended geometry.
 *  \return See brief.
 *
 *  Most commonly, the base degree is equal to the input volume polynomial degree. However, for certain blending
 *  functions, notably that of Lenoir (see \ref Simulation::geom_blending for reference), a sequential correction from
 *  the lowest (degree 2) to desired order is required and this function returns the lowest degree.
 */
static int compute_p_base_min
	(const struct Solver_Volume_T*const s_vol, ///< The current volume.
	 const struct Simulation*const sim         ///< \ref Simulation.
	);

/** \brief Constructor for the values of the blending to be used for each of the volume geometry nodes.
 *  \return See brief.
 *
 *  The blending functions for tensor-product and simplex element types are selected according to their definitions in
 *  \todo [add ref] Zwanenburg2017 (Discrete_Curvature).
 */
static const struct const_Vector_R* constructor_blend_values
	(const int ind_b,                                    ///< The index of the boundary under consideration.
	 const ptrdiff_t n_n,                                ///< The 'n'umber of geometry 'n'odes.
	 const struct Boundary_Comp_Elem_Data_T*const b_ce_d ///< \ref Boundary_Comp_Elem_Data_T.
	);

/** \brief Constructor for a \ref const_Matrix_T\* storing the difference between the corrected and current geometry
 *         nodal values with volume geometry nodes projected to the boundary.
 *  \return See brief. */
static const struct const_Matrix_R* constructor_xyz_diff_T
	(const struct Boundary_Comp_Elem_Data_T*const b_ce_d, ///< \ref Boundary_Comp_Elem_Data_T.
	 const struct const_Matrix_R*const xyz_i,             ///< Input xyz coordinates.
	 const struct Solver_Volume_T*const s_vol,            ///< The current \ref Solver_Volume_T.
	 const char n_type,                                   ///< \ref Blended_Parametric_Data_T::n_type.
	 const struct Simulation*const sim                    ///< \ref Simulation.
	);

/** \brief Find the index of the first node in the matrix of boundary nodes which is also a vertex node.
 *  \return See brief. */
static ptrdiff_t find_boundary_vertex_index
	(const struct const_Matrix_R*const xyz,   ///< The xyz coordinates in which to search.
	 const struct const_Matrix_R*const xyz_ve ///< The xyz coordinates of the vertices.
	);

static const struct const_Multiarray_R* constructor_xyz_blended_ce
	(const char ce_type, const char n_type, const struct const_Multiarray_R* xyz_i,
	 const struct Solver_Volume_T* s_vol, const struct Simulation* sim)
{
	assert(n_type == 'v' || n_type == 'g');
	const int p_geom = ( n_type == 'v' ? 2 : s_vol->p_ref ),
	          p_min = GSL_MIN(compute_p_base_min(s_vol,sim),p_geom);
	struct Boundary_Comp_Elem_Data_T b_ce_d =
		constructor_static_Boundary_Comp_Elem_Data_T(ce_type,n_type,p_geom,s_vol); // destructed

	struct Multiarray_R* xyz = constructor_copy_Multiarray_R((struct Multiarray_R*)xyz_i); // returned
	struct Matrix_R xyz_M = interpret_Multiarray_as_Matrix_R(xyz);

	const ptrdiff_t n_n = xyz_i->extents[0];
	for (int p = p_min; p <= p_geom; ++p) {
	for (int b = 0; b < b_ce_d.bc_boundaries->ext_0; ++b) {
		if (!is_face_bc_curved(b_ce_d.bc_boundaries->data[b]))
			continue;

		set_Boundary_Comp_Elem_operators_T(&b_ce_d,s_vol,ce_type,n_type,p,b);

		const struct const_Vector_R*const blend_values = constructor_blend_values(b,n_n,&b_ce_d); // destructed
		const struct const_Matrix_R*const xyz_diff =
			constructor_xyz_diff_T(&b_ce_d,(struct const_Matrix_R*)&xyz_M,s_vol,n_type,sim); // destructed

		for (int d = 0; d < DIM; ++d) {
			Real*const data_xyz            = get_col_Matrix_R(d,&xyz_M);
			const Real*const data_xyz_diff = get_col_const_Matrix_R(d,xyz_diff);
			for (int n = 0; n < n_n; ++n)
				data_xyz[n] += blend_values->data[n]*data_xyz_diff[n];
		}

		destructor_const_Vector_d(blend_values);
		destructor_const_Matrix_R(xyz_diff);
	}}
	destructor_static_Boundary_Comp_Elem_Data_T(&b_ce_d);

	return (struct const_Multiarray_R*) xyz;
}

static double compute_radius_from_xyz_ve
	(const int n_component_radius, const struct const_Matrix_R*const xyz, const struct const_Matrix_R*const xyz_ve,
	 const struct Blended_Parametric_Data_T*const b_p_d)
{
	const ptrdiff_t ind_b = find_boundary_vertex_index(xyz,xyz_ve);

	double r = 0.0;
	switch (b_p_d->domain_type) {
	case DOM_BLENDED:
		r = norm_d(n_component_radius,get_row_const_Matrix_R(ind_b,xyz),"L2");
		break;
	case DOM_PARAMETRIC: {
		const ptrdiff_t extents[] = { 1, xyz->ext_1, };
		const double*const data = get_row_const_Matrix_R(ind_b,xyz);
		const struct const_Multiarray_R xyz_Ma = { .order = 2, .extents = extents, .layout = 'C', .data = data, };
		const struct const_Multiarray_R* xyz_surf = b_p_d->constructor_xyz(0,&xyz_Ma,NULL,NULL); // destructed

		r = norm_d(n_component_radius,xyz_surf->data,"L2");
		destructor_const_Multiarray_R(xyz_surf);
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",b_p_d->domain_type);
		break;
	}
	return r;
}

// Level 1 ********************************************************************************************************** //

static int compute_p_base_min (const struct Solver_Volume_T*const s_vol, const struct Simulation*const sim)
{
	const struct Volume*const vol      = (struct Volume*) s_vol;
	const struct const_Element*const e = vol->element;

	int p_base_min = 0;

	const int s_type = e->s_type;
	switch (s_type) {
	case ST_TP:
/// \todo make `geom_blending` an array of `int`s such that switch statements can be used here.
		if (strcmp(sim->geom_blending[s_type],"gordon_hall") == 0)
			p_base_min = s_vol->p_ref;
		else
			EXIT_ERROR("Unsupported: %s\n",sim->geom_blending[s_type]);
		break;
	case ST_SI:
		if (strcmp(sim->geom_blending[s_type],"szabo_babuska_gen") == 0 ||
		    strcmp(sim->geom_blending[s_type],"scott") == 0             ||
		    strcmp(sim->geom_blending[s_type],"lenoir_simple") == 0     ||
		    strcmp(sim->geom_blending[s_type],"nielson") == 0           )
			p_base_min = s_vol->p_ref;
		else if (strcmp(sim->geom_blending[s_type],"lenoir") == 0)
			p_base_min = 2;
		else
			EXIT_ERROR("Unsupported: %s\n",sim->geom_blending[s_type]);
		break;
	case ST_PYR:
	case ST_WEDGE:
		EXIT_ADD_SUPPORT; ///< Ensure that all is working as expected.
	default:
		EXIT_ERROR("Unsupported: %d\n",s_type);
		break;
	}
	assert(p_base_min != 0);

	return p_base_min;
}

static const struct const_Vector_R* constructor_blend_values
	(const int ind_b, const ptrdiff_t n_n, const struct Boundary_Comp_Elem_Data_T*const b_ce_d)
{
	assert(DIM >= 2);
	struct Vector_d*const blend_values = constructor_empty_Vector_d(n_n); // returned
	double*const data_blend = blend_values->data;

	const struct Operator*const vv0_vv_vX = b_ce_d->vv0_vv_vX;
	const struct const_Vector_i*const b_ve_b = b_ce_d->b_ve->data[ind_b];

	switch (b_ce_d->s_type) {
	case ST_TP:
		for (int n = 0; n < n_n; ++n) {
			data_blend[n] = 0.0;
			const Real*const data_b_coords = get_row_const_Matrix_R(n,vv0_vv_vX->op_std);
			for (int ve = 0; ve < b_ve_b->ext_0; ++ve)
				data_blend[n] += data_b_coords[b_ve_b->data[ve]];
		}
		break;
	case ST_SI: {
		const struct Operator*const vv0_bv_vX = b_ce_d->vv0_bv_vX->data[ind_b];
		for (int n = 0; n < n_n; ++n) {
			const Real*const data_b_coords_num = get_row_const_Matrix_R(n,vv0_vv_vX->op_std),
			          *const data_b_coords_den = get_row_const_Matrix_R(n,vv0_bv_vX->op_std);

			Real blend_num = 1.0,
			     blend_den = 1.0;
			for (int ve = 0; ve < b_ve_b->ext_0; ++ve) {
				blend_num *= data_b_coords_num[b_ve_b->data[ve]];
				blend_den *= data_b_coords_den[ve];
			}
			data_blend[n] = ( (blend_num < EPS) ? 0.0 : blend_num/blend_den );
		}
		break;
	} default:
		EXIT_ERROR("Unsupported: %d\n",b_ce_d->s_type);
		break;
	}
	return (struct const_Vector_d*) blend_values;
}

static const struct const_Matrix_R* constructor_xyz_diff_T
	(const struct Boundary_Comp_Elem_Data_T*const b_ce_d, const struct const_Matrix_R*const xyz_i,
	 const struct Solver_Volume_T*const s_vol, const char n_type, const struct Simulation*const sim)
{
	const struct const_Matrix_R*const xyz_surf_diff =
		constructor_xyz_surf_diff_T(b_ce_d,xyz_i,s_vol,n_type,sim); // destructed

	const struct const_Matrix_R*const xyz_diff =
		constructor_mm_const_Matrix_R('N','N',1.0,b_ce_d->vv0_bX_vX->op_std,xyz_surf_diff,'C'); // returned
	destructor_const_Matrix_R(xyz_surf_diff);

	return xyz_diff;
}

static ptrdiff_t find_boundary_vertex_index
	(const struct const_Matrix_R*const xyz, const struct const_Matrix_R*const xyz_ve)
{
	bool found = false;
	ptrdiff_t ind_ve = -1;

	const ptrdiff_t n_n  = xyz->ext_0,
	                n_ve = xyz_ve->ext_0,
	                dim  = xyz->ext_1;
	assert(dim == xyz_ve->ext_1);

	const bool transpose = (xyz_ve->layout != 'R');
	if (transpose)
		transpose_Matrix_R((struct Matrix_R*)xyz_ve,true);

	for (int n = 0; !found && n < n_n; ++n) {
		const double*const data_xyz = get_row_const_Matrix_R(n,xyz);
		for (int ve = 0; ve < n_ve; ++ve) {
			const double*const data_xyz_ve = get_row_const_Matrix_R(ve,xyz_ve);
			const double diff = norm_diff_d(dim,data_xyz_ve,data_xyz,"Inf");
			if (diff < EPS) {
				found = true;
				ind_ve = n;
				break;
			}
		}
	}
	assert(found);

	if (transpose)
		transpose_Matrix_R((struct Matrix_R*)xyz_ve,true);
	return ind_ve;
}
