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

#include "macros.h"

#include "def_templates_matrix.h"
#include "def_templates_multiarray.h"
#include "def_templates_vector.h"

#include "def_templates_compute_all_rlhs_dpg.h"
#include "def_templates_compute_volume_rlhs.h"
#include "def_templates_volume_solver.h"
#include "def_templates_volume_solver_dpg.h"

// Static function declarations ************************************************************************************* //

/** \brief Constructor for the H0 norm operator of the input volume.
 *  \return See brief. */
static const struct const_Matrix_T* constructor_norm_op_H0
	(const struct DPG_Solver_Volume_T* dpg_s_vol ///< \ref DPG_Solver_Volume_T.
	);

/** \brief Constructor for the H1 norm operator of the input volume.
 *  \return See brief. */
static const struct const_Matrix_T* constructor_norm_op_H1
	(const struct DPG_Solver_Volume_T* dpg_s_vol ///< \ref DPG_Solver_Volume_T.
	);

// Interface functions ********************************************************************************************** //

void constructor_derived_DPG_Solver_Volume_T (struct Volume* volume_ptr, const struct Simulation* sim)
{
	struct DPG_Solver_Volume_T* dpg_s_vol = (struct DPG_Solver_Volume_T*) volume_ptr;
	UNUSED(sim);

	dpg_s_vol->norm_op_H0 = constructor_norm_op_H0(dpg_s_vol); // destructed
	dpg_s_vol->norm_op_H1 = constructor_norm_op_H1(dpg_s_vol); // destructed
}

void destructor_derived_DPG_Solver_Volume_T (struct Volume* volume_ptr)
{
	struct DPG_Solver_Volume_T* dpg_s_vol = (struct DPG_Solver_Volume_T*) volume_ptr;

	destructor_const_Matrix_T(dpg_s_vol->norm_op_H0);
	destructor_const_Matrix_T(dpg_s_vol->norm_op_H1);
}

// Static functions ************************************************************************************************* //
// Level 0 ********************************************************************************************************** //

static const struct const_Matrix_T* constructor_norm_op_H0 (const struct DPG_Solver_Volume_T* dpg_s_vol)
{
	struct Solver_Volume_T*const s_vol = (struct Solver_Volume_T*) dpg_s_vol;

	const struct Operator*const cv0_vt_vc = get_operator__cv0_vt_vc_T(s_vol);
	const struct const_Vector_R* w_vc = get_operator__w_vc__s_e_T(s_vol);

	const struct const_Vector_T jacobian_det_vc = interpret_const_Multiarray_as_Vector_T(s_vol->jacobian_det_vc);
	const struct const_Vector_T* wJ_vc = constructor_dot_mult_const_Vector_T_RT(1.0,w_vc,&jacobian_det_vc,1); // destructed

	const struct const_Matrix_R* H0_l = cv0_vt_vc->op_std;
	const struct const_Matrix_T* H0_r = constructor_mm_diag_const_Matrix_R_T(1.0,H0_l,wJ_vc,'L',false); // destructed
	destructor_const_Vector_T(wJ_vc);

	const struct const_Matrix_T* H0 = constructor_mm_RT_const_Matrix_T('T','N',1.0,H0_l,H0_r,'R'); // returned
	destructor_const_Matrix_T(H0_r);

	return H0;
}

static const struct const_Matrix_T* constructor_norm_op_H1 (const struct DPG_Solver_Volume_T* dpg_s_vol)
{
	struct Solver_Volume_T* s_vol = (struct Solver_Volume_T*) dpg_s_vol;

	const struct Multiarray_Operator cv1_vt_vc = get_operator__cv1_vt_vc_T(s_vol);
	const struct const_Vector_R* w_vc = get_operator__w_vc__s_e_T(s_vol);
	const struct const_Vector_T J_vc  = interpret_const_Multiarray_as_Vector_T(s_vol->jacobian_det_vc);

	const struct const_Vector_T* J_inv_vc = constructor_inverse_const_Vector_T(&J_vc);                // destructed
	const struct const_Vector_T* wJ_vc    = constructor_dot_mult_const_Vector_T_RT(1.0,w_vc,J_inv_vc,1); // destructed
	destructor_const_Vector_T(J_inv_vc);

	struct Matrix_T* H1 = (struct Matrix_T*) constructor_norm_op_H0(dpg_s_vol); // returned

	for (int d = 0; d < DIM; ++d) {
		const struct const_Matrix_R* H1_l = cv1_vt_vc.data[d]->op_std;
		const struct const_Matrix_T* H1_r =
			constructor_mm_diag_const_Matrix_R_T(1.0,H1_l,wJ_vc,'L',false); // destructed
		mm_RTT('T','N',1.0,1.0,H1_l,H1_r,H1);
		destructor_const_Matrix_T(H1_r);
	}
	destructor_const_Vector_T(wJ_vc);

	return (struct const_Matrix_T*) H1;
}

#include "undef_templates_matrix.h"
#include "undef_templates_multiarray.h"
#include "undef_templates_vector.h"

#include "undef_templates_compute_all_rlhs_dpg.h"
#include "undef_templates_compute_volume_rlhs.h"
#include "undef_templates_volume_solver.h"
#include "undef_templates_volume_solver_dpg.h"
