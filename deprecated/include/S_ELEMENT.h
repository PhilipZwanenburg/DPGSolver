// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__S_ELEMENT_h__INCLUDED
#define DPG__S_ELEMENT_h__INCLUDED

#include <stdbool.h>

#include "S_OpCSR.h"
#include "matrix_structs.h"

#include "simulation.h"


/** \brief Struct holding data related to the base Element. */
struct S_ELEMENT {
	// New
	const bool tensor_product;

	// Mesh
	unsigned int present, type, d, Nve, NveP2, Ne, Nf, Nvref, NvrefSF, Eclass, NEhref,
	             Neve, *Nfve, *VeCGmsh, *VeEcon, *VeFcon, *NrefV, *type_h;

	// Operators
	unsigned int *connect_NE, *NvnP, *Nvve,
	             *NvnGs, *NvnG2, *NvnGc, *NvnCs, *NvnCc, *NvnJs, *NvnJc, *NvnIs, *NvnIc, *NvnS,
	             **NfnG2, **NfnGc, **NfnS, **NfnIs, **NfnIc,
	             *NenG2, *NenGc,
	             Neref, *Nfref, *NfMixed,
	             **connectivity, **connect_types, **connectivityE,
	             ***nOrd_fS, ***nOrd_fIs, ***nOrd_fIc, ****Fmask, ****VeMask;
	double       volume, **VeE, **VeF, **VeV, *nr,
	             **w_vIs, **w_vIc, ***w_fIs, ***w_fIc,
	             ****ChiS_vP, ****ChiS_vS, ****ChiS_vIs, ****ChiS_vIc,
	             ****ChiInvS_vS, ****ChiInvGs_vGs, ****ChiBezInvS_vS,
	             ****IG2, ****IGc, ****ICs, ****ICc,
	             ****TGs, ****TS, ****TS_vB, ****TInvS_vB,
	             ****I_vGs_vP, ****I_vGs_vGs, ****I_vGs_vG2, ****I_vGs_vGc, ****I_vGs_vCs, ****I_vGs_vS, ****I_vGs_vIs, ****I_vGs_vIc,
	             ****I_vGc_vP,                               ****I_vGc_vCc, ****I_vGc_vS, ****I_vGc_vIs, ****I_vGc_vIc,
	             ****I_vCs_vS, ****I_vCs_vIs, ****I_vCs_vIc,
	             ****I_vCc_vS, ****I_vCc_vIs, ****I_vCc_vIc,
	             ****Ihat_vS_vS,
	             *****GradChiS_vS, *****GradChiS_vIs, *****GradChiS_vIc,
	             *****D_vGs_vCs, *****D_vGs_vIs,
	             *****D_vGc_vCc, *****D_vGc_vIc,
	             *****D_vCs_vCs,
	             *****D_vCc_vCc,
	             ****ChiS_fS, ****ChiS_fIs, ****ChiS_fIc,
	             *****GradChiS_fIs, *****GradChiS_fIc,
	             ****I_vGs_fS, ****I_vGs_fIs, ****I_vGs_fIc,
	             ****I_vGc_fGc, ****I_vG2_fG2, ****I_vGc_fS, ****I_vGc_fIs, ****I_vGc_fIc,
	             ****I_vGc_eGc, ****I_vG2_eG2,
	             ****I_vCs_fS, ****I_vCs_fIs, ****I_vCs_fIc,
	             ****I_vCc_fS, ****I_vCc_fIs, ****I_vCc_fIc,
	             *****D_vGs_fIs, *****D_vGs_fIc,
	             *****D_vGc_fIs, *****D_vGc_fIc,
	             ****I_fGs_vGc, ****I_fGc_vGc, ****I_fGs_vG2, ****I_fG2_vG2,
	             ****I_eGs_vGc, ****I_eGc_vGc, ****I_eGs_vG2, ****I_eG2_vG2,
	             ****Is_Weak_VV, ****Ic_Weak_VV, ****Is_Strong_VV, ****Ic_Strong_VV,
	             ****Is_Weak_FV, ****Ic_Weak_FV,
	             *****Ds_Weak_VV, *****Dc_Weak_VV, *****Ds_Strong_VV, *****Dc_Strong_VV,
	             ****L2hat_vS_vS,
	             ****GfS_fIs, ****GfS_fIc;

	// S_OpCSR should be converted to S_MATRIX in the future (ToBeModified)
	struct S_OpCSR ****ChiS_fIs_sp, ****ChiS_fIc_sp,
	               *****Ds_Weak_VV_sp, *****Dc_Weak_VV_sp,
	               ****Is_Weak_FV_sp, ****Ic_Weak_FV_sp;

	struct S_OPS {
		struct S_OPS_SOLVER {
			struct S_OPS_SOLVER_DG {
				// VOLUME
				struct S_MATRIX ****ChiS_vIs, ****ChiS_vIc,
				                *****Ds_Weak_VV, *****Dc_Weak_VV,
				                ****I_vGs_vIs, ****I_vGc_vIc;

				// FACE
				struct S_MATRIX ****ChiS_fIs, ****ChiS_fIc,
				                ****Is_Weak_FV, ****Ic_Weak_FV;
			} DG;

			struct S_OPS_SOLVER_HDG {
				struct S_MATRIX **ChiTRS_vIs, **ChiTRS_vIc,
				                **Is_FF, **Ic_FF;
			} HDG;
		} solver;
	} ops;

	struct S_ELEMENT *next;
	struct S_ELEMENT **ELEMENTclass, **ELEMENT_FACE;
};

struct S_ELEMENT* get_element_by_type2 (struct S_ELEMENT*const element_head, const size_t type);
const struct S_ELEMENT* get_const_element_by_type (const struct S_ELEMENT*const element_head, const size_t type);

///< \brief Set up operators for the base Element.
void set_up_operators
	(const struct Simulation*const simulation ///< Standard.
	);

#endif // DPG__S_ELEMENT_h__INCLUDED
