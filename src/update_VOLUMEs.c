// Copyright 2016 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/master/LICENSE)

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "database.h"
#include "parameters.h"
#include "functions.h"

/*
 *	Purpose:
 *		Update VOLUME related information/operators in ELEMENTs which have undergone hp refinement.
 *
 *	Comments:
 *		May be better not to allow P adaptivity to go down to P = 0, but have the limit at P = 1. Try/THINK (ToBeDeleted).
 *		Sum factorization not currently used here. Profile and re-evaluate. (ToBeDeleted)
 *
 *	Notation:
 *
 *	References:
 */

struct S_OPERATORS {
	unsigned int NvnGs, NvnGc, NvnS, NvnSP, NvnI;
	double       *I_vGs_vGc, **I_vGs_vGs, **Ihat_vS_vS, **Ghat_vS_vS, *w_vI, *ChiS_vI;
};

static void init_ops(struct S_OPERATORS *OPS, const struct S_VOLUME *VOLUME, const unsigned int IndClass)
{
	unsigned int P, PNew, type, curved;
	struct S_ELEMENT *ELEMENT;

	P      = VOLUME->P;
	PNew   = VOLUME->PNew;
	type   = VOLUME->type;
	curved = VOLUME->curved;

	ELEMENT = get_ELEMENT_type(type);

	OPS->NvnGs = ELEMENT->NvnGs[1];
	OPS->NvnGc = ELEMENT->NvnGc[PNew];
	OPS->NvnS  = ELEMENT->NvnS[P];
	OPS->NvnSP = ELEMENT->NvnS[PNew];
	OPS->I_vGs_vGs  = ELEMENT->I_vGs_vGs[1][1];
	OPS->I_vGs_vGc  = ELEMENT->I_vGs_vGc[1][PNew][0];
	OPS->Ihat_vS_vS = ELEMENT->Ihat_vS_vS[P][PNew]; // ToBeDeleted: Remove all instances of Ihat_vS_vS from the code if
	                                                //              not used here.
	OPS->Ghat_vS_vS = ELEMENT->Ghat_vS_vS[P][PNew];
	if (!curved) {
		OPS->NvnI = ELEMENT->NvnIs[P];

		OPS->w_vI    = ELEMENT->w_vIs[P];
		OPS->ChiS_vI = ELEMENT->ChiS_vIs[P][P][0];
	} else {
		OPS->NvnI  = ELEMENT->NvnIc[P];

		OPS->w_vI    = ELEMENT->w_vIc[P];
		OPS->ChiS_vI = ELEMENT->ChiS_vIc[P][P][0];
	}
}

static unsigned int get_VOLUMEc_type(const unsigned int VType, const unsigned int vh)
{
	switch (VType) {
	case TET:
		if (vh < 4)
			return TET;
		else
			return PYR;
		break;
	case PYR:
		if (vh < 4 || vh > 7)
			return PYR;
		else
			return TET;
		break;
	default:
		printf("Error: Unsupported VType in get_VOLUMEc_type.\n"), exit(1);
		break;
	}
}

void update_VOLUME_hp(void)
{
	// Initialize DB Parameters
	unsigned int d         = DB.d,
	             NV        = DB.NV,
	             AC        = DB.AC,
	             PMax      = DB.PMax,
				 LevelsMax = DB.LevelsMax,
	             Nvar      = DB.Nvar;

	char         *MeshType = DB.MeshType;

	// Standard datatypes
	unsigned int P, PNew, f, level, adapt_type, vh, vhMin, vhMax, href_type, VType, Nf,
	             NvnGs, NvnGc, NvnS, NvnSP, NCols, update, maxP;
	double       *I_vGs_vGc, *XYZ_vC, *XYZ_S,
	             **Ihat_vS_vS, **I_vGs_vGs, **Ghat_vS_vS, *What, *RES, *WhatP, *RESP;

	struct S_OPERATORS *OPS;
	struct S_ELEMENT   *ELEMENT;
	struct S_VOLUME    *VOLUME, *VOLUMEc, *VOLUMEp;

	// silence
	I_vGs_vGs = NULL;
	VOLUMEc   = NULL;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
//printf("upV: %d %d %d\n",VOLUME->indexg,VOLUME->Vadapt,VOLUME->update);
		if (VOLUME->Vadapt) {
			P     = VOLUME->P;
			level = VOLUME->level;
			adapt_type = VOLUME->adapt_type;

			switch(adapt_type) {
			case PREFINE:
				if (P < PMax)
					PNew = P+1;
				else
					printf("Error: Should not be entering PREFINE in update_VOLUME_hp for P = %d.\n",P), exit(1);
				VOLUME->PNew = PNew;
				break;
			case PCOARSE:
				if (P >= 1)
					PNew = P-1;
				else
					printf("Error: Should not be entering PCOARSE in update_VOLUME_hp for P = %d.\n",P), exit(1);
				VOLUME->PNew = PNew;
				break;
			case HREFINE:
				if (level >= LevelsMax)
					printf("Error: Should not be entering HREFINE in update_VOLUME_hp for level = %d.\n",level), exit(1);
				VOLUME->PNew = P;
				break;
			case HCOARSE:
				if (level == 0)
					printf("Error: Should not be entering HCOARSE in update_VOLUME_hp for level = %d.\n",level), exit(1);
				VOLUME->PNew = P;
				break;
			default:
				printf("Error: Unsupported adapt_type = %d in update_VOLUME_hp.\n",adapt_type), exit(1);
				break;
			}

			init_ops(OPS,VOLUME,0);

			VOLUME->update = 1;
			switch(adapt_type) {
			default: // PREFINE or PCOARSE
				VOLUME->P = PNew;

				// Update geometry
				if (VOLUME->curved) {
					NvnGs      = OPS->NvnGs;
					NvnGc      = OPS->NvnGc;
					I_vGs_vGc  = OPS->I_vGs_vGc;

					NCols = d;

					XYZ_vC = VOLUME->XYZ_vC;
					XYZ_S  = malloc(NvnGc*NCols * sizeof *XYZ_S); // keep
					mm_CTN_d(NvnGc,NCols,NvnGs,I_vGs_vGc,XYZ_vC,XYZ_S);

					free(VOLUME->XYZ_S);
					VOLUME->XYZ_S = XYZ_S;
					VOLUME->NvnG  = NvnGc;
				}

				free(VOLUME->XYZ);
				if (strstr(MeshType,"ToBeCurved") != NULL)
					setup_ToBeCurved(VOLUME);
				else
					printf("Error: Add in support for MeshType != ToBeCurved\n"), exit(1);

				free(VOLUME->detJV_vI);
				free(VOLUME->C_vI);
				setup_geom_factors(VOLUME);

				// Project What and RES
				NvnS       = OPS->NvnS;
				NvnSP      = OPS->NvnSP;
				Ihat_vS_vS = OPS->Ihat_vS_vS;
				Ghat_vS_vS = OPS->Ghat_vS_vS;

				VOLUME->NvnS = NvnSP;

				What  = VOLUME->What;
				RES   = VOLUME->RES;

				WhatP = malloc(NvnSP*Nvar * sizeof *WhatP); // keep
				RESP  = malloc(NvnSP*Nvar * sizeof *RESP);  // keep

//				mm_CTN_d(NvnSP,Nvar,NvnS,Ihat_vS_vS[0],What,WhatP);
//				mm_CTN_d(NvnSP,Nvar,NvnS,Ihat_vS_vS[0],RES,RESP);
				mm_CTN_d(NvnSP,Nvar,NvnS,Ghat_vS_vS[0],What,WhatP);
				mm_CTN_d(NvnSP,Nvar,NvnS,Ghat_vS_vS[0],RES,RESP);

				free(What);
				free(RES);

				VOLUME->What = WhatP;
				VOLUME->RES  = RESP;
				break;
			case HREFINE:
				VType = VOLUME->type;


				NvnGs     = OPS->NvnGs;
				NvnGc     = OPS->NvnGc;
				I_vGs_vGs = OPS->I_vGs_vGs;
				I_vGs_vGc = OPS->I_vGs_vGc;

				NCols = d;


				href_type = VOLUME->hrefine_type;

				get_vh_range(VOLUME,&vhMin,&vhMax);
				for (vh = vhMin; vh <= vhMax; vh++) {
					if (vh == vhMin) {
						VOLUMEc = New_VOLUME();
						VOLUME->child0 = VOLUMEc;
					} else {
						VOLUMEc->next = New_VOLUME();
						VOLUMEc = VOLUMEc->next;
					}
					VOLUMEc->update = 1;
					VOLUMEc->parent = VOLUME;
					VOLUMEc->indexg = NV++;

					VOLUMEc->P = VOLUME->P;
					VOLUMEc->level = (VOLUME->level)+1;
					switch (VType) {
					default: // LINE, TRI, QUAD, HEX, WEDGE
						VOLUMEc->type = VType;
						break;
					case TET:
					case PYR:
						VOLUMEc->type = get_VOLUMEc_type(VType,vh);
						break;
					}
					ELEMENT = get_ELEMENT_type(VOLUMEc->type);
					Nf = ELEMENT->Nf;
					for (f = 0; f < Nf; f++)
						VOLUMEc->NsubF[f] = 1;

					VOLUMEc->Eclass = get_Eclass(VOLUMEc->type);

					if (AC) {
						VOLUMEc->curved = 1;
					} else if (VOLUME->curved) {
						printf("Error: Add support for h-refinement VOLUMEc->curved.\n"), exit(1);
						// Use VToBC and knowledge of whether the new VOLUME shares the BC.
					} else {
						VOLUMEc->curved = 0;
					}

					// Update geometry
// When updating XYZ_vC, ensure that corners on curved boundaries are placed on the boundary.
					VOLUMEc->XYZ_vC = malloc(NvnGs*d * sizeof *XYZ_vC); // keep
					XYZ_vC = VOLUMEc->XYZ_vC;

					mm_CTN_d(NvnGs,NCols,NvnGs,I_vGs_vGs[vh],VOLUME->XYZ_vC,VOLUMEc->XYZ_vC);
					if (!VOLUMEc->curved) {
						double *XYZ;

						VOLUMEc->NvnG = NvnGs;

						VOLUMEc->XYZ_S = malloc(NvnGs*NCols * sizeof *XYZ_S); // keep
						VOLUMEc->XYZ   = malloc(NvnGs*NCols * sizeof *XYZ);   // keep
						XYZ_S = VOLUMEc->XYZ_S;
						XYZ   = VOLUMEc->XYZ;
						for (unsigned int i = 0, iMax = NCols*NvnGs; i < iMax; i++) {
							XYZ_S[i] = XYZ_vC[i];
							XYZ[i]   = XYZ_S[i];
						}
					} else {
						VOLUMEc->NvnG = NvnGc;

						VOLUMEc->XYZ_S = malloc(NvnGc*NCols * sizeof *XYZ_S); // keep
						mm_CTN_d(NvnGc,NCols,NvnGs,I_vGs_vGc,XYZ_vC,VOLUMEc->XYZ_S);

						if (strstr(MeshType,"ToBeCurved") != NULL) {
							setup_ToBeCurved(VOLUMEc);
						} else {
							printf("Add in support for MeshType != ToBeCurved");
							exit(1);
						}
					}
					setup_geom_factors(VOLUMEc);
				}

				// Fix VOLUME linked list and Vgrp linked list
				// Also update indexg and indexl

				// Project What and RES
				// free VOLUME->What and VOLUME->RES after projection.
				break;
			case HCOARSE:
				VOLUMEp = VOLUME->parent;

				if (VOLUME == VOLUMEp->child0) {
					level = VOLUME->level;

					update = 1;
					maxP = VOLUME->P;

					get_vh_range(VOLUMEp,&vhMin,&vhMax);
					VOLUMEc = VOLUME;
					for (vh = vhMin+1; vh <= vhMax; vh++) {
						VOLUMEc = VOLUMEc->next;
						maxP = max(maxP,VOLUMEc->P);

						if (VOLUMEc->level != level || !VOLUMEc->Vadapt) {
							update = 0;
							break;
						}
					}
					VOLUMEp->update = update;

					if (update) {
						// Most of the information of VOLUMEp has been stored.
						VOLUMEp->update = 1;
						VOLUMEp->indexg = NV++;

						VOLUMEp->P = maxP;
						VOLUMEp->adapt_type = HCOARSE;
						// Project What and RES
						// free all child VOLUMEs when finished.
					} else {
						// Ensure that all children are marked as not to be updated.
						VOLUMEc = VOLUME;
						for (vh = vhMin; vh <= vhMax; vh++) {
							VOLUMEc->Vadapt = 0;
							VOLUMEc->update = 0;
							VOLUMEc = VOLUMEc->next;
						}
					}
				} else {
					VOLUMEc = VOLUMEp->child0;
					if (!(VOLUMEc->adapt_type == HCOARSE && VOLUMEc->Vadapt)) {
//						VOLUME->Vadapt = 0;
						VOLUME->update = 0;
					}
				}


				break;
			}
		}
//printf("upV: %d %d %d\n",VOLUME->indexg,VOLUME->Vadapt,VOLUME->update);
	}
	free(OPS);
}

void update_VOLUME_list(void)
{
	/*
	 *	Comments:
	 *		This is done inelegantly because it was required to keep the pointer to the parent VOLUME.
	 */

	unsigned int adapt_type;

	struct S_VOLUME *VOLUME, *VOLUMEc, *VOLUMEnext;

	// Fix list head if necessary
	VOLUME = DB.VOLUME;

	if (VOLUME->update) {
		adapt_type = VOLUME->adapt_type;
		if (adapt_type == HREFINE) {
			DB.VOLUME = VOLUME->child0;
			for (VOLUMEc = DB.VOLUME; VOLUMEc->next != NULL; VOLUMEc = VOLUMEc->next)
				;
			VOLUMEc->next = VOLUME->next;
		} else if (adapt_type == HCOARSE) {
			DB.VOLUME = VOLUME->parent;
			for (VOLUMEc = VOLUME; VOLUMEc->next->parent == DB.VOLUME; VOLUMEc = VOLUMEc->next)
				VOLUMEc->update = 0;
			DB.VOLUME->next = VOLUMEc->next;
			VOLUMEc->next = NULL;
		}
	}

	// Fix remainder of list
	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		VOLUMEnext = VOLUME->next;
		if (VOLUMEnext && VOLUMEnext->update) {
			adapt_type = VOLUMEnext->adapt_type;
			if (adapt_type == HREFINE) {
				VOLUME->next = VOLUMEnext->child0;
				for (VOLUMEc = VOLUME->next; VOLUMEc->next != NULL; VOLUMEc = VOLUMEc->next)
					;
				VOLUMEc->next = VOLUMEnext->next;
			} else if (adapt_type == HCOARSE) {
				VOLUME->next = VOLUMEnext->parent;
				for (VOLUMEc = VOLUMEnext; VOLUMEc->next && VOLUMEc->next->parent == VOLUME->next; VOLUMEc = VOLUMEc->next)
					VOLUMEc->update = 0;
				VOLUME->next->next = VOLUMEc->next;
				VOLUMEc->next = NULL;
			}
		}
	}

	for (VOLUME = DB.VOLUME; VOLUME; VOLUME = VOLUME->next) {
		VOLUME->Vadapt = 0;
		VOLUME->update = 0;
	}
}

void update_Vgrp(void)
{
	// Initialize DB Parameters
	unsigned int NP     = DB.NP,
	             NTVgrp = DB.NTVgrp,
	             *NVgrp = DB.NVgrp;

	// Standard datatypes
	unsigned int i, P, Eclass, curved, IndVgrp;

	struct S_VOLUME *VOLUME, **VOLUME_prev;

	if (DB.update) {
		DB.update = 0;

		VOLUME_prev = malloc(NTVgrp * sizeof *VOLUME_prev); // free

		for (i = 0; i < NTVgrp; i++)
			NVgrp[i] = 0;

		for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
			Eclass = VOLUME->Eclass;
			P      = VOLUME->P;
			curved = VOLUME->curved;

			IndVgrp = Eclass*NP*2 + P*2 + curved;

			if (!NVgrp[IndVgrp])
				DB.Vgrp[IndVgrp] = VOLUME;
			else
				VOLUME_prev[IndVgrp]->grpnext = VOLUME;

			NVgrp[IndVgrp]++;
			VOLUME_prev[IndVgrp] = VOLUME;
		}
		free(VOLUME_prev);
	}
}

void update_VOLUME_Ops(void)
{
	// Initialize DB Parameters
	char         *SolverType = DB.SolverType;
	unsigned int Collocated  = DB.Collocated;

	// Standard datatypes
	unsigned int iMax, jMax,
	             NvnS, NvnI;
	double       *detJV_vI, *detJV_vI_ptr, *w_vI, *w_vI_ptr, *wdetJV_vI, *wdetJV_vI_ptr, *ChiS_vI, *ChiS_vI_ptr,
	             *wdetJVChiS_vI, *wdetJVChiS_vI_ptr, *IS, *M, *MInv;

	struct S_OPERATORS *OPS;
	struct S_VOLUME    *VOLUME;

	OPS = malloc(sizeof *OPS); // free

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		if (VOLUME->update) {
			VOLUME->update = 0;
//printf("update_VOLUME_Ops: %d\n",VOLUME->indexg);
			if (strstr(SolverType,"Explicit") != NULL) {
				if (!Collocated) {
					init_ops(OPS,VOLUME,0);

					NvnS = OPS->NvnS;
					NvnI = OPS->NvnI;
					w_vI    = OPS->w_vI;
					ChiS_vI = OPS->ChiS_vI;

					detJV_vI = VOLUME->detJV_vI;

					// Compute required portion of MInv
					wdetJV_vI = malloc(NvnI * sizeof *wdetJV_vI); // free
					MInv      = NULL;

					w_vI_ptr      = w_vI;
					detJV_vI_ptr  = detJV_vI;
					wdetJV_vI_ptr = wdetJV_vI;
					for (iMax = NvnI; iMax--; )
						*wdetJV_vI_ptr++ = (*w_vI_ptr++)*(*detJV_vI_ptr++);

					wdetJVChiS_vI = malloc(NvnI*NvnS * sizeof *wdetJVChiS_vI); // free

					ChiS_vI_ptr       = ChiS_vI;
					wdetJVChiS_vI_ptr = wdetJVChiS_vI;
					for (iMax = NvnI*NvnS; iMax--; )
						*wdetJVChiS_vI_ptr++ = *ChiS_vI_ptr++;

					wdetJV_vI_ptr     = wdetJV_vI;
					wdetJVChiS_vI_ptr = wdetJVChiS_vI;
					for (iMax = NvnI; iMax--; ) {
						for (jMax = NvnS; jMax--; )
							*wdetJVChiS_vI_ptr++ *= *wdetJV_vI_ptr;
						wdetJV_vI_ptr++;
					}

					M    = mm_Alloc_d(CBRM,CBT,CBNT,NvnS,NvnS,NvnI,1.0,ChiS_vI,wdetJVChiS_vI); // free
					IS   = identity_d(NvnS); // free
					MInv = inverse_d(NvnS,NvnS,M,IS); // keep
//array_print_d(NvnS,NvnS,M,'R');
//array_print_d(NvnI,NvnS,wdetJVChiS_vI,'R');
//exit(1);

					free(wdetJVChiS_vI);
					free(M);
					free(IS);
					free(wdetJV_vI);

					free(VOLUME->MInv);
					VOLUME->MInv = MInv;
				}
			} else {
				; // Updates for implicit.
			}
		}
	}

	free(OPS);
}

void update_VOLUME_finalize(void)
{
	struct S_VOLUME *VOLUME;

	for (VOLUME = DB.VOLUME; VOLUME != NULL; VOLUME = VOLUME->next) {
		if (VOLUME->Vadapt) {
			VOLUME->Vadapt = 0;
		}
	}
}