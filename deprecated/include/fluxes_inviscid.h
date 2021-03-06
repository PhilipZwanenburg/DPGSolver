// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__fluxes_inviscid_h__INCLUDED
#define DPG__fluxes_inviscid_h__INCLUDED

#include "fluxes_structs.h"

extern void flux_inviscid     (struct S_FLUX *const FLUXDATA);
extern void flux_num_inviscid (struct S_NUMERICALFLUX *const NUMFLUXDATA);
extern void flux_Euler        (struct S_FLUX *const FLUXDATA);

extern void flux_inviscid_MA (struct S_FLUX_MA *const FLUXDATA_MA);

#endif // DPG__fluxes_inviscid_h__INCLUDED
