// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__test_code_fluxes_h__INCLUDED
#define DPG__test_code_fluxes_h__INCLUDED

extern double *initialize_W   (unsigned int *Nn,      unsigned int *Nel,      const unsigned int d);
extern double **initialize_Q  (unsigned int const Nn, unsigned int const Nel, unsigned int const d);
extern double *initialize_n   (const unsigned int Nn, const unsigned int Nel, const unsigned int d);
extern double *initialize_XYZ (const unsigned int Nn, const unsigned int Nel, const unsigned int d);

extern void set_FiTypes                       (unsigned int *NFiTypesOut, char ***FiTypeOut);
extern void set_FNumTypes                     (unsigned int *NFNumTypesOut, char ***FNumTypeOut);
extern void set_parameters_test_flux_Num      (char const *const FNumType, unsigned int const d);
extern void set_parameters_test_flux_inviscid (char const *const FiType, unsigned int const d);
extern void reset_entered_test_num_flux       (char const *const FNumType);

#endif // DPG__test_code_fluxes_h__INCLUDED
