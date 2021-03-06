// Copyright 2017 Philip Zwanenburg
// MIT License (https://github.com/PhilipZwanenburg/DPGSolver/blob/master/LICENSE)

#ifndef DPG__S_OpCSR_h__INCLUDED
#define DPG__S_OpCSR_h__INCLUDED

struct S_OpCSR {
	unsigned int NRows, NVals, *rowIndex, *columns;
	double       *values;
};

#endif // DPG__S_OpCSR_h__INCLUDED
