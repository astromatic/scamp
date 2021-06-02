/*
*				astrsolve.h
*
* Include file astrsolve.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2021 IAP/CNRS/SorbonneU
*
*	License:		GNU General Public License
*
*	SCAMP is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	SCAMP is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		02/06/2021
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _SOLVE_H_
#define _SOLVE_H_

#include "samples.h"
#include "wcs/poly.h"
#include "fgroup.h"
#include "field.h"

//----------------------------- Internal constants --------------------------

#define ASTROM_MAXITER		1000	// Maximum number of solution iterations
#define ASTROM_UPDATEFACTOR	1e-12	// Gradient factor at each iteration
#define ASTREF_WEIGHTFACTOR	1.0	// Fudge factor for reference weights
#define ASTROM_REGULFACTOR	0.001	// Fudge factor for regularization

//-------------------------------- flags ------------------------------------

#define REPROJ_NONE		0x0000	// No extra computation
#define REPROJ_PROPER_MOTION	0x0001	// Compute proper motions
#define REPROJ_JACOBIAN		0x0002	// Compute Jacobians


//--------------------------- structure definitions -------------------------
//---------------------------------- protos --------------------------------
extern int	compute_jacobian(samplestruct *sample);

extern void	astr_orthopoly(polystruct *poly),
		astrsolve_fgroups(fgroupstruct **fgroups,
				fieldstruct **reffields, int nfgroup),
		astrweight_fgroups(fgroupstruct **fgroups, int nfgroup),
		mat_to_wcs(polystruct *poly, polystruct *poly2, double *mat,
				setstruct *set),
		reproj_fgroup(fgroupstruct *fgroup,fieldstruct *field,
				int flags),
		regul_mat(fgroupstruct **fgroups, int nfgroup,
			double *alpha, int ncoefftot),
		shrink_mat(double *alpha, double *beta, int ncoefftot,
				int index, int nmiss);

#endif // _SOLVE_H_
