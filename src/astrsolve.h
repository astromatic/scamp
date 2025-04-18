/*
*				astrsolve.h
*
* Include file astrsolve.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		23/03/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _SOLVE_H_
#define _SOLVE_H_

#include "samples.h"
#include "wcs/poly.h"
#include "fgroup.h"
#include "field.h"

/*----------------------------- Internal constants --------------------------*/

#define	ASTREF_WEIGHTFACTOR	1.0	/* Fudge factor applied to ref.weights*/
#define	ASTROM_REGULFACTOR	0.001	/* Fudge factor applied to regul. */

/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/
extern int	compute_jacobian(samplestruct *sample, double *dprojdred);

extern void	astr_orthopoly(polystruct *poly),
		astrsolve_fgroups(fgroupstruct **fgroups, int nfgroup),
		astrweight_fgroups(fgroupstruct **fgroups, int nfgroup),
		mat_to_wcs(polystruct *poly, polystruct *poly2, double *mat,
				setstruct *set),
		reproj_fgroup(fgroupstruct *fgroup,fieldstruct *field,
				int propflag),
		regul_mat(fgroupstruct **fgroups, int nfgroup,
			double *alpha, int ncoefftot),
		shrink_mat(double *alpha, double *beta, int ncoefftot,
				int index, int nmiss);

#endif // _SOLVE_H_
