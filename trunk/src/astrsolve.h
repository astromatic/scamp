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
*	Last modified:		25/05/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

#ifndef _POLY_H_
#include "wcs/poly.h"
#endif

#ifndef _SAMPLES_H_
#include "samples.h"
#endif

#ifndef _SOLVE_H_
#define _SOLVE_H_

/*----------------------------- Internal constants --------------------------*/

#define	ASTREF_WEIGHTFACTOR	1.0	/* Fudge factor applied to ref.weights*/

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
		shrink_mat(double *alpha, double *beta, int ncoefftot,
				int index, int nmiss);

#endif

