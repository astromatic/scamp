/*
*				photsolve.h
*
* Include file for photsolve.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		07/02/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef _PHOTSOLVE_H_
#define _PHOTSOLVE_H_

#include "field.h"
#include "samples.h"

/*----------------------------- Internal constants --------------------------*/
#define	PHOTOM_MINMAGERR	0.001	/* Mag error floor for weighting */
/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/
int		photclip_fgroup(fgroupstruct *fgroup, int instru,
			double nsigma);

extern void	avermags_fgroup(fgroupstruct *fgroup),
		compmags_fgroup(fgroupstruct *fgroup),
		photsolve_fgroups(fgroupstruct **fgroups, int nfgroup),
		photstats_fgroup(fgroupstruct *fgroup, int instru,
			double hsn_thresh);

#endif // _PHOTSOLVE_H_
