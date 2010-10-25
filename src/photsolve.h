/*
*				photsolve.h
*
* Include file for photsolve.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _SAMPLES_H_
#include "samples.h"
#endif

#ifndef _PHOTSOLVE_H_
#define _PHOTSOLVE_H_

/*----------------------------- Internal constants --------------------------*/
/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/
int		photclip_fgroup(fgroupstruct *fgroup, int instru,
			double nsigma);

extern void	avermags_fgroup(fgroupstruct *fgroup),
		compmags_fgroup(fgroupstruct *fgroup),
		photsolve_fgroups(fgroupstruct **fgroups, int nfgroup),
		photstats_fgroup(fgroupstruct *fgroup, int instru,
			double hsn_thresh);

#endif

