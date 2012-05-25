/*
*				astrstats.h
*
* Include file for astrstats.c.
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
*	Last modified:		25/05/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

#ifndef _ASTRSTAT_H_
#define _ASTRSTAT_H_

/*----------------------------- Internal constants --------------------------*/
/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/
extern int	astrclip_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
				double nsigma);

extern void	astrstats_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
				double hsn_thresh);

#endif

