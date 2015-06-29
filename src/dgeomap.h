/**
* @file		dgeomap.h
* @brief	Include file for dgeomap.c.
* @date		20/01/2015
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2011-2015 IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SCAMP is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SCAMP is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _DGEOMAP_H_
#define _DGEOMAP_H_

//------------------------------- constants ----------------------------------

#define	DGEOMAP_NNEIGHBOURMAX	1025	/// Max number of neighbours used in map
#define	DGEOMAP_RES		32	/// Default diff. geom. map resolution
#define	DGEOMAP_MEMINC		8192	/// dgeopoint memory alloc. increment

//--------------------------------- typedefs ---------------------------------

typedef struct dgeopoint {
   double	pos[NAXIS];			/// Position vector
   double	dpos[NAXIS];			/// Shift vector
   float	weight;				/// Relative weight
}	dgeopointstruct;

/*------------------------------ Prototypes ---------------------------------*/

int	dgeomap_compdx(const void *dgeopoint1, const void *dgeopoint2),
	dgeomap_compdy(const void *dgeopoint1, const void *dgeopoint2),
	dgeomap_instru(fieldstruct **fields, int nfield, int instru,
			char *filename);

#endif

