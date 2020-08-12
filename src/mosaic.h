/*
*				mosaic.h
*
* Include file for mosaic.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		12/08/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef _MOSAIC_H_
#define _MOSAIC_H_

#include "field.h"
#include "fitswcs.h"

/*------------------------------- functions ---------------------------------*/

extern void	adjust_mosaic(fieldstruct **fields, int nfield),
		adjust_set(fieldstruct **fields, int nfield, int s),
		crval_to_crpix(wcsstruct *wcs, double *wcspos);	
#ifdef USE_THREADS
extern void	*pthread_adjust_set(void *arg),
		pthread_adjust_sets(fieldstruct **fields, int nfield, int l);
#endif // USE_THREADS
#endif // _MOSAIC_H_
