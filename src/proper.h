/*
*				proper.h
*
* Include file for proper.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2008-2020 IAP/CNRS/SorbonneU
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

#ifndef _PROPER_H_
#define _PROPER_H_

#include "fitswcs.h"

/*--------------------------------- constants -------------------------------*/

#define	PROPER_MAXNBAD		2	/* Maximum number of bad measurements */
#define	PROPER_MAXFRACBAD	0.2	/* Max. fraction of bad measurements */
#define	PROPER_MAXCHI2		6.0	/* Trigger chi2/d.o.f. threshold for */
					/* source clipping */

/*--------------------------------- typedefs --------------------------------*/


/*--------------------------- structure definitions -------------------------*/


/*-------------------------------- protos -----------------------------------*/
extern int	astrpropclip_fgroup(fgroupstruct *fgroup, double maxpropmod);

extern void	astrcolshift_fgroup(fgroupstruct *fgroup, fieldstruct *reffield),
		astrconnect_fgroup(fgroupstruct *fgroup),
		astrprop_fgroup(fgroupstruct *fgroup);

#endif // _PROPER_H_
