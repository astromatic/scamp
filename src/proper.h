/*
 				proper.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP, Leiden observatory & ESO)
*
*	Contents:	Declarations related to proper motions
*
*	Last modify:	06/03/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _PROPER_H_
#define _PROPER_H_

/*--------------------------------- constants -------------------------------*/

/*--------------------------------- typedefs --------------------------------*/


/*--------------------------- structure definitions -------------------------*/


/*-------------------------------- protos -----------------------------------*/

void	astrcolshift_fgroup(fgroupstruct *fgroup, fieldstruct *reffield),
	astrprop_fgroup(fgroupstruct *fgroup);

#endif
