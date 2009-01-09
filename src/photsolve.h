 /*
 				photsolve.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for photsolve.c.
*
*	Last modify:	25/05/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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

