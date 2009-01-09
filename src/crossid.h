 /*
 				crossid.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	include for crossid.c.
*
*	Last modify:	25/02/2005
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

/*---------------------------- Internal constants ---------------------------*/

/*--------------------------------- typedefs --------------------------------*/
/*------------------------------- functions ---------------------------------*/

extern int	check_fieldoverlap(fieldstruct *field1, fieldstruct *field2),
		check_fieldphotomoverlap(fieldstruct *field, int instru);

extern void	crossid_fgroup(fgroupstruct *fgroup,
				fieldstruct *reffield, double tolerance);

