 /*
 				mosaic.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	include for mosaic.c.
*
*	Last modify:	27/09/2004
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _MOSAIC_H_
#define _MOSAIC_H_

/*------------------------------- functions ---------------------------------*/

extern void	adjust_mosaic(fieldstruct **fields, int nfield),
		adjust_set(fieldstruct **fields, int nfield, int s),
		crval_to_crpix(wcsstruct *wcs, double *wcspos);	
#ifdef USE_THREADS
void		*pthread_adjust_set(void *arg),
		pthread_adjust_sets(fieldstruct **fields, int nfield, int l);
#endif
#endif
