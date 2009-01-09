 /*
 				check.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include file for check.c.
*
*	Last modify:	05/09/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _CHECK_H_
#define _CHECK_H_

/*----------------------------- Internal constants --------------------------*/
#define	MAXCHECK	64	/* max. # of CHECKimages */

/*--------------------------------- typedefs --------------------------------*/
typedef enum {CHECK_NONE, CHECK_ASPAIR, CHECK_ASREFPAIR, CHECK_ASXCORR,
			CHECK_LLPAIR, CHECK_LLREFPAIR, CHECK_LLXCORR}
		checkenum;

/*------------------------------- functions ---------------------------------*/

extern int	check_check(checkenum checktype);

extern void	write_check(char *filename, float *pix, int width, int height);
#endif
