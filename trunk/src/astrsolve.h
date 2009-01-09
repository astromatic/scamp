 /*
 				astrsolve.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for solve.c.
*
*	Last modify:	30/10/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

#ifndef _POLY_H_
#include "wcs/poly.h"
#endif

#ifndef _SAMPLES_H_
#include "samples.h"
#endif

#ifndef _SOLVE_H_
#define _SOLVE_H_

/*----------------------------- Internal constants --------------------------*/

#define	ASTREF_WEIGHTFACTOR	1.0	/* Fudge factor applied to ref.weights*/

/*--------------------------- structure definitions -------------------------*/
/*---------------------------------- protos --------------------------------*/
extern int	astrclip_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
				double nsigma),
		compute_jacobian(samplestruct *sample, double *dprojdred);

extern void	astrsolve_fgroups(fgroupstruct **fgroups, int nfgroup),
		astrstats_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
				double hsn_thresh),
		astrweight_fgroups(fgroupstruct **fgroups, int nfgroup),
		mat_to_wcs(polystruct *poly, double *mat, setstruct *set),
		reproj_fgroup(fgroupstruct *fgroup,fieldstruct *field),
		shrink_mat(double *alpha, double *beta, int ncoefftot,
			int index, int nmiss);

#endif

