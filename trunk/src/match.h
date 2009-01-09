 /*
 				match.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	include for match.c.
*
*	Last modify:	19/02/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _FGROUPS_H_
#include "fgroup.h"
#endif

#ifndef _MATCH_H_
#define _MATCH_H_

/*---------------------------- Internal constants ---------------------------*/
#define		MATCH_CONFPROB	0.01	/* Confusion probability in auto-resol*/
#define		MATCH_MINCONT	2.0	/* Minimum constrast factor */
#define		AS_NSOURCESTART	2000	/* Start source nb for auto AS match */
#define		AS_NSOURCEMAX	8000	/* Max source nb for auto AS match */
#define		SCALE_RANGE	16.0	/* Initial scale search interval */
#define		SCALE_ZERO	1.0	/* Initial scale zero-point */
#define		SCALE_MAXSIZE	4096	/* Max resol. elements in scale */
#define		ANGLE_MAXSIZE	4096	/* Max resol. elements in pos angle */
#define		ASLOPASS_LAMBDA	2.0 	/* Lo-pass filter cutoff (pixels) */
#define		ASHIPASS_LAMBDA	(1/64.0)/* Hi-pass filter cutoff (frame frac.)*/
#define		AS_FLUXEXP	0.125	/* Exponent used for flux in LL match*/
#define		LL_NSOURCEMAX	200000	/* Max source nb/set for LL match */
#define		LL_WRAPPING	1	/* Wrapping factor along each dim. */
#define		LL_MAXSIZE	8192	/* Max resolution elements in LL */
#define		LLLOPASS_LAMBDA	2.0	/* Lo-pass filter cutoff (pixels) */
#define		LLHIPASS_LAMBDA	(1/32.0)/* Hi-pass filter cutoff (frame frac.)*/
#define		LL_FLUXEXP	0.125	/* Exponent used for flux in LL match*/
#define		GAUSS_MAXSIG	4.0	/* Max. Gaussian sigma (for speedup) */
#define		GAUSS_MAXNSIG	3.0	/* Max. Gauss mask size (for speedup) */
#define		PEAKFIND_NITER	1	/* Max. number of peakfind iterations*/

/*------------------------------- functions ---------------------------------*/
extern setstruct	*frame_set(setstruct *setin, wcsstruct *wcs,
				double *centerpos, double radius),
			*new_fieldset(fieldstruct *field);

extern double		findcrosspeak(float *histo,
				int width, int height,
				double xrange, double yrange,
				double xresol, double yresol,
				double *xpeak, double *ypeak),
			match_setas(setstruct *set, setstruct *refset, int nmax,
				double matchresol, double *angle,double *scale),
			match_setll(setstruct *set, setstruct *refset,
				double matchresol, double *dlng, double *dlat),
			mean_rawposvar(setstruct *set);

extern void		compute_rawpos(wcsstruct *wcs, samplestruct *refsample,
				int nrefsample),
			compute_wcsss(wcsstruct *wcs,
				double *sangle, double *shear),
			flipimage(float *histo, int width, int height),
			match_field(fieldstruct *field, fieldstruct *reffield),
			match_refine(setstruct *set, setstruct *refset,
				double matchresol, double *angle, double *scale,
				double *sangle, double *ratio,
				double *dlng,double *dlat),
			match_test(setstruct *set, setstruct *refset,
				double matchresol, double *angle, double *scale,
				double *sangle, double *ratio,
				double *dlng,double *dlat),
			print_matchinfo(fieldstruct *field),
			update_wcsas(wcsstruct *wcs, double angle,
				double scale),
			update_wcscc(wcsstruct *wcs, double drawlng,
				double drawlat),
			update_wcsll(wcsstruct *wcs, double dlng, double dlat),
			update_wcsss(wcsstruct *wcs, double sangle,
				double shear);

#ifdef USE_THREADS
void			*pthread_match_field(void *arg),
			pthread_match_fields(fgroupstruct **fgroups,
				fieldstruct **reffields, int ngroup);
#endif
#endif
