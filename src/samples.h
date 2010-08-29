/*
 				samples.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Type definitions related to samples
*
*	Last modify:	29/08/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _SAMPLES_H_
#define _SAMPLES_H_

/*--------------------------------- constants -------------------------------*/

#define	LSAMPLE_DEFSIZE	1000		/* Sample stacksize at the beginning */
#define MAXCONTEXT	10		/* Maximum number of contexts */
#define MAX_FLUX	1e26		/* Maximum flux allowed */
#define MIN_FLUXERR	1e-15		/* Minimum flux error estimate allowed*/
#define MIN_POSERR	1e-6		/* Minimum pos. error estimate allowed*/

/*-------------------------------- flags ------------------------------------*/

#define		OBJ_CROWDED     0x0001
#define		OBJ_MERGED      0x0002
#define		OBJ_SATUR       0x0004
#define		OBJ_TRUNC       0x0008

/*--------------------------------- typedefs --------------------------------*/

typedef enum {UNION_RAW, UNION_PROJ, UNION_WCS}	unionmodenum;

/*--------------------------- structure definitions -------------------------*/

typedef struct sample
  {
  struct sample	*prevsamp;		/* Link to previous sample */
  struct sample	*nextsamp;		/* Link to next sample */
  struct set	*set;			/* Link to parent set */
  double	*context;		/* Context vector */
  double	projpos[NAXIS];		/* Projected coordinates */
  double	rawpos[NAXIS];		/* Raw coordinates */
  double	wcspos[NAXIS];		/* World Coordinate positions */
  float		rawposerr[NAXIS];	/* Uncertainty on raw coordinates */
  float		wcsposerr[NAXIS];	/* Errors on WCS positions */
  float		wcsprop[NAXIS];		/* Proper motion vectors in the WCS */
  float		wcsproperr[NAXIS];	/* P. motion vector errors in the WCS */
  float		wcsparal;		/* Trigonometric parallax */
  float		wcsparalerr;		/* Trigonometric parallax error */
  float		weight;			/* Relative weight in solutions */
  float		flux;			/* Flux */
  float		fluxerr;		/* Flux uncertainty (1-sigma) */
  float		mag;			/* Magnitude */
  float		magerr;			/* Magnitude uncertainty (1-sigma) */
  float		colour;			/* A colour index */
  float		fwhm;			/* Full Width at Half Maximum */
  int		flags;			/* Source flags */
  }	samplestruct;

typedef struct set
  {
  struct sample	*sample;		/* Array of samples */
  int		nsample;		/* Number of samples in stack */
  int		nsamplemax;		/* Max number of samples in stack */
  tabstruct	*imatab;		/* Image-related table structure */
  wcsstruct	*wcs;			/* WCS information */
/* ---- matching parameters */
  float		match_dscale;		/* Computed scale correction */
  float		match_dangle;		/* Computed angle correction */
  float		match_shear;		/* Computed shear amplitude */
  float		match_sangle;		/* Computed direction of the shear */
  float		match_dlng, match_dlat;	/* Positional corrections */
  float		match_asig, match_sig;	/* Contrasts for angle/scale and pos */
/* ---- astrometric parameters */
  double	wcspos[NAXIS];		/* Central pixel coordinate */
  double	wcsscale[NAXIS];	/* Central pixel scale */
  double	radius;			/* Approximate radius of set (deg)*/
  double	projposmin[NAXIS];	/* Minimum projected position in set */
  double	projposmax[NAXIS];	/* Maximum projected position in set */
  int		naxis;			/* Number of axes */
  int		lng,lat;		/* Longitude and latitude indices */
  int		ncontext;		/* Number of contexts */
  char		**contextname;		/* List of context keywords used */
  double	*contextoffset;		/* Offset to apply to context data */
  double	*contextscale;		/* Scaling to apply to context data */
  int		contextx;		/* Context associated to x (-1=none)*/
  int		contexty;		/* Context associated to y (-1=none)*/
  int		index;			/* Set index for CONTEXTs */
  int		index2;			/* Set index for field-dependent prm */
  double	weightfac;		/* Weight factor for astrom. refs */
/* ---- photometric parameters */
  double	airmass;		/* Air mass */
  double	expotime;		/* Exposure time */
  double	extcoeff;		/* Extinction coefficient */
  double	magzero;		/* Magnitude zero-point */
  double	dmagzero;		/* Computed magnitude z.-p. offset*/
  double	fluxscale;		/* Relative flux scale */
  int		nsaturated;		/* Number of saturated detections */
  int		ncoeff;			/* Number of fit coefficients kept */
  int		nconst;			/* Number of constraints */
  struct field	*field;			/* Link to parent field */
  }	setstruct;

/*-------------------------------- protos -----------------------------------*/

samplestruct	*remove_sample(setstruct *set, int isample);

setstruct	*init_set(void),
		*load_samples(char **filename, int ncat),
		*read_samples(setstruct *set, tabstruct *tab, char *rfilename);

void		copy_samples(samplestruct *samplein, setstruct *set,
			int nsample),
		end_set(setstruct *set),
		free_samples(setstruct *set),
		locate_set(setstruct *set),
 		malloc_samples(setstruct *set, int nsample),
		make_weights(setstruct *set, samplestruct *sample),
		mix_samples(setstruct *set),
		realloc_samples(setstruct *set, int nsample),
		sort_samples(setstruct *set),
		union_samples(samplestruct *samplein, setstruct *set,
			int nsamplein, double radius, unionmodenum mode),
		unlink_samples(setstruct *set),
		update_retina(setstruct *set, samplestruct *sample,
			float pixstep);


#endif
