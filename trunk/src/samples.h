/*
*				samples.h
*
* Include file for samples.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		25/02/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _SAMPLES_H_
#define _SAMPLES_H_

/*--------------------------------- constants -------------------------------*/

#define	LSAMPLE_DEFSIZE	1000	/* Sample stacksize at the beginning */
#define MAXCONTEXT	10	/* Maximum number of contexts */
#define MAX_FLUX	1e26	/* Maximum flux allowed [ADU] */
#define MIN_FLUXERR	1e-15	/* Minimum flux uncertainty allowed [ADU] */
#define FLOOR_FLUXERR	1e-3	/* Floor of fractional flux uncertainties */
#define MIN_POSERR	1e-6	/* Minimum position uncertainties [pix] */
#define FLOOR_POSERR	1e-2	/* Floor of position uncertainties [pix] */

/*-------------------------------- flags ------------------------------------*/

#define		OBJ_CROWDED		0x0001
#define		OBJ_MERGED		0x0002
#define		OBJ_SATUR		0x0004
#define		OBJ_TRUNC		0x0008

#define		SCAMP_ASTRCLIPPED	0x0001
#define		SCAMP_PHOTCLIPPED	0x0010
#define		SCAMP_PHOTNOCOLOR	0x0020

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
  short		sexflags;		/* Source extraction flags */
  short		scampflags;		/* SCAMP flags */		
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
