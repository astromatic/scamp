/*
*				catout.h
*
* Include file for catout.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		23/08/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

#ifndef _CATOUT_H_
#define _CATOUT_H_
/*--------------------------------- constants -------------------------------*/
/*--------------------------------- typedefs --------------------------------*/

typedef enum {CAT_NONE, CAT_ASCII_HEAD, CAT_ASCII, CAT_ASCII_SKYCAT,
		CAT_ASCII_VOTABLE, CAT_FITS_LDAC} cattypenum;

/*--------------------------- structure definitions -------------------------*/
typedef struct mergedsample
  {
  int		sourceindex;		/* Object index */
  double	wcspos[NAXIS];		/* Mean World Coordinate positions */
  int		npos_tot;		/* Total number of available positions*/
  int		npos_ok;		/* Number of available positions OK */
  float		wcsposerr[NAXIS];	/* Errors on mean WCS positions */
  float		wcspostheta;		/* WCS error position angle */
  float		wcsposdisp[NAXIS];	/* Dispersion on mean WCS positions */
  float		wcsprop[NAXIS];		/* Proper motion vectors in WCS */
  float		wcsproperr[NAXIS];	/* Proper motion vector errors in WCS */
  float		wcsparal;		/* Parallax in mas */
  float		wcsparalerr;		/* Parallax error mas */
  float		wcschi2;		/* P. motion/parallax fit Chi2/d.o.f.*/
  double	epochmin;		/* Min epoch for observations */
  double	epoch;			/* Mean epoch for observations */
  double	epochmax;		/* Max epoch for observations */
  float		flux[MAXPHOTINSTRU];	/* Mean flux */
  float		fluxerr[MAXPHOTINSTRU];	/* Mean flux uncertainty (1-sigma) */
  float		mag[MAXPHOTINSTRU];	/* "Mean" magnitude */
  int		nmag[MAXPHOTINSTRU];	/* Number of available measurements */
  float		magerr[MAXPHOTINSTRU];	/* Mean mag. uncertainty (1-sigma) */
  float		magdisp[MAXPHOTINSTRU];	/* Mean mag. dispersion (1-sigma) */
  int		nband;			/* Number of available bands */
  float		colour;			/* Colour index */
  float		spread;			/* SPREAD_MODEL weighted average*/
  float		spreaderr;		/* SPREAD_MODEL uncertainty */
  short		sexflags;		/* Merged SExtractor flags */
  short		scampflags;		/* Merged SCAMP flags */
  }	mergedsamplestruct;

typedef struct fullsample
  {
  int		sourceindex;		/* Source index */
  int		fieldindex;		/* Field index */
  short		setindex;		/* Set index */
  short		astrinstruindex;	/* Astrometric instrument index */
  short		photinstruindex;	/* Photometric instrument index */
  double	rawpos[NAXIS];		/* Mean World Coordinate positions */
  float		rawposerr[NAXIS];	/* Errors on mean pixel positions */
  float		rawpostheta;		/* Pixel error position angle */
  double	wcspos[NAXIS];		/* World Coordinate positions */
  float		wcsposerr[NAXIS];	/* Errors on WCS positions */
  float		wcspostheta;		/* WCS error position angle */
  double	epoch;			/* Epoch for observations */
  float		mag;			/* Magnitude */
  float		magerr;			/* Mag. uncertainty (1-sigma) */
  float		spread;			/* SExtractor's SPREAD_MODEL */
  float		spreaderr;		/* SPREAD_MODEL uncertainty */
  short		sexflags;		/* Merged SExtractor flags */
  short		scampflags;		/* Merged SCAMP flags */
  }	fullsamplestruct;

/*-------------------------------- protos -----------------------------------*/

void		writefullcat_fgroup(char *filename, fgroupstruct *fgroup),
		writemergedcat_fgroup(char *filename, fgroupstruct *fgroup),
		write_vo_fields(FILE *file, tabstruct *objtab);

#endif
