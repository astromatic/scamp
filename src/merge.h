/*
*				merge.h
*
* Include file for merge.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		28/06/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _MERGE_H_
#define _MERGE_H_

#include "fgroup.h"

/*--------------------------------- constants -------------------------------*/
/*--------------------------------- typedefs --------------------------------*/
/*--------------------------- structure definitions -------------------------*/
typedef struct msample
  {
  int		sourceindex;		/* Object index */
  struct sample	*samp;			/* Pointer to last sample (detection) */
  double	wcspos[NAXIS];		/* Mean World Coordinate positions */
  int		npos_tot;		/* Total number of available positions*/
  int		npos_ok;		/* Number of available positions OK */
  float		wcsposerr[NAXIS];	/* Errors on mean WCS positions */
  float		wcspostheta;		/* WCS error position angle */
  float		wcsposdisp[NAXIS];	/* Dispersion on mean WCS positions */
  float		wcsprop[NAXIS];		/* Proper motion vectors in mas/yr */
  float		wcsproperr[NAXIS];	/* Proper motion vector errors in mas/yr */
  float		wcsparal;		/* Parallax in mas */
  float		wcsparalerr;		/* Parallax error mas */
  float		wcschi2;		/* P. motion/parallax fit Chi2/d.o.f.*/
  double	epochmin;		/* Min epoch for observations */
  double	epoch;			/* Mean epoch for observations */
  double	epochmax;		/* Max epoch for observations */
  float		colour;			/* Colour index */
  float		spread;			/* SPREAD_MODEL weighted average*/
  float		spreaderr;		/* SPREAD_MODEL uncertainty */
  unsigned short sexflags;		/* Merged SExtractor flags */
  unsigned short scampflags;		/* Merged SCAMP flags */
  unsigned int	imaflags;		/* Merged image flags */
  }	msamplestruct;

/*-------------------------------- protos -----------------------------------*/

msamplestruct	*merge_fgroup(fgroupstruct *fgroup, fieldstruct *reffield);
;

#endif // _MERGE_H_
