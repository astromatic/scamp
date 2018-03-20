/*
*				astrefcat.h
*
* Include file for astrefcat.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2018 IAP/CNRS/UPMC
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
*	Last modified:		20/03/2018
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _ASTREFCAT_H_
#define _ASTREFCAT_H_

#include "field.h"
#include "samples.h"

/*----------------------------- Internal constants --------------------------*/
#define		MAX_SERVER	16
#define		MAX_BAND	16	/* Maximum number of bands */
#define		MAX_COLUMN	32	/* Maximum number of columns */
#define		COLUMN_SIZE	32	/* Number of characters per column */
#define		DENIS3_POSERR	(0.20*ARCSEC/DEG)
#define		USNOA2_POSERR	(0.25*ARCSEC/DEG)
#define		USNOB1_POSERR	(0.25*ARCSEC/DEG)	/* if not given */
#define		CMC15_POSERR	(0.10*ARCSEC/DEG)	/* if not given */
#define		XPM_PROPERR	(10*MAS/DEG)
#define		USNOA2_BMAGERR	0.40
#define		USNOB1_BMAGERR	0.40
#define         NOMAD1_MAGERR   0.30
#define         GSC_MAGERR	0.20
#define         TWOMASS_MAGERR  0.1		/* Just a default value */
#define         UCAC_MAGERR     0.12		/* Just a default value */
#define         DEFAULT_MAGERR  0.1		/* Just a default value */

#define		ASTREF_ASSOCRADIUS	(0.2*ARCSEC/DEG)

/*--------------------------------- typedefs --------------------------------*/
typedef enum {ASTREFCAT_NONE, ASTREFCAT_FILE,
		ASTREFCAT_USNOA2, ASTREFCAT_USNOB1, ASTREFCAT_GSC23,
		ASTREFCAT_2MASS, ASTREFCAT_DENIS3, ASTREFCAT_UCAC4,
		ASTREFCAT_URAT1, ASTREFCAT_SDSSR9, ASTREFCAT_NOMAD1,
		ASTREFCAT_PPMX, ASTREFCAT_CMC15, ASTREFCAT_TYCHO2,
		ASTREFCAT_IGSL, ASTREFCAT_ALLWISE, ASTREFCAT_GAIADR1,
		ASTREFCAT_PANSTARRS1}	astrefenum;

typedef struct
  {
  char		name[16];		/* Catalog name */
  char		viziername[32];		/* Vizier catalog name */
  char		viziercolumns[MAX_COLUMN][COLUMN_SIZE];
					/* List of Vizier column names */
  char		vizierbandnames[MAX_BAND][32];
					/* Vizier names of available bands */
  char		bandnames[MAX_BAND][32];/* Real names of available bands */
  int		nband;			/* Number of available bands */
  int		defband;		/* Default band */
  int		band;			/* Chosen band */
  char		*bandname;		/* Name of chosen band */
  }	astrefstruct;

extern astrefstruct   astrefcats[];

/*------------------------------- functions ---------------------------------*/

extern fieldstruct	*get_astreffield(astrefenum refcat, double *wcspos,
                                int lng, int lat, int naxis, double maxradius),
			*load_astreffield(char *filename, double *wcspos,
				int lng, int lat,
				int naxis, double maxradius, int band,
				double *maglim);

extern setstruct	*read_astrefsamples(setstruct *set, tabstruct *tab,
				char *rfilename,
				double *wcspos,
				int lng, int lat,
				int naxis, double maxradius, int band,
				double *maglim);

extern void		save_astreffield(char *filename, fieldstruct *reffield);

#endif // _ASTREFCAT_H_
