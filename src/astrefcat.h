/*
*				astrefcat.h
*
* Include file for astrefcat.c.
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
*	Last modified:		20/04/2012
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _SAMPLES_H_
#include "samples.h"
#endif

#ifndef _ASTREFCAT_H_
#define _ASTREFCAT_H_

/*----------------------------- Internal constants --------------------------*/
#define		MAX_SERVER	16
#define		MAX_BAND	16	/* Maximum number of bands */
#define		DENIS3_POSERR	(0.20*ARCSEC/DEG)
#define		SDSSR3_POSERR	(0.10*ARCSEC/DEG)
#define		USNOA1_POSERR	(0.25*ARCSEC/DEG)
#define		USNOA2_POSERR	(0.25*ARCSEC/DEG)
#define		USNOB1_POSERR	(0.25*ARCSEC/DEG)	/* if not given */
#define		USNOA1_BMAGERR	0.40
#define		USNOA2_BMAGERR	0.40
#define		USNOB1_BMAGERR	0.40
#define         NOMAD1_MAGERR   0.30
#define         GSC_MAGERR	0.20

#define		ASTREF_ASSOCRADIUS	(0.2*ARCSEC/DEG)

/*--------------------------------- typedefs --------------------------------*/
typedef enum {ASTREFCAT_NONE, ASTREFCAT_FILE,
		ASTREFCAT_USNOA1, ASTREFCAT_USNOA2, ASTREFCAT_USNOB1,
		ASTREFCAT_GSC1, ASTREFCAT_GSC22, ASTREFCAT_GSC23,
		ASTREFCAT_2MASS, ASTREFCAT_DENIS3,
		ASTREFCAT_UCAC1, ASTREFCAT_UCAC2, ASTREFCAT_UCAC3,
		ASTREFCAT_SDSSR3, ASTREFCAT_SDSSR5, ASTREFCAT_SDSSR6,
		ASTREFCAT_SDSSR7, ASTREFCAT_SDSSR8,
		ASTREFCAT_NOMAD1, ASTREFCAT_PPMX}
			astrefenum;

typedef struct
  {
  char		name[16];		/* Catalog name */
  int		nband;			/* Number of available bands */
  int		defband;		/* Default band */
  char		bandnames[MAX_BAND][32];/* Real names of available bands */
  char		cdsbandnames[MAX_BAND][32];/* CDS names of available bands */
  int		band;			/* Chosen band */
  char		*bandname;		/* Name of chosen band */
  }	astrefstruct;

extern astrefstruct   astrefcat[];

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

#endif
