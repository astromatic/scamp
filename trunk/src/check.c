/*
*				check.c
*
* Create and manage "check images".
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2010 IAP/CNRS/UPMC
*
*	Author:			Emmanuel Bertin (IAP)
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
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "define.h"
#include "globals.h"
#include "check.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "prefs.h"

/****** write_check ***********************************************************
PROTO	void write_check(char *filename, float *pix, int width, int height);
PURPOSE	Write a check-image.
INPUT	Check-image filename,
	ptr to the image array to be written,
	image width in pixels,
	image height in pixels.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2002
 ***/
void	write_check(char *filename, float *pix, int width, int height)

  {
   catstruct	*cat;
   tabstruct	*tab;
   char		str[MAXCHAR],
		*rfilename;

/* Display an abridged version of the filename */
  if (!(rfilename = strrchr(filename, '/')))
    rfilename = filename;
  else
    rfilename++;
  sprintf(str,"Writing check-image %s", rfilename);
  NFPRINTF(OUTPUT, str);


/* Create the new cat (well it is not a "cat", but simply a FITS table) */
  cat = new_cat(1);
  init_cat(cat);
  tab = cat->tab;
  tab->naxis = 2;       /* This is an image */
  QMALLOC(tab->naxisn, int, tab->naxis);
  tab->bitpix =  BP_FLOAT;
  tab->bytepix = t_size[T_FLOAT];
  tab->naxisn[0] = width;
  tab->naxisn[1] = height;
  tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1];
  tab->bodybuf = (void *)pix;
  save_cat(cat, filename);
  tab->bodybuf = NULL;
  free_cat(&cat, 1);

  return;
  }


/****** check_check ***********************************************************
PROTO	int check_check(checkenum checktype)
PURPOSE	Check that the specified check-image type has been requested by user,
	by returning its index in CHECKIMAGE_TYPE keyword list, or RETURN_ERROR
	(-1) otherwise.
INPUT	Check-image type.
OUTPUT	Index in CHECKIMAGE_TYPE keyword list, or RETURN_ERROR (-1) otherwise.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2002
 ***/
int	 check_check(checkenum checktype)

  {
   int		i;

  for (i=0; i<prefs.ncheck_type; i++)
    if (checktype == prefs.check_type[i])
      return i;

  return RETURN_ERROR;
  }

