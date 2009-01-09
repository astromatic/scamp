/*
                                  check.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        SCAMP
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       "check-image" handling routines
*
*       Last modify:    26/07/2002
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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

