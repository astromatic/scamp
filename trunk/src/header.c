/*
*				header.c
*
* Manage FITS or ASCII headers.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		16/11/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "field.h"
#include "fitswcs.h"
#include "fgroup.h"
#include "prefs.h"
#include "samples.h"

/*-------------------------- Worthy WCS keywords ----------------------------*/
char	wcskey[][12] = {"EQUINOX", "RADECSYS", "CTYPE???", "CUNIT???",
		"CRVAL???", "CRPIX???", "CDELT???", "CD?_?", "PV?_????",
		"FGROUPNO", "ASTIRMS?", "ASTRRMS?", "ASTINST ",
		"FLXSCALE", "MAGZEROP", "PHOTIRMS", "PHOTINST", "PHOTLINK",
		""};


/****** read_aschead ********************************************************
PROTO	int read_aschead(char *filename, int frameno, tabstruct *tab)
PURPOSE	Read an ASCII header file and update the current field's tab
INPUT	Name of the ASCII file,
	Frame number (if extensions),
	Tab structure.
OUTPUT	RETURN_OK if the file was found, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/09/2004
 ***/
int	read_aschead(char *filename, int frameno, tabstruct *tab)
  {
   char		keyword[88],data[88],comment[88], str[88];
   FILE         *file;
   h_type       htype;
   t_type       ttype;
   int          i;

  if ((file=fopen(filename, "r")))
    {
/*- Skip previous ENDs in multi-FITS extension headers */
    for (i=frameno; i--;)
      while (fgets(str, MAXCHAR, file)
                && strncmp(str,"END ",4)
                && strncmp(str,"END\n",4));
    memset(str, ' ', 80);
    while (fgets(str, 81, file) && strncmp(str,"END ",4)
                                && strncmp(str,"END\n",4))
      {
      fitspick(str, keyword, data, &htype, &ttype, comment);
/*---- Block critical keywords */
      if (!wstrncmp(keyword, "SIMPLE  ", 8)
        ||!wstrncmp(keyword, "BITPIX  ", 8)
        ||!wstrncmp(keyword, "NAXIS   ", 8)
        ||!wstrncmp(keyword, "BSCALE  ", 8)
        ||!wstrncmp(keyword, "BZERO   ", 8))
        continue;
      addkeywordto_head(tab, keyword, comment);
      fitswrite(tab->headbuf, keyword, data, htype, ttype);
      memset(str, ' ', 80);
      }
    fclose(file);
/*-- Update the tab data */
    readbasic_head(tab);
    return RETURN_OK;
    }
  else
    return RETURN_ERROR;
  }


/****** write_aschead ********************************************************
PROTO	int write_aschead(char *filename, fieldstruct *field)
PURPOSE	Write an ASCII header file from a series of tabs
INPUT	Name of the ASCII file,
	Field structure.
OUTPUT	RETURN_OK if the file was found, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/11/2010
 ***/
int	write_aschead(char *filename, fieldstruct *field)
  {
   fgroupstruct	*fgroup;
   tabstruct	*tab;
   setstruct	**set;
   FILE         *file;
   char		str[88],
		*ptr;
   double	val;
   int          d,i,n,s, instru, naxis;

  if (!(file=fopen(filename, "w")))
    return RETURN_ERROR;

  set = field->set;
  fgroup = field->fgroup;
  if (!fgroup)
    error(EXIT_FAILURE, "*Internal Error*: no field group found for ",
	field->filename);
  naxis = fgroup->naxis;

  for (s=0; s<field->nset; s++)
    {
    tab = new_tab("");
    update_head(tab);
    write_wcs(tab, set[s]->wcs);
    switch(prefs.header_type)
      {
      case NORMAL:
/*------ Update header with the latest astrometric and photometric info */
        fitswrite(tab->headbuf,"FGROUPNO", &fgroup->no, H_INT, T_LONG);
        addkeywordto_head(tab, "FGROUPNO","SCAMP field group label");
        fitswrite(tab->headbuf,"FGROUPNO", &fgroup->no, H_INT, T_LONG);
        for (d=0; d<naxis; d++)
          {
          sprintf(str, "ASTIRMS%1d", d+1);
          addkeywordto_head(tab, str,
		"Astrom. dispersion RMS (intern., high S/N)");
          fitswrite(tab->headbuf, str,&fgroup->sig_interr_hsn[d],
		H_EXPO, T_DOUBLE);
          sprintf(str, "ASTRRMS%1d", d+1);
          val = (double)field->sig_referr_hsn[d];
          addkeywordto_head(tab, str,
		"Astrom. dispersion RMS (ref., high S/N)");
          fitswrite(tab->headbuf, str, &val, H_EXPO, T_DOUBLE);
          }
        instru = field->astromlabel+1;
        addkeywordto_head(tab,"ASTINST ","SCAMP astrometric instrument label");
        fitswrite(tab->headbuf,"ASTINST ", &instru, H_INT, T_LONG);
        addkeywordto_head(tab, "FLXSCALE", "SCAMP relative flux scale");
        fitswrite(tab->headbuf, "FLXSCALE", &set[s]->fluxscale,
		H_EXPO, T_DOUBLE);
        addkeywordto_head(tab, "MAGZEROP", "SCAMP zero-point");
        fitswrite(tab->headbuf, "MAGZEROP",
		&prefs.magzero_out[field->photomlabel], H_FLOAT, T_DOUBLE);
        addkeywordto_head(tab, "PHOTIRMS",
		"mag dispersion RMS (internal, high S/N)");
        fitswrite(tab->headbuf,"PHOTIRMS",
		&fgroup->sig_intmagerr_hsn[field->photomlabel],
		H_FLOAT,T_DOUBLE);
        addkeywordto_head(tab, "PHOTRRMS",
		"mag dispersion RMS (ref., high S/N)");
        fitswrite(tab->headbuf,"PHOTRRMS",
		&fgroup->sig_refmagerr_hsn[field->photomlabel],
		H_FLOAT,T_DOUBLE);
        instru = field->photomlabel+1;
        addkeywordto_head(tab,"PHOTINST","SCAMP photometric instrument label");
        fitswrite(tab->headbuf,"PHOTINST", &instru, H_INT, T_LONG);
        addkeywordto_head(tab,"PHOTLINK",
		"True if linked to a photometric field");
        fitswrite(tab->headbuf,"PHOTLINK", &field->photomlink, H_BOOL, T_LONG);
/*------ We are only interested in astrometric and photometric informations */
        sprintf(str, "HISTORY   Astrometric solution by %s version %s (%s)",
		BANNER,MYVERSION,DATE);
        fprintf(file, "%.79s\n", str); 
        sprintf(str, "COMMENT   (c) %s", COPYRIGHT);
        fprintf(file, "%.79s\n", str); 
        fprintf(file, "COMMENT   \n");
        for (i=0; *wcskey[i]; i++)
          for (ptr=tab->headbuf; (n=fitsfind(ptr, wcskey[i])) != RETURN_ERROR;
		ptr+=80)
            fprintf(file, "%.79s\n", ptr += n*80);          
        break;
      case FOCAL_PLANE:
        for (ptr=tab->headbuf; (n=fitsfind(ptr, "CRPIX???")) != RETURN_ERROR;
		ptr+=80)
          fprintf(file, "%.79s\n", ptr += n*80);          
        for (ptr=tab->headbuf; (n=fitsfind(ptr, "CD?_?")) != RETURN_ERROR;
		ptr+=80)
          fprintf(file, "%.79s\n", ptr += n*80);          
        break;
      default:
        error(EXIT_FAILURE, "*Internal Error*: Unknown option in ",
	    "write_aschead()");
      }

    fprintf(file, "END     \n");
    free_tab(tab);
    }

  fclose(file);

  return RETURN_OK;
  }

