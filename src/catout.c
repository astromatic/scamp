/*
*				catout.c
*
* Produce and write merged catalogs.
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
*	Last modified:		09/02/2011
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
#include "cathead.h"
#include "catout.h"
#include "fgroup.h"
#include "field.h"
#include "prefs.h"
#include "samples.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include "xml.h"


/****** writemergedcat_fgroup *************************************************
PROTO	void writemergedcat_fgroup(char *filename, fgroupstruct *fgroup)
PURPOSE	Save a SExtractor-like catalog containing merged detections calibrated
	by SCAMP.
INPUT	File name,
	pointer to the fgroup structure.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 09/02/2011
*/
void	writemergedcat_fgroup(char *filename, fgroupstruct *fgroup)

  {
   static char  imtabtemplate[][80] = {
"SIMPLE  =                    T / This is a FITS file",
"BITPIX  =                    8 / ",
"NAXIS   =                    2 / 2D data",
"NAXIS1  =                    1 / Number of rows",
"NAXIS2  =                    1 / Number of columns",
"EXTEND  =                    T / This file may contain FITS extensions",
"END                            "};
   catstruct		*cat;
   tabstruct		*asctab, *imtab, *objtab;
   keystruct		*key, *objkeys;
   mergedsamplestruct	msample;
   fieldstruct		*field;
   setstruct		*set;
   samplestruct		*samp,*samp2;
   FILE			*ascfile;
   double		mag[MAXPHOTINSTRU], magerr[MAXPHOTINSTRU],
			magdisp[MAXPHOTINSTRU], magchi2[MAXPHOTINSTRU],
			magref[MAXPHOTINSTRU],
			wcspos[NAXIS], wcsposerr[NAXIS], wcsposdisp[NAXIS],
			wcsposref[NAXIS], wcsprop[NAXIS],wcsproperr[NAXIS],
			wcsparal,wcsparalerr,
			epoch,epochmin,epochmax, err2, colour, dummy;
   char			str[80],
			*buf, *rfilename;
   long			dptr;
   int			nmag[MAXPHOTINSTRU],
			d,f,i,k,n,p,s,nm, npinstru, naxis, N;

  if (prefs.mergedcat_type == CAT_NONE)
    return;

  naxis = fgroup->naxis;
  refmergedsample.nband = npinstru = prefs.nphotinstrustr;

/* LDAC Object header */
  objtab = new_tab("LDAC_OBJECTS");
/* Set key pointers */
  QCALLOC(objkeys, keystruct, (sizeof(refmergedkey) / sizeof(keystruct)));
  dptr = (long)((char *)&msample - (char *)&refmergedsample);
  for (k=0; refmergedkey[k].name[0]; k++)
    {
    objkeys[k] = refmergedkey[k];
    key = objkeys+k;
/*-- A trick to access the fields of the dynamic mergedsample structure */
    key->ptr = (void *)((char *)key->ptr + dptr);
    key->nbytes = t_size[key->ttype]*(key->naxis? *key->naxisn : 1);
    add_key(key,objtab, 0);
    }
/* Create a new output catalog */
  if (prefs.mergedcat_type == CAT_ASCII_HEAD
	|| prefs.mergedcat_type == CAT_ASCII
	|| prefs.mergedcat_type == CAT_ASCII_SKYCAT
	|| prefs.mergedcat_type == CAT_ASCII_VOTABLE)
    {
    cat = NULL;
    if (prefs.mergedcatpipe_flag)
      ascfile = stdout;
    else
      if (!(ascfile = fopen(filename, "w+")))
        error(EXIT_FAILURE,"*Error*: cannot open ", filename);
    if (prefs.mergedcat_type == CAT_ASCII_HEAD && (key = objtab->key))
      for (i=0,n=1; i++<objtab->nkey; key=key->nextkey)
        {
        if (*key->unit)
          fprintf(ascfile, "# %3d %-22.22s %-58.58s [%s]\n",
                n, key->name,key->comment, key->unit);
        else
          fprintf(ascfile, "# %3d %-22.22s %-58.58s\n",
                n, key->name,key->comment);
        n += key->naxis? *key->naxisn : 1;
        }
    else if (prefs.mergedcat_type == CAT_ASCII_SKYCAT && (key = objtab->key))
      {
      if (objtab->nkey<3)
        error(EXIT_FAILURE,"The SkyCat format requires at least 4 parameters:",
              " Id Ra Dec Mag");
/*--- We add a tab between rows, as required by Skycat */
      fprintf(ascfile, skycathead, 8.0);
      for (i=1,key=key->nextkey; i++<objtab->nkey; key=key->nextkey)
        {
        if (i>4)
          fprintf(ascfile, "\t%s", key->name);
        sprintf(str, "\t%s", key->printf);
        strcpy(key->printf, str);
        }
      fprintf(ascfile, "\n------------------\n");
      }
    else if (prefs.mergedcat_type == CAT_ASCII_VOTABLE && objtab->key) 
      {
/*---- A short, "relative" version of the filename */
      if (!(rfilename = strrchr(filename, '/')))
        rfilename = filename;
      else
        rfilename++;
      write_xml_header(ascfile);
      fprintf(ascfile,
	" <TABLE ID=\"Merged_List\" name=\"%s/out\">\n", rfilename);
      fprintf(ascfile,
        "  <DESCRIPTION>Table of detections merged by %s</DESCRIPTION>\n",
	BANNER);
      fprintf(ascfile,
        "  <!-- Now comes the definition of each %s parameter -->\n", BANNER);
      write_vo_fields(ascfile, objtab);
      fprintf(ascfile, "   <DATA><TABLEDATA>\n");
      }
    }
  else
    {
    cat = new_cat(1);
    init_cat(cat);
    strcpy(cat->filename, filename);
    if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
      error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);

/*-- Primary header */
    save_tab(cat, cat->tab);

/*-- We create a dummy table (only used through its header) */
    QCALLOC(asctab, tabstruct, 1);
    asctab->headnblock = 1 + (sizeof(imtabtemplate)-1)/FBSIZE;
    QCALLOC(asctab->headbuf, char, asctab->headnblock*FBSIZE);
    memcpy(asctab->headbuf, imtabtemplate, sizeof(imtabtemplate));
    for (buf = asctab->headbuf, i=FBSIZE*asctab->headnblock; i--; buf++)
      if (!*buf)
        *buf = ' ';
    write_wcs(asctab, fgroup->wcs);
/*-- (dummy) LDAC Image header */

    imtab = new_tab("LDAC_IMHEAD");
    key = new_key("Field Header Card");
    key->ptr = asctab->headbuf;
    asctab->headbuf = NULL;
    free_tab(asctab);
    key->naxis = 2;
    QMALLOC(key->naxisn, int, key->naxis);
    key->naxisn[0] = 80;
    key->naxisn[1] = 36;
    key->htype = H_STRING;
    key->ttype = T_STRING;
    key->nobj = 1;
    key->nbytes = 80*(fitsfind(key->ptr, "END     ")+1);
    add_key(key, imtab, 0);
    save_tab(cat, imtab);
    free_tab(imtab);
    objtab->cat = cat;
    init_writeobj(cat, objtab, &buf);
    }
N=0;
  for (f=0; f<fgroup->nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      samp = set->sample;
      for (n=set->nsample; n--; samp++)
        if (!samp->nextsamp && samp->prevsamp)
          {
          memset(&msample, 0, sizeof(msample));
/*-------- Photometry */
          for (p=0; p<npinstru; p++)
            {
            mag[p] = magerr[p] = magdisp[p] = magchi2[p] = magref[p] = 0.0;
            nmag[p] = 0;
            }
          nm = 0;
          colour = 0.0;
          for (samp2 = samp;
		samp2 && (p=samp2->set->field->photomlabel)>=0;
                samp2=samp2->prevsamp)
            {
            colour += samp2->colour;
            if (samp2->flux > 0.0 && (err2 = samp2->magerr*samp2->magerr)>0.0)
              {
              magerr[p] += 1.0 / err2;
              mag[p] += samp2->mag / err2;
              if (!nmag[p])
                magref[p] = samp2->mag;
              magdisp[p] += (samp2->mag-magref[p]) * (samp2->mag - magref[p]);
              nmag[p]++;
              }
            nm++;
            }
          msample.colour = nm? colour/nm : 0.0;
          for (p=0; p<npinstru; p++)
            {
            if ((nm=nmag[p]))
              {
              msample.mag[p] = mag[p] / magerr[p];
              msample.magerr[p] = sqrt(1.0 / magerr[p]);
              msample.magdisp[p] = nm > 1?
		sqrt(fabs(magdisp[p]
		 - nm*(msample.mag[p] - magref[p])
		 *(msample.mag[p] - magref[p]))/(nm-1.0))
		: 0.0;
              }
            else
              msample.mag[p] = msample.magerr[p] = msample.magdisp[p] = 99.0;
            msample.nmag[p] = nmag[p];
            }
/*-------- Astrometry */
          nm = 0;
          for (d=0; d<naxis; d++)
            wcspos[d] = wcsposerr[d] = wcsposdisp[d] = wcsposref[d]
		= wcsprop[d] = wcsproperr[d] = 0.0;
          wcsparal = wcsparalerr = 0.0;
          epoch = 0.0;
          epochmin = BIG;
          epochmax = -BIG;
          for (samp2 = samp;
		samp2 && (p=samp2->set->field->photomlabel)>=0;
                samp2=samp2->prevsamp)
            {
            for (d=0; d<naxis; d++)
              {
              err2 = samp2->wcsposerr[d]*samp2->wcsposerr[d];
              wcspos[d] += samp2->wcspos[d];
              if (err2 <= 0.0)
                err2 += 1*MAS / DEG;
              wcsposerr[d] += 1.0 / err2;
              wcspos[d] += samp2->wcspos[d] / err2;
              if (!nm)
                wcsposref[d] = samp2->wcspos[d];
              wcsposdisp[d] += (samp2->wcspos[d] - wcsposref[d])
				* (samp2->wcspos[d] - wcsposref[d]);
              wcsprop[d] += samp2->wcsprop[d];
              wcsproperr[d] += samp2->wcsproperr[d]*samp2->wcsproperr[d];
              }
/*
raw_to_wcs(fgroup->wcs, samp2->toto, samp2->toto);
printf("%d %.7f %.7f %.5f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %d\n", N,
samp2->wcspos[0],
samp2->wcspos[1],
samp2->set->field->epoch,
(samp2->wcspos[0]-wcsposref[0])*cos(samp2->wcspos[0]*DEG)*DEG/MAS,
(samp2->wcspos[1]-wcsposref[1])*DEG/MAS,
samp2->wcsposerr[0]*DEG/MAS,
samp2->wcsposerr[1]*DEG/MAS,
samp2->wcsprop[0]*DEG/MAS,
samp2->wcsprop[1]*DEG/MAS,
samp2->wcsproperr[0]*DEG/MAS,
samp2->wcsproperr[1]*DEG/MAS,
samp2->wcspos[0]-
-(samp2->wcspos[0]-samp2->toto[0])*cos(samp2->wcspos[0]*DEG)*DEG/MAS,
-(samp2->wcspos[1]-samp2->toto[1])*DEG/MAS,
samp2->flags);
*/
            wcsparal += samp2->wcsparal;
            wcsparalerr += samp2->wcsparalerr*samp2->wcsparalerr;
            epoch += samp2->set->field->epoch;
            if (samp2->set->field->epoch < epochmin)
              epochmin = samp2->set->field->epoch;
            if (samp2->set->field->epoch > epochmax)
              epochmax = samp2->set->field->epoch;
            msample.sexflags |= samp2->sexflags;
            msample.scampflags |= samp2->scampflags;
            nm++;
            }
          if (nm)
            {
            for (d=0; d<naxis; d++)
              {
              msample.wcspos[d] = wcspos[d] / wcsposerr[d];
              msample.wcsposerr[d] = sqrt(1.0/wcsposerr[d]);
              msample.wcsposdisp[d] = nm > 1?
		sqrt(fabs(wcsposdisp[d]
		 - nm*(msample.wcspos[d] - wcsposref[d])
			*(msample.wcspos[d] - wcsposref[d]))/(nm-1.0))
		: 0.0;
              msample.wcsprop[d] = (wcsprop[d]/nm) * DEG/MAS;
              msample.wcsproperr[d] = sqrt(wcsproperr[d]/nm) * DEG/MAS;
              }
            if (msample.wcsposerr[0] < msample.wcsposerr[1])
              {
              dummy = msample.wcsposerr[0];
              msample.wcsposerr[0] = msample.wcsposerr[1];
              msample.wcsposerr[1] = dummy;
              msample.wcspostheta = 90.0;
              }
            else
              msample.wcspostheta = 0.0;
            msample.wcsparal = (wcsparal/nm) * DEG/MAS;
            msample.wcsparalerr = sqrt(wcsparalerr/nm) * DEG/MAS;
            msample.epoch = epoch / nm;
            msample.epochmin = epochmin;
            msample.epochmax = epochmax;
            msample.npos = nm;
            }
          if (prefs.mergedcat_type == CAT_ASCII_HEAD
		|| prefs.mergedcat_type == CAT_ASCII
		|| prefs.mergedcat_type == CAT_ASCII_SKYCAT)
            print_obj(ascfile, objtab);
          else if (prefs.mergedcat_type == CAT_ASCII_VOTABLE)
            voprint_obj(ascfile, objtab);
          else
            write_obj(objtab, buf);
N++;
          }
      }
    }

  if (prefs.mergedcat_type == CAT_ASCII_HEAD
	|| prefs.mergedcat_type == CAT_ASCII
	|| prefs.mergedcat_type == CAT_ASCII_SKYCAT)
    {
    if (prefs.mergedcat_type == CAT_ASCII_SKYCAT)
      fprintf(ascfile, skycattail);
    if (!prefs.mergedcatpipe_flag)
      fclose(ascfile);
    }
  else if (prefs.mergedcat_type == CAT_ASCII_VOTABLE)
    {
    fprintf(ascfile, "    </TABLEDATA></DATA>\n");
    fprintf(ascfile, "  </TABLE>\n");
/*-- Add configuration file meta-data */
    write_xml_meta(ascfile, NULL);
    fprintf(ascfile, "</RESOURCE>\n");
    fprintf(ascfile, "</VOTABLE>\n");
    }
  else
    end_writeobj(cat, objtab, buf);

  objtab->key = NULL;
  objtab->nkey = 0;
  free_tab(objtab);
  free(objkeys);

  if (cat)
    free_cat(&cat, 1);

  return;
  }

/****** write_vo_fields *******************************************************
PROTO	void	write_vo_fields(FILE *file, tabstruct *objtab)
PURPOSE	Write the list of columns to an XML-VOTable file or stream
INPUT	Pointer to the output file (or stream),
	Pointer to the object table.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	19/10/2009
 ***/
void	write_vo_fields(FILE *file, tabstruct *objtab)
  {
   keystruct	*key;
   char		datatype[40], arraysize[40], str[40];
   int		i, d;

  if (!objtab || !objtab->key)
    return;
  key=objtab->key;
  for (i=0; i++<objtab->nkey; key=key->nextkey)
    {
/*--- indicate datatype, arraysize, width and precision attributes */
/*--- Handle multidimensional arrays */
    arraysize[0] = '\0';
    if (key->naxis>1)
      {
      for (d=0; d<key->naxis; d++)
        {
        sprintf(str, "%s%d", d?"x":" arraysize=\"", key->naxisn[d]);
        strcat(arraysize, str);
        }
      strcat(arraysize, "\"");
      }
    switch(key->ttype)
      {
      case T_BYTE:	strcpy(datatype, "unsignedByte"); break;
      case T_SHORT:	strcpy(datatype, "short"); break;
      case T_LONG:	strcpy(datatype, "int"); break;
      case T_FLOAT:	strcpy(datatype, "float"); break;
      case T_DOUBLE:	strcpy(datatype, "double"); break;
      default:		error(EXIT_FAILURE,
			"*Internal Error*: Unknown datatype in ",
			"initcat()");
      }
    fprintf(file,
	"  <FIELD name=\"%s\" ucd=\"%s\" datatype=\"%s\" unit=\"%s\"%s>\n",
	key->name, key->voucd, datatype,key->vounit, arraysize);
    fprintf(file, "   <DESCRIPTION>%s</DESCRIPTION>\n", key->comment);
    fprintf(file, "  </FIELD>\n");
    }

  return;
  }


