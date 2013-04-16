/*
*				catout.c
*
* Produce and write merged and full catalogs.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		07/04/2013
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
#include "key.h"
#include "merge.h"
#include "prefs.h"
#include "samples.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include "xml.h"


/****** writemergedcat_fgroup *************************************************
PROTO	void writemergedcat_fgroup(char *filename, fgroupstruct *fgroup)
PURPOSE	Save a SExtractor-like catalog containing merged detections as
	calibrated by SCAMP.
INPUT	File name,
	pointer to the fgroup structure.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 07/04/2013
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
   msamplestruct	*msamp;
   mergedsamplestruct	mergedsample;
   fieldstruct		*field;
   setstruct		*set;
   samplestruct		*samp,*samp2;
   FILE			*ascfile;
   double		mag[MAXPHOTINSTRU], magerr[MAXPHOTINSTRU],
			magdisp[MAXPHOTINSTRU], magchi2[MAXPHOTINSTRU],
			magref[MAXPHOTINSTRU],
			wcspos[NAXIS], wcsposerr[NAXIS], wcsposdisp[NAXIS],
			wcsposref[NAXIS], wcsprop[NAXIS],wcsproperr[NAXIS],
			wcsparal,wcsparalerr, wcschi2,
			epoch,epochmin,epochmax, err2,
			weight,weights, dummy;
   char			str[80],
			*buf, *rfilename;
   long			dptr;
   short		astrflagmask,photflagmask;
   int			nmag[MAXPHOTINSTRU],
			d,f,i,k,m,n,p,s, nall,nphotok,nposok, npinstru, naxis,
			index, refflag;

  if (prefs.mergedcat_type == CAT_NONE)
    return;

  astrflagmask = (short)prefs.astr_flagsmask;
  photflagmask = (short)prefs.phot_flagsmask;
  naxis = fgroup->naxis;
  refmergedsample.nband = npinstru = prefs.nphotinstrustr;
  refflag = prefs.astrefinprop_flag;

/* LDAC Object header */
  objtab = new_tab("LDAC_OBJECTS");
/* Set key pointers */
  QCALLOC(objkeys, keystruct, (sizeof(refmergedkey) / sizeof(keystruct)));
  dptr = (long)((char *)&mergedsample - (char *)&refmergedsample);
  for (k=0; refmergedkey[k].name[0]; k++)
    {
    if (!prefs.spread_flag
	&& !cistrcmp(refmergedkey[k].name, "SPREAD", FIND_NOSTRICT))
      continue;
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
    key->naxisn[1] = fitsfind(key->ptr, "END     ")+1;
    key->htype = H_STRING;
    key->ttype = T_STRING;
    key->nobj = 1;
    key->nbytes = key->naxisn[0]*key->naxisn[1];
    add_key(key, imtab, 0);
    save_tab(cat, imtab);
    free_tab(imtab);
    objtab->cat = cat;
    init_writeobj(cat, objtab, &buf);
    }

  index=0;
  msamp = fgroup->msample;
  for (m=fgroup->nmsample; m--; msamp++)
    {
    samp = msamp->samp;
    memset(&mergedsample, 0, sizeof(mergedsample));
    mergedsample.sourceindex = msamp->sourceindex;
/*-- Photometry */
    for (p=0; p<npinstru; p++)
      {
      mag[p] = magerr[p] = magdisp[p] = magchi2[p] = magref[p] = 0.0;
      nmag[p] = 0;
      }
    nphotok = 0;
    for (samp2 = samp; samp2 && (p=samp2->set->field->photomlabel)>=0;
                samp2=samp2->prevsamp)
      {
      if (samp2->sexflags & photflagmask)
        continue;
      if (samp2->flux > 0.0 && (err2 = samp2->magerr*samp2->magerr)>0.0)
        {
        magerr[p] += 1.0 / err2;
        mag[p] += samp2->mag / err2;
        if (!nmag[p])
          magref[p] = samp2->mag;
        magdisp[p] += (samp2->mag-magref[p]) * (samp2->mag - magref[p]);
        nmag[p]++;
        }
      nphotok++;
      }
    for (p=0; p<npinstru; p++)
      {
      if ((nphotok=nmag[p]))
        {
        mergedsample.mag[p] = mag[p] / magerr[p];
        mergedsample.magerr[p] = sqrt(1.0 / magerr[p]);
        mergedsample.magdisp[p] = nphotok > 1? sqrt(fabs(magdisp[p]
			 - nphotok*(mergedsample.mag[p] - magref[p])
		 	*(mergedsample.mag[p] - magref[p]))/(nphotok-1.0))
		: 0.0;
        }
      else
        mergedsample.mag[p] = mergedsample.magerr[p]
		= mergedsample.magdisp[p] = 99.0;
        mergedsample.nmag[p] = nmag[p];
      }

    mergedsample.colour = msamp->colour;

/*-- Astrometry */
    for (d=0; d<naxis; d++)
      {
      mergedsample.wcspos[d] = msamp->wcspos[d];
      mergedsample.wcsposerr[d] = msamp->wcsposerr[d];
      mergedsample.wcsprop[d] = msamp->wcsprop[d] * DEG/MAS;
      mergedsample.wcsproperr[d] = msamp->wcsproperr[d] * DEG/MAS;
      }
    mergedsample.wcschi2 = msamp->wcschi2;

/*-- Epochs */
    mergedsample.epoch = msamp->epoch;
    mergedsample.epochmin = msamp->epochmin;
    mergedsample.epochmax = msamp->epochmax;
    mergedsample.npos_tot = msamp->npos_tot;
    mergedsample.npos_ok = msamp->npos_ok;

/*-- Morphometry */
    mergedsample.spread = msamp->spread;
    mergedsample.spreaderr = msamp->spreaderr;

/*-- Flags */
    mergedsample.sexflags = msamp->sexflags;
    mergedsample.scampflags = msamp->scampflags;

/*-- Write to the catalog */
    if (prefs.mergedcat_type == CAT_ASCII_HEAD
		|| prefs.mergedcat_type == CAT_ASCII
		|| prefs.mergedcat_type == CAT_ASCII_SKYCAT)
      print_obj(ascfile, objtab);
    else if (prefs.mergedcat_type == CAT_ASCII_VOTABLE)
      voprint_obj(ascfile, objtab);
    else
      write_obj(objtab, buf);
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


/****** writefullcat_fgroup **************************************************
PROTO	void writefullcat_fgroup(char *filename, fgroupstruct *fgroup)
PURPOSE	Save a SExtractor-like catalog containing all individual detections as
	calibrated by SCAMP.
INPUT	File name,
	pointer to the fgroup structure.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 27/08/2012
*/
void	writefullcat_fgroup(char *filename, fgroupstruct *fgroup)

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
   fullsamplestruct	fsample;
   fieldstruct		*field;
   setstruct		*set;
   samplestruct		*samp,*samp2;
   FILE			*ascfile;
   double		wcspos[NAXIS], wcsposerr[NAXIS],
			mag, magerr,epoch;
   char			str[80],
			*buf, *rfilename;
   long			dptr;
   int			nmag[MAXPHOTINSTRU],
			d,f,i,k,n,p,s,nm, npinstru, naxis, index;

  if (prefs.fullcat_type == CAT_NONE)
    return;

  naxis = fgroup->naxis;

/* LDAC Object header */
  objtab = new_tab("LDAC_OBJECTS");
/* Set key pointers */
  QCALLOC(objkeys, keystruct, (sizeof(reffullkey) / sizeof(keystruct)));
  dptr = (long)((char *)&fsample - (char *)&reffullsample);
  for (k=0; reffullkey[k].name[0]; k++)
    {
    if (!prefs.spread_flag
	&& !cistrcmp(reffullkey[k].name, "SPREAD", FIND_NOSTRICT))
      continue;
    objkeys[k] = reffullkey[k];
    key = objkeys+k;
/*-- A trick to access the fields of the dynamic fullsample structure */
    key->ptr = (void *)((char *)key->ptr + dptr);
    key->nbytes = t_size[key->ttype]*(key->naxis? *key->naxisn : 1);
    add_key(key,objtab, 0);
    }
/* Create a new output catalog */
  if (prefs.fullcat_type == CAT_ASCII_HEAD
	|| prefs.fullcat_type == CAT_ASCII
	|| prefs.fullcat_type == CAT_ASCII_SKYCAT
	|| prefs.fullcat_type == CAT_ASCII_VOTABLE)
    {
    cat = NULL;
    if (prefs.fullcatpipe_flag)
      ascfile = stdout;
    else
      if (!(ascfile = fopen(filename, "w+")))
        error(EXIT_FAILURE,"*Error*: cannot open ", filename);
    if (prefs.fullcat_type == CAT_ASCII_HEAD && (key = objtab->key))
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
    else if (prefs.fullcat_type == CAT_ASCII_SKYCAT && (key = objtab->key))
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
    else if (prefs.fullcat_type == CAT_ASCII_VOTABLE && objtab->key) 
      {
/*---- A short, "relative" version of the filename */
      if (!(rfilename = strrchr(filename, '/')))
        rfilename = filename;
      else
        rfilename++;
      write_xml_header(ascfile);
      fprintf(ascfile,
	" <TABLE ID=\"Full_List\" name=\"%s/out\">\n", rfilename);
      fprintf(ascfile,
        "  <DESCRIPTION>Table of all detections by %s</DESCRIPTION>\n",
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
    key->naxisn[1] = fitsfind(key->ptr, "END     ")+1;
    key->htype = H_STRING;
    key->ttype = T_STRING;
    key->nobj = 1;
    key->nbytes = key->naxisn[0]*key->naxisn[1];
    add_key(key, imtab, 0);
    save_tab(cat, imtab);
    free_tab(imtab);
    objtab->cat = cat;
    init_writeobj(cat, objtab, &buf);
    }

  index = 0;
  for (f=0; f<fgroup->nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      samp = set->sample;
      for (n=set->nsample; n--; samp++)
        if (!samp->nextsamp)
          {
          ++index;
          for (samp2 = samp; samp2; samp2=samp2->prevsamp)
            {
/*---------- Indices */
            fsample.sourceindex = index;
            fsample.fieldindex = samp2->set->field->astromlabel>=0?
					samp2->set->field->fieldindex+1 : 0;
            fsample.setindex = samp2->set->setindex+1;
            fsample.astrinstruindex = samp2->set->field->astromlabel>=0?
					samp2->set->field->astromlabel+1 : 0;
            fsample.photinstruindex = samp2->set->field->photomlabel>=0?
					samp2->set->field->photomlabel+1 : 0;
/*---------- Astrometry */
            for (d=0; d<naxis; d++)
              {
              fsample.rawpos[d] = samp2->rawpos[d];
              fsample.rawposerr[d] = samp2->rawposerr[d];
              fsample.rawpostheta = 0.0;
              fsample.wcspos[d] = samp2->wcspos[d];
              fsample.wcsposerr[d] = samp2->wcsposerr[d];
              fsample.wcspostheta = 0.0;
              fsample.epoch = samp2->epoch;
              }
/*---------- Photometry */
            fsample.mag = samp2->mag;
            fsample.magerr = samp2->magerr;
/*---------- Morphometry */
            fsample.spread = samp2->spread;
            fsample.spreaderr = samp2->spreaderr;
/*---------- Flags */
            fsample.sexflags = samp2->sexflags;
            fsample.scampflags = samp2->scampflags;
/*-------- Write to the catalog */
            if (prefs.fullcat_type == CAT_ASCII_HEAD
		|| prefs.fullcat_type == CAT_ASCII
		|| prefs.fullcat_type == CAT_ASCII_SKYCAT)
              print_obj(ascfile, objtab);
            else if (prefs.fullcat_type == CAT_ASCII_VOTABLE)
              voprint_obj(ascfile, objtab);
            else
              write_obj(objtab, buf);
            }
          }
      }
    }

  if (prefs.fullcat_type == CAT_ASCII_HEAD
	|| prefs.fullcat_type == CAT_ASCII
	|| prefs.fullcat_type == CAT_ASCII_SKYCAT)
    {
    if (prefs.fullcat_type == CAT_ASCII_SKYCAT)
      fprintf(ascfile, skycattail);
    if (!prefs.fullcatpipe_flag)
      fclose(ascfile);
    }
  else if (prefs.fullcat_type == CAT_ASCII_VOTABLE)
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


