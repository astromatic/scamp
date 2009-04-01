 /*
				catout.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Produce and write merged catalogs.
*
*	Last modify:	28/09/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
#include "fgroup.h"
#include "field.h"
#include "prefs.h"
#include "samples.h"
#ifdef USE_THREADS
#include "threads.h"
#endif


mergedsamplestruct    refmergedsample;
keystruct       refmergedkey[] = {
  {"X_WORLD", "Barycenter position along world x axis",
        &refmergedsample.wcspos[0], H_FLOAT, T_DOUBLE,
	"%15e", "deg", "pos.eq.ra;stat.mean", "deg"},
  {"Y_WORLD", "Barycenter position along world y axis",
        &refmergedsample.wcspos[1], H_FLOAT, T_DOUBLE,
	"%15e", "deg", "pos.eq.de;stat.mean", "deg"},
  {"ERRA_WORLD", "RMS position error along major world axis",
        &refmergedsample.wcsposerr[0], H_FLOAT, T_FLOAT,
	"%12e", "deg", "stat.error;stat.max;pos.errorEllipse;meta.main", "deg"},
  {"ERRB_WORLD", "RMS position error along minor world axis",
        &refmergedsample.wcsposerr[1], H_FLOAT, T_FLOAT,
	"%12e", "deg", "stat.error;stat.min;pos.errorEllipse;meta.main", "deg"},
  {"DISPA_WORLD", "RMS dispersion of pos along major world axis",
        &refmergedsample.wcsposdisp[0], H_FLOAT, T_FLOAT,
	"%12e", "deg", "stat.stdev;stat.max;pos.errorEllipse;meta.main", "deg"},
  {"DISPB_WORLD", "RMS dispersion of pos along minor world axis",
        &refmergedsample.wcsposdisp[1], H_FLOAT, T_FLOAT,
	"%12e", "deg", "stat.stdev;stat.min;pos.errorEllipse;meta.main", "deg"},
  {"PMX_WORLD", "Proper motion along world x axis",
        &refmergedsample.wcsprop[0], H_FLOAT, T_FLOAT,
	"%12e", "deg", "pos.pm;pos.eq.ra;stat.fit", "mas/yr"},
  {"PMY_WORLD", "Proper motion along world y axis",
        &refmergedsample.wcsprop[1], H_FLOAT, T_FLOAT,
	"%12e", "deg", "pos.pm;pos.eq.de;stat.fit", "mas/yr"},
  {"PMXERR_WORLD", "P.motion uncertainty along world x axis",
        &refmergedsample.wcsproperr[0], H_FLOAT, T_FLOAT,
	"%12e", "deg", "stat.error; pos.pm;pos.eq.ra;stat.fit", "mas/yr"},
  {"PMYERR_WORLD", "P.motion uncertainty along world y axis",
        &refmergedsample.wcsproperr[1], H_FLOAT, T_FLOAT,
	"%12e", "deg", "stat.error;pos.pm;pos.eq.de;stat.fit", "mas/yr"},
  {"EPOCH", "Mean epoch",
        &refmergedsample.epoch, H_FLOAT, T_FLOAT,
	"%12e", "time.epoch;stat.mean", ""},
  {"NMAG", "Generic magnitude vector",
        &refmergedsample.nmag, H_INT, T_LONG,
	"%4d", "", "meta.number", "",
	1, &refmergedsample.nband},
  {"MAG_VEC", "Generic magnitude vector",
        &refmergedsample.mag, H_FLOAT, T_FLOAT,
	"%8.4f", "mag", "phot.mag", "mag",
	1, &refmergedsample.nband},
  {"MAGERR_VEC", "Generic magnitude RMS error estimate vector",
        &refmergedsample.magerr, H_FLOAT, T_FLOAT,
	"%8.4f", "mag", "stat.error;phot.mag", "mag",
	1, &refmergedsample.nband},
  {"MAG_DISP", "Generic magnitude RMS dispersion",
        &refmergedsample.magdisp, H_FLOAT, T_FLOAT,
	"%8.4f", "mag", "stat.stdev;phot.mag", "mag",
	1, &refmergedsample.nband},
  {"FLAGS", "SCAMP flags",
        &refmergedsample.flags, H_INT, T_SHORT,
	"%3d", "", "meta.code.qual", ""},
  {""},
  };


/****** writemergedcat_fgroup *************************************************
PROTO	void writemergedcat_fgroup(char *filename, fgroupstruct *fgroup)
PURPOSE	Save a SExtractor-like catalog containing merged detections calibrated
	by SCAMP.
INPUT	File name,
	pointer to the fgroup structure.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 05/10/2008
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
			wcspos[NAXIS], wcsprop[NAXIS],wcsproperr[NAXIS],
			epoch, err2;
   char			*buf;
   int			nmag[MAXPHOTINSTRU],
			d,f,i,k,n,p,s,nm, npinstru, naxis;

  if (prefs.mergedcat_type == CAT_NONE)
    return;

  naxis = fgroup->naxis;
  refmergedsample.nband = npinstru = prefs.nphotinstrustr;

/* Create a new output catalog */
  if (prefs.mergedcat_type == CAT_ASCII_HEAD
	|| prefs.mergedcat_type == CAT_ASCII
	|| prefs.mergedcat_type == CAT_ASCII_SKYCAT)
    {
    if (prefs.mergedcatpipe_flag)
      ascfile = stdout;
    else
      if (!(ascfile = fopen(filename, "w+")))
        error(EXIT_FAILURE,"*Error*: cannot open ", filename);
    fclose(ascfile);
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

/*-- LDAC Object header */
    objtab = new_tab("LDAC_OBJECTS");
    objtab->cat = cat;
/*-- Set key pointers */
    QCALLOC(objkeys, keystruct, (sizeof(refmergedkey) / sizeof(keystruct)));
    for (k=0; refmergedkey[k].name[0]; k++)
      {
      objkeys[k] = refmergedkey[k];
/*---- A trick to access the fields of the dynamic mergedsample structure */
      objkeys[k].ptr += (void *)&msample - (void *)&refmergedsample;
      add_key(&objkeys[k],objtab, 0);
      }
    init_writeobj(cat, objtab, &buf);
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
            for (p=0; p<npinstru; p++)
              {
              mag[p] = magerr[p] = magdisp[p] = magchi2[p] = 0.0;
              nmag[p] = 0;
              }
            for (samp2 = samp;
		samp2 && (p=samp2->set->field->photomlabel)>=0;
                samp2=samp2->prevsamp)
              {
              magerr[p] += 1.0/(err2 = samp2->magerr*samp2->magerr);
              mag[p] += samp2->mag / err2;
              magdisp[p] += (samp2->mag*samp2->mag);
              nmag[p]++;
              }
            for (p=0; p<npinstru; p++)
              {
              if ((nm=nmag[p]))
                {
                msample.mag[p] = mag[p] / magerr[p];
                msample.magerr[p] = sqrt(nm / magerr[p]);
                msample.magdisp[p] = sqrt(fabs(magdisp[p] / nm
					- msample.mag[p]*msample.mag[p]));
                }
              else
                msample.mag[p] = msample.magerr[p] = msample.magdisp[p] = 99.0;
              msample.nmag[p] = nmag[p];
              }
            nm = 0;
            for (d=0; d<naxis; d++)
                wcspos[d] = wcsprop[d] = wcsproperr[d] =  0.0;
            epoch = 0.0;
            for (samp2 = samp;
		samp2 && (p=samp2->set->field->photomlabel)>=0;
                samp2=samp2->prevsamp)
              {
              for (d=0; d<naxis; d++)
                {
                wcspos[d] += samp2->wcspos[d];
//			- colshiftscale[d][p+npinstru*refplabel];
                wcsprop[d] += samp2->wcsprop[d];
                wcsproperr[d] += samp2->wcsproperr[d]*samp2->wcsproperr[d];
                }
              epoch += samp2->set->field->epoch;
              nm++;
              }
            if (nm)
              {
              for (d=0; d<naxis; d++)
                {
                msample.wcspos[d] = wcspos[d] / nm;
                msample.wcsprop[d] = (wcsprop[d]/nm) * DEG/MAS;
                msample.wcsproperr[d] = sqrt(wcsproperr[d]/nm) * DEG/MAS;
                }
              msample.epoch = epoch / nm;
              }
            write_obj(objtab, buf);
            }
        }
      }

    end_writeobj(cat, objtab, buf);
    objtab->key = NULL;
    objtab->nkey = 0;
    free_tab(objtab);
    free(objkeys);
    free_cat(&cat, 1);
    }

  return;
  }
