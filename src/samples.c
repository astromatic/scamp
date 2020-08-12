/*
 *    samples.c
 *
 * Read and filter input detections from catalogs.
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 * This file part of: SCAMP
 *
 * Copyright:  (C) 2002-2020 IAP/CNRS/SorbonneU
 *
 * License:  GNU General Public License
 *
 * SCAMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * SCAMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
 *
 * Last modified:  27/04/2020
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "field.h"
#include "prefs.h"
#include "samples.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include "wcs/wcs.h"

//#define WFCAM_FIX 1

/*------------------- global variables for multithreading -------------------*/
#ifdef USE_THREADS
extern pthread_mutex_t sortmutex;
#endif

int compraw_samples(const void *, const void *);
int compproj_samples(const void *, const void *);
int compwcs_samples(const void *, const void *);
int sort_coord;

/****** read_samples *********************************************************
  PROTO setstruct *read_samples(setstruct *set, tabstruct *tab,char *rfilename)
  PURPOSE Read a set of samples.
  INPUT set structure pointer,
  pointer to the tab that contains the catalog,
  reduced filename.
  OUTPUT  setstruct pointer (allocated if the input setstruct pointer is NULL).
  NOTES   The filename is used for error messages only. Global preferences are
  used.
  AUTHOR  E. Bertin (IAP)
  VERSION  29/03/2019
 */
setstruct *read_samples(setstruct *set, tabstruct *tab, char *rfilename)

{
    tabstruct  *keytab;
    keystruct  *key;
    samplestruct  *sample;
    wcsstruct  *wcs;
    t_type  contexttyp[MAXCONTEXT];
    void   *context[MAXCONTEXT];
    char   str[MAXCHAR],
    **kstr,
    *buf,*head,*head0;
    double  contextval[MAXCONTEXT], *cmin,*cmax, *dfluxrad, *delong,
            *dxm, *dym, *derra,*derrb, *dflux,*dfluxerr, *dspread,
            *dspreaderr,
            dminrad,dmaxrad, dmaxelong, dval,
            x,y, ea,eb, f, ferr, xmax,ymax, max;
    float  *xm,*ym, *flux, *fluxerr, *erra,*errb, *fluxrad, *elong,
           *spread,*spreaderr,
           minrad,maxrad, maxelong;
    unsigned int  *imaflags;
    int   *lxm,*lym,
          i, n, nsample,nsamplemax, nnobjflag,
          nobj, headflag, head0flag, errorflag, fflag;
    unsigned short *flags, *wflags,
		sexflags;
    short  *sxm,*sym;

#ifdef WFCAM_FIX
    int wfcamflag;
#endif

    head = head0 = NULL; /* to avoid gcc -Wall warnings */
    cmin = cmax = NULL; /* to avoid gcc -Wall warnings */
    dxm = dym = derra = derrb = dflux = dfluxerr = dfluxrad = delong = NULL;
    xm = ym = erra = errb = flux = fluxerr = fluxrad = elong = NULL;
    lxm = lym = NULL;
    sxm = sym = NULL;
    minrad = (float)(dminrad = prefs.fwhm_thresh[0]/2.0);
    maxrad = (float)(dmaxrad = prefs.fwhm_thresh[1]/2.0);
    maxelong = (float)(dmaxelong = prefs.maxellip < 1.0?
            (prefs.maxellip + 1.0)/(1.0 - prefs.maxellip)
            : 100.0);

    /* If a NULL pointer is provided, we allocate a new set */
    if (!set)
    {
        set = init_set();
        nsample = nsamplemax = 0;
    }
    else
        nsample = nsamplemax = set->nsample;

    if (set->ncontext)
    {
        QMALLOC(cmin, double, set->ncontext);
        QMALLOC(cmax, double, set->ncontext);
        for (i=0; i<set->ncontext; i++)
        {
            cmin[i] = nsample? set->contextoffset[i]-set->contextscale[i]/2.0 : BIG;
            cmax[i] = nsample? cmin[i] + set->contextscale[i] : -BIG;
        }
    }

    /* Init the single-row tab */
    keytab = init_readobj(tab, &buf);

    if (!(key = name_to_key(keytab, prefs.centroid_key[0])))
    {
        sprintf(str, "*Error*: %s parameter not found in catalog ",
                prefs.centroid_key[0]);
        error(EXIT_FAILURE, str, rfilename);
    }
    if (key->ttype == T_DOUBLE)
        dxm = (double *)key->ptr;
    else if (key->ttype == T_FLOAT)
        xm = (float *)key->ptr;
    else if (key->ttype == T_LONG)
        lxm = (int *)key->ptr;
    else
        sxm = (short *)key->ptr;

    nobj = key->nobj;

    if (!(key = name_to_key(keytab, prefs.centroid_key[1])))
    {
        sprintf(str, "*Error*: %s parameter not found in catalog ",
                prefs.centroid_key[1]);
        error(EXIT_FAILURE, str, rfilename);
    }
    if (key->ttype == T_DOUBLE)
        dym = (double *)key->ptr;
    else if (key->ttype == T_FLOAT)
        ym = (float *)key->ptr;
    else if (key->ttype == T_LONG)
        lym = (int *)key->ptr;
    else
        sym = (short *)key->ptr;

    if (!(key = name_to_key(keytab, prefs.centroiderr_key[0])))
    {
        sprintf(str, "*Error*: %s parameter not found in catalog ",
                prefs.centroiderr_key[0]);
        error(EXIT_FAILURE, str, rfilename);
    }
    if (key->ttype == T_DOUBLE)
        derra = (double *)key->ptr;
    else
        erra = (float *)key->ptr;

    if (!(key = name_to_key(keytab, prefs.centroiderr_key[1])))
    {
        sprintf(str, "*Error*: %s parameter not found in catalog ",
                prefs.centroiderr_key[1]);
        error(EXIT_FAILURE, str, rfilename);
    }
    if (key->ttype == T_DOUBLE)
        derrb = (double *)key->ptr;
    else
        errb = (float *)key->ptr;

    if (!(key = name_to_key(keytab, prefs.photflux_rkey)))
    {
        sprintf(str, "*Error*: %s parameter not found in catalog ",
                prefs.photflux_rkey);
        error(EXIT_FAILURE, str, rfilename);
    }
    if (key->ttype == T_DOUBLE)
        dflux = (double *)key->ptr;
    else
        flux = (float *)key->ptr;
    n = prefs.photflux_num - 1;
    if (n)
    {
        if (key->naxis==1 && n<key->naxisn[0])
        {
            if (key->ttype == T_DOUBLE)
                dflux += n;
            else
                flux += n;
        }
        else
        {
            sprintf(str, "Not enough apertures for %s in catalog %s: ",
                    prefs.photflux_rkey, rfilename);
            warning(str, "using first aperture");
        }
    }

    if (!(key = name_to_key(keytab, prefs.photfluxerr_rkey)))
    {
        sprintf(str, "*Error*: %s parameter not found in catalog ",
                prefs.photfluxerr_rkey);
        error(EXIT_FAILURE, str, rfilename);
    }
    if (key->ttype == T_DOUBLE)
        dfluxerr = (double *)key->ptr;
    else
        fluxerr = (float *)key->ptr;
    n = prefs.photfluxerr_num - 1;
    if (n)
    {
        if (key->naxis==1 && n<key->naxisn[0])
        {
            if (key->ttype == T_DOUBLE)
                dfluxerr += n;
            else
                fluxerr += n;
        }
        else
        {
            sprintf(str, "Not enough apertures for %s in catalog %s: ",
                    prefs.photfluxerr_rkey, rfilename);
            warning(str, "using first aperture");
        }
    }

    /* Load optional SExtractor FLAGS parameter */
    if (!(key = name_to_key(keytab, "FLAGS")))
        warning("FLAGS parameter not found in catalog ", rfilename);
    flags = key? (unsigned short *)key->ptr : NULL;

    /* Load optional SExtractor IMAFLAGS_ISO parameter */
    if ((key = name_to_key(keytab, "IMAFLAGS_ISO")))
        imaflags = (unsigned int *)key->ptr;
    else
        imaflags = NULL;

    /* Load optional SExtractor FLAGS_WEIGHT parameter */
    if ((key = name_to_key(keytab, "FLAGS_WEIGHT")))
        wflags = (unsigned short *)key->ptr;
    else
        wflags = NULL;

    /* Load optional SExtractor FLUX_RADIUS parameter */
    if ((key = name_to_key(keytab, "FLUX_RADIUS")))
    {
        if (key->ttype == T_DOUBLE)
            dfluxrad = (double *)key->ptr;
        else
            fluxrad = (float *)key->ptr;
    }
    else
    {
        dfluxrad = NULL;
        fluxrad = NULL;
    }

    /* Load optional SExtractor ELONGATION parameter */
    if ((key = name_to_key(keytab, "ELONGATION")))
    {
        if (key->ttype == T_DOUBLE)
            delong = (double *)key->ptr;
        else
            elong = (float *)key->ptr;
    }
    else
    {
        delong = NULL;
        elong = NULL;
    }

    /* Load optional SExtractor SPREAD_MODEL parameter */
    if ((key = name_to_key(keytab, "SPREAD_MODEL")))
    {
        if (key->ttype == T_DOUBLE)
            dspread = (double *)key->ptr;
        else
            spread = (float *)key->ptr;
        prefs.spread_flag = 1;
    }
    else
    {
        dspread = NULL;
        spread = NULL;
    }

    /* Load optional SExtractor SPREADERR_MODEL parameter */
    if ((key = name_to_key(keytab, "SPREADERR_MODEL")))
    {
        if (key->ttype == T_DOUBLE)
            dspreaderr = (double *)key->ptr;
        else
            spreaderr = (float *)key->ptr;
    }
    else
    {
        dspreaderr = NULL;
        spreaderr = NULL;
    }

    /* Try to load the set of context keys */
    head0flag = (set->imatab && set->imatab->cat && set->imatab->cat->tab
            && (head0=set->imatab->cat->tab->headbuf));
    headflag = (set->imatab && (head = set->imatab->headbuf));

    /* Try to load the set of context keys */
    kstr = prefs.context_name;
    for (i=0; i<set->ncontext; i++, kstr++)
        if (**kstr==(char)':' && (headflag|head0flag))
        {
            context[i] = &contextval[i];
            contexttyp[i] = T_DOUBLE;
            if (headflag)
            {
                errorflag = fitsread(head, *kstr+1, context[i], H_FLOAT,T_DOUBLE);
                if (errorflag==RETURN_ERROR && head0flag)
                    errorflag = fitsread(head0, *kstr+1, context[i],H_FLOAT,T_DOUBLE);
                if (errorflag==RETURN_ERROR)
                {
                    sprintf(str, "*Error*: %s parameter not found in the header of ",
                            *kstr+1);
                    error(EXIT_FAILURE, str, rfilename);
                }
            }
            strcpy(set->contextname[i], *kstr);
        }
        else
        {
            if (!(key = name_to_key(keytab, *kstr)))
            {
                sprintf(str, "*Error*: %s parameter not found in catalog ", *kstr);
                error(EXIT_FAILURE, str, rfilename);
            }
            context[i] = key->ptr;
            contexttyp[i] = key->ttype;
            strcpy(set->contextname[i], key->name);
        }
#ifdef WFCAM_FIX
    wfcamflag = strncmp(rfilename,"w2", 2)? 0 : 1;
#endif

    /* Read photometric parameters */
    if (headflag)
    {
        if (fitsread(head, prefs.astraccuracy_key, &set->astraccuracy,
                    H_FLOAT, T_DOUBLE) != RETURN_OK)
            set->astraccuracy = prefs.astraccuracy;
        if (fitsread(head, prefs.airmass_key, &set->airmass,
                    H_FLOAT, T_DOUBLE) != RETURN_OK)
            set->airmass = 1.0;
        set->airmass = fabs(set->airmass);
        if (fitsread(head, prefs.expotime_key, &set->expotime,
                    H_FLOAT, T_DOUBLE) != RETURN_OK)
            set->expotime = 1.0;
        set->epochmin = set->wcs->obsdate;
        set->epochmax = set->epochmin==0.0? 0.0 : set->epochmin+set->expotime/YEAR;
        set->epoch = set->epochmin==0.0? 0.0 : set->epochmin+set->expotime/2.0/YEAR;
        if (fitsread(head, prefs.extcoeff_key, &set->extcoeff,
                    H_FLOAT, T_DOUBLE) != RETURN_OK)
            set->extcoeff = 0.0;
        set->extcoeff = fabs(set->extcoeff);
        if (fitsread(head, prefs.magzero_key, &set->magzero,
                    H_FLOAT, T_DOUBLE) != RETURN_OK)
            set->magzero = 0.0;
        /*-- Photometric flag must not be unset by other sets of the same field */
        fitsread(head, prefs.photomflag_key,&set->field->photomflag,H_BOOL,T_LONG);
    }

    nnobjflag = 0;

    /* Set maximum x and y */
    xmax = ymax = 1.5;
    if (set->wcs && set->wcs->naxis>0)
        xmax = set->wcs->naxisn[0]+0.5;
    if (set->wcs && set->wcs->naxis>1)
        ymax = set->wcs->naxisn[1]+0.5;
    max = xmax>ymax ? xmax : ymax;

    /* Now examine each vector of the shipment */
    for (n=0; nobj--; n++)
    {
        sexflags = 0;
        read_obj(keytab,tab, buf);
#ifdef USE_THREADS
        if (!(n%10000) && prefs.nthreads<2)
#else
            if (!(n%10000))
#endif
            {
                sprintf(str,"%-.36s : %d / %d detections stored",
                        rfilename,nsample,n);
                NFPRINTF(OUTPUT, str);
            }
        /*---- Apply some selection over flags, fluxes... */
        /*---- No saturated or cropped detection */
        if (flags)
        {
            if (*flags & prefs.sexflags_mask)
                continue;
            /*---- Mapping from SExtractor flags is straightforward */
            sexflags = *flags;
            if (sexflags & OBJ_SATUR)  /* A saturated object */
                set->nsaturated++;
            else if ((fluxrad && (*fluxrad < minrad || *fluxrad > maxrad))
                    || (dfluxrad && (*dfluxrad < dminrad || *dfluxrad > dmaxrad))
                    || (elong && *elong > maxelong)
                    || (delong && *delong > dmaxelong))
                continue;
        }
        else if ((fluxrad && (*fluxrad < minrad || *fluxrad > maxrad))
                || (dfluxrad && (*dfluxrad < dminrad || *dfluxrad > dmaxrad))
                || (elong && *elong > maxelong)
                || (delong && *delong > dmaxelong))
            continue;
        if (wflags)
        {
            if ((*wflags & prefs.wflags_mask))
                continue;
            if (*wflags)
                sexflags |= OBJ_TRUNC;
        }

        if (imaflags && (*imaflags & prefs.imaflags_mask))
            continue;

        if (dxm)
            x = *dxm;
        else if (xm)
            x = *xm;
        else if (lxm)
            x = *lxm;
        else
            x = *sxm;
        if (dym)
            y = *dym;
        else if (ym)
            y = *ym;
        else if (lym)
            y = *lym;
        else
            y = *sym;
        if (x<0.5 || x>xmax || y<0.5 || y>ymax)
            continue;

#ifdef WFCAM_FIX
        if (wfcamflag)
        {
            double dx,dy;
            int xflag,yflag;

            xflag = x<1024.5;
            yflag = y<1024.5;
            if (yflag)
            {
                if (xflag)
                    /*-------- lower left quadrant */
                {
                    dx = -0.015;
                    dy = -0.004;
                }
                else
                {
                    /*-------- lower right quadrant */
                    dx = -0.014;
                    dy = -0.019;
                }
            }
            else
            {
                if (xflag)
                    /*-------- upper left quadrant */
                {
                    dx = 0.013;
                    dy = 0.026;
                }
                else
                {
                    /*-------- upper right quadrant */
                    dx = 0.021;
                    dy = -0.022;
                }
            }
            x -= dx;
            if (dxm)
                *dxm -= dx;
            else if (xm)
                *xm -= dx;
            y -= dy;
            if (dym)
                *dym -= dy;
            else if (ym)
                *ym -= dy;
        }
#endif

        f = flux? *flux:*dflux;
#ifdef HAVE_ISNAN2
        fflag = isnan(f);
#endif
        ferr = fluxerr? *fluxerr:*dfluxerr;
        if (ferr<MIN_FLUXERR || fflag || f>MAX_FLUX
                || (ferr>0.0 && f/ferr < prefs.sn_thresh[0]))
            continue;
        ea= erra? *erra:*derra;
        eb= errb? *errb:*derrb;

        if (ea<MIN_POSERR || eb<MIN_POSERR || ea>max || eb>max)
            continue;

        ferr = sqrt(ferr*ferr+f*f*prefs.photaccuracy*prefs.photaccuracy);

        /*-- ... and check the integrity of the sample */
        /*-- Allocate memory for the first shipment */
        if (!set->sample)
        {
            nsample = 0;
            nsamplemax = LSAMPLE_DEFSIZE;
            malloc_samples(set, nsamplemax);
        }

        /*-- Increase storage space to receive new candidates if needed */
        if (nsample>=nsamplemax)
        {
            int nadd=(int)(1.62*nsamplemax);
            nsamplemax = nadd>nsamplemax?nadd:nsamplemax+1;
            realloc_samples(set, nsamplemax);
        }

        if (!sexflags)
            nnobjflag++;

        sample = set->sample + nsample;
        sample->set = set;
        sample->sexflags = sexflags;
        sample->scampflags = 0;
        sample->flux = f;
        sample->fluxerr = ferr;
        sample->mag = (sample->flux>0.0)? -2.5*log10(sample->flux) : 99.0;
        sample->rawpos[0] = x;
        sample->rawpos[1] = y;
        if (raw_to_wcs(set->wcs, sample->rawpos, sample->wcspos) != RETURN_OK) {
            warning("Invalid deprojected coordinates in " , rfilename);
        }

        sample->epoch = set->epoch;

        if (fluxrad)
            sample->fwhm = 2.0**fluxrad;
        else if (dfluxrad)
            sample->fwhm = 2.0**dfluxrad;
        else
            sample->fwhm = 0.0;

        if (spread)
            sample->spread = *spread;
        else if (dspread)
            sample->spread = *dspread;
        else
            sample->spread = 0.0;

        if (spreaderr)
            sample->spreaderr = *spreaderr;
        else if (dspreaderr)
            sample->spreaderr = *dspreaderr;
        else
            sample->spreaderr = 0.0;

        if (imaflags)
            sample->imaflags = *imaflags;

        sample->rawposerr[0] = sample->rawposerr[1] = sqrt(0.5*(ea*ea+eb*eb));

        /*-- In case of a contamination, position errors are strongly degraded */
        if (flags)
        {
            if (*flags&2)
                sample->rawposerr[0] = (sample->rawposerr[1]
                        = sqrt(sample->rawposerr[0]*sample->rawposerr[0]+0.01));
            /*
               else if (*flags&1)
               sample->rawposerr[0] = (sample->rawposerr[1] *= 3.0);
             */
        }

        for (i=0; i<set->ncontext; i++)
        {
            ttypeconv(context[i], &dval, contexttyp[i], T_DOUBLE);
            sample->context[i] = dval;
            /*---- Update min and max */
            if (dval<cmin[i])
                cmin[i] = dval;
            if (dval>cmax[i])
                cmax[i] = dval;
        }

        nsample++;
    }

    /* Update the scaling */
    if (set->ncontext && nsample)
    {
        for (i=0; i<set->ncontext; i++)
        {
            set->contextscale[i] = cmax[i] - cmin[i];
            set->contextoffset[i] = (cmin[i] + cmax[i])/2.0;
        }
        free(cmin);
        free(cmax);
    }

    end_readobj(keytab, tab, buf);
    if (nsample && !nnobjflag)
        warning("All sources have non-zero flags in ", rfilename);
    if (!nsample)
        warning("No valid source found in ", rfilename);

    set->nsample = nsample;

    /* Don't waste memory! */
    if (nsample)
        realloc_samples(set, nsample);

    return set;
}


/****** malloc_samples *******************************************************
  PROTO   void malloc_samples(setstruct *set, int nsample)
  PURPOSE Allocate memory for a set of samples.
  INPUT   set structure pointer,
  desired number of samples.
  OUTPUT  -.
  NOTES   -.
  AUTHOR  E. Bertin (IAP)
  VERSION 27/04/2020
 */
void malloc_samples(setstruct *set, int nsample)

{
    samplestruct *sample;
    int  n;

    QCALLOC(set->sample, samplestruct, nsample);
    if (set->ncontext) {
        sample = set->sample;
        for (n=nsample; n--; sample++) {
            QMALLOC(sample->context, double, set->ncontext);
        }
    }

    set->nsamplemax = nsample;

    return;
}


/****** realloc_samples ******************************************************
  PROTO   void realloc_samples(setstruct *set, int nsample)
  PURPOSE Re-allocate memory for a set of samples.
  INPUT   set structure pointer,
  desired number of samples.
  OUTPUT  -.
  NOTES   -.
  AUTHOR  E. Bertin (IAP), S.Serre (LAB)
  VERSION 15/04/2020
 */
void realloc_samples(setstruct *set, int nsample)

{
    samplestruct *sample;
    int  n;

    /* If we want to reallocate 0 samples, better free the whole thing! */
    if (!nsample)
        free_samples(set);

    /* Two cases: either more samples are required, or the opposite! */
    if (nsample>set->nsamplemax)
    {
        QREALLOC(set->sample, samplestruct, nsample);
        sample = set->sample + set->nsamplemax;
        memset(sample, '\0', (nsample - set->nsamplemax)*sizeof(samplestruct));
        if (set->ncontext)
            for (n = nsample - set->nsamplemax; n--; sample++)
                QMALLOC(sample->context, double, set->ncontext);
    }
    else if (nsample<set->nsamplemax)
    {
        sample = set->sample + nsample;
        if (set->ncontext)
            for (n = set->nsamplemax - nsample; n--; sample++)
                free(sample->context);
        QREALLOC(set->sample, samplestruct, nsample);
    }

    set->nsamplemax = nsample;

    return;
}


/****** free_samples *********************************************************
  PROTO   void free_samples(setstruct *set, int nsample)
  PURPOSE free memory for a set of samples.
  INPUT   set structure pointer,
  desired number of samples.
  OUTPUT  -.
  NOTES   -.
  AUTHOR  E. Bertin (IAP)
  VERSION 27/04/2020
 */
void free_samples(setstruct *set)

{
    samplestruct *sample;
    int  n;

    sample = set->sample;
    if (set->ncontext)
        for (n = set->nsamplemax; n--; sample++)
            free(sample->context);

    free(set->sample);
    set->sample = NULL;
    set->nsample = set->nsamplemax = 0;

    return;
}


/****** remove_sample ********************************************************
  PROTO   samplestruct *remove_sample(setstruct *set, int isample)
  PURPOSE Remove an element from a set of samples.
  INPUT   set structure pointer,
  sample number.
  OUTPUT  The new pointer for the element that replaced the removed one.
  NOTES   -.
  AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
  VERSION 10/09/2009
 */
samplestruct *remove_sample(setstruct *set, int isample)

{
    samplestruct  exsample;
    samplestruct  *sample;
    int   nsample;

    /* If we want to reallocate 0 samples, better free the whole thing! */
    nsample = set->nsample-1;
    if (nsample>0)
    {
        sample = set->sample + isample;
        exsample = *(set->sample+nsample);
        *(set->sample+nsample) = *sample;
        *sample = exsample;
    }
    else
        nsample=0;
    realloc_samples(set, nsample);
    set->nsample = nsample;

    return set->sample+isample;
}


/****** init_set *************************************************************
  PROTO   setstruct *init_set(void)
  PURPOSE Generate a new set.
  INPUT   -.
  OUTPUT  Pointer to the new mallocated set structure.
  NOTES   -.
  AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
  VERSION 22/07/2002
 */
setstruct *init_set(void)

{
    setstruct *set;
    int i;

    QCALLOC(set, setstruct, 1);
    set->ncontext = prefs.ncontext_name;
    if (set->ncontext)
    {
        QMALLOC(set->contextoffset, double, set->ncontext);
        QMALLOC(set->contextscale, double, set->ncontext);
        QMALLOC(set->contextname, char *, set->ncontext);
        for (i=0; i<set->ncontext; i++)
            QMALLOC(set->contextname[i], char, 80);
    }

    return set;
}


/****** end_set *************************************************************
  PROTO   void end_set(setstruct *set)
  PURPOSE free memory allocated by a complete set structure.
  INPUT   set structure pointer,
  OUTPUT  -.
  NOTES   -.
  AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
  VERSION 08/08/2002
 */
void end_set(setstruct *set)

{
    int i;

    free_samples(set);
    if (set->ncontext)
    {
        for (i=0; i<set->ncontext; i++)
            free(set->contextname[i]);
        free(set->contextname);
        free(set->contextoffset);
        free(set->contextscale);
    }

    if (set->imatab)
        free_tab(set->imatab);
    if (set->wcs)
        end_wcs(set->wcs);
    free(set);

    return;
}


/****** sort_samples **********************************************************
  PROTO   void sort_samples(setstruct *set)
  PURPOSE Sort samples in y, and compute min/max raw pos. to speed up cross-id.
  INPUT   set structure pointer.
  OUTPUT  -.
  NOTES   Not reentrant (yet).
  AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
  VERSION 28/12/2004
 */
void sort_samples(setstruct *set)

{
    samplestruct *samp;
    double *projmin, *projmax;
    int  e,i, naxis;

    naxis = set->naxis;
    /* Sort sources in latitude */
#ifdef USE_THREADS
    QPTHREAD_MUTEX_LOCK(&sortmutex);
#endif
    sort_coord = (set->lat<0)? set->naxis-1 : set->lat;
    qsort(set->sample, set->nsample, sizeof(samplestruct), compproj_samples);
#ifdef USE_THREADS
    QPTHREAD_MUTEX_UNLOCK(&sortmutex);
#endif
    projmin = set->projposmin;
    projmax = set->projposmax;
    for (i=0; i<naxis; i++)
        projmin[i] = -(projmax[i] = -BIG);
    samp = set->sample;
    for (e=set->nsample; e--; samp++)
    {
        for (i=0; i<naxis; i++)
        {
            if (samp->projpos[i]<projmin[i])
                projmin[i] = samp->projpos[i];
            if (samp->projpos[i]>projmax[i])
                projmax[i] = samp->projpos[i];
        }
    }

    return;
}


/****** unlink_samples *******************************************************
  PROTO   void unlink_samples(setstruct *set)
  PURPOSE Reset links between (possibly) overlapping samples.
  INPUT   set structure pointer.
  OUTPUT  -.
  NOTES   -.
  AUTHOR  E. Bertin (IAP)
  VERSION 18/03/2004
 */
void unlink_samples(setstruct *set)

{
    samplestruct *samp;
    int  e;

    samp = set->sample;
    for (e=set->nsample; e--; samp++)
        samp->prevsamp = samp->nextsamp = NULL;

    return;
}


/****** compproj_samples ******************************************************
  PROTO   int (*compproj_samples(const void *sample1, const void *sample2))
  PURPOSE Provide a sample comparison function for sort_samples() based on
  projected coordinates.
  INPUT   pointer to 1st sample,
  pointer to 2nd sample.
  OUTPUT  <0 if y1<y2, >0 if y1>y2, 0 otherwise .
  NOTES   Not fully reentrant because of the sort_coord global variable.
  AUTHOR  E. Bertin (IAP, Leiden observatory & ESO)
  VERSION 27/12/2004
 */
int compproj_samples(const void *sample1, const void *sample2)
{
    double dy;

    dy = ((samplestruct *)sample1)->projpos[sort_coord]
        - ((samplestruct *)sample2)->projpos[sort_coord];

    return dy>0.0 ? 1: (dy<0.0? -1 : 0);
}


/****** compraw_samples ******************************************************
  PROTO   int (*compraw_samples(const void *sample1, const void *sample2))
  PURPOSE Provide a sample comparison function for sort_samples() based on
  raw coordinates.
  INPUT   pointer to 1st sample,
  pointer to 2nd sample.
  OUTPUT  <0 if y1<y2, >0 if y1>y2, 0 otherwise .
  NOTES   Not fully reentrant because of the sort_coord global variable.
  AUTHOR  E. Bertin (IAP)
  VERSION 27/12/2004
 */
int compraw_samples(const void *sample1, const void *sample2)
{
    double dy;

    dy = ((samplestruct *)sample1)->rawpos[sort_coord]
        - ((samplestruct *)sample2)->rawpos[sort_coord];

    return dy>0.0 ? 1: (dy<0.0? -1 : 0);
}


/****** compwcs_samples ******************************************************
  PROTO   int (*compwcs_samples(const void *sample1, const void *sample2))
  PURPOSE Provide a sample comparison function for sort_samples() based on
  raw coordinates.
  INPUT   pointer to 1st sample,
  pointer to 2nd sample.
  OUTPUT  <0 if y1<y2, >0 if y1>y2, 0 otherwise .
  NOTES   Not fully reentrant because of the sort_coord global variable.
  AUTHOR  E. Bertin (IAP)
  VERSION 27/12/2004
 */
int compwcs_samples(const void *sample1, const void *sample2)
{
    double dy;

    dy = ((samplestruct *)sample1)->wcspos[sort_coord]
        - ((samplestruct *)sample2)->wcspos[sort_coord];

    return dy>0.0 ? 1: (dy<0.0? -1 : 0);
}


/****** union_samples *********************************************************
  PROTO   void union_samples(samplestruct *samplein, setstruct *set,
  int nsample)
  PURPOSE Copy a series of samples to a set, if they are not already there.
  INPUT   pointer to 1st sample,
  pointer to output set,
  number of samples,
  association radius,
  association mode (UNION_RAW, UNION_PROJ or UNION_WCS).
  OUTPUT  -.
  NOTES   Memory for the new samples is reallocated if needed. All input data are
  reordered.
  AUTHOR  E. Bertin (IAP)
  VERSION 27/04/2020
 */
void union_samples(samplestruct *samplein, setstruct *set,
        int nsamplein, double radius, unionmodenum mode)

{
    samplestruct *sample,*sample2,*sample3;
    double *context,
           radius2,minin,min,max, dx,dy, dfac;
    int  n,n2,n3, nin, nsample,nsamplemax, matchflag, lng,lat;

    /* Allocate memory for new samples if necessary */
    nsample = set->nsample;
    if (!set->nsamplemax)
        malloc_samples(set, nsamplein);
    else if ((nsamplemax=nsample+nsamplein) > set->nsamplemax)
        realloc_samples(set, (int)nsamplemax);
    radius2 = radius*radius;
    lng = set->lng;
    lat = set->lat;
    sort_coord = lat<0? set->naxis-1 : lat;
    switch(mode)
    {
        case UNION_RAW:
            sample = set->sample;
            qsort(sample, nsample, sizeof(samplestruct), compraw_samples);
            qsort(samplein, nsamplein, sizeof(samplestruct), compraw_samples);
            /*---- Find the starting element */
            sample3 = sample + nsample;
            n = n3 = nsample;
            minin = sample->rawpos[sort_coord] - radius;
            /*---- Go! */
            for (nin=nsamplein; nin--; samplein++)
            {
                matchflag = 0;
                if (samplein->rawpos[sort_coord] >= minin && n)
                {
                    min = samplein->rawpos[sort_coord] - radius;
                    max = min + 2.0*radius;
                    for (; n && (sample++)->rawpos[sort_coord] < min; n--);
                    n2 = n;
                    for (sample2=sample; n2-- && sample2->rawpos[sort_coord]<max;
                            sample2++)
                    {
                        dx = sample2->rawpos[lng] - samplein->rawpos[lng];
                        dy = sample2->rawpos[lat] - samplein->rawpos[lat];
                        if (dx*dx+dy*dy < radius2)
                        {
                            matchflag++;
                            break;
                        }
                    }
                }
                if (!matchflag)
                {
                    context = sample3->context;
                    if (set->ncontext)
                        memcpy(sample3->context, samplein->context,
                            set->ncontext*sizeof(double));
                    *sample3 = *samplein;
                    sample3->context = context;
                    sample3++;
                    n3++;
                }
            }
            break;
        case UNION_PROJ:
            sample = set->sample;
            qsort(sample, nsample, sizeof(samplestruct), compproj_samples);
            qsort(samplein, nsamplein, sizeof(samplestruct), compproj_samples);
            /*---- Find the starting element */
            sample3 = sample + nsample;
            n = n3 = nsample;
            minin = sample->projpos[sort_coord] - radius;
            /*---- Go! */
            for (nin=nsamplein; nin--; samplein++)
            {
                matchflag = 0;
                if (samplein->projpos[sort_coord] >= minin && n)
                {
                    min = samplein->projpos[sort_coord] - radius;
                    max = min + 2.0*radius;
                    for (; n && (sample++)->projpos[sort_coord] < min; n--);
                    n2 = n;
                    for (sample2=sample; n2-- && sample2->projpos[sort_coord]<max;
                            sample2++)
                    {
                        dx = sample2->projpos[lng] - samplein->projpos[lng];
                        dy = sample2->projpos[lat] - samplein->projpos[lat];
                        if (dx*dx+dy*dy < radius2)
                        {
                            matchflag++;
                            break;
                        }
                    }
                }
                if (!matchflag)
                {
                    context = sample3->context;
                    if (set->ncontext)
                        memcpy(sample3->context, samplein->context,
                            set->ncontext*sizeof(double));
                    *sample3 = *samplein;
                    sample3->context = context;
                    sample3++;
                    n3++;
                }
            }
            break;
        case UNION_WCS:
            sample = set->sample;
            qsort(sample, nsample, sizeof(samplestruct), compwcs_samples);
            qsort(samplein, nsamplein, sizeof(samplestruct), compwcs_samples);
            /*---- Find the starting element */
            sample3 = sample + nsample;
            n = n3 = nsample;
            minin = sample->wcspos[sort_coord] - radius;
            /*---- Go! */
            for (nin=nsamplein; nin--; samplein++)
            {
                matchflag = 0;
                if (samplein->wcspos[sort_coord] >= minin && n)
                {
                    min = samplein->wcspos[sort_coord] - radius;
                    max = min + 2.0*radius;
                    for (; n && (sample++)->wcspos[sort_coord] < min; n--);
                    dfac = (lng!=lat)? 1.0: cos(samplein->wcspos[lat]*DEG);
                    n2 = n;
                    for (sample2=sample; n2-- && sample2->wcspos[sort_coord]<max;
                            sample2++)
                    {
                        dx = (sample2->wcspos[lng] - samplein->wcspos[lng])*dfac;
                        dy = sample2->wcspos[lat] - samplein->wcspos[lat];
                        if (dx*dx+dy*dy < radius2)
                        {
                            matchflag++;
                            break;
                        }
                    }
                }
                if (!matchflag)
                {
                    context = sample3->context;
                    if (set->ncontext)
                        memcpy(sample3->context, samplein->context,
                            set->ncontext*sizeof(double));
                    *sample3 = *samplein;
                    sample3->context = context;
                    sample3->set = set;
                    sample3++;
                    n3++;
                }
            }
            break;
        default:
            n3 = 0; /* To avoid gcc -Wall warnings */
            error(EXIT_FAILURE, "*Internal Error*: Unknown sample union mode in ",
                    "union_samples()");
            break;
    }

    set->nsample = n3;

    return;
}


/****** copy_samples **********************************************************
  PROTO   void copy_samples(samplestruct *samplein, setstruct *set,
  int nsample)
  PURPOSE Copy a series of samples to a set.
  INPUT   pointer to 1st sample,
  pointer to output set,
  number of samples.
  OUTPUT  -.
  NOTES   Memory for the new samples is reallocated if needed.
  AUTHOR  E. Bertin (IAP)
  VERSION 27/04/2020
 */
void copy_samples(samplestruct *samplein, setstruct *set,
        int nsample)

{
    samplestruct *sample;
    double *context;
    int  i, nsamplemax;

    /* Allocate memory for new samples if necessary */
    if (!set->nsamplemax)
        malloc_samples(set, 3);
    if ((nsamplemax=set->nsample+nsample) > set->nsamplemax)
        realloc_samples(set, (int)(nsamplemax*1.3));
    sample = set->sample+set->nsample;
    for (i=0; i<nsample; i++)
    {
        context = sample->context;
        if (set->ncontext)
            memcpy(sample->context, samplein->context, set->ncontext*sizeof(double));
        *sample = *(samplein++);
        sample->context = context;
        sample++;
    }
    set->nsample += nsample;

    return;
}


/****** update_samples *******************************************************
  PROTO   void update_samples(setstruct *set, double *radius)
  PURPOSE Update "effective" positions and position uncertainties.
  INPUT   set structure pointer,
  field (or set) radius (for computing atm. turbulence contribution).
  OUTPUT  -.
  NOTES   Global preferences are used.
  AUTHOR  E. Bertin (IAP)
  VERSION 10/10/2012
 */
void update_samples(setstruct *set, double radius)

{
    samplestruct *samp;
    double dwcsastacc2;
    float wcsscale[NAXIS],wcsscale2[NAXIS], wcsastacc2[NAXIS];
    int  e,i, naxis;

    naxis = set->naxis = set->wcs->naxis;

    for (i=0; i<naxis; i++)
    {
        wcsscale[i] = set->wcsscale[i];
        wcsscale2[i] = wcsscale[i]*wcsscale[i];
    }

    /* Update virtual raw coordinates */
    samp = set->sample;
    for (e=set->nsample; e--; samp++)
        for (i=0; i<naxis; i++)
            samp->vrawpos[i] = samp->rawpos[i];

    if (set->astraccuracy > 0.0)
    {
        /*-- Update astrometric accuracy floor */
        if (prefs.astraccuracy_type==ASTACC_SIGMA_PIXEL)
            /*---- accuracy floor in pixels */
            for (i=0; i<naxis; i++)
                wcsastacc2[i] = set->astraccuracy*set->astraccuracy*wcsscale2[i];
        else if (prefs.astraccuracy_type==ASTACC_SIGMA_ARCSEC)
            /*---- accuracy floor in arcsec */
            for (i=0; i<naxis; i++)
                wcsastacc2[i] = set->astraccuracy*set->astraccuracy;
        else
            /*---- Atmospheric turbulence scaled from seeing (e.g. Han & Gatewood 1995) */
        {
            /*---- 2/3 factor because integ of (dx2+dy2)^1/6 ~ ((w^2+h^2)^(1/2)/3)^(1/3) */
            dwcsastacc2 = set->astraccuracy * sqrt(1/2.0)
                * pow(radius*2.0/3.0*DEG/(10.0*ARCMIN), 1/3.0)
                / sqrt(set->expotime) * ARCSEC/DEG;
            dwcsastacc2 *= dwcsastacc2;
            for (i=0; i<naxis; i++)
                wcsastacc2[i] = dwcsastacc2;
        }

        /*-- Update sample positional errors */
        samp = set->sample;
#ifdef __INTEL_COMPILER
// Turn off auto-vectorization because of an INTEL bug.
#pragma novector
#endif
        for (e=set->nsample; e--; samp++)
            for (i=0; i<naxis; i++)
                samp->wcsposerr[i] = sqrtf(wcsastacc2[i]
                        + samp->rawposerr[i]*samp->rawposerr[i]*wcsscale2[i]);
    }
    else
    {
        /*-- Update sample positional errors */
        samp = set->sample;
        for (e=set->nsample; e--; samp++)
            for (i=0; i<naxis; i++)
                samp->wcsposerr[i] = samp->rawposerr[i]*wcsscale[i];
    }

    return;
}


/****** locate_set **********************************************************
  PROTO   void locate_set(setstruct *set)
  PURPOSE Find set center coordinates and radius
  INPUT   set structure pointer.
  OUTPUT  -.
  NOTES   Global preferences are used.
  AUTHOR  E. Bertin (IAP)
  VERSION 29/09/2012
 */
void locate_set(setstruct *set)

{
    wcsstruct *wcs;
    double pixpos[NAXIS], wcspos[NAXIS],
    dwcsscale;
    int  i, naxis, lng,lat;

    wcs = set->wcs;
    dwcsscale = 0.0; /* to avoid gcc -Wall warnings */
    lng = set->lng = wcs->lng;
    lat = set->lat = wcs->lat;
    naxis = set->naxis = wcs->naxis;

    /* Find center of set*/
    for (i=0; i<naxis; i++)
        pixpos[i] = (wcs->naxisn[i]+1.0)/2.0;
    raw_to_wcs(wcs, pixpos, set->wcspos);

    /* Find center pixel scale */
    if (lat != lng)
        dwcsscale = sqrt(wcs_scale(wcs, pixpos));
    for (i=0; i<naxis; i++)
    {
        if (lat==lng || (i!=lng && i!=lat))
            set->wcsscale[i] = fabs(wcs->cd[i*(naxis+1)]);
        else
            set->wcsscale[i] = dwcsscale;
    }

    /* Find set "radius" */
    for (i=0; i<naxis; i++)
        pixpos[i] = 0.0;
    raw_to_wcs(wcs, pixpos, wcspos);
    set->radius = wcs_dist(wcs, set->wcspos, wcspos);

    return;
}


