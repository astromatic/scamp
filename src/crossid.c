/*
 *    crossid.c
 *
 * Manage source cross-identifications.
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 * This file part of: SCAMP
 *
 * Copyright:  (C) 2002-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
 * Last modified:  19/02/2018
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "globals.h"
#include "crossid.h"
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "match.h"
#include "misc.h"
#include "prefs.h"
#include "samples.h"

/****** crossid_fgroup *******************************************************
  PROTO void crossid_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
  double tolerance)
  PURPOSE Perform source cross-identifications in a group of fields.
  INPUT ptr to the group of fields,
  ptr to the reference field,
  Tolerance (in deg. if angular coordinates).
  OUTPUT -.
  NOTES Uses the global preferences.
  AUTHOR E. Bertin (IAP)
  VERSION 19/02/2018
 ***/
void crossid_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
        double tolerance)
{
    fieldstruct **field, *field1, *field2;
    wcsstruct *wcs;
    setstruct **pset1,**pset2, **pset,
              *set1,*set2, *set;
    samplestruct *samp1,*nexsamp1, *samp2,*samp2b,*samp2min,
                 *prevsamp2,*nextsamp1,*nextsamp2;
    double projmin2[NAXIS], projmax2[NAXIS],
    *proj1,
    lng1,lat1, latmin1,latmax1, lngmin2,lngmax2,latmin2,latmax2,
    dlng,dlat, dx, rlim,rlimmin,r2,r2n,r2p,r2min;
    float fmax;
    int  i, f,f1,f2, s1,s2, nset1,nset2, nsamp, nsamp2,nsamp2b,
         s, nfield, naxis, lng,lat, yaxis;

    field1 = NULL; /* to avoid gcc -Wall warnings */
    proj1 = NULL; /* to avoid gcc -Wall warnings */
    lng1 = lngmin2 = lngmax2 = latmin2 = latmax2 = 0.0;
    field = fgroup->field;
    nfield = fgroup->nfield;
    naxis = fgroup->naxis;
    lng = fgroup->lng;
    lat = fgroup->lat;
    wcs = fgroup->wcs;

    /* Compute the largest possible error in pixels allowed in previous matching */
    rlimmin = 0.0;
    for (i=0; i<naxis; i++)
        if ((rlim=tolerance/fgroup->meanwcsscale[i])>rlimmin)
            rlimmin = rlim;
    rlim = rlimmin;

    /* Sort samples to accelerate further processing and reset pointers */
    for (f=0; f<nfield; f++)
    {
        pset = field[f]->set;
        field[f]->prevfield = field[f]->nextfield = NULL;
        for (s=field[f]->nset; s--;)
        {
            set = *(pset++);
            //sort_samples(set);
            unlink_samples(set);
        }
    }

    /* Now start the real cross-id loop */
    for (f1=1; f1<nfield; f1++)
    {
        field1 = field[f1];
        pset1 = field1->set;
        nset1 = field1->nset;
        for (s1=nset1; s1--;)
        {
            set1 = *(pset1++);
            for (f2=0; f2<f1; f2++)
            {
                field2 = field[f2];
                pset2 = field2->set;
                nset2 = field2->nset;
                for (s2=nset2; s2--;)
                {
                    set2 = *(pset2++);
                    /*-------- Exclude non-overlapping frames */
                    if (lng != lat)
                    {
                        if (set1->projposmin[lng] > (lngmax2=set2->projposmax[lng]+rlim)
                                || (lngmin2=set2->projposmin[lng]-rlim) > set1->projposmax[lng]
                                || set1->projposmin[lat] > (latmax2=set2->projposmax[lat]+rlim)
                                || (latmin2=set2->projposmin[lat]-rlim)> set1->projposmax[lat])
                            continue;
                    }
                    else
                        for (i=0; i<naxis; i++)
                            if (set1->projposmin[i] > (projmax2[i]=set2->projposmax[i]+rlim)
                                    || (projmin2[i]=set2->projposmin[i]-rlim)> set1->projposmax[i])
                                continue;
                    samp1 = set1->sample;
                    samp2b = set2->sample;
                    nsamp2b = set2->nsample;
                    for (nsamp=set1->nsample; nsamp--; samp1++)
                    {
                        if (lat!=lng)
                        {
                            lng1 = samp1->projpos[lng];
                            lat1 = samp1->projpos[lat];
                            yaxis = lat;
                            /*------------ Jump over sources in the non-overlapping region */
                            if (lat1<latmin2 || lat1>latmax2 || lng1<lngmin2 || lng1>lngmax2)
                                continue;
                        }
                        else
                        {
                            proj1 = samp1->projpos;
                            for (i=0; i<naxis; i++)
                                if (proj1[i] < projmin2[i] || proj1[i]>projmax2[i])
                                    continue;
                            lat1 = (naxis<2) ? proj1[yaxis=0] : proj1[yaxis=1];
                        }
                        latmin1 = lat1-rlim;
                        latmax1 = lat1+rlim;
                        r2min = rlim*rlim;
                        samp2min = NULL;
                        samp2 = samp2b;
                        /*---------- Jump over sources that can't match in y */
                        for (nsamp2=nsamp2b; nsamp2-- && samp2->projpos[yaxis]<latmin1;
                                samp2++);
                        samp2b = samp2;
                        nsamp2b = ++nsamp2;
                        for (; nsamp2-- && samp2->projpos[yaxis]<latmax1; samp2++)
                        {
                            if (samp2->nextsamp && samp2->nextsamp->set->field!=field1)
                                continue;
                            if (lat!=lng)
                            {
                                dlng = lng1 - samp2->projpos[lng];
                                dlat = lat1 - samp2->projpos[lat];
                                r2 = dlng*dlng + dlat*dlat;
                            }
                            else
                            {
                                r2 = 0.0;
                                for (i=0; i<naxis; i++)
                                {
                                    dx = proj1[i] - samp2->projpos[i];
                                    r2 += dx*dx;
                                }
                            }
                            /*------------ Finally select the closest source within the search disk */
                            if (r2<r2min)
                            {
                                r2min = r2;
                                samp2min = samp2;
                            }
                        }
                        if (samp2min)
                        {
                            r2p = r2n = BIG;
                            if ((prevsamp2=samp1->prevsamp))
                            {
                                /*-------------- Check if it is a better match than the previous one */
                                if (lat!=lng)
                                {
                                    dlng = prevsamp2->projpos[lng] - samp1->projpos[lng];
                                    dlat = prevsamp2->projpos[lat] - samp1->projpos[lat];
                                    r2p = dlng*dlng + dlat*dlat;
                                }
                                else
                                {
                                    r2p = 0.0;
                                    for (i=0; i<naxis; i++)
                                    {
                                        dx = prevsamp2->projpos[i] - samp1->projpos[i];
                                        r2p += dx*dx;
                                    }
                                }
                            }
                            if ((nextsamp1=samp2min->nextsamp))
                            {
                                /*-------------- Check if it is a better match than the previous one */
                                if (lat!=lng)
                                {
                                    dlng = samp2min->projpos[lng] - nextsamp1->projpos[lng];
                                    dlat = samp2min->projpos[lat] - nextsamp1->projpos[lat];
                                    r2n = dlng*dlng + dlat*dlat;
                                }
                                else
                                {
                                    r2n = 0.0;
                                    for (i=0; i<naxis; i++)
                                    {
                                        dx = samp2min->projpos[i] - nextsamp1->projpos[i];
                                        r2n += dx*dx;
                                    }
                                }
                            }
                            if (r2min<r2p && r2min<r2n)
                                /*------------ unlink from previous match if this is a better match */
                            {
                                if (prevsamp2)
                                    prevsamp2->nextsamp = NULL;
                                if (nextsamp1)
                                    nextsamp1->prevsamp = NULL;
                                samp1->prevsamp = samp2min;
                                samp2min->nextsamp = samp1;
                            }
                        }
                    }
                }
            }
        }
    }

    /* Now bring also the reference field samples to the common projection */
    /* Sort samples to accelerate further processing and reset pointers */
    if (reffield)
    {
        set1 = reffield->set[0];
        sort_samples(set1);
        unlink_samples(set1);

        for (f2=0; f2<nfield; f2++)
        {
            field2 = field[f2];
            pset2 = field2->set;
            nset2 = field2->nset;
            for (s2=nset2; s2--;)
            {
                set2 = *(pset2++);
                /*---------- Exclude non-overlapping frames */
                if (lng != lat)
                {
                    if (set1->projposmin[lng] > (lngmax2=set2->projposmax[lng]+rlim)
                            || (lngmin2=set2->projposmin[lng]-rlim) > set1->projposmax[lng]
                            || set1->projposmin[lat] > (latmax2=set2->projposmax[lat]+rlim)
                            || (latmin2=set2->projposmin[lat]-rlim)> set1->projposmax[lat])
                        continue;
                }
                else
                    for (i=0; i<naxis; i++)
                        if (set1->projposmin[i] > (projmax2[i]=set2->projposmax[i]+rlim)
                                || (projmin2[i]=set2->projposmin[i]-rlim) >set1->projposmax[i])
                            continue;
                samp1 = set1->sample;
                samp2b = set2->sample;
                nsamp2b = set2->nsample;
                for (nsamp=set1->nsample; nsamp--; samp1++)
                {
                    if (lat!=lng)
                    {
                        lng1 = samp1->projpos[lng];
                        lat1 = samp1->projpos[lat];
                        yaxis = lat;
                        /*---------- Jump over sources in the non-overlapping region */
                        if (lat1<latmin2 || lat1>latmax2 || lng1<lngmin2 || lng1>lngmax2)
                            continue;
                    }
                    else
                    {
                        proj1 = samp1->projpos;
                        for (i=0; i<naxis; i++)
                            if (proj1[i] < projmin2[i] || proj1[i]>projmax2[i])
                                continue;
                        lat1 = (naxis<2) ? proj1[yaxis=0] : proj1[yaxis=1];
                    }
                    latmin1 = lat1-rlim;
                    latmax1 = lat1+rlim;
                    r2min = rlim*rlim;
                    fmax = 0.0;
                    samp2min = NULL;
                    samp2 = samp2b;
                    /*-------- Jump over sources that can't match in y */
                    for (nsamp2=nsamp2b; nsamp2-- && samp2->projpos[yaxis]<latmin1;
                            samp2++);
                    samp2b = samp2;
                    nsamp2b = ++nsamp2;
                    for (; nsamp2-- && samp2->projpos[yaxis]<latmax1; samp2++)
                    {
                        if (samp2->prevsamp)
                            continue;
                        if (lat!=lng)
                        {
                            dlng = lng1 - samp2->projpos[lng];
                            dlat = lat1 - samp2->projpos[lat];
                            r2 = dlng*dlng + dlat*dlat;
                        }
                        else
                        {
                            r2 = 0.0;
                            for (i=0; i<naxis; i++)
                            {
                                dx = proj1[i] - samp2->projpos[i];
                                r2 += dx*dx;
                            }
                        }
                        /*---------- Finally select the closest source within the search disk */
                        if (r2<r2min)
                        {
                            r2min = r2;
                            samp2min = samp2;
                        }
                    }
                    if (samp2min)
                    {
                        r2n = BIG;
                        if ((nextsamp2=samp1->nextsamp))
                        {
                            /*------------ Check if it is a better match than the previous one */
                            if (lat!=lng)
                            {
                                dlng = nextsamp2->projpos[lng] - samp1->projpos[lng];
                                dlat = nextsamp2->projpos[lat] - samp1->projpos[lat];
                                r2n = dlng*dlng + dlat*dlat;
                            }
                            else
                            {
                                r2n = 0.0;
                                for (i=0; i<naxis; i++)
                                {
                                    dx = nextsamp2->projpos[i] - samp1->projpos[i];
                                    r2n += dx*dx;
                                }
                            }
                        }
                        if (r2min<r2n)
                            /*---------- unlink from previous match if this is a better match */
                        {
                            if (nextsamp2)
                                nextsamp2->prevsamp = NULL;
                            samp1->nextsamp = samp2min;
                            samp2min->prevsamp = samp1;
                        }
                    }
                }
            }
        }
    }

    return;
}


/****** recenter_fgroup *******************************************************
  PROTO void recenter_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
  PURPOSE Perform field recentering with respect to a reference catalog in a
  group of fields.
  INPUT ptr to the group of fields,
  ptr to the reference field.
  OUTPUT -.
  NOTES Uses the global preferences.
  AUTHOR E. Bertin (IAP)
  VERSION 09/06/2011
 ***/
void recenter_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
{
    fieldstruct *field;
    setstruct **sets, *set;
    samplestruct *samp, *samp2;
    double *offsetbuf[NAXIS],
    offset[NAXIS], rawpos[NAXIS], wcspos[NAXIS], dwcspos[NAXIS];  
    int  d,f,s, naxis, nsamp, o,omax;

    NFPRINTF(OUTPUT, "Re-centering fields...");

    set = reffield->set[0];
    naxis = fgroup->naxis;
    omax = 0;
    for (f=0; f<fgroup->nfield; f++)
    {
        o = 0;
        field = fgroup->field[f];
        samp = set->sample;
        for (nsamp=set->nsample; nsamp--; samp++)
        {
            for (samp2=samp; (samp2=samp2->nextsamp);)
            {
                if (samp2->set && samp2->set->field == field)
                {
                    if (o>=omax)
                    {
                        omax += 8192;
                        if (o)
                            for (d=0; d<naxis; d++)
                            {
                                QREALLOC(offsetbuf[d], double, omax);
                            }
                        else
                            for (d=0; d<naxis; d++)
                            {
                                QMALLOC(offsetbuf[d], double, omax);
                            }
                    }
                    for (d=0; d<naxis; d++)
                        offsetbuf[d][o] = samp2->projpos[d] - samp->projpos[d];
                    o++;
                }
            }
        }
        /*-- Compute the median reprojected shift in each dimension */
        for (d=0; d<naxis; d++)
            offset[d] = fast_median(offsetbuf[d], o);
        /*-- Convert it to a shift in world coordinates */
        for (d=0; d<naxis; d++)
            rawpos[d] = field->set[0]->wcs->crpix[d] - offset[d];
        raw_to_wcs(field->set[0]->wcs, rawpos, wcspos);
        for (d=0; d<naxis; d++)
            dwcspos[d] = wcspos[d] - field->set[0]->wcs->crval[d];
        sets = field->set;
        for (s=0; s<field->nset; s++)
            update_wcsll(sets[s]->wcs, dwcspos[set->lng], dwcspos[set->lat]);
    }

    for (d=0; d<naxis; d++)
        free(offsetbuf[d]);

    return;
}


/****** check_fieldoverlap ****************************************************
  PROTO int check_fieldoverlap(fieldstruct *field1, fieldstruct *field2)
  PURPOSE Check if two fields overlap or not.
  INPUT ptr to the first field,
  ptr to the second field.
  OUTPUT 1 if they overlap, 0 otherwise.
  NOTES -.
  AUTHOR E. Bertin (IAP)
  VERSION 07/02/2005
 ***/
int check_fieldoverlap(fieldstruct *field1, fieldstruct *field2)

{
    setstruct **pset,
              *set;
    samplestruct *samp,*samp2;
    int  n,s;

    pset = field1->set;
    set = *(pset++);
    for (s=field1->nset; s--; set=*(pset++))
    {
        samp = set->sample;
        for (n=set->nsample; n--; samp++)
        {
            if (samp->nextsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->nextsamp))
                    if (samp2->set->field == field2)
                        return 1;
            }
            if (samp->prevsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->prevsamp))
                    if (samp2->set->field == field2)
                        return 1;
            }
        }
    }

    /* No link found between both fields */
    return 0;
}


/****** check_fieldphotomoverlap **********************************************
  PROTO int check_fieldphotomoverlap(fieldstruct *field, int instru)
  PURPOSE Check if a field overlaps a photometric field or not.
  INPUT ptr to the field to check,
  photometric instrument index.
  OUTPUT Photometric code (1 for genuine, 2 for dummy) if it overlaps, 0
  otherwise.
  NOTES -.
  AUTHOR E. Bertin (IAP)
  VERSION 19/02/2018
 ***/
int check_fieldphotomoverlap(fieldstruct *field, int instru)

{
    setstruct **pset,
              *set;
    samplestruct *samp,*samp2;
    int  n,s;

    pset = field->set;
    for (s=field->nset; s--;)
    {
        set = *(pset++);
        samp = set->sample;
        for (n=set->nsample; n--; samp++)
        {
            if (samp->nextsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->nextsamp))
                    if (samp2->set->field->photomflag
                            && samp2->set->field->photomlabel==instru)
                        return samp2->set->field->photomflag;
            }
            if (samp->prevsamp)
            {
                samp2 = samp;
                while ((samp2=samp2->prevsamp))
                    if (samp2->set->field->photomflag
                            && samp2->set->field->photomlabel==instru)
                        return samp2->set->field->photomflag;
            }
        }
    }

    /* No photometric field found */
    return 0;
}

