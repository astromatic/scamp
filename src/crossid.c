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

