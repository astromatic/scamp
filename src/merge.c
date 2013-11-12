/*
*				merge.c
*
* merge detections into sources.
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
*	Last modified:		12/11/2013
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
#include "fgroup.h"
#include "field.h"
#include "merge.h"
#include "prefs.h"
#include "samples.h"

#ifdef USE_THREADS
#include "threads.h"
#endif


/****** merge_fgroup ********************************************************
PROTO	msamplestruct *merge_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
PURPOSE	Create a global, merged catalogue.
INPUT	Pointer to the fgroup structure,
	Pointer to the reference field.
OUTPUT  Pointer to an allocated array of merged samples (sources).
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 12/11/2013
*/
msamplestruct	*merge_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)

  {
   fieldstruct		*field;
   setstruct		*set;
   msamplestruct	*msamp;
   samplestruct		*samp,*samp2;
   double		wcspos[NAXIS], wcsposerr[NAXIS], wcsposdisp[NAXIS],
			wcsposref[NAXIS],
			epoch,epochmin,epochmax, err2, spread, wspread,
			weight,weights, dummy;
   long			dptr;
   short		sexflagmask;
   unsigned int		imaflagmask;
   int			d,f,i,k,n,p,s, nall,nphotok,nposok, npinstru, naxis,
			nmsample, index, refflag;

  sexflagmask = (short)prefs.astr_sexflagsmask;
  imaflagmask = prefs.astr_imaflagsmask;
  naxis = fgroup->naxis;
  npinstru = prefs.nphotinstrustr;
  refflag = prefs.astrefinprop_flag;

/* Clear any existing msample array */
  free(fgroup->msample);

/* Clear out msamp pointers in reference samples */
  for (s=0; s<reffield->nset; s++)
    {
    set = reffield->set[s];
    samp = set->sample;
    for (n=set->nsample; n--; samp++)
      samp->msamp = NULL;
    }

/* Count the total number of sources */
  nmsample=0;
  for (f=0; f<fgroup->nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      samp = set->sample;
      for (n=set->nsample; n--; samp++)
        {
        samp->msamp = NULL;
        if (!samp->nextsamp)
          nmsample++;
        }
      }
    }


  fgroup->nmsample = nmsample;
  QCALLOC(fgroup->msample, msamplestruct, fgroup->nmsample);
  msamp = fgroup->msample;

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
          msamp->sourceindex = ++index;
          msamp->samp = samp;

/*-------- Astrometry */
          nall = nposok = 0;
          for (d=0; d<naxis; d++)
            wcspos[d] = wcsposerr[d] = wcsposdisp[d] = wcsposref[d] = 0.0;
          weights = 0.0;
          epoch = spread = wspread = 0.0;
          epochmin = BIG;
          epochmax = -BIG;
          for (samp2 = samp;
		samp2 && ((p=samp2->set->field->photomlabel)>=0 || refflag);
                samp2=samp2->prevsamp)
            {
            nall++;
            if ((samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask)
		|| (samp2->scampflags & SCAMP_BADPROPER))
              continue;
            for (d=0; d<naxis; d++)
              {
              err2 = samp2->wcsposerr[d]*samp2->wcsposerr[d];
              wcspos[d] += samp2->wcspos[d];
              if (err2 <= 0.0)
                err2 += 1.0;
              weight = 1.0/err2;
              weights += weight;
              epoch += weight*samp2->epoch;
              wcsposerr[d] += weight;
              wcspos[d] += weight*samp2->wcspos[d];
              if (!nposok)
                wcsposref[d] = samp2->wcspos[d];
              wcsposdisp[d] += (samp2->wcspos[d] - wcsposref[d])
				* (samp2->wcspos[d] - wcsposref[d]);
              }

/*---------- Epochs */
            if (p>=0)
              {
              if (samp2->set->epochmin < epochmin)
                epochmin = samp2->set->epochmin;
              if (samp2->set->epochmax > epochmax)
                epochmax = samp2->set->epochmax;
              }
            else	/* Special treatment for astrometric ref. catalog */
              {
              if (samp2->epoch < epochmin)
                epochmin = samp2->epoch;
              if (samp2->epoch > epochmax)
                epochmax = samp2->epoch;
              }
/*---------- Morphometry */
            if (prefs.spread_flag && p>=0)	/* Exclude ref. from stats */
              {
              weight = samp2->spreaderr>TINY?
		1.0/(samp2->spreaderr*samp2->spreaderr) : 1.0;
              spread += weight*samp2->spread;
              wspread += weight;
              }
/*---------- Flags */
            msamp->sexflags |= samp2->sexflags;
            msamp->scampflags |= samp2->scampflags;
            msamp->imaflags |= samp2->imaflags;
            nposok++;
            }
          if (nposok)
            {
            for (d=0; d<naxis; d++)
              {
              msamp->wcspos[d] = wcspos[d] / wcsposerr[d];
              msamp->wcsposerr[d] = sqrt(1.0/wcsposerr[d]);
              msamp->wcsposdisp[d] = nposok > 1?
		sqrt(fabs(wcsposdisp[d]
		 - nposok*(msamp->wcspos[d] - wcsposref[d])
			*(msamp->wcspos[d] - wcsposref[d]))/(nposok-1.0))
		: 0.0;
              }
            if (msamp->wcsposerr[0] < msamp->wcsposerr[1])
              {
              dummy = msamp->wcsposerr[0];
              msamp->wcsposerr[0] = msamp->wcsposerr[1];
              msamp->wcsposerr[1] = dummy;
              msamp->wcspostheta = 90.0;
              }
            else
              msamp->wcspostheta = 0.0;
            msamp->epoch = epoch / weights;
            msamp->epochmin = epochmin;
            msamp->epochmax = epochmax;
            if (prefs.spread_flag)
              {
              msamp->spread = spread / wspread;
              msamp->spreaderr = sqrt(1.0 / wspread);
              }
            }
          msamp->npos_tot = nall;
          msamp->npos_ok = nposok;
          for (samp2 = samp; samp2; samp2=samp2->prevsamp)
            samp2->msamp = msamp;
          msamp++;
          }
      }
    }

  return fgroup->msample;
  }



