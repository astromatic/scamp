/*
*				astrstats.c
*
* Compute astrometric statistics.
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

#ifdef USE_THREADS
#include <pthread.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "astrstats.h"
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "prefs.h"
#include "samples.h"


/****** astrstats_fgroup *****************************************************
PROTO	void astrstats_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
			double hsn_thresh)
PURPOSE	Compute stats of a global astrometric solution among a group of fields.
INPUT	ptr to a group of fields pointers,
	ptr to a reference catalog,
	S/N threshold for the high-S/N sample.
OUTPUT	-.
NOTES	Input structures must have gone through crossid_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	12/11/2013
 ***/
void	astrstats_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
			double hsn_thresh)
  {
   fieldstruct	**fields,
		*field,*field2,*field3;
   setstruct	*set;
   samplestruct	*samp,*samp2,*samp3;
   long double	lsig[NAXIS], lsig_hsn[NAXIS],
		gcorr_int, gcorr_int_hsn, lchi2, lchi2_hsn;
   double	mean[NAXIS],mean_hsn[NAXIS], chi2[NAXIS],
		goffset_ref[NAXIS], goffset_ref_hsn[NAXIS],
		gcorr_ref,gcorr_ref_hsn, gcorr_sub,gcorr_sub_hsn,
		corr,corr_hsn,
		dx, sn2;
   long long	ndeg, ndeg_hsn;
   short	sexflagmask;
   unsigned int	imaflagmask;
   int		i,f,n,s, naxis,nfield,nsamp, nsource,nsource_hsn,
		nmean,nmean_hsn, nmatch,nmatch_hsn;

  sexflagmask = (short)prefs.astr_sexflagsmask;
  imaflagmask = prefs.astr_imaflagsmask;
  naxis = fgroup->naxis;
  nfield = fgroup->nfield;
  fields = fgroup->field;

/* Reset all counters and integrators */
  for (f=0; f<nfield; f++)
    {
    field = fields[f];
    field->chi2_int = field->chi2_int_hsn
	= field->chi2_ref = field->chi2_ref_hsn = 0.0;
    field->nchi2_int = field->nchi2_int_hsn
	= field->nchi2_ref = field->nchi2_ref_hsn
	= field->sig_corr_ref = field->sig_corr_ref_hsn = 0.0;
    for (i=0; i<naxis; i++)
      field->offset_ref[i] = field->offset_ref_hsn[i]
	= field->sig_referr[i] = field->sig_referr_hsn[i] = 0.0;
    }
  for (i=0; i<naxis; i++)
    lsig[i] = lsig_hsn[i] = 0.0;
  lchi2 = lchi2_hsn = 0.0;
  sn2 = 0.0;			/* to avoid gcc -Wall warnings */
  nsource = nsource_hsn = nmatch = nmatch_hsn = 0;
  ndeg = ndeg_hsn = 0;
  gcorr_int = gcorr_int_hsn = 0.0;

/* Internal stats */
  for (f=0; f<nfield; f++)
    {
    field=fields[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      samp = set->sample;
      for (n=nsamp; n--; samp++)
        if (!samp->nextsamp && samp->prevsamp)
	  {
          for (i=0; i<naxis; i++)
            mean[i] = mean_hsn[i] = 0.0;
          nmean = nmean_hsn = 0;
          for (samp2 = samp; samp2 && samp2->set->field->astromlabel>=0;
		samp2=samp2->prevsamp)
	    {
            if ((samp2->sexflags & sexflagmask)
            	|| (samp2->imaflags & imaflagmask))
              continue;
            for (i=0; i<naxis; i++)
              mean[i] += samp2->projpos[i];
            nmean++;
            sn2 = samp2->fluxerr>0.0? samp2->flux/samp2->fluxerr: 0.0;
            if (sn2 >= hsn_thresh)
	      {
              for (i=0; i<naxis; i++)
                mean_hsn[i] += samp2->projpos[i];
              nmean_hsn++;
              }
            }
          if (nmean>1)
	    {
            nsource++;
            for (i=0; i<naxis; i++)
              mean[i] /= (double)nmean;
            if (nmean_hsn>1)
              {
              nsource_hsn++;
              for (i=0; i<naxis; i++)
                mean_hsn[i] /= (double)nmean_hsn;
              }
            for (samp2 = samp; samp2 && samp2->set->field->astromlabel>=0;
		samp2=samp2->prevsamp)
	      {
              if ((samp2->sexflags & sexflagmask)
			|| (samp2->imaflags & imaflagmask))
                continue;
              corr = corr_hsn = 1.0;
              for (i=0; i<naxis; i++)
                {
                dx = samp2->projpos[i] - mean[i];
                corr *= dx;
                lsig[i] += dx*dx / (nmean-1.0);                
                }
              gcorr_int += corr / (nmean - 1.0);
              if (nmean_hsn>1)
	        {
                sn2 = samp2->fluxerr>0.0? samp2->flux/samp2->fluxerr: 0.0;
                if (sn2 >= hsn_thresh)
                  {
                  for (i=0; i<naxis; i++)
                    {
                    dx = samp2->projpos[i] - mean_hsn[i];
                    corr_hsn *= dx;
                    lsig_hsn[i] += dx*dx / (nmean_hsn-1.0);
                    }
                  gcorr_int_hsn += corr / (nmean_hsn - 1.0);
                  }
                }
              field2 = samp2->set->field;
              for (samp3=samp2; (samp3=samp3->prevsamp)
			&& samp3->set->field->astromlabel>=0;)
                {
                if ((samp3->sexflags & sexflagmask)
			|| (samp3->imaflags & imaflagmask))
                  continue;
                field3 = samp3->set->field;
                for (i=0; i<naxis; i++)
                  {
                  dx = samp2->projpos[i] - samp3->projpos[i];
                  lchi2 += (chi2[i] = dx*dx
			*fgroup->meanwcsscale[i]*fgroup->meanwcsscale[i]
			/(samp2->wcsposerr[i]*samp2->wcsposerr[i]
			  + samp3->wcsposerr[i]*samp3->wcsposerr[i]));
                  field2->chi2_int += chi2[i];
                  field3->chi2_int += chi2[i];
                  }
                field2->nchi2_int++;
                field3->nchi2_int++;
                ndeg++;
                if (nmean_hsn>1 && sn2>hsn_thresh && samp3->fluxerr>0.0
			&& samp3->flux/samp3->fluxerr>hsn_thresh)
		  {
                  for (i=0; i<naxis; i++)
                    {
                    lchi2_hsn += chi2[i];
                    field2->chi2_int_hsn += chi2[i];
                    field3->chi2_int_hsn += chi2[i];
                    }
                  field2->nchi2_int_hsn++;
                  field3->nchi2_int_hsn++;
                  ndeg_hsn++;
                  }
                }
              }
            }
	  nmatch += nmean;
          nmatch_hsn += nmean_hsn;
          }
      }
    }

  for (i=0; i<naxis; i++)
    {
    fgroup->sig_interr[i] = nsource?
	sqrt(lsig[i]/nsource)*fgroup->meanwcsscale[i] : 0.0;
    fgroup->sig_interr_hsn[i] = nsource_hsn?
	sqrt(lsig_hsn[i]/nsource_hsn)*fgroup->meanwcsscale[i] : 0.0;
    lsig[i] = lsig_hsn[i] = 0.0;
    goffset_ref[i] = goffset_ref_hsn[i] = 0.0;
    if (fgroup->sig_interr[i] > TINY)
      gcorr_int *= fgroup->meanwcsscale[i] / fgroup->sig_interr[i];
    if (fgroup->sig_interr_hsn[i] > TINY)
      gcorr_int_hsn *= fgroup->meanwcsscale[i] / fgroup->sig_interr_hsn[i];
    }

  fgroup->sig_corr_int = nsource? gcorr_int / nsource : 0.0;
  fgroup->sig_corr_int_hsn = nsource_hsn? gcorr_int / nsource_hsn : 0.0;
  fgroup->chi2_int = ndeg? lchi2 / ((double)naxis*ndeg) : 0.0;
  fgroup->chi2_int_hsn = ndeg_hsn?
		lchi2_hsn / ((double)naxis*ndeg_hsn) : 0.0;
  fgroup->nintmatch = nmatch;
  fgroup->nintmatch_hsn = nmatch_hsn;

/* Now the stats with respect to the reference catalog */
  lchi2 = lchi2_hsn = 0.0;
  gcorr_ref = gcorr_ref_hsn = gcorr_sub = gcorr_sub_hsn = 0.0;
  nmatch = nmatch_hsn = 0;
  ndeg = ndeg_hsn = 0;
  set = reffield->set[0];
  nsamp = set->nsample;
  samp = set->sample;
  for (n=nsamp; n--; samp++)
    {
    if ((samp->sexflags & sexflagmask)
	|| (samp->imaflags & imaflagmask))
      continue;
    for (i=0; i<naxis; i++)
      mean[i] = mean_hsn[i] = 0.0;
    nmean = nmean_hsn = 0;
    for (samp2 = samp; (samp2=samp2->nextsamp);)
      {
      if ((samp2->sexflags & sexflagmask)
	|| (samp2->imaflags & imaflagmask))
        continue;
      field2 = samp2->set->field;
      sn2 = samp2->fluxerr>0.0? samp2->flux/samp2->fluxerr: 0.0;
      corr = corr_hsn = 1.0;
      for (i=0; i<naxis; i++)
        {
        mean[i] += samp2->projpos[i];
        dx = samp2->projpos[i] - samp->projpos[i];
        field2->offset_ref[i] += dx;
        field2->sig_referr[i] += dx*dx;
        corr *= dx;
        if (sn2 >= hsn_thresh)
          {
          field2->offset_ref_hsn[i] += dx;
          field2->sig_referr_hsn[i] += dx*dx;
          corr_hsn *= dx;
          }
        lchi2 += (chi2[i] = dx*dx
			*fgroup->meanwcsscale[i]*fgroup->meanwcsscale[i]
			/(samp2->wcsposerr[i]*samp2->wcsposerr[i]
			  + samp->wcsposerr[i]*samp->wcsposerr[i]));
        field2->chi2_ref += chi2[i];
        }
      field2->sig_corr_ref += corr;
      field2->nchi2_ref++;
      nmean++;
      if (sn2 >= hsn_thresh)
        {
        for (i=0; i<naxis; i++)
	  {
          mean_hsn[i] += samp2->projpos[i];
	  lchi2_hsn += chi2[i];
          field2->chi2_ref_hsn += chi2[i];
          }
        field2->sig_corr_ref_hsn += corr;
        field2->nchi2_ref_hsn++;
        nmean_hsn++;
        }
      }
    if (nmean)
      {
      if (nmean>1)
        {
        for (i=0; i<naxis; i++)
          mean[i] /= (double)nmean;
        if (nmean_hsn>1)
          for (i=0; i<naxis; i++)
            mean_hsn[i] /= (double)nmean_hsn;
        }
      corr = 1.0;
      for (i=0; i<naxis; i++)
        {
        dx = mean[i] - samp->projpos[i];
        goffset_ref[i] += dx;
        lsig[i] += dx * dx;
        corr *= dx;
        }
      gcorr_ref += corr;
      nmatch++;
      if (nmean_hsn)
        {
        corr_hsn = 1.0;
        for (i=0; i<naxis; i++)
          {
          dx = mean_hsn[i] - samp->projpos[i];
          goffset_ref_hsn[i] += dx;
          lsig_hsn[i] += dx * dx;
          corr_hsn *= dx;
          }
        gcorr_ref_hsn += corr_hsn;
        nmatch_hsn++;
        }
      ndeg += nmean;
      ndeg_hsn += nmean_hsn;
      }
    }

  for (i=0; i<naxis; i++)
    {
    if (nmatch)
      goffset_ref[i] /= (double)nmatch;
    fgroup->sig_referr[i] = nmatch>1?
	sqrt((lsig[i]-goffset_ref[i]*goffset_ref[i]/nmatch/nmatch)/(nmatch-1))
	*fgroup->meanwcsscale[i] : 0.0;
    gcorr_sub *= goffset_ref[i];
    fgroup->offset_ref[i] = goffset_ref[i]*fgroup->meanwcsscale[i];
    if (nmatch_hsn)
      goffset_ref_hsn[i] /= nmatch_hsn;
    fgroup->sig_referr_hsn[i] = nmatch_hsn>1?
	sqrt((lsig_hsn[i]-goffset_ref_hsn[i]*goffset_ref_hsn[i])/(nmatch_hsn-1))
	*fgroup->meanwcsscale[i] : 0.0;
    gcorr_sub_hsn *= goffset_ref_hsn[i];
    fgroup->offset_ref_hsn[i] = goffset_ref_hsn[i] * fgroup->meanwcsscale[i];
    }
  gcorr_ref -= gcorr_sub;
  gcorr_ref_hsn -= gcorr_sub_hsn;
  for (i=0; i<naxis; i++)
    {
    if (fgroup->sig_referr[i] > TINY)
      gcorr_ref *= fgroup->meanwcsscale[i] / fgroup->sig_referr[i];
    if (fgroup->sig_referr_hsn[i] > TINY)
      gcorr_ref_hsn *= fgroup->meanwcsscale[i] / fgroup->sig_referr_hsn[i];
    }

  fgroup->sig_corr_ref = nmatch? gcorr_ref / nmatch : 0.0;
  fgroup->sig_corr_ref_hsn = nmatch_hsn? gcorr_ref / nmatch_hsn : 0.0;
  fgroup->chi2_ref = ndeg? lchi2 / ((double)naxis*ndeg) : 0.0;
  fgroup->chi2_ref_hsn = ndeg_hsn?
			lchi2_hsn / ((double)naxis*ndeg_hsn) : 0.0;
  fgroup->nrefmatch = nmatch;
  fgroup->nrefmatch_hsn = nmatch_hsn;

  for (f=0; f<nfield; f++)
    {
    gcorr_sub = gcorr_sub_hsn = 0.0;
    field = fields[f];
    if (field->nchi2_int)
      field->chi2_int /= (long double)field->nchi2_int*naxis;
    if (field->nchi2_int_hsn)
      field->chi2_int_hsn /= (long double)field->nchi2_int_hsn*naxis;
    if (field->nchi2_ref)
      {
      field->chi2_ref /= (long double)field->nchi2_ref*naxis;
      for (i=0; i<naxis; i++)
        {
        field->offset_ref[i] /= (double)field->nchi2_ref;
        field->sig_referr[i] = field->nchi2_ref > 1?
		(long double)sqrt((field->sig_referr[i]
		- field->offset_ref[i]*field->offset_ref[i])
		/ (field->nchi2_ref-1)) * fgroup->meanwcsscale[i] : 0.0;
        gcorr_sub *= field->offset_ref[i];
        field->offset_ref[i] *= fgroup->meanwcsscale[i];
        }
      field->sig_corr_ref = (field->sig_corr_ref - gcorr_sub);
      for (i=0; i<naxis; i++)
        if (field->sig_referr[i] > TINY)
          field->sig_corr_ref *= fgroup->meanwcsscale[i]/field->sig_referr[i];
      if (field->nchi2_ref > 1)
        field->sig_corr_ref /= (long double)field->nchi2_ref;
      }
    if (field->nchi2_ref_hsn)
      {
      field->chi2_ref_hsn /= (long double)field->nchi2_ref_hsn*naxis;
      for (i=0; i<naxis; i++)
        {
        field->offset_ref_hsn[i] /= (double)field->nchi2_ref_hsn;
        field->sig_referr_hsn[i] = field->nchi2_ref_hsn > 1?
		(long double)sqrt((field->sig_referr_hsn[i]
		- field->offset_ref_hsn[i]*field->offset_ref_hsn[i])
		/ (field->nchi2_ref_hsn-1)) * fgroup->meanwcsscale[i] : 0.0;
        gcorr_sub_hsn *= field->offset_ref_hsn[i];
        field->offset_ref_hsn[i] *= fgroup->meanwcsscale[i];
        }
      field->sig_corr_ref_hsn = (field->sig_corr_ref_hsn - gcorr_sub_hsn);
      for (i=0; i<naxis; i++)
        if (field->sig_referr_hsn[i] > TINY)
          field->sig_corr_ref_hsn *= fgroup->meanwcsscale[i]
				/ field->sig_referr_hsn[i];
      if (field->nchi2_ref_hsn>1)
        field->sig_corr_ref_hsn /= (long double)field->nchi2_ref_hsn;
      }
    }

  return;
  }


/****** astrclip_fgroup *****************************************************
PROTO	void astrclip_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
			double nsigma)
PURPOSE	Removes outliers from cross-identification linked lists.
INPUT	ptr to a group of fields pointers,
	ptr to a reference catalog,
	threshold in number of sigma
OUTPUT	-.
NOTES	Input structures must have gone through crossid_fgroup() and
	astrstats_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	12/11/2013
 ***/
int	astrclip_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
				double nsigma)
  {
   fieldstruct	**fields,
		*field;
   setstruct	*set;
   samplestruct	*samp,*sampr,*samp2, *prevsamp2;
   double	clipi[NAXIS],clipr[NAXIS], meani[NAXIS],meanr[NAXIS],
		dx,dr2, ksig2;
   short	sexflagmask;
   unsigned int	imaflagmask;
   int		i,f,n,s, naxis,nfield,nsamp, flag, nmeani,nmeanr, nclipi,nclipr;

  sexflagmask = (short)prefs.astr_sexflagsmask;
  imaflagmask = prefs.astr_imaflagsmask;
  naxis = fgroup->naxis;

  ksig2 = nsigma;
  ksig2 *= ksig2;
  for (i=0; i<naxis; i++)
    {
    clipi[i] = fgroup->sig_interr[i] > 0.0 ?
			ksig2*fgroup->sig_interr[i]*fgroup->sig_interr[i]
			/(fgroup->meanwcsscale[i]*fgroup->meanwcsscale[i])
			: BIG;
    clipr[i] = fgroup->sig_referr[i] > 0.0 ?
			ksig2*fgroup->sig_referr[i]*fgroup->sig_referr[i]
			/(fgroup->meanwcsscale[i]*fgroup->meanwcsscale[i])
			: BIG;
    }

  nclipi = nclipr = 0;
  nfield = fgroup->nfield;
  fields = fgroup->field;
  for (f=0; f<nfield; f++)
    {
    field=fields[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      samp = set->sample;
      for (n=nsamp; n--; samp++)
        if (!samp->nextsamp && samp->prevsamp)
          {
          for (i=0; i<naxis; i++)
            meani[i] = meanr[i] = 0.0;
          nmeani = nmeanr = 0;
          for (samp2 = samp; samp2 && samp2->set->field->astromlabel>=0;
		samp2=samp2->prevsamp)
            {
            if ((samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
              continue;
            for (i=0; i<naxis; i++)
              meani[i] += samp2->projpos[i];
            nmeani++;
            }
          if (!(nmeani))
            continue;
          sampr = samp2;
          if (nmeani>1)
            {
            for (i=0; i<naxis; i++)
              meani[i] /= (double)nmeani;
            for (samp2 = samp; samp2 && samp2->set->field->astromlabel>=0;
		samp2 = prevsamp2)
              {
              flag = 0;
              dr2 = 0.0;
              for (i=0; i<naxis; i++)
                {
                dx = samp2->projpos[i] - meani[i];
                if ((dr2 += dx*dx) > clipi[i])
                  {
                  flag = 1;
                  break;
                  }
                }
              prevsamp2 = samp2->prevsamp;
              if ((flag))
                {
/*-------------- Remove (unlink) outlier */
                if (samp2->nextsamp)
                  samp2->nextsamp->prevsamp = samp2->prevsamp;
                if (samp2->prevsamp)
                  samp2->prevsamp->nextsamp = samp2->nextsamp;
                samp2->prevsamp = samp2->nextsamp = NULL;
                nclipi++;
                }
              else
                {
                for (i=0; i<naxis; i++)
                  meanr[i] += samp2->projpos[i];
                nmeanr++;
                }
              }
            }
          else
            {
            for (i=0; i<naxis; i++)
              meanr[i] = meani[i];
            nmeanr = 1;
            }
          if ((sampr) && (sampr->nextsamp))
            {
            flag = 0;
            dr2 = 0.0;
            for (i=0; i<naxis; i++)
              {
              dx = sampr->projpos[i] - meanr[i]/nmeanr;
              if ((dr2 += dx*dx) > clipr[i])
                {
                flag = 1;
                break;
                }
              }
            if ((flag))
              {
/*------------ Remove (unlink) outlier */
              sampr->nextsamp->prevsamp = NULL;
              sampr->nextsamp = NULL;
              nclipr++;
              }
            }
          }
      }
    }

  return nclipr+nclipi;
  }


