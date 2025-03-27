/*
*				photsolve.c
*
* Compute the "global" photometric solution.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		28/06/2020
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
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "photsolve.h"
#include "prefs.h"
#include "samples.h"
#include "wcs/poly.h"

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_ACCELERATE
#include ACCELERATE_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
#endif

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

/****** photsolve_fgroups ****************************************************
PROTO	void photsolve_fgroups(fgroupstruct **fgroups, int nfgroup)
PURPOSE	Compute a global photometric solution among a group of fields.
INPUT	ptr to an array of group of fields pointers,
	number of groups.
OUTPUT	-.
NOTES	Uses the global preferences. Input structures must have gone through
	crossid_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	28/06/2020
 ***/
void	photsolve_fgroups(fgroupstruct **fgroups, int nfgroup)
  {
   fieldstruct	**fields,
		*field,*field2,*field3;
   setstruct	*set,*set2,*set3;
   samplestruct	*samp,*samp2,*samp3;
   double	*alpha,*beta,*cv, *cvo,*cvoa,
		*cv2,*cvo2,*cvoa2,
		*cc,*cco,*ccoa, *cc2,*cco2,*ccoa2,
		*coeffval,*coeffconst, *coeffval2,*coeffconst2, 
		sigma2,sigma3, sigmaoverlap, sigmaref,
		weight, weightref, mag;
   unsigned short sexflagmask;
   unsigned int	imaflagmask;
   int		*ci,*cio,*cioa,
		*ci2,*cio2,*cioa2,
		*coeffindex, *coeffindex2,
		f,f0,g,s, nfield,ninstru,nifields,ntfields,
		instru, naxis, ncoefftot, nlsamp,
		nsamp, nscoeffmax, nscoeffmax2,
		nlsampmax, offlag, photcode, nsamp2, nsamp3;


  sexflagmask = prefs.phot_sexflagsmask;
  imaflagmask = prefs.phot_imaflagsmask;
  alpha = beta = NULL;	/* To avoid gcc -Wall warnings*/
  naxis = fgroups[0]->naxis;
  ninstru = prefs.nphotinstrustr;
  if (!ninstru)
    ninstru = 1;

/* Solve a different system for each instrument */
  for (instru=0; instru<ninstru; instru++)
    {
    NFPRINTF(OUTPUT, "Initializing the global photometry matrix...");
    nifields = 0;
    offlag = 0;		/* For now until we sort out some issues */
    for (g=0 ; g<nfgroup; g++)
      {
      fields = fgroups[g]->field;
      nfield = fgroups[g]->nfield;
      ntfields = 0;
/*---- Check that there is at least one photometric field */
      f0 = -1;
      for (f=0; f<nfield; f++)
        {
        field=fields[f];
        if (field->photomlabel!=instru)
          continue;
        if (f0<0)
          f0 = f;
/*------ Check that the field is linked to at least one photometric field */
        photcode = (field->photomflag==1)?
		1 : check_fieldphotomoverlap(field, instru);
/*------ If not, make it a dummy photometric field */
        if (!photcode)
          field->photomflag = 2;
/*------ Flag if it is linked to a genuine photometric field */
        field->photomlink = (( photcode==1)? 1 : 0);
        ntfields++;
	}
      for (f=0; f<nfield; f++)
        {
        field=fields[f];
        if (field->photomlabel!=instru)
          continue;
/*------ Index current field */
/*
        if (field->photomflag)
          field->index = -1;
        else
*/
        field->index = nifields++;
/*------ Index sets */
        for (s=0; s<field->nset; s++)
          {
          set = field->set[s];
          set->index = s;
/*-------- Set flux scales to 1 by default */
          set->fluxscale = 1.0;
          }
        }
      }
/*-- Ignore Instruments with no calibrable fields */
    if (nifields>=1)
      {
      sigmaoverlap = prefs.magzero_interr[instru]*prefs.magzero_interr[instru];
      sigmaref = prefs.magzero_referr[instru]*prefs.magzero_referr[instru];
      ncoefftot = nifields;
      if (offlag)
        ncoefftot += fields[f0]->nset-1;

      QCALLOC(alpha, double, (size_t)ncoefftot*ncoefftot);
      QCALLOC(beta, double, ncoefftot);
      coeffval = coeffval2 = coeffconst = coeffconst2 = (double *)NULL;
      coeffindex = coeffindex2 = (int *)NULL;  /* To avoid gcc -Wall warnings*/

/*     Matrix parameters will be ordered as :               */
/*     Local reduced coordinate (NAXISn)                    */
/*     Fields                                               */
/*       (if prefs.phot_type[instru] == PHOT_ADJUSTOFFSETS) */
/*     Sets                                                 */
/*       (endif)                                            */

/*     Go through all groups, all fields, all sets, all samples */
      nlsampmax = 0;
      for (g=0; g<nfgroup; g++)
        {
        nfield = fgroups[g]->nfield;
        fields = fgroups[g]->field;
        for (f=0; f<nfield; f++)
          {
          field = fields[f];
          for (s=0; s<field->nset; s++)
            {
            set = field->set[s];
            samp = set->sample;
            for (nsamp=set->nsample; nsamp--; samp++)
              {
/*------------ Care only about one end of the series of linked detections */
              if (samp->prevsamp && !samp->nextsamp)
                {
                samp2 = samp;
/*-------------- Count the samples in the pair-list */
                for (nlsamp=1;
		samp2->set->field->photomlabel>=0 && (samp2=samp2->prevsamp);)
                  if (samp2->set->field->photomlabel==instru
			&& !(samp2->sexflags & sexflagmask)
			&& !(samp2->imaflags & imaflagmask))
                    nlsamp++;
/*-------------- Allocate memory for storing list sample information */
                if (nlsamp>nlsampmax)
                  {
                  if (nlsampmax)
                    {
                    free(coeffindex);
                    free(coeffval);
                    free(coeffconst);
                    free(coeffindex2);
                    free(coeffval2);
                    free(coeffconst2);
                    }
                  nlsampmax = nlsamp;
                  nscoeffmax = nlsamp;
                  QMALLOC(coeffindex, int, nscoeffmax);
                  QMALLOC(coeffval, double, nscoeffmax);
                  QMALLOC(coeffconst, double, nscoeffmax);
                  nscoeffmax2 = nlsamp;
                  QMALLOC(coeffindex2, int, nscoeffmax2);
                  QMALLOC(coeffval2, double, nscoeffmax2);
                  QMALLOC(coeffconst2, double, nscoeffmax2);
                  }
/*-------------- Fill list-sample coefficients */
                ci = coeffindex;
                ci2 = coeffindex2;
                cv = coeffval;
                cv2 = coeffval2;
                cc = coeffconst;
                cc2 = coeffconst2;
                nsamp2 = 0;
                for (samp2 = samp; samp2 && samp2->set->field->photomlabel>=0;
			samp2=samp2->prevsamp)
                  {
                  nsamp2++;
                  if (samp2->set->field->photomlabel != instru
			|| (samp2->sexflags & sexflagmask)
			|| (samp2->imaflags & imaflagmask))
                    continue;
                  cio = ci;
                  cio2 = ci2;
                  cvo = cv;
                  cvo2 = cv2;
                  cco = cc;
                  cco2 = cc2;
                  set2 = samp2->set;
                  field2 = set2->field;
/*---------------- Pre-compute numbers that will be used in the beta matrix  */
                  if (samp2->flux>0.0)
                    {
                    mag = -2.5*log10(samp2->flux/set2->expotime)
			- set2->extcoeff*(set2->airmass-1.0) + set2->magzero;
/*------------------ Compute the magnitude variance for this source */
                    sigma2 = (2.5/log(10.0))*samp2->fluxerr/samp2->flux;
                    }
                  else
                    {
                    mag = 99.0;
                    sigma2 = mag*mag;
                    }
/*---------------- Fill something equivalent to a row of a design matrix */
                  *(ci++) = field2->index;
                  *(cv++) = 0.0;
                  *(cc++) = mag;
                  if (offlag)
                    {
                    *(ci2++) = set2->index - 1;
                    *(cv2++) = 0.0;
                    *(cc2++) = mag;
                    }
/*---------------- Now fill the matrices */
                  cvoa = coeffval;
                  cioa = coeffindex;
                  ccoa = coeffconst;
                  cvoa2 = coeffval2;
                  cioa2 = coeffindex2;
                  ccoa2 = coeffconst2;
                  nsamp3 = nsamp2;
                  for (samp3 = samp; nsamp3--; samp3=samp3->prevsamp)
                    {
                    if (samp3->set->field->photomlabel != instru
			|| (samp3->sexflags & sexflagmask)
			|| (samp3->imaflags & imaflagmask))
                      continue;
                    set3 = samp3->set;
                    field3 = set3->field;
/*------------------ Compute the relative weight of the pair */
                    sigma3 = samp3->flux>0.0?
			(2.5/log(10.0))*samp3->fluxerr/samp3->flux
			: 99.0*99.0;
                    weight = 1.0/(sigma3+sigma2+sigmaoverlap);
                    weightref = 1.0/(sigma3+sigma2+sigmaref+sigmaoverlap);
/*------------------ Fill the matrices */
                    if (field2 != field3)
                      {
                      alpha[*cio+(size_t)ncoefftot**cio] += weight;
                      alpha[*cioa+(size_t)ncoefftot**cioa] += weight;
                      alpha[*cio+(size_t)ncoefftot**cioa] -= weight;
                      alpha[*cioa+(size_t)ncoefftot**cio] -= weight;
                      beta[*cio] += weight*(*ccoa - *cco);
                      beta[*cioa] += weight*(*cco - *ccoa);
                      if (field3->photomflag)
                        {
                        alpha[*cio+(size_t)ncoefftot**cio] += weightref;
                        beta[*cio] += weightref*(*ccoa+*cvoa - *cco);
                        }
                      }
                    if (field2->photomflag)
                      {
                      alpha[*cioa+(size_t)ncoefftot**cioa] += weightref;
                      beta[*cioa] += weightref*(*cco+*cvo - *ccoa);
                      }
                    cioa++;
                    cvoa++;
                    ccoa++;
                    }
                  }
                }
              }
            }
          }
        }

      free(coeffindex);
      free(coeffval);
      free(coeffconst);
      free(coeffindex2);
      free(coeffval2);
      free(coeffconst2);

/*---- Solve! */
      NFPRINTF(OUTPUT, "Solving the global photometry matrix...");
#if defined(HAVE_LAPACKE)
      LAPACKE_dposv(LAPACK_COL_MAJOR, 'L',
		ncoefftot, 1, alpha, ncoefftot, beta, ncoefftot);
#elif defined(HAVE_ACCELERATE)
      int info,n_beta=1;
      char* upper_or_lower="U";
      dposv_(upper_or_lower, &ncoefftot, &n_beta, alpha, &ncoefftot, beta, &ncoefftot, &info);
 #else
      clapack_dposv(CblasRowMajor, CblasUpper,
		ncoefftot, 1, alpha, ncoefftot, beta, ncoefftot);
#endif
      }

/*-- Update the field structures with the derived parameters */
    for (g=0; g<nfgroup; g++)
      {    
      nfield = fgroups[g]->nfield;
      fields = fgroups[g]->field;
      for (f=0; f<nfield; f++)
        {
        field=fields[f];
        if (field->photomlabel != instru)
          continue;
        field->dmagzero = 0.0;
        for (s=0; s<field->nset; s++)
          {
          set = field->set[s];
          field->dmagzero += (set->dmagzero = beta[field->index]);
          set->fluxscale = DEXP(0.4*(prefs.magzero_out[instru]
				- set->magzero - set->dmagzero
				+ set->extcoeff*(set->airmass-1.0)))
				/set->expotime;
          }
        field->dmagzero = fabs(field->dmagzero)>1e-6?
				field->dmagzero / field->nset : 0.0;
	}
      }

    if (nifields>=1)
      {
      free(alpha);
      free(beta);
      }
    }

  return;
  }


/****** photstats_fgroup *****************************************************
PROTO	void photstats_fgroups(fgroupstruct *fgroup, int instru,
				double hsn_thresh)
PURPOSE	Compute statistics about a global photometric solution among a group
	of fields.
INPUT	ptr to a group of fields,
	photometric instrument index,
	S/N threshold for the high-S/N sample.
OUTPUT	-.
NOTES	Input structures must have gone through crossid_fgroup() and
	compmags_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	28/06/2020
 ***/
void	photstats_fgroup(fgroupstruct *fgroup, int instru, double hsn_thresh)
  {
   fieldstruct	**fields,
		*field,*field2,*field3;
   setstruct	*set;
   samplestruct	*samp,*samp1,*samp2,*samp3;
   long double	lsig, lsig_hsn, lchi2, lchi2_hsn;
   double	mean,mean_hsn, dm, sn2, chi2, meanr, sig2;
   long long	ndeg,ndeg_hsn;
   unsigned short sexflagmask;
   unsigned int	imaflagmask;
   int		f,n,s, nfield,nsamp, nsource,nsource_hsn,
		nmean,nmean_hsn, nmatch,nmatch_hsn,
		nmeanr;

  sexflagmask = prefs.phot_sexflagsmask;
  imaflagmask = prefs.phot_imaflagsmask;
  lsig = lsig_hsn = 0.0;
  lchi2 = lchi2_hsn = 0.0;
  sn2 = 0.0;
  nsource = nsource_hsn = nmatch = nmatch_hsn = 0;
  ndeg = ndeg_hsn = 0;
  nfield = fgroup->nfield;
  fields = fgroup->field;

/* Internal stats */
  for (f=0; f<nfield; f++)
    {
    field=fields[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      samp1 = set->sample;
      for (n=nsamp; n--; samp1++)
        if (!samp1->nextsamp && samp1->prevsamp)
	  {
          mean = mean_hsn = 0.0;
          nmean = nmean_hsn = 0;
          samp2 = NULL;
/*-------- Look for a counterpart from the right photometric instrument */
          for (samp = samp1; samp && samp->set->field->photomlabel>=0;
		samp=samp->prevsamp)
            if (samp->set->field->photomlabel == instru
		&& !(samp->sexflags & sexflagmask)
		&& !(samp->imaflags & imaflagmask))
              {
              samp2 = samp;
              break;
              }
/*-------- No right instrument found: skip */
          if (!samp2)
            continue;
          for (; samp2 && samp2->set->field->photomlabel>=0;
		samp2=samp2->prevsamp)
	    {
            field2 = samp2->set->field;
/*---------- Don't bother if field is a different instru or photometric ref */
/*---------- or the flux is negative, or the object is saturated */
            if (field2->photomlabel != instru || samp2->flux <= 0.0
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
              continue;
            mean += samp2->mag;
            nmean++;
            sn2 = samp2->fluxerr>0.0? samp2->flux/samp2->fluxerr: 0.0;
            if (sn2 >= hsn_thresh)
	      {
              mean_hsn += samp2->mag;
              nmean_hsn++;
              }
            }
          if (nmean>1)
	    {
            nsource++;
            mean /= (double)nmean;
            if (nmean_hsn>1)
              {
              nsource_hsn++;
              mean_hsn /= (double)nmean_hsn;
              }
            for (samp2 = samp; samp2 && samp2->set->field->photomlabel>=0;
		samp2=samp2->prevsamp)
	      {
              field2 = samp2->set->field;
              if (field2->photomlabel != instru || samp2->flux <= 0.0
			|| (samp2->sexflags & sexflagmask)
			|| (samp2->imaflags & imaflagmask))
                continue;
              dm = samp2->mag - mean;
              lsig += dm*dm / (nmean-1.0);
              if (nmean_hsn>1)
	        {
                sn2 = samp2->fluxerr>0.0? samp2->flux/samp2->fluxerr: 0.0;
                if (sn2 >= hsn_thresh)
                  {
                  dm = samp2->mag - mean_hsn;
                  lsig_hsn += dm*dm / (nmean_hsn-1.0);
                  }
                }
              for (samp3=samp2; (samp3=samp3->prevsamp)
			&& samp3->set->field->photomlabel>=0;)
                {
                field3 = samp3->set->field;
                if (field3->photomlabel != instru || samp3->flux <= 0.0
			|| (samp3->sexflags & sexflagmask)
			|| (samp3->imaflags & imaflagmask))
                  continue;
                dm = samp2->mag - samp3->mag;
                lchi2 += (chi2 = dm*dm / (samp2->magerr*samp2->magerr
					  + samp3->magerr*samp3->magerr));
                ndeg++;
                if (nmean_hsn>1 && sn2>hsn_thresh && samp3->fluxerr>0.0
			&& samp3->flux/samp3->fluxerr>hsn_thresh)
		  {
                  lchi2_hsn += chi2;
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

  fgroup->sig_intmagerr[instru] = nsource? sqrt(lsig/nsource) : 0.0;
  fgroup->sig_intmagerr_hsn[instru] = nsource_hsn?
		sqrt(lsig_hsn/nsource_hsn) : 0.0;

  fgroup->chi2_intmag[instru] = ndeg? lchi2 / ((double)ndeg) : 0.0;
  fgroup->chi2_intmag_hsn[instru] = ndeg_hsn?
		lchi2_hsn / ((double)ndeg_hsn) : 0.0;
  fgroup->nintmagmatch[instru] = nmatch;
  fgroup->nintmagmatch_hsn[instru] = nmatch_hsn;

/* Now the stats with respect to the photometric reference fields */
  lsig = lsig_hsn = 0.0;
  lchi2 = lchi2_hsn = 0.0;
  nmatch = nmatch_hsn = 0;
  ndeg = ndeg_hsn = 0;
  for (f=0; f<nfield; f++)
    {
    field=fields[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      samp1 = set->sample;
      for (n=nsamp; n--; samp1++)
        if (!samp1->nextsamp && samp1->prevsamp)
	  {
          sig2 = meanr = mean = mean_hsn = 0.0;
          nmeanr = nmean = nmean_hsn = 0;
/*-------- Look for a counterpart from the right photometric instrument */
          for (samp = samp1; samp && samp->set->field->photomlabel>=0;
		samp=samp->prevsamp)
            {
            field = samp->set->field;
            if (field->photomlabel == instru && field->photomflag == 1
		&& samp->flux > 0.0
		&& !(samp->sexflags & sexflagmask)
		&& !(samp->imaflags & imaflagmask))
              {
              meanr += samp->mag;
              sig2 += samp->magerr*samp->magerr;
              nmeanr++;
              }
            }
/*-------- No right instrument found: skip */
          if (!nmeanr)
            continue;
          meanr /= (double)nmeanr;
          sig2 /= (double)nmeanr;
          for (samp2=samp1; samp2 && samp2->set->field->photomlabel>=0;
		samp2=samp2->prevsamp)
            {
            field2 = samp2->set->field;
            if (field2->photomlabel != instru || field2->photomflag == 1
		|| samp2->flux <= 0.0
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
              continue;
            mean += samp2->mag;
            dm = samp2->mag - meanr;
            lchi2 += (chi2 = dm*dm / (sig2 + samp2->magerr*samp2->magerr));
            nmean++;
            sn2 = samp2->fluxerr>0.0? samp2->flux/samp2->fluxerr: 0.0;
            if (sn2 >= hsn_thresh)
              {
              mean_hsn += samp2->mag;
              lchi2_hsn += chi2;
              nmean_hsn++;
              }
            }
          if (nmean)
            {
            if (nmean>1)
              {
              mean /= (double)nmean;
              if (nmean_hsn>1)
                mean_hsn /= (double)nmean_hsn;
              }
            dm = mean - meanr;
            lsig += dm * dm;
            nmatch++;
            if (nmean_hsn)
              {
              dm = mean_hsn - meanr;
              lsig_hsn += dm * dm;
              nmatch_hsn++;
              ndeg_hsn += nmean_hsn;
              }
            ndeg += nmean;
            }
          }
      }
    }

  fgroup->sig_refmagerr[instru] = nmatch? sqrt(lsig/nmatch) : 0.0;
  fgroup->sig_refmagerr_hsn[instru] = nmatch_hsn?
			sqrt(lsig_hsn/nmatch_hsn) : 0.0;
  fgroup->chi2_refmag[instru] = ndeg? lchi2 / (double)ndeg : 0.0;
  fgroup->chi2_refmag_hsn[instru] = ndeg_hsn?
			lchi2_hsn / (double)ndeg_hsn : 0.0;
  fgroup->nrefmagmatch[instru] = nmatch;
  fgroup->nrefmagmatch_hsn[instru] = nmatch_hsn;

  return;
  }


/****** photclip_fgroup *****************************************************
PROTO	void photclip_fgroup(fgroupstruct *fgroup, int instru, double nsigma)
PURPOSE	Removes outliers from cross-identification linked lists.
INPUT	ptr to a group of fields pointers,
	threshold in number of sigma
OUTPUT	-.
NOTES	Input structures must have gone through crossid_fgroup() and
	photstats_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	28/06/2020
 ***/
int	photclip_fgroup(fgroupstruct *fgroup, int instru, double nsigma)
  {
   fieldstruct	**fields,
		*field,*field2;
   setstruct	*set;
   samplestruct	*samp,*samp2;
   double	dm, clip, mean;
   unsigned short sexflagmask;
   unsigned int	imaflagmask;
   int		f,n,s, nfield,nsamp, nmean, nclip;

  sexflagmask = prefs.phot_sexflagsmask;
  imaflagmask = prefs.phot_imaflagsmask;
  clip = fgroup->sig_intmagerr[instru] > 0.0 ? nsigma*nsigma
		*fgroup->sig_intmagerr[instru]*fgroup->sig_intmagerr[instru]
			: BIG;
  nclip = 0;
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
          mean = 0.0;
          nmean = 0;
          for (samp2 = samp; samp2 && samp2->set->field->photomlabel>=0;
		samp2=samp2->prevsamp)
            {
            if (samp2->set->field->photomlabel != instru
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
              continue;
            mean += samp2->mag;
            nmean++;
            }
          if (nmean>1)
            {
            for (samp2 = samp; samp2 && samp2->set->field->photomlabel>=0;
		samp2=samp2->prevsamp)
              {
              field2 = samp2->set->field;
              if (field2->photomlabel != instru)
                continue;
              dm = samp2->mag - mean/nmean;
              if (dm*dm > clip)
                {
/*-------------- Remove outlier */
                if (samp2->nextsamp)
                  samp2->nextsamp->prevsamp = samp2->prevsamp;
                if (samp2->prevsamp)
                  samp2->prevsamp->nextsamp = samp2->nextsamp;
                samp2->prevsamp = samp2->nextsamp = NULL;
                nclip++;
                }
              }
            }
          }
      }
    }

  return nclip;
  }


/****** compmags_fgroup *******************************************************
PROTO	void compmags_fgroup(fgroupstruct *fgroup)
PURPOSE	Computes magnitudes from fluxes and zero-points in a group of fields.
INPUT	ptr to the group of fields.
OUTPUT	-.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	19/02/2018
 ***/
void	compmags_fgroup(fgroupstruct *fgroup)
  {
   fieldstruct	**fields;
   setstruct	**pset,
		*set;
   samplestruct	*samp;
   int		f, nsamp, s, nfield;

  fields = fgroup->field;
  nfield = fgroup->nfield;

  for (f=0; f<nfield; f++)
    {
    pset = fields[f]->set;
    for (s=fields[f]->nset; s--;)
      {
      set = *(pset++);
      samp = set->sample;
      for (nsamp=set->nsample; nsamp--; samp++)
        if (samp->flux>0.0)
	  {
          samp->mag = set->magzero
			- 2.5*log10(samp->flux/set->expotime)
			- set->extcoeff*(set->airmass-1.0) + set->dmagzero;
          samp->magerr = 1.0857*samp->fluxerr/samp->flux;
          }
        else
          samp->mag = samp->magerr = 99.0;
      }
    }

  return;
  }


/****** avermags_fgroup *******************************************************
PROTO	void avermags_fgroup(fgroupstruct *fgroup)
PURPOSE	Average magnitudes for a given photom instrument in a group of fields.
INPUT	ptr to the group of fields.
OUTPUT	-.
NOTES	Must have gone through compmags_fgroup() first.
	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	25/05/2005
 ***/
void	avermags_fgroup(fgroupstruct *fgroup)
  {
   fieldstruct	**fields;
   setstruct	**pset,
		*set;
   samplestruct	*samp,*samp1,*samp2;
   double	mag;
   int		f, s, instru, ninstru, nsamp, nfield, nmag;

  ninstru = prefs.nphotinstrustr;
  fields = fgroup->field;
  nfield = fgroup->nfield;

  for (f=0; f<nfield; f++)
    {
    pset = fields[f]->set;
    set = *(pset++);
    for (s=fields[f]->nset; s--; set=*(pset++))
      {
      samp1 = set->sample;
      for (nsamp=set->nsample; nsamp--; samp1++)
        if (!samp1->nextsamp && samp1->prevsamp)
          {
/*--------- Look for a counterpart from the right photometric instrument */
          for (instru = 0; instru<ninstru; instru++)
            {
            mag = 0.0;
            nmag = 0;
            for (samp = samp1; samp && samp->set->field->photomlabel>=0;
		samp=samp->prevsamp)
              if (samp->set->field->photomlabel == instru)
                {
                for (samp2=samp1; samp2 && samp2->set->field->photomlabel>=0;
			samp2=samp2->prevsamp)
                  {
/*---------------- Don't bother if field is a different instru or photom. ref*/
/*---------------- or if the flux is negative */
                  if (samp2->set->field->photomlabel != instru
			|| samp2->flux <= 0.0)
                    continue;
                  mag += samp2->mag;
                  nmag++;
                  }
                break;
                }
            if (nmag>1)
              {
              mag /= (double)nmag;
              for (samp = samp1; samp && samp->set->field->photomlabel>=0;
			samp=samp->prevsamp)
                if (samp->set->field->photomlabel == instru)
                  samp->mag = mag;
              }
            }
          }
      }
    }

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
 **
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

    No link found between both fields 
    return 0;
}
*/

/****** selecphotom_fgroup *****************************************************
PROTO	void selecphotom_fgroup(fgroupstruct *fgroup)
PURPOSE	Select fields according to their photometric properties and position in
	an fgroup.
INPUT	ptr to the group of fields.
OUTPUT	-.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	07/02/2005
 ***
void	selecphotom_fgroup(fgroupstruct *fgroup, int instru)
  {
   fieldstruct	**fields;
   setstruct	**pset,
		*set;
   samplestruct	*samp;
   int		*flist,
		f, nsamp, s, nfield;

  fields = fgroup->field;
  nfield = fgroup->nfield;
* Count the number of fields supposed to be photometric *
  for (f=0; f<nfield; f++)
    if (fields[f]->photomlabel==instru && fields[f]->photomflag==1)
      nphotom++;

* Identify connected groups of photometric exposures *
  QMALLOC(pfield, int, nfield);
  p = 1;
  for (f=0; f<nfield; f++)
    if (fields[f]->photomlabel==instru && fields[f]->photomflag==1)
      {
      if (p>1)
        {
        for (f2=f-1; f2--;)
          if (fields[f2]->photomlabel==instru && fields[f2]->photomflag==1)
            if (check_fieldoverlap(fields[f2], fields[f]))
              {
              pfield[f] = pfield[f2];
              break;
              }
        if (f2<0)
          pfield[f] = p++;
        }
      else
        pfield[f] = p++;
      p++;
      }

* Intra-group selection *
  np = p;
  if (prefs.photomselec_type[0] != PHOTOMSELEC_ALL)
    {
    for (p=0; p<np; p++)
      switch(prefs.photomselec_type[0])
        {
        case PHOTOMSELEC_AIRMASSMIN:
          value = BIG;
          fmin = -1;
          for (f=0; f<nfield; f++)
            if (pfield[f]==p && fields[f]->airmass < value)
              {
              value = fields[f]->airmass;
              fmin = f;
              }
          break;
        case PHOTOMSELEC_AIRMASSMAX:
          value = -BIG;
          fmin = -1;
          for (f=0; f<nfield; f++)
            if (pfield[f]==p && fields[f]->airmass > value)
              {
              value = fields[f]->airmass;
              fmin = f;
              }
          break;
        case PHOTOMSELEC_AIRMASSMEDIAN:
          QMALLOC(values, double, nphotom);
          nvalue = 0
          for (f=0; f<nfield; f++)
            if (pfield[f]==p)
              values[nvalues++] = field->airmass;
          value = fast_median(values, nvalue);
          free(values);
          fmin = -1;
          for (f=0; f<nfield; f++)
            if (pfield[f]==p && field->airmass==value)
              {
              fmin = f;
              break;
              }
          break;
        case PHOTOMSELEC_ZPCORRMIN:
          value = BIG;
          fmin = -1;
          for (f=0; f<nfield; f++)
            if (pfield[f]==p && field->dmagzero < value)
              {
              value = fields[f]->dmagzero;
              fmin = f;
              }
          break;
        case PHOTOMSELEC_ZPCORRMAX:
          value = -BIG;
          fmin = -1;
          for (f=0; f<nfield; f++)
            if (pfield[f]==p && field->dmagzero > value)
              {
              value = fields[f]->dmagzero;
              fmin = f;
              }
          break;
        case PHOTOMSELEC_ZPCORRMEDIAN:
          QMALLOC(values, double, nphotom);
          nvalue = 0
          for (f=0; f<nfield; f++)
            if (pfield[f]==p)
              values[nvalues++] = field->dmagzero;
          value = fast_median(values, nvalue);
          free(values);
          fmin = -1;
          for (f=0; f<nfield; f++)
            if (pfield[f]==p && field->dmagzero==value)
              {
              fmin = f;
              break;
              }
          break;
        case PHOTOMSELEC_DISTMIN:
          value = BIG;
          fmin = -1;
          for (f=0; f<nfield; f++)
            if (pfield[f]==p && field->dmagzero < value)
              {
              value = fields[f]->dmagzero;
              fmin = f;
              }
          break;
        default: 
          error(EXIT_FAILURE, "*Internal ERROR*: Type Unknown",
				" in selecphotom_type()");
        }
    }

* Select within overlapping photometric fields *
  QMALLOC(flist, int, nfield);

  for (l = 0; l<nfield; l++)
    {
    if (
    for (f=0; f<nfield; f++)
      {      
      pset = fields[f]->set;
      set = *(pset++);
      for (s=fields[f]->nset; s--; set=*(pset++))
        {
        samp = set->sample;
        for (nsamp=set->nsample; nsamp--; samp++)
        if (samp->flux>0.0)
	  {
          samp->mag = set->magzero
			- 2.5*log10(samp->flux/set->expotime)
			- set->extcoeff*(set->airmass-1.0) + set->dmagzero;
          samp->magerr = 1.0857*samp->fluxerr/samp->flux;
          }
        else
          samp->mag = samp->magerr = 99.0;
      }
    }

  return;
  }

*/
