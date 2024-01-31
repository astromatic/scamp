/*
*				proper.c
*
* Compute differential chromatic refraction, proper motions and parallaxes.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2008-2020 IAP/CNRS/SorbonneU
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

#ifdef USE_THREADS
#include <pthread.h>
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
#include "merge.h"
#include "prefs.h"
#include "proper.h"
#include "samples.h"

#ifdef USE_THREADS
#include "threads.h"
#endif

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_ACCELERATE
#include ACCELERATE_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
#endif

/*------------------- global variables for multithreading -------------------*/
#ifdef USE_THREADS
#endif

/*-------------------------------- protos -----------------------------------*/
static int	astrprop_solve(fgroupstruct *fgroup, samplestruct *samp,
				wcsstruct *wcsec, double *alpha,double *beta,
				double wis, double *chi2);

/****** astrcolshift_fgroup ***************************************************
PROTO	void astrcolshift_fgroup(fgroupstruct *fgroup)
PURPOSE	Compute colour shift terms for sources in a group of fields.
INPUT	ptr to a group of field.
OUTPUT	-.
NOTES	Uses the global preferences. Input structures must have gone through
	reproj_fgroup() and crossid_fgroup() first, and preferably through
	astrsolve_fgroups and photsolve_fgroups() too.
AUTHOR	E. Bertin (IAP)
VERSION	28/06/2020
 ***/
void	astrcolshift_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
  {
   fieldstruct	*field1, *field2;
   setstruct	*set;
   samplestruct	*samp,*samp2;
   double	*ss[NAXIS],*sx[NAXIS],*sxx[NAXIS],*sy[NAXIS],*sxy[NAXIS],
		sigma[NAXIS],
		sig, a,b, wi,xi,yi, delta;
   unsigned short sexflagmask;
   unsigned int	imaflagmask;
   int		*ncolshift[NAXIS],
		d,f1,f2,n,s, naxis, nfield, instru1,instru2, ninstru;

  sexflagmask = prefs.astr_sexflagsmask | prefs.phot_sexflagsmask;
  imaflagmask = prefs.astr_imaflagsmask | prefs.phot_imaflagsmask;
  nfield = fgroup->nfield;
  naxis = fgroup->naxis;
  ninstru = prefs.nphotinstrustr;

/* Index fields in a very simple way */
  for (f2=0; f2<fgroup->nfield; f2++)
    fgroup->field[f2]->index = f2;

/* Allocate memory for the relative colour-shift matrices */
  for (d=0; d<naxis; d++)
    {
    QCALLOC(fgroup->intcolshiftscale[d], double, nfield*nfield);
    QCALLOC(fgroup->intcolshiftzero[d], double, nfield*nfield);
    QCALLOC(fgroup->refcolshiftscale[d], double, nfield);
    QCALLOC(fgroup->refcolshiftzero[d], double, nfield);
    QCALLOC(fgroup->colshiftscale[d], double, ninstru*ninstru);
    QCALLOC(fgroup->colshiftzero[d], double, ninstru*ninstru);
    QCALLOC(ncolshift[d], int, ninstru*ninstru);
    QMALLOC(ss[d], double, nfield+1);
    QMALLOC(sx[d], double, nfield+1);
    QMALLOC(sxx[d], double, nfield+1);
    QMALLOC(sy[d], double, nfield+1);
    QMALLOC(sxy[d], double, nfield+1);
    }

  for (f1=0; f1<nfield; f1++)
    {
/*-- Reset accumulators */
    for (d=0; d<naxis; d++)
      {
      memset(ss[d], 0, (nfield+1)*sizeof(double));
      memset(sx[d], 0, (nfield+1)*sizeof(double));
      memset(sxx[d], 0, (nfield+1)*sizeof(double));
      memset(sy[d], 0, (nfield+1)*sizeof(double));
      memset(sxy[d], 0, (nfield+1)*sizeof(double));
      }
    field1 = fgroup->field[f1];
    instru1 = field1->photomlabel;
    for (s=0; s<field1->nset; s++)
      {
      set = field1->set[s];
      samp = set->sample;
      for (n=set->nsample; n--; samp++)
        {
        for (d=0; d<naxis; d++)
          {
/*-------- Another simplistic treatment of position uncertainties */
          sig = samp->wcsposerr[d];
          sigma[d] = sig*sig;
          }
        xi = samp->msamp ? samp->msamp->colour : 0.0;
/*------ Explore the forward direction */
        for (samp2=samp; (samp2=samp2->nextsamp);)
          {
          if (samp2->flux <= 0.0
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
            continue;
          field2 = samp2->set->field;
          f2 = (field2==reffield ? nfield : field2->index);
          for (d=0; d<naxis; d++)
            {
            sig = samp2->wcsposerr[d];
            sig = 1.0 /*sig*sig + sigma[d]*/;
            if (sig>0.0)
              {
              wi = 1.0/sig;
              yi = samp2->projpos[d] - samp->projpos[d];
              ss[d][f2] += wi;
              sx[d][f2] += wi*xi;
              sy[d][f2] += wi*yi;
              sxx[d][f2] += wi*xi*xi;
              sxy[d][f2] += wi*xi*yi;
              }
            }
          }
/*------ Explore the backward direction */
        for (samp2=samp; (samp2=samp2->prevsamp);)
          {
          if (samp2->flux <= 0.0
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
            continue;
          field2 = samp2->set->field;
          f2 = (field2==reffield ? nfield : field2->index);
          for (d=0; d<naxis; d++)
            {
            sig = samp2->wcsposerr[d];
            sig = 1.0 /*sig*sig + sigma[d]*/;
            if (sig>0.0)
              {
              wi = 1.0/sig;
              yi = samp2->projpos[d] - samp->projpos[d];
              ss[d][f2] += wi;
              sx[d][f2] += wi*xi;
              sy[d][f2] += wi*yi;
              sxx[d][f2] += wi*xi*xi;
              sxy[d][f2] += wi*xi*yi;
              }
            }
          }
        }
      }
    for (f2=0; f2<nfield; f2++)
      {
      instru2 = fgroup->field[f2]->photomlabel;
      for (d=0; d<naxis; d++)
        if (ss[d][f2] > 0.0)
          {
          delta = ss[d][f2]*sxx[d][f2] - sx[d][f2]*sx[d][f2];
          if (delta == 0.0)
            continue;
          a = (ss[d][f2]*sxy[d][f2] - sx[d][f2]*sy[d][f2]) / delta;
          b = (sxx[d][f2]*sy[d][f2] - sx[d][f2]*sxy[d][f2]) / delta;
          fgroup->intcolshiftscale[d][f1*nfield+f2] = a;
          fgroup->intcolshiftzero[d][f1*nfield+f2] = b;
          fgroup->colshiftscale[d][instru1*ninstru+instru2] += a;
          fgroup->colshiftzero[d][instru1*ninstru+instru2] += b;
          ++ncolshift[d][instru1*ninstru+instru2];
          }
      }
    }

  for (n=0; n<ninstru*ninstru; n++)
    if (ncolshift[n])
      for (d=0; d<naxis; d++)
        {
        fgroup->colshiftscale[d][n] /= ncolshift[d][n];
        fgroup->colshiftzero[d][n] /= ncolshift[d][n];
        }
  for (d=0; d<naxis; d++)
    {
    free(ncolshift[d]);
    free(ss[d]);
    free(sx[d]);
    free(sxx[d]);
    free(sy[d]);
    free(sxy[d]);
    }

  return;
  }


/****** astrprop_fgroup *****************************************************
PROTO	void astrprop_fgroup(fgroupstruct *fgroup)
PURPOSE	Compute proper motions for sources in a group of fields.
INPUT	ptr to a group of field.
OUTPUT	-.
NOTES	Uses the global preferences. Input structures must have gone through
	crossid_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	28/06/2020
 ***/
void	astrprop_fgroup(fgroupstruct *fgroup)
  {
   fieldstruct	*field;
   setstruct	*set;
   msamplestruct	*msamp;
   samplestruct	*samp,*samp2, *sampmin, *samppropmod;
   wcsstruct	*wcs, *wcsec;
   char		*wcsectype[NAXIS];
   double	alpha[25], beta[5],
		coord[NAXIS], coorderr[NAXIS], propmean[NAXIS],
		wis, paral,paralerr, chi2,chi2min, propmod,propmodmin;
   unsigned short sexflagmask;
   unsigned int	imaflagmask;
   int		d,f,s,m,n, naxis, nfield, nsamp,nsamp2, lng,lat, celflag,
		propflag, paralflag, refflag, ncoeff,ncoeffp1,
		nfree, nbad,nbadmax, bad, nfreemin, npropmean;

  sexflagmask = prefs.astr_sexflagsmask;
  imaflagmask = prefs.astr_imaflagsmask;
  paralflag = prefs.parallax_flag;
  propflag = prefs.propmotion_flag;
  refflag = prefs.astrefinprop_flag;
  naxis = fgroup->naxis;
  nfield = fgroup->nfield;
  wcs = fgroup->wcs;
  lng = wcs->lng;
  lat = wcs->lat;
  celflag = (lat>=0 && lng>=0);
  paral = paralerr = 0.0;
  ncoeff = paralflag? 5 : 4;
  ncoeffp1 = ncoeff+1;
  for (d=0; d<naxis; d++)
    propmean[d] = 0.0;
  npropmean = 0;

/* Set up a WCS structure to handle ecliptic coordinates */
  if (paralflag)
    {
    for (d=0; d<naxis; d++)
      QCALLOC(wcsectype[d], char, 16); 
    strcat(wcsectype[lng], "ELON-AIT");
    strcat(wcsectype[lat], "ELAT-AIT");
    wcsec = create_wcs(wcsectype, NULL, NULL, NULL, NULL, 2);
    for (d=0; d<naxis; d++)
      free(wcsectype[d]); 
    }
  else
    wcsec = NULL;


/* Index fields within the group for local usage */
  for (f=0; f<nfield; f++)
    fgroup->field[f]->index = f;


  msamp = fgroup->msample;
  for (m=fgroup->nmsample; m--; msamp++)
    {
    samp = msamp->samp;
    if (!samp->prevsamp)
      {
      for (d=0; d<naxis; d++)
         msamp->wcsprop[d] = msamp->wcsproperr[d] = 0.0;
      msamp->wcsparal = msamp->wcsparalerr = msamp->wcschi2 = 0.0;
      continue;
      }

/*-- Switch off SCAMP_BADPROPER flag */
    for (samp2=samp; samp2; samp2 = samp2->prevsamp)
      samp2->scampflags ^= (samp2->scampflags&SCAMP_BADPROPER);
/*-- Count detections */
    nsamp2 = msamp->npos_ok;
    if ((nbadmax=(int)(PROPER_MAXFRACBAD*nsamp2)) < PROPER_MAXNBAD)
      nbadmax = PROPER_MAXNBAD;
    wis = wcs_scale(wcs, samp->projpos);
    nfreemin = astrprop_solve(fgroup, samp, wcsec, alpha,beta, wis, &chi2min);
    sampmin = NULL;
    propmodmin = BIG;
    for (nbad=0; nbad<nbadmax; nbad++)
      {
/*----- Exit if chi2 is OK or if less than 2 degrees of freedom remain */
      if (chi2min<nfreemin*PROPER_MAXCHI2 || nfreemin<2)
        break;
      sampmin = samppropmod = NULL;
      for (samp2=samp; samp2; samp2 = samp2->prevsamp)
        {
        if (!(refflag || samp2->set->field->astromlabel>=0)
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask)
		|| (samp2->scampflags & SCAMP_BADPROPER))
          continue;
/*------ Switch detection to "bad" state */
        samp2->scampflags |= SCAMP_BADPROPER;
        nfree = astrprop_solve(fgroup,samp, wcsec,alpha,beta, wis, &chi2);
        if (chi2*nfreemin<chi2min*nfree)
          {
          chi2min = chi2;
          sampmin = samp2;
          nfreemin = nfree;
          }
        propmod = sqrt(beta[0]*beta[0]+beta[1]*beta[1]);
        if (propmod<propmodmin)
          {
          propmodmin = propmod;
          samppropmod = samp2;
          }
/*------ Revert detection to "good" state */
        samp2->scampflags ^= SCAMP_BADPROPER;
        }
      if (sampmin)
        {
        if (nfreemin<=0 && samppropmod)
          samppropmod->scampflags |= SCAMP_BADPROPER;
        else
          sampmin->scampflags |= SCAMP_BADPROPER;
        msamp->scampflags |= SCAMP_BADPROPER;
        }
      }

    if (nfreemin<0)
      {
/*---- Not enough "good" detections to derive a solution */
      for (d=0; d<naxis; d++)
        msamp->wcsprop[d] = msamp->wcsproperr[d] = 0.0;
      msamp->wcsparal = msamp->wcsparalerr = 0.0;
      msamp->wcschi2 = 0.0;
      continue;
      }

    if (sampmin)
      nfreemin = astrprop_solve(fgroup,samp,wcsec,alpha,beta,wis, &chi2min);
#if defined(HAVE_LAPACKE)
    LAPACKE_dpotri(LAPACK_COL_MAJOR, 'L',ncoeff, alpha, ncoeff);
#elif defined(HAVE_ACCELERATE)
    int info,n_beta=1;
    char* upper_or_lower="U";
    dpotri_(upper_or_lower, &ncoeff, alpha, &ncoeff, &info);
#else
    clapack_dpotri(CblasRowMajor, CblasUpper, ncoeff, alpha, ncoeff);
#endif
/*-- Prepare a normalization factor to avoid problems with large motions */
    propmod = 0.0;
    for (d=0; d<naxis; d++)
       propmod += beta[d]*beta[d];
    propmod = sqrt(propmod + 1.0);
    for (d=0; d<naxis; d++)
      {
      coord[d] = samp->projpos[d] + beta[d] / propmod;
      coorderr[d] = sqrt(alpha[d*ncoeffp1]*wis);
      }
    if (paralflag)
      {
      paral = beta[4];
      paralerr = sqrt(alpha[24]);
      }
/*-- Project shifted coordinates onto the sky... */
    raw_to_wcs(wcs, coord, coord);
/*-- ... and recover the proper motion vector in celestial coords */
    for (d=0; d<naxis; d++)
      coord[d] = (coord[d] - samp->wcspos[d])*propmod;
    if (celflag)
      coord[lng] *= cos(samp->wcspos[lat]*DEG);
    propmod = 0.0;
    for (d=0; d<naxis; d++)
      propmod += coord[d]*coord[d];
    if (sqrt(propmod)*DEG/MAS<20.0)
      {
      for (d=0; d<naxis; d++)
        propmean[d] += coord[d];
      npropmean++;
      }

    for (d=0; d<naxis; d++)
      {
      msamp->wcsprop[d] = coord[d];
      msamp->wcsproperr[d] = fabs(coorderr[d]);
      }
    msamp->wcsparal = paral;
    msamp->wcsparalerr = paralerr;
    msamp->wcschi2 = nfreemin? chi2min / nfreemin : 0.0;
    }

 for (d=0; d<naxis; d++)
   fgroup->meanwcsprop[d] = npropmean? propmean[d]/npropmean : 0.0;

  if (wcsec)
    end_wcs(wcsec);

  return;
  }


/****** astrprop_solve *****************************************************
PROTO	int astrprop_solve(samplestruct *samp, double *alpha, double *beta,
		double wis, double *chi2)
PURPOSE	Fill system alpha and beta matrices for solving proper motion equations.
INPUT	Ptr to the field group,
	ptr to last sample in matching detection chain,
	ptr to the ecliptic WCS structure (or NULL if parallax is turned off),
	ptr to alpha normal equation matrix,
	ptr to beta vector,
	projected pixel scale,
	ptr to output chi2 value (or chi2 if not needed).
OUTPUT	Number of "good" detections in the chain.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	18/08/2018
 ***/
static int	astrprop_solve(fgroupstruct *fgroup, samplestruct *samp,
			wcsstruct *wcsec, double *alpha, double *beta,
			double wis, double *chi2)
  {
   fieldstruct	*field, *field1,*field2;
   setstruct	*set, *set1;
   samplestruct	*samp1,*samp2;
   wcsstruct	*wcs;
   double 	mpos[NAXIS], mposw[NAXIS],  projpos[NAXIS],
		ecpos[NAXIS], pfac[NAXIS],a[25],b[5],
		**csscale, **cszero,
		sig, wi,yi, sy2, dt, bec,cbec,lec, dval, colour;
   unsigned short sexflagmask;
   unsigned int	imaflagmask;
   int		d,j,k, f1,ff, naxis, nfield, ncoeff,nfree, ncoeffp1, lng,lat,
		colcorflag,paralflag,meanflag,refflag;

  sexflagmask = prefs.astr_sexflagsmask;
  imaflagmask = prefs.astr_imaflagsmask;
  nfield = fgroup->nfield;
  naxis = fgroup->naxis;
  wcs = fgroup->wcs;
  lng = wcs->lng;
  lat = wcs->lat;
  csscale = fgroup->intcolshiftscale;
  cszero = fgroup->intcolshiftzero;
  colour = samp->msamp ? samp->msamp->colour : 0.0;
  refflag = prefs.astrefinprop_flag;
  colcorflag = prefs.colourshiftcorr_flag;
  paralflag = prefs.parallax_flag;
  meanflag = 1;	/* Always set for the 1st detection in chain */
  ncoeff = paralflag? 5 : 4;
  ncoeffp1 = ncoeff+1;
/* Initialize matrices and accumulators */
  memset(alpha, 0, ncoeff*ncoeff*sizeof(double));
  memset(beta, 0, ncoeff*sizeof(double));
  sy2 = 0.0;
  nfree = -ncoeff;

  set = samp->set;
  field = set->field;

  for (samp1=samp; samp1; samp1 = samp1->prevsamp)
    {
    if (!(refflag || samp1->set->field->astromlabel>=0)
		|| (samp1->sexflags & sexflagmask)
		|| (samp1->imaflags & imaflagmask)
		|| (samp1->scampflags & SCAMP_BADPROPER))
      continue;
    set1 = samp1->set;
    field1 = set1->field;
    f1 = field1->index;
    if (meanflag)
      {
      for (d=0; d<naxis; d++)
        mpos[d] = mposw[d] = 0.0;

/*---- First pass through overlapping detections to compute averages*/
      for (samp2=samp; samp2; samp2 = samp2->prevsamp)
        {
        if (!(refflag || samp2->set->field->astromlabel)
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask)
		|| (samp2->scampflags & SCAMP_BADPROPER))
          continue;
        field2 = samp2->set->field;
        ff = f1*nfield + field2->index;
        for (d=0; d<naxis; d++)
          {
/*-------- Simplistic treatment of position uncertainties */
          sig = samp2->wcsposerr[d];
          if (sig>0.0)
            wi = wis/(sig*sig);
          else
            wi = 0.0;
/*-------- DCR correction */
          mpos[d] += colcorflag?
		wi*(samp2->projpos[d] + colour*csscale[d][ff] + cszero[d][ff])
		: wi*samp2->projpos[d];
          mposw[d] += wi;
          }
        }
      meanflag = colcorflag;	/* Always set if colcorflag is set */
      }

    if (paralflag)
      {
/*---- Compute parallax factors from ecliptic coords */
/*---- (still quick and dirty, needs more proper determination) */
      ecpos[lng] = samp1->wcspos[lng];
      ecpos[lat] = samp1->wcspos[lat];
      eq_to_celsys(wcsec, ecpos);
      bec = ecpos[lat]*DEG;
      lec = samp1->epoch - 0.246; /* l_sun ~ 0 when frac. part=0 */
      lec = ecpos[lng]*DEG - (lec - floor(lec))*2*PI;	/* l - l_sun */
      ecpos[lng] -= (fabs(cbec=cos(bec)) > 1e-12?
			(ARCSEC/DEG)*sin(lec)/cbec : 0.0);
      ecpos[lat] -= (ARCSEC/DEG)*cos(lec)*sin(bec);
      celsys_to_eq(wcsec, ecpos);
      wcs_to_raw(wcs, ecpos, projpos);
      pfac[lng] = projpos[lng] - samp1->projpos[lng];
      pfac[lat] = projpos[lat] - samp1->projpos[lat];
      }

    dt = samp1->epoch - samp->epoch;

/* Fill the normal equation matrix */
/* Basis functions are              */
/* X0_i = t_i - t_ref * is_lng(d_i) */
/* X1_i = t_i - t_ref * is_lat(d_i) */
/* X2_i = is_lng(d_i)               */
/* X3_i = is_lat(d_i)               */
/* X4_i = dproj_d_i(sin(l_i-l_sun)/cos(b_i),cos(l_i-l_sun)*sin(b_i)), */
/* where d_i is the axis index (lng or lat) of measurement i          */

    for (d=0; d<naxis; d++)
      {
      if (mposw[d]<=0.0)
        continue;
      nfree++;
      sig = samp1->wcsposerr[d];
      wi = sig>0.0? wis/(sig*sig) : 0.0;
      alpha[d*ncoeffp1] += wi*dt*dt;
      alpha[(d+2)*ncoeffp1] += wi;
      alpha[2*ncoeff+ncoeffp1*d] = (alpha[d*ncoeffp1+2] += wi*dt);
      yi = samp1->projpos[d] - mpos[d] / mposw[d];
      sy2 += wi*yi*yi;
      beta[d] += wi*yi*dt;
      beta[d+2] += wi*yi;
      if (paralflag)
        {
        alpha[d+20] = (alpha[d*5+4] += wi*dt*pfac[d]);
        alpha[d+22] = (alpha[d*5+14] += wi*pfac[d]);
        alpha[24] += wi*pfac[d]*pfac[d];
        beta[4] += wi*yi*pfac[d];
        }
      }
    }

  if (nfree<0)
    {
    if (chi2)
      *chi2 = 0.0;
    return nfree;
    }

/* Make a copy of beta before it is overwritten */
  memcpy(b,beta,ncoeff*sizeof(double));
  memcpy(a,alpha,ncoeff*ncoeff*sizeof(double));
#if defined(HAVE_LAPACKE)
  LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', ncoeff, 1, alpha, ncoeff, beta, ncoeff);
#elif defined(HAVE_ACCELERATE)
  int info,n_beta=1;
  char* upper_or_lower="U";
  dposv_(upper_or_lower, &ncoeff, &n_beta, alpha, &ncoeff, beta, &ncoeff, &info);
#else
  clapack_dposv(CblasRowMajor,CblasUpper,ncoeff,1, alpha,ncoeff,beta,ncoeff);
#endif

  if (chi2)
    {
/*-- Compute chi2 */
    dval = 0.0;
    for (k=0; k<ncoeff; k++)
      {
      for (j=0; j<ncoeff; j++)
        dval += a[k*ncoeff+j]*beta[j]*beta[k];
      dval -= 2.0*beta[k]*b[k];
      }
    *chi2 = dval + sy2;
    }

  return nfree;
  }


/****** astrconnect_fgroup ****************************************************
PROTO	void astrconnect_fgroup(fgroupstruct *fgroup)
PURPOSE	Connect different "sources" (matched detections) in a group of fields,
	based on their proper motions.
INPUT	ptr to a group of field.
OUTPUT	-.
NOTES	Not complete yet. Input structures must have gone through
	astrprop_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	26/07/2012
 ***/
void	astrconnect_fgroup(fgroupstruct *fgroup)
  {
   fieldstruct	*field;
   setstruct	*set;
   samplestruct	*samp;
   double	depoch;
   int		f,s,n, nfield, nsamp;

  depoch = fgroup->epochmax - fgroup->epochmin;
  if (depoch<=0.0)
    return;

  nfield = fgroup->nfield;
  for (f=0; f<nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      samp = set->sample;
      for (n=nsamp; n--; samp++)
        {
        if (samp->nextsamp)
          continue;
        }
      }
    }

  return;
  }


