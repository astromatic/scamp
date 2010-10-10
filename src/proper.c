/*
*				proper.c
*
* Compute differential chromatic refraction, proper motions and parallaxes.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2008-2010 IAP/CNRS/UPMC
*
*	Author:			Emmanuel Bertin (IAP)
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
*	Last modified:		10/10/2010
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
#include "prefs.h"
#include "samples.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include ATLAS_LAPACK_H

/*------------------- global variables for multithreading -------------------*/
#ifdef USE_THREADS
#endif

/****** astrcolshift_fgroup ***************************************************
PROTO	void astrcolshift_fgroup(fgroupstruct *fgroup)
PURPOSE	Compute colour shift terms for sources in a group of fields.
INPUT	ptr to a group of field.
OUTPUT	-.
NOTES	Uses the global preferences. Input structures must have gone through
	reproj_fgroup() and crossid_fgroup() first, and preferably through
	astrsolve_fgroups and photsolve_fgroups() too.
AUTHOR	E. Bertin (IAP)
VERSION	29/06/2009
 ***/
void	astrcolshift_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
  {
   fieldstruct	*field1, *field2;
   setstruct	*set;
   samplestruct	*samp,*samp2;
   double	*ss[NAXIS],*sx[NAXIS],*sxx[NAXIS],*sy[NAXIS],*sxy[NAXIS],
		sigma[NAXIS],
		sig, a,b, wi,xi,yi, delta;
   int		*ncolshift[NAXIS],
		d,f1,f2,n,s, naxis, nfield, instru1,instru2, ninstru;

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
        if (samp->flux <= 0.0 || (samp->flags & (OBJ_SATUR|OBJ_TRUNC)))
          continue;
/*------ Explore the forward direction */
        if (samp->nextsamp)
	  {
          samp2 = samp;
          while ((samp2=samp2->nextsamp))
            {
            if (samp2->flux <= 0.0 || (samp2->flags & (OBJ_SATUR|OBJ_TRUNC)))
              continue;
            field2 = samp2->set->field;
            f2 = (field2==reffield ? nfield : field2->index);
            for (d=0; d<naxis; d++)
              {
              sig = samp2->wcsposerr[d];
              sig = sig*sig + sigma[d];
              if (sig>0.0)
                {
                wi = 1.0/sig;
                xi = samp2->mag - samp->mag;
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
/*------ Explore the backward direction */
        if (samp->prevsamp)
	  {
          samp2 = samp;
          while ((samp2=samp2->prevsamp))
            {
            if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
              continue;
            field2 = samp2->set->field;
            f2 = (field2==reffield ? nfield : field2->index);
            for (d=0; d<naxis; d++)
              {
              sig = samp2->wcsposerr[d];
              sig = sig*sig + sigma[d];
              if (sig>0.0)
                {
                wi = 1.0/sig;
                xi = samp2->mag - samp->mag;
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
      }

    for (f2=f1+1; f2<=nfield; f2++)
      {
      instru2 = f2<nfield? fgroup->field[f2]->photomlabel : ninstru;
      for (d=0; d<naxis; d++)
        if (ss[d][f2] > 0.0)
          {
          delta = ss[d][f2]*sxx[d][f2] - sx[d][f2]*sx[d][f2];
          if (delta == 0.0)
            continue;
          a = (ss[d][f2]*sxy[d][f2] - sx[d][f2]*sy[d][f2]) / delta;
          b = (sxx[d][f2]*sy[d][f2] - sx[d][f2]*sxy[d][f2]) / delta;
          if (f2==nfield)
            {
            fgroup->refcolshiftscale[d][f1] = a;
            fgroup->refcolshiftzero[d][f1] = b;
            }
          else
            {
            fgroup->intcolshiftscale[d][f2*nfield+f1]
		= fgroup->intcolshiftscale[d][f1*nfield+f2] = a;
            fgroup->intcolshiftzero[d][f2*nfield+f1]
		= -(fgroup->intcolshiftzero[d][f1*nfield+f2] = b);
            fgroup->colshiftscale[d][instru2*ninstru+instru1]
		= (fgroup->colshiftscale[d][instru1*ninstru+instru2] += a);
            fgroup->colshiftzero[d][instru2*ninstru+instru1]
		= -(fgroup->colshiftzero[d][instru1*ninstru+instru2] += b);
            ncolshift[d][instru1*ninstru+instru2]
		= ++ncolshift[d][instru2*ninstru+instru1];
            }
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
VERSION	29/08/2010
 ***/
void	astrprop_fgroup(fgroupstruct *fgroup)
  {
   fieldstruct	*field,*field1,*field2;
   setstruct	*set;
   samplestruct	*samp,*samp1,*samp2;
   wcsstruct	*wcs, *wcsec;
   char		*wcsectype[NAXIS];
   double	alpha[25], beta[5],
		sigma[NAXIS], coord[NAXIS], coorderr[NAXIS],
		projpos[NAXIS], ecpos[NAXIS], pfac1[NAXIS],pfac2[NAXIS],
		**csscale, **cszero,
		*wcsscale,
		sig, wi,wis,yi, dt,epoch, dmag, bec,cbec,lec,
		pfac,paral,paralerr;
   int		d,f,f1,ff,s,n, naxis, nfield, nsamp, lng,lat, celflag,
		colcorflag, propflag, paralflag, ncoeff,ncoeffp1;

  NFPRINTF(OUTPUT, "Computing proper motions...");

  colcorflag = prefs.colourshift_flag;
  paralflag = prefs.parallax_flag;
  naxis = fgroup->naxis;
  nfield = fgroup->nfield;
  wcs = fgroup->wcs;
  lng = wcs->lng;
  lat = wcs->lat;
  celflag = (lat>=0 && lng>=0);
  wcsscale = fgroup->meanwcsscale;
  csscale = fgroup->intcolshiftscale;
  cszero = fgroup->intcolshiftzero;
  paral = paralerr = 0.0;
  ncoeff = paralflag? 5 : 4;
  ncoeffp1 = ncoeff+1;

/* Set up a WCS structure to handle ecliptic coordinates */
  for (d=0; d<naxis; d++)
    QMALLOC(wcsectype[d], char, 16); 
  strcat(wcsectype[lng], "ELON-AIT");
  strcat(wcsectype[lat], "ELAT-AIT");
  wcsec = create_wcs(wcsectype, NULL, NULL, NULL, NULL, 2);
  for (d=0; d<naxis; d++)
    free(wcsectype[d]); 

/* Index fields within the group for local usage */
  for (f=0; f<nfield; f++)
    fgroup->field[f]->index = f;

  for (f=0; f<nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      samp = set->sample;
      for (n=nsamp; n--; samp++)
        if (!samp->nextsamp && samp->prevsamp)
          {
/*-------- Initialize matrices */
          memset(alpha, 0, ncoeff*ncoeff*sizeof(double));
          memset(beta, 0, ncoeff*sizeof(double));
          samp1 = samp;
          field1 = samp1->set->field;
          f1 = field1->index;
          epoch = field1->epoch;
          if (paralflag)
            {
/*---------- Compute parallax factors for detection #1 from ecliptic coords */
/*---------- (still quick and dirty, needs more proper determination) */
            ecpos[lng] = samp1->wcspos[lng];
            ecpos[lat] = samp1->wcspos[lat];
            eq_to_celsys(wcsec, ecpos);
            bec = ecpos[lat]*DEG;
            lec = field1->epoch - 0.246; /* l_sun ~ 0 when fractional part=0 */
            lec = ecpos[lng]*DEG - (lec - floor(lec))*2*PI;	/* l - l_sun */
            ecpos[lng] -= (fabs(cbec=cos(bec)) > 1e-12?
			(ARCSEC/DEG)*sin(lec)/cbec : 0.0);
            ecpos[lat] -= (ARCSEC/DEG)*cos(lec)*sin(bec);
            celsys_to_eq(wcsec, ecpos);
            wcs_to_raw(wcs, ecpos, projpos);
            pfac1[lng] = projpos[lng] - samp1->projpos[lng];
            pfac1[lat] = projpos[lat] - samp1->projpos[lat];
            }
          wis = wcs_scale(wcs, samp1->projpos);
          for (d=0; d<naxis; d++)
            {
/*---------- Another simplistic treatment of position uncertainties */
            sig = samp1->wcsposerr[d];
            sig = sig*sig;
             if (sig>0.0)
              {
              wi = wis/sig;
              alpha[(d+2)*ncoeffp1] += wi;
              }
            }
          samp2 = samp;
          while ((samp2=samp2->prevsamp) && samp2->set->field->astromlabel>=0)
            {
            field2 = samp2->set->field;
            ff = f1*nfield + field2->index;
            dt = field2->epoch - epoch;
            dmag = samp2->mag - samp1->mag;
            if (paralflag)
              {
/*------------ Compute parallax factors for detection #n */
/*------------ (still quick and dirty, needs more proper determination) */
              ecpos[lng] = samp2->wcspos[lng];
              ecpos[lat] = samp2->wcspos[lat];
              eq_to_celsys(wcsec, ecpos);
              bec = ecpos[lat]*DEG;
              lec = field2->epoch - 0.246;
              lec = ecpos[lng]*DEG - (lec - floor(lec))*2*PI;
              ecpos[lng] -= (fabs(cbec=cos(bec)) > 1e-12?
			(ARCSEC/DEG)*sin(lec)/cbec : 0.0);
              ecpos[lat] -= (ARCSEC/DEG)*cos(lec)*sin(bec);
              celsys_to_eq(wcsec, ecpos);
              wcs_to_raw(wcs, ecpos, projpos);
              pfac2[lng] = projpos[lng] - samp2->projpos[lng];
              pfac2[lat] = projpos[lat] - samp2->projpos[lat];
              }
            for (d=0; d<naxis; d++)
              {
              sig = samp2->wcsposerr[d];
              sig *= sig;
              if (sig>0.0)
                {
                yi = samp2->projpos[d] - samp1->projpos[d];
                if (colcorflag)
/*-------------- DCR correction */
                  yi -= dmag*csscale[d][ff] + cszero[d][ff];
                wi = wis/sig;
                alpha[d*ncoeffp1] += wi*dt*dt;
                alpha[(d+2)*ncoeffp1] += wi;
                alpha[2*ncoeff+ncoeffp1*d] = (alpha[d*ncoeffp1+2] += wi*dt);
                beta[d] += wi*yi*dt;
                beta[d+2] += wi*yi;
                if (paralflag)
                  {
                  pfac = pfac2[d] - pfac1[d];
                  alpha[d+20] = (alpha[d*5+4] += wi*dt*pfac);
                  alpha[d+22] = (alpha[d*5+14] += wi*pfac);
                  alpha[24] += wi*(pfac*pfac+1.0);
                  beta[4] += wi*yi*pfac;
                  }
                }
              }
            }
          clapack_dposv(CblasRowMajor, CblasUpper, ncoeff, 1, alpha, ncoeff,
			beta, ncoeff);
          clapack_dpotri(CblasRowMajor, CblasUpper, ncoeff, alpha, ncoeff);
          propflag = 1;
          for (d=0; d<naxis; d++)
            {
            coord[d] = samp->projpos[d] + beta[d];
            coorderr[d] = sqrt(alpha[d*ncoeffp1]*wis);
            }
          if (paralflag)
            {
            paral = beta[4];
            paralerr = sqrt(alpha[24]);
            }
          samp2 = samp;
          if (propflag)
            {
/*---------- Project shifted coordinates onto the sky... */
            raw_to_wcs(wcs, coord, coord);
/*---------- ... and recover the proper motion vector in celestial coords */
            for (d=0; d<naxis; d++)
              coord[d] -= samp->wcspos[d];
            if (celflag)
              coord[lng] *= cos(samp->wcspos[lat]*DEG);
            while ((samp2=samp2->prevsamp)
		&& samp2->set->field->astromlabel>=0)
              {
              for (d=0; d<naxis; d++)
                {
                samp2->wcsprop[d] = coord[d];
                samp2->wcsproperr[d] = fabs(coorderr[d]);
                }
              samp2->wcsparal = paral;
              samp2->wcsparalerr = paralerr;
              }
            }
          else
            while ((samp2=samp2->prevsamp)
		&& samp2->set->field->astromlabel>=0)
              {
              for (d=0; d<naxis; d++)
                samp2->wcsprop[d] = samp2->wcsproperr[d] = 0.0;
              samp2->wcsparal = samp2->wcsparalerr = 0.0;
              }
          }
      }
    }

  return;
  }


