/*
 *				astrsolve.c
 *
 * Compute the "global" astrometric solution.
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *	This file part of:	SCAMP
 *
 *	Copyright:		(C) 2002-2023 IAP/CNRS/SorbonneU
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
 *	Last modified:		31/03/2021
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

#ifdef HAVE_MKL
#include MKL_H
#endif

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
#define MATSTORAGE_PACKED 1
#endif

#ifdef HAVE_OPENBLASP
#include BLAS_H
#endif

#include <sys/mman.h>

#include "define.h"
#include "globals.h"
#include "astrsolve.h"
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "merge.h"
#include "prefs.h"
#include "samples.h"
#include "wcs/poly.h"

#ifdef USE_THREADS
#include "threads.h"
#endif

//------------------- global variables for multithreading -------------------
#ifdef USE_THREADS
#define	NMAT_MUTEX	32
#define	NMAT_MUTEX2	32

pthread_t		*thread;
pthread_mutex_t		matmutex[NMAT_MUTEX], matmutex2[NMAT_MUTEX2],
			fillastromutex, nconstastromutex;

int			pthread_endflag,
			pthread_gindex, pthread_findex, pthread_sindex;

void			pthread_astrom_eval(fgroupstruct **fgroups,
				int ngroup, int ncoefftot, int *nconst,
				double *params),
			*pthread_astrom_eval_thread(void *arg);
#endif


//--------------------------- global variables  ---------------------------

fgroupstruct		**astrom_fgroups;
polystruct		*astrom_poly, *astrom_poly2;

double			*astrom_coeffs,
			*astrom_dcoeffs;

int			*astrom_nconst,
			astrom_nfgroup, astrom_ncoefftot;


static double	astrom_eval(void *extra, const double *coeffs, double *dcoeffs,
			const int ncoefftot, const double step),
		astrom_eval_update(setstruct *set, const double *coeffs,
			double *dcoeffs, int ncoefftot,
			polystruct *poly, polystruct *poly2, int *nconst);

static int	add_params(double *params, int naxis, double weight,
			int *ci, double *cv, double *cc, int ncoeff);

static void	astrom_add_dcoeffs(double *dcoeffs, int naxis,
        		double weight, int *ci, double *cv, double *cc,
        		int ncoeff);

/****** astrsolve_fgroups *****************************************************
PROTO	void astrsolve_fgroups(fgroupstruct **fgroups, int nfgroup)
PURPOSE	Compute a global astrometric solution among a group of fields.
INPUT	ptr to an array of group of fields pointers,
	number of groups.
OUTPUT	-.
NOTES	Uses the global preferences. Input structures must have gone through
	crossid_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	24/02/2021
 ***/
void	astrsolve_fgroups(fgroupstruct **fgroups, fieldstruct **reffields,
			 int nfgroup) {

   fieldstruct	**fields,**fields2,
		*field,*field2;
   setstruct	*set,*set2;
   samplestruct	*samp;
   polystruct	*poly, *poly2;
   char		str[64],
		**contextname,
		lap_equed;
   double	cmin[MAXCONTEXT],cmax[MAXCONTEXT],
		cscale[MAXCONTEXT], czero[MAXCONTEXT],
		*coeffs, *dcoeffs,
		dval;
   size_t	size;
   int		group2[NAXIS],
		*findex, *findex2, *nsetmax, *nconst,*nc,*nc2,
		c,cm,f,f2,g,g2,i,n,s,s2,sm, nfield,nfield2,
		ninstru, instru, naxis, ncoeff, ncoeff2, ncoefftot,
		npcoeff, npcoeff2,
		d, index, index2, minindex2,
		ncontext, cx,cy, nicoeff, nicoeff2, groupdeg2,
		startindex2, nmiss, it;

// Compute weight factors for each set based on the relative number of
// astrometric references and detections
  NFPRINTF(OUTPUT, "Initializing detection weight factors...");
  astrweight_fgroups(fgroups, nfgroup);

  NFPRINTF(OUTPUT, "Initializing the global astrometry matrix...");

  naxis = fgroups[0]->naxis;

// CONTEXT polynomial
  poly = poly_init(prefs.context_group, prefs.ncontext_name, prefs.group_deg,
            prefs.ngroup_deg);

// Linear field-dependent polynomial
  for (d=0; d<naxis; d++)
    group2[d] = 1;
  groupdeg2 = prefs.focal_deg;
  poly2 = poly_init(group2, naxis, &groupdeg2, 1);

// Use a different index for each instrument
  ninstru = prefs.nastrinstrustr;
  if (!ninstru)
    ninstru = 1;
  QCALLOC(findex, int, ninstru+1);
  QCALLOC(findex2, int, ninstru+1);
  QCALLOC(nsetmax, int, ninstru+1);

// Compute the total number of fields and the max number of sets per field
  ncontext = fgroups[0]->field[0]->set[0]->ncontext;
  contextname = fgroups[0]->field[0]->set[0]->contextname;

  for (c=0; c<ncontext; c++) {
    cmin[c] = BIG;
    cmax[c] = -BIG;
  }

  for (g=0 ; g<nfgroup; g++) {
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f=0; f<nfield; f++) {
      field=fields[f];
      instru = field->astromlabel;
//---- Index current field
      if (field->nset > nsetmax[instru])
        nsetmax[instru] = field->nset;
      switch(prefs.stability_type[instru]) {
        case STABILITY_INSTRUMENT:
          field->index = 0;
          field->index2 = findex2[instru+1]++ - 1;
          findex[instru+1] = nsetmax[instru];
          break;
        case STABILITY_EXPOSURE:
          field->index = findex[instru+1];
          findex[instru+1] += field->nset;
          field->index2 = -1;
          findex2[instru+1] = -1;
          break;
        case STABILITY_PREDISTORTED:
           field->index = -1;
           findex[instru+1] = 0;
           field->index2 = findex2[instru+1]++;
           break;
        default:
          break;
      }
//---- Index sets
      for (s=0; s<field->nset; s++) {
        set = field->set[s];
        set->index = s;
//------ Compute normalization factors for context values
        samp = set->sample;
        for (n=set->nsample; n--;) {
          for (c=0; c<ncontext; c++) {
            dval = samp->context[c];
//---------- Update min and max
            if (dval<cmin[c])
              cmin[c] = dval;
            if (dval>cmax[c])
              cmax[c] = dval;
          }
        }
      }
    }
  }

// Compute context rescaling coefficients (to avoid dynamic range problems)
  cx = cy = -1;
  for (c=0; c<ncontext; c++) {
    czero[c] = cmin[c];
    cscale[c] =  (dval = cmax[c] - cmin[c]) != 0.0 ? 1.0/dval : 1.0;
    if (contextname) {
//---- Identify X_IMAGE and Y_IMAGE parameters (will be processed separately)
      if (!strcmp(contextname[c], "X_IMAGE")
		|| !strcmp(contextname[c], "XWIN_IMAGE")
		|| !strcmp(contextname[c], "XPSF_IMAGE")
		|| !strcmp(contextname[c], "XMODEL_IMAGE")
		|| !strcmp(contextname[c], "XPEAK_IMAGE")
		|| !strcmp(set->contextname[c], prefs.centroid_key[0]))
        cx = c;
      else if (!strcmp(contextname[c], "Y_IMAGE")
		|| !strcmp(contextname[c], "YWIN_IMAGE")
		|| !strcmp(contextname[c], "YPSF_IMAGE")
		|| !strcmp(contextname[c], "YMODEL_IMAGE")
		|| !strcmp(contextname[c], "YPEAK_IMAGE")
		|| !strcmp(set->contextname[c], prefs.centroid_key[1]))
        cy = c;
    }
  }

// Compute the required number of coefficients
  npcoeff = poly->ncoeff;
  ncoeff = npcoeff*naxis;
  nicoeff = ncoeff*naxis;
  npcoeff2 = poly2->ncoeff;
  ncoeff2 = npcoeff2*naxis;
  nicoeff2 = ncoeff2*naxis;

  ncoefftot = startindex2 = 0;

  for (i=0; i<ninstru; i++) {
    if (prefs.stability_type[instru]!=STABILITY_PREDISTORTED)
	findex2[i+1]--;	// Remove first field (too many degrees of freedom)
//------ Set the extra exposure-dependent parameter flag
    if (findex[i+1]>0)
      ncoefftot += ncoeff*findex[i+1];
    startindex2 += ncoeff*findex[i+1];
    if (findex2[i+1]>0)
      ncoefftot += ncoeff2*findex2[i+1];
    if (i) {
      findex[i] += findex[i-1];
      findex2[i] += findex2[i-1];
    }

  }

// Compute the "definitive" index for each set
  for (g=0 ; g<nfgroup; g++) {
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f=0; f<nfield; f++) {
      field=fields[f];
      instru = field->astromlabel;
      index = ncoeff*(field->index+findex[instru]);
      index2 = (field->index2 <0)? -1
		: startindex2 + ncoeff2*(field->index2+findex2[instru]);
      field->index2 = index2;
      field->ncoeff2 = ncoeff2;
      for (s=0; s<field->nset; s++) {
        set = field->set[s];
        set->index = field->index<0? -1 : index + set->index*ncoeff;
        set->index2 = index2;
        set->nconst = 0;
        set->ncoeff = ncoeff;
//------ Take the opportunity to set the context scales and offset
        for (c=0; c<ncontext; c++) {
          set->contextoffset[c] = czero[c];
          set->contextscale[c] = cscale[c];
        }
        set->contextx = cx;
        set->contexty = cy;
      }
    }
 }

// Distortion coefficients ordered as :
// Local reduced coordinate (NAXISn)
//   Instru               (ninstru at most)
//   (if stability_type == EXPOSURE)
//     Field
//       Set
//         Polynomial term
//           Projected coordinate (NAXISn)
//   (else)
//     Set
//       Polynomial term
//         Projected coordinate (NAXISn)
//     Extra field
//       Field polynom term
//         Projected coordinate (NAXISn)
//   (endif)

  QCALLOC(nconst, int, ncoefftot);
  QCALLOC(coeffs, double, ncoefftot);
  QCALLOC(dcoeffs, double, ncoefftot);

  astrom_fgroups = fgroups;
  astrom_nfgroup = nfgroup;
  astrom_ncoefftot = ncoefftot;
  astrom_nconst = nconst;
  astrom_coeffs = coeffs;
  astrom_dcoeffs = dcoeffs;
  astrom_poly = poly;
  astrom_poly2 = poly2;


  for (it=0; it<ASTROM_MAXITER; it++) {

    sprintf(str, "Solving astrometry: iteration #%d", it + 1);
    NFPRINTF(OUTPUT, str);
    for (g=0; g<nfgroup; g++)
      reproj_fgroup(fgroups[g], reffields[g], REPROJ_NONE);

#ifdef USE_THREADS
    pthread_astrom_eval(fgroups, nfgroup, ncoefftot, nconst, dcoeffs);
#else
//-- Go through all groups, all fields, all sets, all samples
    for (g=0; g<nfgroup; g++) {
      nfield = fgroups[g]->nfield;
      fields = fgroups[g]->field;
      for (f=0; f<nfield; f++) {
        field = fields[f];
        for (s=0; s<field->nset; s++) {
          set = field->set[s];
          astrom_eval(set, dcoeffs, ncoefftot, poly,poly2, nconst);
        }
      }
    }
#endif
//-- Gradient descent
    for (c=0; c<ncoefftot; c++)
      coeffs[c] = ASTROM_UPDATEFACTOR * dcoeffs[c];  

//-- Fill the astrom structures with the derived parameters
    for (g=0; g<nfgroup; g++) {
      nfield = fgroups[g]->nfield;
      fields = fgroups[g]->field;
      for (f=0; f<nfield; f++) {
        field = fields[f];
        for (s=0; s<field->nset; s++) {
          set = field->set[s];
          mat_to_wcs(poly, poly2, coeffs, set);
        }
      }
    }
  }

  poly_end(poly);
  poly_end(poly2);
  free(coeffs);
  free(nconst);
  free(findex);
  free(findex2);
  free(nsetmax);

  return;
}

/****** astrom_eval ***************************************************
PROTO	double astrom_eval(void *extra, const double *coeffs, double *dcoeffs,
		const int ncoefftot, const double step)
PURPOSE Compute current chi2 and chi2 gradient for astrometry.
INPUT	Ptr to data,
	ptr to the vector of coefficients,
	ptr to the computed gradient (on return)
	total number of coefficients,
	input step from the line search routine.
OUTPUT	Current chi2 value.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION 31/03/2021
 ***/
static double astrom_eval(void *extra, const double *coeffs, double *dcoeffs,
		const int ncoefftot, const double step) {

   fgroupstruct	**fgroups;
   fieldstruct	**fields,
   		*field;
   setstruct	*set;
   polystruct	*poly, *poly2;
   double	chi2;
   int		*nconst,
   		nfgroup, nfield, naxis, f, g, s;


  naxis = astrom_fgroups[0]->naxis;

  nfgroup = astrom_nfgroup;
  poly = astrom_poly;
  poly2 = astrom_poly2;
  nconst = astrom_nconst;

  chi2 = 0.0;
// Go through all groups, all fields, all sets, all samples
  for (g=0; g<nfgroup; g++) {
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f=0; f<nfield; f++) {
      field = fields[f];
      for (s=0; s<field->nset; s++) {
        set = field->set[s];
        chi2 += astrom_eval_update(set, coeffs, dcoeffs, ncoefftot, poly, poly2,
        	nconst);
      }
    }
  }

  return chi2;    
}


/****** astrom_eval_update ***************************************************
PROTO	double astrom_eval_update(setstruct *set, double *coeffs,
		double *dcoeffs, int ncoefftot,
		polystruct *poly, polystruct *poly2, int *nconst)
PURPOSE Provide contributions to the astrometry chi2 and chi2 gradient from
	detections of a given set.
INPUT	Ptr to the set structure,
	ptr to the array of coefficients,
	ptr to the chi2 gradient (content modified in output),
	total number of coefficients,
	pointer to the 1st polynomial (steady, instrument-dependent),
	pointer to the 2nd polynomial (variable, exposure-dependent),
	pointer to the number of constraints per chip/exposure.
OUTPUT	Contribution to the total chi2.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION 31/03/2021
 ***/
static double	astrom_eval_update(setstruct *set, const double *coeffs,
		double *dcoeffs,
		int ncoefftot, polystruct *poly, polystruct *poly2,
		int *nconst) {

   fieldstruct		*field2;
   setstruct		*set2;
   samplestruct	*samp,*samp2,*samp3;
   double		dprojdred[NAXIS*NAXIS], context[MAXCONTEXT],
			redpos[NAXIS],
			*coeffval, *coeffval2, *coeffconst,*coeffconst2,
			*cv, *cva,*cvb, *cv2,*cva2,*cvb2,
			*cc,*cca,*ccb, *cc2,*cca2,*ccb2, *ccat, *ccbt,
			*co, *co2, *x,
			*basis,*basis2, *czero, *cscale,
			dx, sigma2,sigma3, weight, weightfac, val, chi2;
   float		*jac,*jact;
   unsigned short	sexflagmask;
   unsigned int	imaflagmask;
   int			*coeffindex,*coeffindex2,
			*ci,*cia,*cib, *ci2,*cia2,*cib2,
			nlsamp,nlsampmax, nscoeffmax, nscoeffmax2, instru,
			index,indext,index2,indext2, redflag, naxis, ncontext,
			npcoeff,ncoeff,nicoeff, npcoeff2,ncoeff2,nicoeff2,
			nsamp,
			cx,cy, c, d, d1,d2, i, p, lng,lat;

#ifdef USE_THREADS
   int			imut=0;
#endif

    sexflagmask = prefs.astr_sexflagsmask;
    imaflagmask = prefs.astr_imaflagsmask;
    naxis = set->naxis;
    npcoeff = poly->ncoeff;
    ncoeff = npcoeff*naxis;
    nicoeff = ncoeff*naxis;
    npcoeff2 = poly2->ncoeff;
    ncoeff2 = npcoeff2*naxis;
    nicoeff2 = ncoeff2*naxis;

    chi2 = 0.0;
    nlsampmax = 0;
    ncontext = set->ncontext;
    coeffval = coeffval2 = coeffconst = coeffconst2 = (double *)NULL;
    coeffindex = coeffindex2 = (int *)NULL;      // To avoid gcc -Wall warnings
    czero = set->contextoffset;
    cscale = set->contextscale;
    cx = set->contextx;
    cy = set->contexty;
    lng = set->lng;
    lat = set->lat;
    samp = set->sample;
    for (nsamp=set->nsample; nsamp--; samp++) {
//---- Care only about one end of the series of linked detections
      if (!samp->prevsamp || samp->nextsamp)
        continue;
      samp2 = samp;
//---- Count the samples in the pair-list
      for (nlsamp=1; (samp2=samp2->prevsamp);)
        if (!(samp2->sexflags & sexflagmask)
        	&& !(samp2->imaflags & imaflagmask))
           nlsamp++;
//---- Allocate memory for storing list sample information
      if (nlsamp > nlsampmax) {
        if ((nlsampmax)) {
          free(coeffindex);
          free(coeffval);
          free(coeffconst);
          free(coeffindex2);
          free(coeffval2);
          free(coeffconst2);
        }
        nlsampmax = nlsamp;
        nscoeffmax = nicoeff*nlsamp;
        QMALLOC(coeffindex, int, nscoeffmax);
        QMALLOC(coeffval, double, nscoeffmax);
        QMALLOC(coeffconst, double, nscoeffmax);
        nscoeffmax2 = nicoeff2*nlsamp;
        QMALLOC(coeffindex2, int, nscoeffmax2);
        QMALLOC(coeffval2, double, nscoeffmax2);
        QMALLOC(coeffconst2, double, nscoeffmax2);
      }
//---- Fill list-sample coefficients
      ci = coeffindex;
      ci2 = coeffindex2;
      cv = coeffval;
      cv2 = coeffval2;
      cc = coeffconst;
      cc2 = coeffconst2;
      for (samp2 = samp; samp2; samp2=samp2->prevsamp) {
//------ Drop it if the object is saturated or truncated
        if ((samp2->sexflags & sexflagmask) || (samp2->imaflags & imaflagmask))
            continue;
        cia = ci;
        cia2 = ci2;
        cva = cv;
        cva2 = cv2;
        cca = cc;
        cca2 = cc2;
        set2 = samp2->set;
        field2 = set2->field;
        instru = field2->astromlabel;
//------ Negative indices are for the reference field
        if (instru>=0) {
          index = set2->index;
          index2 = set2->index2;
#ifdef USE_THREADS
          QPTHREAD_MUTEX_LOCK(&nconstastromutex);
#endif
          set2->nconst += naxis;
          if (index>=0)
            nconst[index]+=naxis;
          if (index2>=0)
            nconst[index2]+=naxis;
#ifdef USE_THREADS
          QPTHREAD_MUTEX_UNLOCK(&nconstastromutex);
#endif
        } else
          index = index2 = -1;

//------ Compute the value of the basis functions at each sample
        if (instru >= 0) {
          redflag = 0;
          for (c=0; c<ncontext; c++) {
            if (c==cx || c==cy) {
//------------ Image coordinates receive a special treatment
              if (!redflag) {
                raw_to_red(set2->wcs, samp2->vrawpos, redpos);
                redflag = 1;
              }
              context[c] = redpos[c==cx?lng:lat];
            } else
              context[c] = (samp2->context[c]-czero[c])*cscale[c];
          }
//-------- The basis functions are computed from reduced context values
          if (index>=0)
            poly_func(poly, context);
          if (index2>=0) {
            if (!redflag)
              raw_to_red(set2->wcs, samp2->vrawpos, redpos);
            poly_func(poly2, redpos);
          }
        }

//------ Pre-compute numbers that will be used in the dcoeffs vector
//------ The innermost index of the jac array is the numerator of the
//------ Jacobian operator
        x = samp2->projpos;
//------ Fill something equivalent to a row of a design matrix
        for (d2=0; d2<naxis; d2++) {
          indext = index;
          dx = *(x++);
          co = (double *)coeffs + index;
          basis = poly->basis;
          jac = samp2->dprojdred + d2;
          for (p=npcoeff; p--; basis++) {
            jact = jac;
            for (d1=naxis; d1--; jact+=naxis) {
              *(ci++) = indext;
              *(cv++) = (val = *jact * *basis);
              if ((++indext))
                dx += val * *(co++);
              else
                indext--;
              }
          }
//-------- Now the extra field-dependent parameters
          indext2 = index2;
          co2 = (double *)coeffs + index2;
          basis2 = poly2->basis;
          for (p=npcoeff2; p--; basis2++) {
            jact = jac;
            for (d1=naxis; d1--; jact+=naxis) {
              *(ci2++) = indext2;
              *(cv2++) = (val = *jact * *basis2);
              if ((++indext2))
                dx += val * *(co2++);
              else
                indext2--;
            }
          }
          for (i=nicoeff; i--;)
            *(cc++) = dx;
          for (i=nicoeff2; i--;)
            *(cc2++) = dx;
        }
//------ Compute the relative positional variance for this source
        sigma2 = 0.0;
        for (d=0; d<naxis; d++)
          sigma2 += samp2->wcsposerr[d]*samp2->wcsposerr[d];
//------ Add extra-weight to the reference to compensate for all the overlaps
        weightfac = nlsamp*prefs.astref_weight;
//------ Now fill the array
        cvb = coeffval;
        cib = coeffindex;
        ccb = coeffconst;
        cvb2 = coeffval2;
        cib2 = coeffindex2;
        ccb2 = coeffconst2;
        for (samp3 = samp; samp3 != samp2; samp3=samp3->prevsamp) {
//-------- Drop it if the object is saturated or truncated
          if ((samp3->sexflags & sexflagmask)
          	|| (samp3->imaflags & imaflagmask))
            continue;

//-------- Compute the relative weight of the pair
          sigma3 = 0.0;
          for (d=0; d<naxis; d++)
            sigma3 += samp3->wcsposerr[d]*samp3->wcsposerr[d];
          weight = (instru<0 ? samp3->set->weightfac*weightfac : 1.0)
			/(sigma3+sigma2);

//-------- Compute contribution from the pair to the chi2
          ccat = cca;
          ccbt = ccb;
          for (d=naxis; d--; ccat += nicoeff, ccbt += nicoeff)
            chi2 += weight * (*ccat - *ccbt)  * (*ccat - *ccbt);
          ccat = cca2;
          ccbt = ccb2;
          for (d=naxis; d--; ccat += nicoeff2, ccbt += nicoeff2)
            chi2 += weight * (*ccat - *ccbt)  * (*ccat - *ccbt);

//-------- Fill the matrices
//-------- First, the contributions involving Jacobians for star A
#ifdef USE_THREADS
          if (*cia >= 0) {
            imut = (*cia/ncoeff)%NMAT_MUTEX;
             QPTHREAD_MUTEX_LOCK(&matmutex[imut]);
          }
#endif
          astrom_add_dcoeffs(dcoeffs,naxis,weight, cia,cva,cca,ncoeff);
          astrom_add_dcoeffs(dcoeffs,naxis,-weight, cia,cva,ccb,ncoeff);
#ifdef USE_THREADS
          if (*cia >= 0) {
            QPTHREAD_MUTEX_UNLOCK(&matmutex[imut]);
          }
          if (*cia2 >= 0) {
            imut = (*cia2/ncoeff)%NMAT_MUTEX2;
          }
#endif
            astrom_add_dcoeffs(dcoeffs,naxis,weight, cia2,cva2,cca2,ncoeff2);
            astrom_add_dcoeffs(dcoeffs,naxis,-weight, cia2,cva2,ccb2,ncoeff2);

//-------- Second, the contributions involving Jacobians for star B
#ifdef USE_THREADS
            if (*cia2 >= 0) {
              QPTHREAD_MUTEX_UNLOCK(&matmutex2[imut]);
            }
            if (*cib >= 0) {
              imut = (*cib/ncoeff)%NMAT_MUTEX;
              QPTHREAD_MUTEX_LOCK(&matmutex[imut]);
            }
#endif
          astrom_add_dcoeffs(dcoeffs,naxis,-weight, cib,cvb,cca,ncoeff);
          astrom_add_dcoeffs(dcoeffs,naxis,weight, cib,cvb,ccb,ncoeff);
#ifdef USE_THREADS
          if (*cib >= 0) {
            QPTHREAD_MUTEX_UNLOCK(&matmutex[imut]);
          }
          if (*cib2 >= 0) {
            imut = (*cib2/ncoeff)%NMAT_MUTEX2;
            QPTHREAD_MUTEX_LOCK(&matmutex2[imut]);
          }
#endif
          astrom_add_dcoeffs(dcoeffs,naxis,-weight,cib2,cvb2,cca2,ncoeff2);
          astrom_add_dcoeffs(dcoeffs,naxis,weight,cib2,cvb2,ccb2,ncoeff2);

#ifdef USE_THREADS
          if (*cib2 >= 0) {
             QPTHREAD_MUTEX_UNLOCK(&matmutex2[imut]);
          }
#endif
          cvb += nicoeff;
          cib += nicoeff;
          ccb += nicoeff;
          cvb2 += nicoeff2;
          cib2 += nicoeff2;
          ccb2 += nicoeff2;
        }
      }
    }

    if (nlsampmax) {
      free(coeffindex);
      free(coeffval);
      free(coeffconst);
      free(coeffindex2);
      free(coeffval2);
      free(coeffconst2);
    }

  return chi2;
}


/****** astrom_add_dcoeffs ***************************************************
PROTO	void astrom_add_dcoeffs(double *dcoeffs, int naxis,
		double weight, int *ci, double *cv, double *cc, int ncoeff)
PURPOSE	Add a product contribution to the chi2 gradient.
INPUT	ptr to the gradient (modified in output),
	number of dimensions,
	weight of the contribution (can be negative for cross-products),
	index array,
	basis array,
	constraint array,
	size of the first component arrays.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	31/03/2021
 ***/
static void astrom_add_dcoeffs(double *dcoeffs, int naxis,
        double weight, int *ci, double *cv, double *cc, int ncoeff) {

   int		d, p;

  if (*ci < 0)
    return;

  for (d=naxis; d--;)
#pragma ivdep
    for (p=ncoeff; p--;)
      dcoeffs[*(ci++)] += 2.0 * weight * *(cc++) * *(cv++);

  return;
}


#ifdef USE_THREADS

/****** pthread_astrom_eval_thread ***************************************
  PROTO   void *pthread_astrom_eval_thread(void *arg)
  PURPOSE thread that takes care of fill the astrometry matrix.
  INPUT   Pointer to the thread number.
  OUTPUT  -.
  NOTES   Relies on global variables.
  AUTHOR  E. Bertin (IAP)
  VERSION 23/03/2021
 ***/
void    *pthread_astrom_eval_thread(void *arg) {

    polystruct	*poly,*poly2;
    char	str[128];
    int		group2[NAXIS],
		findex, gindex, sindex, proc, naxis, groupdeg2, d;

    gindex = findex = -1;
    proc = *((int *)arg);
    naxis = astrom_fgroups[0]->naxis;
    for (d=0; d<naxis; d++)
        group2[d] = 1;
    groupdeg2 = prefs.focal_deg;
    poly = astrom_poly;
    poly2 = astrom_poly2;
    threads_gate_sync(pthread_startgate);
    while (!pthread_endflag)
    {
        QPTHREAD_MUTEX_LOCK(&fillastromutex);
        if (pthread_gindex < astrom_nfgroup)
        {
            gindex = pthread_gindex;
            findex = pthread_findex;
            sindex = pthread_sindex++;
            if (pthread_sindex>=astrom_fgroups[gindex]->field[findex]->nset)
            {
                pthread_sindex = 0;
                pthread_findex++;
            }
            if (pthread_findex>=astrom_fgroups[gindex]->nfield)
            {
                pthread_findex = 0;
                pthread_gindex++;
            }
/*
            if (!sindex)
            {
                sprintf(str, "Filling the global astrometry matrix: "
                        "group %d/%d, field %d/%d",
                        gindex+1, astrom_ngroup,
                        findex+1, astrom_fgroups[gindex]->nfield);
                NFPRINTF(OUTPUT, str);
            }
*/
            QPTHREAD_MUTEX_UNLOCK(&fillastromutex);

            /*---- Fetch field data */
            astrom_eval_update(astrom_fgroups[gindex]->field[findex]->set[sindex],
                    astrom_coeffs, astrom_dcoeffs, astrom_ncoefftot,
                    poly, poly2, astrom_nconst);
        }
        else
        {
            QPTHREAD_MUTEX_UNLOCK(&fillastromutex);
            /*---- Wait for the input buffer to be updated */
            threads_gate_sync(pthread_stopgate);
            /* ( Master thread process loads and saves new data here ) */
            threads_gate_sync(pthread_startgate);
        }
    }

    poly_end(poly);
    poly_end(poly2);

    pthread_exit(NULL);

    return (void *)NULL;
}


/****** pthread_astrom_eval **********************************************
  PROTO   double pthread_astrom_eval(fgroupstruct **fgroups, int ngroup,
  int ncoefftot, int *nconst, double *coeffs)
  PURPOSE	Start multithreaded filling of the astrometry matrix
  INPUT	ptr to the array of groups of fields,
  number of groups involved,
  matrix size,
  pointer to an array of constraint counts,
  pointer of pointer to the beta vector matrix (output).
  OUTPUT	-.
  NOTES	Uses the global preferences.
  AUTHOR	E. Bertin (IAP)
  VERSION	08/03/2007
 ***/
void	pthread_astrom_eval(fgroupstruct **fgroups, int ngroup,
        int ncoefftot, int *nconst, double *params)
{
    static pthread_attr_t	pthread_attr;
    int				*proc,
                    p, i;

    /* Number of active threads */
    nproc = prefs.nthreads;
    /* Set up multi-threading stuff */
    QMALLOC(proc, int, nproc);
    QMALLOC(thread, pthread_t, nproc);
    QPTHREAD_MUTEX_INIT(&fillastromutex, NULL);
    QPTHREAD_MUTEX_INIT(&nconstastromutex, NULL);
    for (i=0; i<NMAT_MUTEX; i++)
    {
        QPTHREAD_MUTEX_INIT(&matmutex[i], NULL);
    }
    for (i=0; i<NMAT_MUTEX2; i++)
    {
        QPTHREAD_MUTEX_INIT(&matmutex2[i], NULL);
    }
    QPTHREAD_ATTR_INIT(&pthread_attr);
    QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
    pthread_startgate = threads_gate_init(nproc+1, NULL);
    pthread_stopgate = threads_gate_init(nproc+1, NULL);

    /* Start the reading threads */
    for (p=0; p<nproc; p++)
    {
        proc[p] = p;
        QPTHREAD_CREATE(&thread[p], &pthread_attr,
                &pthread_astrom_eval_thread, &proc[p]);
    }

    QPTHREAD_MUTEX_LOCK(&fillastromutex);
    pthread_gindex = pthread_findex = pthread_sindex = 0;
    pthread_endflag = 0;
    QPTHREAD_MUTEX_UNLOCK(&fillastromutex);
    /* Release threads!! */
    threads_gate_sync(pthread_startgate);
    /* ( Slave threads process the current buffer data here ) */
    threads_gate_sync(pthread_stopgate);
    pthread_endflag = 1;
    /* (Re-)activate existing threads... */
    threads_gate_sync(pthread_startgate);
    /* ... and shutdown all threads */
    for (p=0; p<nproc; p++)
        QPTHREAD_JOIN(thread[p], NULL);

    /* Clean up multi-threading stuff */
    threads_gate_end(pthread_startgate);
    threads_gate_end(pthread_stopgate);
    pthread_mutex_destroy(&fillastromutex);
    pthread_mutex_destroy(&nconstastromutex);
    for (i=0; i<NMAT_MUTEX; i++)
        pthread_mutex_destroy(&matmutex[i]);
    for (i=0; i<NMAT_MUTEX2; i++)
        pthread_mutex_destroy(&matmutex2[i]);
    QPTHREAD_ATTR_DESTROY(&pthread_attr);

    free(proc);
    free(thread);

    return;
}


#endif


/****** mat_to_wcs *********************************************************
  PROTO	void mat_to_wcs(polystruct *poly, polystruct *poly2, double *mat,
  setstruct *set)
  PURPOSE	Fill a WCS structure according to the matrix solution.
  INPUT	ptr to the instrument-dependent polynom structure,
  ptr to the exposure-dependent polynom structure,
  ptr to the solution vector,
  ptr to the WCS structure to be modified.
  OUTPUT	-.
  NOTES	Unfortunately the present WCS description of distortions is valid
  in 2D only.
  AUTHOR	E. Bertin (IAP)
  VERSION	04/04/2013
 ***/
void mat_to_wcs(polystruct *poly, polystruct *poly2, double *mat,
        setstruct *set)
{
    static const int pvindex[8][8] =	{
        {0,1,4,7,12,17,24,31},
        {2,5,8,13,18,25,32,0},
        {6,9,14,19,26,33,0,0},
        {10,15,20,27,34,0,0,0},
        {16,21,28,35,0,0,0,0},
        {22,29,36,0,0,0,0,0},
        {30,37,0,0,0,0,0,0},
        {38,0,0,0,0,0,0,0}
        };

    wcsstruct	*wcs;
    samplestruct	*sample;
    double	cd[NAXIS*NAXIS],
    *context, *cscale,*czero, *projp, *mat1, *mat2,
    dval, det;
    int		*powers, *powerst,
            c,cx,cy, d, e,ex,ey, p,s, lng,lat, naxis, ncoeff,ncoeff2,
            ncontext, ival;

    wcs = set->wcs;
    naxis = wcs->naxis;
    lng = wcs->lng;
    lat = wcs->lat;
    ncoeff = set->ncoeff/naxis;
    ncoeff2 = poly2->ncoeff;
    ncontext = poly->ndim;
    cscale = set->contextscale;
    czero = set->contextoffset;
    projp = wcs->projp;

    /* Precompute FITS or average parameters */
    QMALLOC(context, double, ncontext);
    cx = cy = -1;
    for (c=0; c<ncontext; c++)
    {
        if (!strcmp(set->contextname[c], "X_IMAGE")
                || !strcmp(set->contextname[c], "XWIN_IMAGE")
                || !strcmp(set->contextname[c], "XPSF_IMAGE")
                || !strcmp(set->contextname[c], "XMODEL_IMAGE")
                || !strcmp(set->contextname[c], "XPEAK_IMAGE")
                || !strcmp(set->contextname[c], prefs.centroid_key[0]))
            cx = c;
        else if (!strcmp(set->contextname[c], "Y_IMAGE")
                || !strcmp(set->contextname[c], "YWIN_IMAGE")
                || !strcmp(set->contextname[c], "YPSF_IMAGE")
                || !strcmp(set->contextname[c], "YMODEL_IMAGE")
                || !strcmp(set->contextname[c], "YPEAK_IMAGE")
                || !strcmp(set->contextname[c], prefs.centroid_key[1]))
            cy = c;
        else if (set->contextname[c][0] == ':')
            /*---- It is a FITS keyword parameter */
        {
            if (set->nsample)
                /*------ Samples are available: lets grab the value here */
                context[c] = (set->sample[0].context[c] - czero[c])*cscale[c];
            /*------ Go and dig the FITS keyword (to be completed) */
        }
        else
            /*---- Average over all samples of this set */
        {
            dval = 0.0;
            sample = set->sample;
            for (s=set->nsample; s--; sample++)
                dval += sample->context[c];
            context[c] = (dval/set->nsample - czero[c])*cscale[c];
        }
    }

    if (set->index>=0 && (cx>=0 || cy>=0))
    {
        for (d=0; d<naxis; d++)
            if (!wcs->nprojp)
                projp[1+100*d] += 1.0;	// Default to the linear term

        /*-- Compute the array of exponents of polynom terms */
        powerst = powers = poly_powers(poly);
        mat1 = mat + set->index;	/* Instrument-dependent section in solution */
        for (p=0; p<ncoeff; p++)
        {
            dval = 1.0;
            ex = ey = 0;
            for (c=0; c<ncontext; c++)
                if (c==cx)
                    ex = *(powerst++);
                else if (c==cy)
                    ey = *(powerst++);
                else
                    for (e=*(powerst++); e--;)
                        dval *= context[c];
            for (d=0; d<naxis; d++, mat1++)
                /*------ Correct only celestial coordinates */
                if (d==lng)
                {
                    ival = pvindex[ey][ex];
                    projp[ival+100*d] += dval**mat1;
                }
                else if (d==lat)
                {
                    ival = pvindex[ex][ey];
                    projp[ival+100*d] += dval**mat1;
                }
        }
        free(powers);
    }

    free(context);

    /* Now take care of the field-dependent parameters */
    if (set->index2>=0)
    {
        /*-- Compute the array of exponents of polynom terms */
        powerst = powers = poly_powers(poly2);
        mat2 = mat + set->index2;	/* Field-dependent section in solution */
        for (p=0; p<ncoeff2; p++)
        {
            dval = 1.0;
            ex = *(powerst++);
            ey = *(powerst++);
            for (d=0; d<naxis; d++, mat2++)
                /*------ Correct only celestial coordinates */
                if (d==lng)
                {
                    ival = pvindex[ey][ex];
                    projp[ival+100*d] += *mat2;
                }
                else if (d==lat)
                {
                    ival = pvindex[ex][ey];
                    projp[ival+100*d] += *mat2;
                }
        }
        free(powers);
    }

    if (ncoeff <= naxis+1 && poly2->degree[0]<2)
        /* linear case: ditch the PV parameters and use only the CD/CRPIX parameters */
    {
        cd[lng*naxis + lng] = projp[1+100*lng]*wcs->cd[lng*naxis + lng]
            + projp[2+100*lng]*wcs->cd[lat*naxis + lng];
        cd[lng*naxis + lat] = projp[1+100*lng]*wcs->cd[lng*naxis + lat]
            + projp[2+100*lng]*wcs->cd[lat*naxis + lat];
        cd[lat*naxis + lng] = projp[1+100*lat]*wcs->cd[lat*naxis + lng]
            + projp[2+100*lat]*wcs->cd[lng*naxis + lng];
        cd[lat*naxis + lat] = projp[1+100*lat]*wcs->cd[lat*naxis + lat]
            + projp[2+100*lat]*wcs->cd[lng*naxis + lat];
        det = cd[lng*naxis + lng]*cd[lat*naxis + lat]
            - cd[lng*naxis + lat]*cd[lat*naxis + lng];
        if (fabs(det) > 1.0/BIG)
        {
            wcs->crpix[lng] += (cd[lng*naxis + lat]*projp[0+100*lat]
                    - cd[lat*naxis + lat]*projp[0+100*lng]) / det;
            wcs->crpix[lat] += (cd[lat*naxis + lng]*projp[0+100*lng]
                    - cd[lng*naxis + lng]*projp[0+100*lat]) / det;
            wcs->cd[lng*naxis + lng] = cd[lng*naxis + lng];
            wcs->cd[lng*naxis + lat] = cd[lng*naxis + lat];
            wcs->cd[lat*naxis + lng] = cd[lat*naxis + lng];
            wcs->cd[lat*naxis + lat] = cd[lat*naxis + lat];
#pragma ivdep
            for (p=100; p--;)
                projp[p+lng*100] = projp[p+lat*100] = 0.0;
        }
        else
            warning("determinant of the computed CD matrix is null: ",
                    "reverting to PV parameters");
    }

    /* Initialize other WCS structures */
    init_wcs(wcs);
    /* Find the range of coordinates */
    range_wcs(wcs);
    /* Invert projection corrections */
    invert_wcs(wcs);

    return;
}


/****** shrink_mat ************************************************************
  PROTO	void	shrink_mat(double *alpha, double *beta, int ncoefftot,
  int index, int nmiss)
  PURPOSE	Remove selected lines and columns from a matrix
  INPUT	Ptr to the alpha matrix,
  ptr to the beta matrix,
  matrix size,
  starting index,
  number of lines/columns to remove.
  OUTPUT	-.
  NOTES	Matrices are not reallocated.
  AUTHOR	E. Bertin (IAP)
  VERSION	09/01/2012
 ***/
void	shrink_mat(double *alpha, double *beta, int ncoefftot,
        int index, int nmiss)
{
    double	*a,*a2, *b,*b2;
    size_t	i,l, nelem, nelem2;

    if (nmiss >= ncoefftot)
        return;

#ifdef MATSTORAGE_PACKED
    /* First, remove lines */
    a = alpha + ((size_t)(2*ncoefftot-index+1)*index)/2;
    l = ((size_t)(2*ncoefftot-index-nmiss+1)*(index+nmiss))/2;
    a2 = alpha + l;
    nelem = ((size_t)(ncoefftot+1)*ncoefftot)/2-l;
    for (i=nelem; i--; )
        *(a++) = *(a2++);
    /* Remove columns */
    a2 = a = alpha + index;
    nelem2 = ncoefftot-nmiss-1;
    for (l=index; l--; nelem2--)
        for (i=nelem2, a2+=nmiss; i--;)
            *(a++) = *(a2++);
    if (index)
        for (i=nelem; i--;)
            *(a++) = *(a2++);
#else
    /* First, remove lines */
    a = alpha + index*(size_t)ncoefftot;
    a2 = a + nmiss*(size_t)ncoefftot;
    nelem = (ncoefftot-nmiss-index)*(size_t)ncoefftot;
    for (i=nelem; i--;)
        *(a++) = *(a2++);

    /* Remove columns */
    a = alpha + index;
    a2 = a + nmiss;
    nelem2 = (ncoefftot-nmiss);
    for (l=nelem2-1; l--; a2+=nmiss)
        for (i=nelem2; i--;)
            *(a++) = *(a2++);
    for (i=nelem2-index; i--;)
        *(a++) = *(a2++);
#endif


    /* The beta matrix */
    b = beta+index;
    b2 = b+nmiss;
    nelem = (ncoefftot-nmiss-index);
    for (i=nelem; i--;)
        *(b++) = *(b2++);

    return;
}


/****** regul_mat ************************************************************
PROTO	void	regul_mat(double *alpha, int ncoefftot)
PURPOSE	Apply Thikonov regularization to a normal equation matrix
INPUT	Ptr to an array of field groups,
	number of groups,
	ptr to the normal equation matrix,
	matrix size.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/04/2020
 ***/
void	regul_mat(fgroupstruct **fgroups, int nfgroup,
		double *alpha, int ncoefftot) {

   fieldstruct	**fields;
   double	alphamin, alphamin2, alphaminc, alphaminc2;
   size_t	stncoefftot, cm;
   int		c, f, g, index2, minindex2, nfield;

  stncoefftot = (size_t)ncoefftot;

// Identify starting index for field dependent parameters
  minindex2 = ncoefftot;
  for (g = 0; g < nfgroup; g++) {    
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f = 0; f < nfield; f++) {
      index2 = fields[f]->index2;
      if (index2 >= 0 && index2 < minindex2)
        minindex2 = index2;
    }
  }

  alphamin = alphamin2 = BIG;

#ifdef MATSTORAGE_PACKED

// Identify smallest constant non-zero diagonal terms
  cm = 0;
  for (c = 0; c < minindex2; c++) {
    if (alpha[cm] < alphamin && alpha[cm] > TINY)
      alphamin = alpha[cm];
    cm += stncoefftot - c;
  }
  for (; c < ncoefftot; c++) {
    if (alpha[cm] < alphamin2 && alpha[cm] > TINY)
      alphamin2 = alpha[cm];
    cm += stncoefftot - c;
  }

alphaminc = alphamin * ASTROM_REGULFACTOR;
alphaminc2 = alphamin2 * ASTROM_REGULFACTOR;

// Crude Tikhonov regularization
  if (alphamin < BIG) {
    cm = 0;
    for (c = 0; c < minindex2; c++) {
      alpha[cm] += alphaminc;
      cm += stncoefftot - c;
    }
  }
  if (alphamin2 < BIG) {
    cm = (minindex2 * (2 * stncoefftot + 1 - minindex2)) / 2;
    for (c = minindex2; c < ncoefftot; c++) {
      alpha[cm] += alphaminc2;
      cm += stncoefftot - c;
    }
  }

#else

// Identify smallest constant non-zero diagonal terms
  for (c = 0; c < minindex2; c++) {
    cm = c + c*stncoefftot;
    if (alpha[cm] < alphamin && alpha[cm] > TINY)
      alphamin = alpha[cm];
  }
  for (c = minindex2; c < ncoefftot; c++) {
    cm = c + c*stncoefftot;
    if (alpha[cm] < alphamin2 && alpha[cm] > TINY)
      alphamin2 = alpha[cm];
  }

alphaminc = alphamin * ASTROM_REGULFACTOR;
alphaminc2 = alphamin2 * ASTROM_REGULFACTOR;

// Crude Tikhonov regularization
  if (alphamin < BIG)
    for (c = 0; c < minindex2; c++)
      alpha[c + c*stncoefftot] += alphaminc;
  if (alphamin2 < BIG)
    for (c = minindex2; c < ncoefftot; c++)
      alpha[c + c*stncoefftot] += alphaminc2;

#endif

  return;
}

/****** compute jacobian ******************************************************
PROTO	void compute_jacobian(samplestruct *sample)
PURPOSE Compute the Jacobian matrices required for solving the astrometric
	equations.
INPUT	ptr to the sample structure,
	ptr to the dproj/dred Jacobian (or NULL if no computation needed).
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	Memory must have been allocated (naxis*naxis*sizeof(double)) for the
	Jacobian array.
AUTHOR	E. Bertin (IAP)
VERSION 03/03/2021
 ***/
int	compute_jacobian(samplestruct *sample) {

   wcsstruct	*wcs;
   fgroupstruct *fgroup;
   fieldstruct	*field;
   setstruct	*set;
   double	rawpos[NAXIS], rawpos2[NAXIS], redpos[NAXIS], wcspos[NAXIS],
		projpos1[NAXIS], projpos2[NAXIS];
   float	*dprojdred;
   int		i,j, naxis;


    if (!(set=sample->set) || !(field=set->field) || !(fgroup=field->fgroup)
            || !(wcs=fgroup->wcs))
        return RETURN_ERROR;

    naxis = wcs->naxis;
    for (i=0; i<naxis; i++)
        rawpos[i] = sample->vrawpos[i];
    dprojdred = sample->dprojdred;
    for (i=0; i<naxis; i++) {
//---- Compute terms of the Jacobian dproj/dred matrix
//---- Things are a bit tortuous because it is computed with respect
//---- to reduced coordinates ("in projected degrees")
      raw_to_red(set->wcs, rawpos, redpos);
      redpos[i] -= 0.5 * ARCSEC/DEG;
      red_to_raw(set->wcs, redpos, rawpos2);
      raw_to_wcs(set->wcs, rawpos2, wcspos);
      wcs_to_raw(wcs, wcspos, projpos1);
      redpos[i] += 1.0 * ARCSEC/DEG;
      red_to_raw(set->wcs, redpos, rawpos2);
      raw_to_wcs(set->wcs, rawpos2, wcspos);
      wcs_to_raw(wcs, wcspos, projpos2);
      for (j=0; j<naxis; j++)
        *(dprojdred++) = 1.0 * (projpos2[j] - projpos1[j])/(1.0 * ARCSEC/DEG);
    }

  return RETURN_OK;
}


/****** reproj_fgroup *******************************************************
PROTO	void reproj_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
		int flags)
PURPOSE Perform re-projection of sources in a group of fields.
INPUT	pointer to the group of fields,
	pointer to the reference field,
	proper motion correction flag.
OUTPUT	-.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION 03/03/2021
 ***/
void	reproj_fgroup(fgroupstruct *fgroup, fieldstruct *reffield, int flags) {

   fieldstruct		**field;
   wcsstruct		*wcs;
   setstruct		**pset,
			*set;
   msamplestruct	*msamp;
   samplestruct	*samp;
   double		wcspos[NAXIS],
			*projposmin,*projposmax,
			prop2, properr2, rprop;
   int			i,f,s, lng,lat, nsamp, nfield, naxis,
   			propflag, jacflag;


  propflag = flags & REPROJ_PROPER_MOTION;
  jacflag = flags & REPROJ_JACOBIAN;
  field = fgroup->field;
  nfield = fgroup->nfield;
  naxis = fgroup->naxis;
  wcs = fgroup->wcs;
  projposmin = fgroup->projposmin;
  projposmax = fgroup->projposmax;
  for (i=0; i<naxis; i++) {
    projposmin[i] = projposmin[i] = BIG;
    projposmax[i] = projposmax[i] = -BIG;
  }

  for (f=0; f<nfield; f++) {
    pset = field[f]->set;
    for (s=field[f]->nset; s--; pset++) {
      set = *pset;
      lng = set->lng;
      lat = set->lat;
      samp = set->sample;
      for (nsamp=set->nsample; nsamp--; samp++) {
        raw_to_wcs(set->wcs, samp->rawpos, samp->wcspos);
        if (propflag && samp->msamp) {
//-------- Correct for proper motions if asked to
          msamp = samp->msamp;
          for (i=0; i<naxis; i++) {
            wcspos[i] = samp->wcspos[i];
            if ((properr2 =
             	msamp->wcsprop[i]*msamp->wcsprop[i]*msamp->wcschi2) > TINY
		&& (prop2=msamp->wcsprop[i]*msamp->wcsprop[i]) > TINY) {
              rprop = msamp->wcsprop[i] * prop2 / (prop2 + properr2);
              wcspos[i] += (msamp->epoch - samp->epoch)
			* (i==lng?rprop/cos(samp->wcspos[lat]*DEG) : rprop);
            }
          }
          wcs_to_raw(set->wcs, wcspos, samp->vrawpos);
//-------- Compute new projection
          wcs_to_raw(wcs, wcspos, samp->projpos);
        } else {
          for (i=0; i<naxis; i++)
            samp->vrawpos[i] = samp->rawpos[i];
//---------- Compute new projection
          wcs_to_raw(wcs, samp->wcspos, samp->projpos);
        }
        if (jacflag)
          compute_jacobian(samp);
      }
//---- Compute the boundaries of projected positions for the whole group
      sort_samples(set);
      for (i=0; i<naxis; i++) {
        if (projposmin[i] > set->projposmin[i])
          projposmin[i] = set->projposmin[i];
        if (projposmax[i] < set->projposmax[i])
          projposmax[i] = set->projposmax[i];
      }
    }
  }

  if (reffield) {
    set = reffield->set[0];
    samp = set->sample;
    for (nsamp=set->nsample; nsamp--; samp++) {
      if (propflag && samp->msamp) {
//------ Correct for proper motions if asked to */
        msamp = samp->msamp;
        for (i=0; i<naxis; i++) {
          wcspos[i] = samp->wcspos[i];
          if ((properr2 =
          	msamp->wcsprop[i]*msamp->wcsprop[i]*msamp->wcschi2) > TINY
			&& (prop2=msamp->wcsprop[i]*msamp->wcsprop[i]) > TINY) {
            rprop = msamp->wcsprop[i] * prop2 / (prop2 + properr2);
            wcspos[i] += (msamp->epoch - samp->epoch)
		* (i==lng?rprop/cos(samp->wcspos[lat]*DEG) : rprop);
          }
        }
        wcs_to_raw(wcs, wcspos, samp->projpos);
      } else
        wcs_to_raw(wcs, samp->wcspos, samp->projpos);
    if (jacflag)
      compute_jacobian(samp);
    }
  }

  return;
}


/****** astrweight_fgroups ****************************************************
  PROTO	void astrweight_fgroups(fgroupstruct **fgroups, int nfgroup)
  PURPOSE	Compute appropriate, relative weights of detection sets prior to solving
  the astrometry.
  INPUT	ptr to an array of group of fields pointers,
  number of groups.
  OUTPUT	-.
  NOTES	Uses the global preferences. Input structures must have gone through
  crossid_fgroup() first.
  AUTHOR	E. Bertin (IAP)
  VERSION	30/10/2006
 ***/
void	astrweight_fgroups(fgroupstruct **fgroups, int nfgroup)
{
    fieldstruct	**fields,
                *field;
    setstruct	*set;
    samplestruct	*samp, *samp2;
    int		f,g,n,s, nfield,nref;

    /* Go through all groups, all fields, all sets, all samples */
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
                nref = 0;
                for (n=set->nsample; n--; samp++) {
                    for (samp2=samp;(samp2=samp2->prevsamp);) {
                        if (samp2->set->field->astromlabel<0)
                            nref++;
                    }
                }
                set->weightfac = nref? (double)set->nsample / nref : 1.0;
            }
        }
    }

    return;
}


/****** astr_orthopoly ********************************************************
  PROTO	void astr_orthopoly(polystruct *poly)
  PURPOSE	Orthonormalize the polynomial basis over the range of possible contexts.
  INPUT	pointer to poly structure.
  OUTPUT  -.
  NOTES   -.
  AUTHOR  E. Bertin (IAP)
  VERSION 16/12/2013
 ***/
void	astr_orthopoly(polystruct *poly)
{
    samplestruct	*sample;
    double	pos[POLY_MAXDIM], posmin[POLY_MAXDIM],
    *basis,*data,*datat,
    norm, step;
    int		count[POLY_MAXDIM],
    c,d,i,n, ndim, ncoeff, ndata, nrange;

    ncoeff = poly->ncoeff;
    ndim = poly->ndim;
    ndata = 1024;
    nrange = (int)pow((double)ndata, 1.0/ndim);
    step = 2.0 / (nrange - 1.0);
    norm = 1.0/sqrt(ndata);

    /* Go through each sample */
    ndata = 1; 
    for (d=0; d<ndim; d++)
    {
        pos[d] = posmin[d] = -1.0;
        ndata *= nrange;
        count[d] = 0;
    }

    QMALLOC(data, double, ndata*ncoeff);

    for (n=0; n<ndata; n++)
    {
        /*-- Get the local context coordinates */
        poly_func(poly, pos);
        basis = poly->basis;
        datat = data + n;
        /*-- Fill basis matrix as a series of row vectors */
#pragma ivdep
        for (c=ncoeff; c--; datat+=ndata)
            *datat = *(basis++)*norm;
        for (d=0; d<ndim; d++)
        {
            pos[d] += step;
            if (++count[d]<nrange)
                break;
            else
            {
                count[d] = 0;       /* No need to initialize it to 0! */
                pos[d] = posmin[d];
            }
        }
    }

    poly_initortho(poly, data, ndata);
    free(data);

    return;
}



