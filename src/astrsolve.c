/*
*				astrsolve.c
*
* Compute the "global" astrometric solution.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		27/03/2012
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
#include "astrsolve.h"
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "prefs.h"
#include "samples.h"
#include "wcs/poly.h"

#ifdef USE_THREADS
#include "threads.h"
#endif

#ifdef HAVE_ATLAS
#include ATLAS_LAPACK_H
#endif

#ifdef HAVE_LAPACKE
#include LAPACKE_H
#define MATSTORAGE_PACKED 1
#endif

/*------------------- global variables for multithreading -------------------*/
#ifdef USE_THREADS
#define	NMAT_MUTEX	32
#define	NMAT_MUTEX2	32
 pthread_t		*thread;
 pthread_mutex_t	matmutex[NMAT_MUTEX],matmutex2[NMAT_MUTEX2],
			fillastromutex, nconstastromutex;
 threads_gate_t		*pthread_startgate, *pthread_stopgate;
 fgroupstruct		**pthread_fgroups;
 double			*pthread_alpha, *pthread_beta;

 int			*pthread_nconst,
			pthread_endflag, pthread_ngroup, pthread_ncoefftot,
			pthread_gindex, pthread_findex, pthread_sindex;

 void		pthread_fill_astromatrix(fgroupstruct **fgroups, int ngroup,
			int ncoefftot, int *nconst,
			double *alpha, double *beta),
		*pthread_fillastromatrix_thread(void *arg);
#endif

 static void	fill_astromatrix(setstruct *set, double *alpha, double *beta,
			int ncoefftot,
			polystruct *poly, polystruct *poly2,
			int *nconst),
		add_alphamat(double *alpha, int naxis, int ncoefftot,
			double weight, int *ci, double *cv, int ncoeff,
			int *ci2, double *cv2, int ncoeff2),
		add_betamat(double *beta, int naxis,
		double weight, int *ci, double *cv, double *cc, int ncoeff);

/****** astrsolve_fgroups *****************************************************
PROTO	void astrsolve_fgroups(fgroupstruct **fgroups, int nfgroup)
PURPOSE	Compute a global astrometric solution among a group of fields.
INPUT	ptr to an array of group of fields pointers,
	number of groups.
OUTPUT	-.
NOTES	Uses the global preferences. Input structures must have gone through
	crossid_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	19/12/2011
 ***/
void	astrsolve_fgroups(fgroupstruct **fgroups, int nfgroup)
  {
   fieldstruct		**fields,**fields2,
			*field,*field2;
   setstruct		*set,*set2;
   samplestruct		*samp;
   polystruct		*poly, *poly2;
   char			str[64],
			**contextname;
   double		cmin[MAXCONTEXT],cmax[MAXCONTEXT],
			cscale[MAXCONTEXT], czero[MAXCONTEXT],
			*alpha, *beta,
			dval;
   int			group2[NAXIS],
			*findex, *findex2, *nsetmax, *nconst,*nc,*nc2,
			c,cm,f,f2,g,g2,i,n,s,s2,sm, nfield,nfield2,
			ninstru, instru, naxis, ncoeff, ncoeff2, ncoefftot,
			npcoeff, npcoeff2,
			d, index, index2,
			ncontext, cx,cy, nicoeff, nicoeff2, groupdeg2,
			startindex2, nmiss;


/* Compute weight factors for each set based on the relative number of */
/* astrometric references and detections */
  NFPRINTF(OUTPUT, "Initializing detection weight factors...");
  astrweight_fgroups(fgroups, nfgroup);

  NFPRINTF(OUTPUT, "Initializing the global astrometry matrix...");

  naxis = fgroups[0]->naxis;

/* CONTEXT polynomial */
  poly = poly_init(prefs.context_group, prefs.ncontext_name, prefs.group_deg,
	prefs.ngroup_deg);

/* Linear field-dependent polynomial */
  for (d=0; d<naxis; d++)
    group2[d] = 1;
  groupdeg2 = 1;
  poly2 = poly_init(group2, naxis, &groupdeg2, 1);

/* Use a different index for each instrument */
  ninstru = prefs.nastrinstrustr;
  if (!ninstru)
    ninstru = 1;
  QCALLOC(findex, int, ninstru+1);
  QCALLOC(findex2, int, ninstru+1);
  QCALLOC(nsetmax, int, ninstru+1);

/* Compute the total number of fields and the max number of sets per field */
  ncontext = fgroups[0]->field[0]->set[0]->ncontext;
  contextname = fgroups[0]->field[0]->set[0]->contextname;

  for (c=0; c<ncontext; c++)
    {
    cmin[c] = BIG;
    cmax[c] = -BIG;
    }

  for (g=0 ; g<nfgroup; g++)
    {
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f=0; f<nfield; f++)
      {
      field=fields[f];
      instru = field->astromlabel;
/*---- Index current field */
      if (field->nset > nsetmax[instru])
        nsetmax[instru] = field->nset;
      if (prefs.stability_type[instru] == STABILITY_EXPOSURE)
        {
        field->index = findex[instru+1];
        findex[instru+1] += field->nset;
        field->index2 = -1;
        findex2[instru+1] = 0;
        }
      else
        {
        field->index2 = findex2[instru+1]++ - 1;
/*
field->index2 = -1;
*/
        findex[instru+1] = nsetmax[instru];
        }
/*---- Index sets */
      for (s=0; s<field->nset; s++)
        {
        set = field->set[s];
        set->index = s;
/*------ Compute normalization factors for context values */
        samp = set->sample;
        for (n=set->nsample; n--;)
          {
          for (c=0; c<ncontext; c++)
            {
            dval = samp->context[c];
/*---------- Update min and max */
            if (dval<cmin[c])
              cmin[c] = dval;
            if (dval>cmax[c])
              cmax[c] = dval;
            }
          }
        }
      }
    }

/* Compute context rescaling coefficients (to avoid dynamic range problems) */
  cx = cy = -1;
  for (c=0; c<ncontext; c++)
    {
    czero[c] = cmin[c];
    cscale[c] =  (dval = cmax[c] - cmin[c]) != 0.0 ? 1.0/dval : 1.0;
    if (contextname)
      {
/*---- Identify X_IMAGE and Y_IMAGE parameters (will be processed separately)*/
      if (!strcmp(contextname[c], "X_IMAGE")
	|| !strcmp(contextname[c], "XWIN_IMAGE")
	|| !strcmp(contextname[c], "XPSF_IMAGE")
	|| !strcmp(contextname[c], "XPEAK_IMAGE"))
        cx = c;
      else if (!strcmp(contextname[c], "Y_IMAGE")
	|| !strcmp(contextname[c], "YWIN_IMAGE")
	|| !strcmp(contextname[c], "YPSF_IMAGE")
	|| !strcmp(contextname[c], "YPEAK_IMAGE"))
        cy = c;
      }
    }

/* Compute the required number of coefficients */
  npcoeff = poly->ncoeff;
  ncoeff = npcoeff*naxis;
  nicoeff = ncoeff*naxis;
  npcoeff2 = poly2->ncoeff;
  ncoeff2 = npcoeff2*naxis;
  nicoeff2 = ncoeff2*naxis;

  ncoefftot = startindex2 = 0;
  for (i=0; i<ninstru; i++)
    {
    findex2[i+1]--;	/* Remove first field (too many degrees of freedom) */
/*-- Set the extra exposure-dependent parameter flag */
    ncoefftot += ncoeff*findex[i+1];
    startindex2 += ncoeff*findex[i+1];
    if (findex2[i+1]>0)
      ncoefftot += ncoeff2*findex2[i+1];
    if (i)
      {
      findex[i] += findex[i-1];
      findex2[i] += findex2[i-1];
      }
    }

/* Compute the "definitive" index for each set */
  for (g=0 ; g<nfgroup; g++)
    {
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f=0; f<nfield; f++)
      {
      field=fields[f];
      instru = field->astromlabel;
      index = ncoeff*(field->index+findex[instru]);
      index2 = (field->index2 <0)? -1
		: startindex2 + ncoeff2*(field->index2+findex2[instru]);
      field->index2 = index2;
      for (s=0; s<field->nset; s++)
        {
        set = field->set[s];
        set->index = index + set->index*ncoeff;
        set->index2 = index2;
        set->nconst = 0;
        set->ncoeff = ncoeff;
/*------ Take the opportunity to set the context scales and offset */
        for (c=0; c<ncontext; c++)
          {
          set->contextoffset[c] = czero[c];
          set->contextscale[c] = cscale[c];
          }
        set->contextx = cx;
        set->contexty = cy;
        }
      }
    }

  QCALLOC(nconst, int, ncoefftot);

#ifdef MATSTORAGE_PACKED
  QCALLOC(alpha, double, ((size_t)(ncoefftot+1)*ncoefftot)/2);
#else
  QCALLOC(alpha, double, (size_t)ncoefftot*ncoefftot);
#endif
  QCALLOC(beta, double, ncoefftot);

#ifdef USE_THREADS
  pthread_fill_astromatrix(fgroups, nfgroup, ncoefftot, nconst, alpha, beta);
#else

/* Matrix parameters will be ordered as :   */
/* Local reduced coordinate (NAXISn)        */
/*   Instru               (ninstru at most) */
/*   (if stability_type == EXPOSURE)        */
/*     Field                                */
/*       Set                                */
/*         Polynom term                     */
/*           Projected coordinate (NAXISn)  */
/*   (else)                                 */
/*     Set                                  */
/*       Polynom term                       */
/*         Projected coordinate (NAXISn)    */
/*     Extra field                          */
/*       Field polynom term                 */
/*         Projected coordinate (NAXISn)    */
/*   (endif)                                */

/* Go through all groups, all fields, all sets, all samples */
  for (g=0; g<nfgroup; g++)
    {
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f=0; f<nfield; f++)
      {
      sprintf(str, "Filling the global astrometry matrix: "
		"group %d/%d, field %d/%d",
	g+1, nfgroup, f+1, nfield);
      NFPRINTF(OUTPUT, str);
      field = fields[f];
      for (s=0; s<field->nset; s++)
        {
        set = field->set[s];
        fill_astromatrix(set, alpha, beta, ncoefftot, poly,poly2, nconst);
        }
      }
    }

#endif

/* Check that there are enough constraints for each set... */
/* And remove unconstrained degrees of freedom if necessary */
/* Warning: ordering of instructions important to avoid overwriting indices! */
/* First: exposure-independent coefficients */
  NFPRINTF(OUTPUT, "Fixing the degrees of freedom");
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
        index = set->index;
        if (nconst[index] < set->ncoeff)
	  {
          if (field->nset > 1)
	    sprintf(str,"set #%d of instrument A%d ",s+1,field->astromlabel+1);
          else
            sprintf(str, "instrument A%d ", field->astromlabel+1);
          warning("Not enough matched detections in ", str);
          nmiss = set->ncoeff - nconst[index];
/*-------- Trim the matrices */
          shrink_mat(alpha,beta, ncoefftot, index+nconst[index], nmiss);
/*-------- Shift all other indexes above */
          nc = nconst+index+nconst[index];
          nc2 = nc + nmiss;
          for (i=ncoefftot-index-set->ncoeff; i--;)
	    *(nc++) = *(nc2++);
          for (g2=0; g2<nfgroup; g2++)
            {    
            nfield2 = fgroups[g2]->nfield;
            fields2 = fgroups[g2]->field;
            for (f2=0; f2<nfield2; f2++)
              {
              field2 = fields2[f2];
              if (field2->index2>=0)
                field2->index2 -= nmiss;
              for (s2=0; s2<field2->nset; s2++)
	        {
                set2 = field2->set[s2];
                if (set2->index2>=0)
                  set2->index2 -= nmiss;
                if (set2->ncoeff && set2->index == index)
		  set2->ncoeff -= nmiss;
                if (set2->index > index)
                  set2->index -= nmiss;
                }
              }
            }
          ncoefftot -= nmiss;
          }
        }
      }
    }

/* Second: exposure dependent coefficients */
  for (g=0; g<nfgroup; g++)
    {    
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f=0; f<nfield; f++)
      {
      field = fields[f];
      index2 = field->index2;
      if (index2 >= 0)
        {
        for (c=s=0; s<field->nset; s++)
          {
          set = field->set[s];
          if (set->ncoeff)
            c += nconst[set->index] - set->ncoeff;
          }
        if (c < ncoeff2)
          {
          nmiss = ncoeff2;
/*-------- Trim the matrices */
          shrink_mat(alpha,beta, ncoefftot, index2, nmiss);
/*-------- Shift all other indexes above */
          nc = nconst+index2;
          nc2 = nc + nmiss;
          for (i=ncoefftot-nmiss-index2; i--;)
            *(nc++) = *(nc2++);
          for (g2=0; g2<nfgroup; g2++)
            {    
            nfield2 = fgroups[g2]->nfield;
            fields2 = fgroups[g2]->field;
            for (f2=0; f2<nfield2; f2++)
              {
              field2 = fields2[f2];
              if (field2->index2>=0 && field2->index2 > index2)
                field2->index2 -= nmiss;
              for (s2=0; s2<field2->nset; s2++)
                {
                set2 = field2->set[s2];
                if (set2->index2 == index2)
                  set2->index2 = -1;
                else if (set2->index2 > index2)
                  set2->index2 -= nmiss;
                }
              }
            }
          ncoefftot -= nmiss;
          field->index2 = -1;
          }
        else
	  {
/*-------- Remove ncoeff2 constraints in a "balanced" way */
          for (n=ncoeff2; n-=naxis;)
            {
            for (sm=-1,cm=s=0; s<field->nset; s++)
              {
              set = field->set[s];
              if (set->ncoeff && (c = nconst[set->index] - set->ncoeff) > cm)
	        {
                cm = c;
                sm = s;
                }
              }
            if (sm<0)
              error(EXIT_FAILURE, "*Internal Error*: mismatched d.o.f. in ",
		"astrsolve_fgroups()");
            nconst[field->set[sm]->index] -= naxis;
            }
          }
        }
      }
    }

/* Solve! */
  NFPRINTF(OUTPUT, "Solving the global astrometry matrix...");

  if (ncoefftot)
    {
#if defined(HAVE_LAPACKE)
  #ifdef MATSTORAGE_PACKED
    if (LAPACKE_dppsv(LAPACK_COL_MAJOR, 'L', ncoefftot, 1, alpha, beta, ncoefftot) !=0)
  #else
    if (LAPACKE_dposv(LAPACK_COL_MAJOR, 'L', ncoefftot, 1, alpha, ncoefftot,
		beta, ncoefftot) !=0)
  #endif
#else
    if (clapack_dposv(CblasRowMajor, CblasUpper,
		ncoefftot, 1, alpha, ncoefftot, beta, ncoefftot) != 0)
#endif
      warning("Not a positive definite matrix", " in astrometry solver");
    }    

/* Fill the astrom structures with the derived parameters */
  for (g=0; g<nfgroup; g++)
    {    
    nfield = fgroups[g]->nfield;
    fields = fgroups[g]->field;
    for (f=0; f<nfield; f++)
      {
      field=fields[f];
      for (s=0; s<field->nset; s++)
        {
        set = field->set[s];
        mat_to_wcs(poly, beta, set);
        }
      }
   }

  poly_end(poly);
  poly_end(poly2);
  free(alpha);
  free(beta);
  free(nconst);
  free(findex);
  free(findex2);
  free(nsetmax);

  return;
  }


/****** fill_astromatrix ******************************************************
PROTO	void fill_astromatrix(setstruct *set, double *alpha, double *beta,
			int ncoefftot,
			polystruct *poly, polystruct *poly2, int *nconst)
PURPOSE	Fill the normal equation matrix for astrometry with elements from a
	given set.
INPUT	Ptr to the set structure,
	ptr to the alpha matrix,
	ptr to the beta matrix.
	matrix size,
	pointer to the polynomial structure (steady ;
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	31/01/2011
 ***/
void	fill_astromatrix(setstruct *set, double *alpha, double *beta,
			int ncoefftot,
			polystruct *poly, polystruct *poly2, int *nconst)
  {
   fieldstruct		*field2;
   setstruct		*set2;
   samplestruct		*samp,*samp2,*samp3;
   double		dprojdred[NAXIS*NAXIS], context[MAXCONTEXT],
			redpos[NAXIS],
			*coeffval, *coeffval2, *coeffconst,*coeffconst2,
			*cv, *cvo,*cvoa, *cv2,*cvo2,*cvoa2,
			*cc,*cco,*ccoa, *cc2,*cco2,*ccoa2, *x, *jac,*jact,
			*basis,*basis2, *czero, *cscale,
			sigma2,sigma3, weight, weightfac;
   int			*coeffindex,*coeffindex2,
			*ci,*cio,*cioa, *ci2,*cio2,*cioa2,
			nlsamp,nlsampmax, nscoeffmax, nscoeffmax2, instru,
			index,indext,index2,indext2, redflag, naxis, ncontext,
			npcoeff,ncoeff,nicoeff, npcoeff2,ncoeff2,nicoeff2,
			nsamp,
			cx,cy, c, d, d1,d2, p, lng,lat;

#ifdef USE_THREADS
   int			imut=0;
#endif

  naxis = set->naxis;
  npcoeff = poly->ncoeff;
  ncoeff = npcoeff*naxis;
  nicoeff = ncoeff*naxis;
  npcoeff2 = poly2->ncoeff;
  ncoeff2 = npcoeff2*naxis;
  nicoeff2 = ncoeff2*naxis;
  
  nlsampmax = 0;
  ncontext = set->ncontext;
  coeffval = coeffval2 = coeffconst = coeffconst2 = (double *)NULL;
  coeffindex = coeffindex2 = (int *)NULL;      /* To avoid gcc -Wall warnings*/
  czero = set->contextoffset;
  cscale = set->contextscale;
  cx = set->contextx;
  cy = set->contexty;
  lng = set->lng;
  lat = set->lat;
  samp = set->sample;
  for (nsamp=set->nsample; nsamp--; samp++)
    {
/*-- Care only about one end of the series of linked detections */
    if (samp->prevsamp && !samp->nextsamp)
      {
      samp2 = samp;
/*---- Count the samples in the pair-list */
      for (nlsamp=1; (samp2=samp2->prevsamp); nlsamp++);
/*---- Allocate memory for storing list sample information */
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
        nscoeffmax = nicoeff*nlsamp;
        QMALLOC(coeffindex, int, nscoeffmax);
        QMALLOC(coeffval, double, nscoeffmax);
        QMALLOC(coeffconst, double, nscoeffmax);
        nscoeffmax2 = nicoeff2*nlsamp;
        QMALLOC(coeffindex2, int, nscoeffmax2);
        QMALLOC(coeffval2, double, nscoeffmax2);
        QMALLOC(coeffconst2, double, nscoeffmax2);
        }
/*---- Fill list-sample coefficients */
      ci = coeffindex;
      ci2 = coeffindex2;
      cv = coeffval;
      cv2 = coeffval2;
      cc = coeffconst;
      cc2 = coeffconst2;
      for (samp2 = samp; samp2; samp2=samp2->prevsamp)
        {
        cio = ci;
        cio2 = ci2;
        cvo = cv;
        cvo2 = cv2;
        cco = cc;
        cco2 = cc2;
        set2 = samp2->set;
        field2 = set2->field;
        instru = field2->astromlabel;
/*------ Negative indices are for the reference field */
        if (instru>=0)
	  {
          index = set2->index;
          index2 = set2->index2;
#ifdef USE_THREADS
          QPTHREAD_MUTEX_LOCK(&nconstastromutex);
#endif
          set2->nconst += naxis;
          nconst[index]+=naxis;
          if (index2>=0)
            nconst[index2]+=naxis;
#ifdef USE_THREADS
          QPTHREAD_MUTEX_UNLOCK(&nconstastromutex);
#endif
          }
        else
          index = index2 = -1;
/*------ Compute the value of the basis functions at each sample */
        if (instru >= 0)
	  {
          redflag = 0;
          for (c=0; c<ncontext; c++)
            {
            if (c==cx || c==cy)
              {
/*------------ Image coordinates receive a special treatment */
              if (!redflag)
                {
                raw_to_red(set2->wcs, samp2->rawpos, redpos);
                redflag = 1;
                }
              context[c] = redpos[c==cx?lng:lat];
              }
            else
              context[c] = (samp2->context[c]-czero[c])*cscale[c];
            }
/*-------- The basis functions are computed from reduced context values*/
          poly_func(poly, context);
/*-------- Compute the local Jacobian matrix */
          compute_jacobian(samp2,dprojdred);
          if (index2>=0)
            {
            if (!redflag)
              raw_to_red(set2->wcs, samp2->rawpos, redpos);
              poly_func(poly2, redpos);
            }
          }

/*------ Pre-compute numbers that will be used in the beta matrix  */
/*------ The innermost index of the jac array is the numerator of the */
/*------ Jacobian operator */
        x = samp2->projpos;
/*------ Fill something equivalent to a row of a design matrix */
        for (d2=0; d2<naxis; x++, d2++)
	  {
          indext = index;
          basis = poly->basis;
          jac = dprojdred+d2;
          for (p=npcoeff; p--; basis++)
            {
            jact = jac;
            for (d1=naxis; d1--; jact+=naxis)
              {
              *(ci++) = indext;
              *(cv++) = *jact**basis;
              *(cc++) = *x;
              if (!(++indext))
                indext--;
              }
            }
/*-------- Now the extra field-dependent parameters */
          indext2 = index2;
          basis2 = poly2->basis;
          for (p=npcoeff2; p--; basis2++)
            {
            jact = jac;
            for (d1=naxis; d1--; jact+=naxis)
              {
              *(ci2++) = indext2;
              *(cv2++) = *jact**basis2;
              *(cc2++) = *x;
              if (!(++indext2))
                indext2--;
              }
            }
          }
/*------ Drop it if the object is saturated or truncated */
        if ((samp2->sexflags & (OBJ_SATUR|OBJ_TRUNC)))
          continue;
/*------ Compute the relative positional variance for this source */
        sigma2 = 0.0;
        for (d=0; d<naxis; d++)
          sigma2 += samp2->wcsposerr[d]*samp2->wcsposerr[d];
/*------ Add extra-weight to the reference to compensate for all the overlaps */
        weightfac = nlsamp*prefs.astref_weight;
/*------ Now fill the matrices */
        cvoa = coeffval;
        cioa = coeffindex;
        ccoa = coeffconst;
        cvoa2 = coeffval2;
        cioa2 = coeffindex2;
        ccoa2 = coeffconst2;
        for (samp3 = samp; samp3 != samp2; samp3=samp3->prevsamp)
          {
/*------ Drop it if the object is saturated or truncated */
          if (!(samp3->sexflags & (OBJ_SATUR|OBJ_TRUNC)))
            {
/*---------- Compute the relative weight of the pair */
            sigma3 = 0.0;
            for (d=0; d<naxis; d++)
              sigma3 += samp3->wcsposerr[d]*samp3->wcsposerr[d];
            weight = (instru<0 ? samp3->set->weightfac*weightfac : 1.0)
			/(sigma3+sigma2);
/*
weight = sqrt(weight);
*/
/*---------- Fill the matrices */
/*---------- First, the contribution from the 1st detection itself */
#ifdef USE_THREADS
            if (*cio>=0)
              {
              imut = (*cio/ncoeff)%NMAT_MUTEX;
              QPTHREAD_MUTEX_LOCK(&matmutex[imut]);
              }
#endif
            add_alphamat(alpha,naxis,ncoefftot,weight,
		cio,cvo,ncoeff, cio,cvo,ncoeff);
            add_alphamat(alpha,naxis,ncoefftot,weight,
		cio,cvo,ncoeff, cio2,cvo2,ncoeff2);
            if (*cioa>=*cio)
              add_alphamat(alpha,naxis,ncoefftot,-weight,
		cio,cvo,ncoeff, cioa,cvoa,ncoeff);
            add_alphamat(alpha,naxis,ncoefftot,-weight,
		cio,cvo,ncoeff, cioa2,cvoa2,ncoeff2);
            add_betamat(beta,naxis,weight, cio,cvo,cco,ncoeff);
            add_betamat(beta,naxis,-weight, cio,cvo,ccoa,ncoeff);
#ifdef USE_THREADS
            if (*cio>=0)
              {
              QPTHREAD_MUTEX_UNLOCK(&matmutex[imut]);
              }
            if (*cio2>=0)
              {
              imut = (*cio2/ncoeff)%NMAT_MUTEX2;
              QPTHREAD_MUTEX_LOCK(&matmutex2[imut]);
              }
#endif
/*
            add_alphamat(alpha,naxis,ncoefftot,weight,
		cio2,cvo2,ncoeff2, cio,cvo,ncoeff);
*/
            add_alphamat(alpha,naxis,ncoefftot,weight,
		cio2,cvo2,ncoeff2, cio2,cvo2,ncoeff2);
/*
            add_alphamat(alpha,naxis,ncoefftot,-weight,
		cio2,cvo2,ncoeff2,  cioa,cvoa,ncoeff);
*/
            if (*cioa2>=*cio2)
              add_alphamat(alpha,naxis,ncoefftot,-weight,
		cio2,cvo2,ncoeff2,  cioa2,cvoa2,ncoeff2);
            add_betamat(beta,naxis,weight, cio2,cvo2,cco2,ncoeff2);
            add_betamat(beta,naxis,-weight, cio2,cvo2,ccoa2,ncoeff2);
/*---------- Then, the contribution from the 2nd detection itself */
#ifdef USE_THREADS
            if (*cio2>=0)
              {
              QPTHREAD_MUTEX_UNLOCK(&matmutex2[imut]);
              }
            if (*cioa>=0)
              {
              imut = (*cioa/ncoeff)%NMAT_MUTEX;
              QPTHREAD_MUTEX_LOCK(&matmutex[imut]);
              }
#endif
            add_alphamat(alpha,naxis,ncoefftot,weight,
		cioa,cvoa,ncoeff, cioa,cvoa,ncoeff);
            add_alphamat(alpha,naxis,ncoefftot,weight,
		cioa,cvoa,ncoeff, cioa2,cvoa2,ncoeff2);
            if (*cio>=*cioa)
              add_alphamat(alpha,naxis,ncoefftot,-weight,
		cioa,cvoa,ncoeff, cio,cvo,ncoeff);
            add_alphamat(alpha,naxis,ncoefftot,-weight,
		cioa,cvoa,ncoeff, cio2,cvo2,ncoeff2);
            add_betamat(beta,naxis,weight, cioa,cvoa,ccoa,ncoeff);
            add_betamat(beta,naxis,-weight, cioa,cvoa,cco,ncoeff);
#ifdef USE_THREADS
            if (*cioa>=0)
              {
              QPTHREAD_MUTEX_UNLOCK(&matmutex[imut]);
              }
            if (*cioa2>=0)
              {
              imut = (*cioa2/ncoeff)%NMAT_MUTEX2;
              QPTHREAD_MUTEX_LOCK(&matmutex2[imut]);
              }
#endif
/*
            add_alphamat(alpha,naxis,ncoefftot,weight,
		cioa2,cvoa2,ncoeff2, cioa,cvoa,ncoeff);
*/
            add_alphamat(alpha,naxis,ncoefftot,weight,
		cioa2,cvoa2,ncoeff2, cioa2,cvoa2,ncoeff2);
            add_betamat(beta,naxis,weight, cioa2,cvoa2,ccoa2,ncoeff2);
/*---------- Then, the crossed contributions */
/*
            add_alphamat(alpha,naxis,ncoefftot,-weight,
		cioa2,cvoa2,ncoeff2, cio,cvo,ncoeff);
*/
            if (*cio2>=*cioa2)
              add_alphamat(alpha,naxis,ncoefftot,-weight,
		cioa2,cvoa2,ncoeff2, cio2,cvo2,ncoeff2);
            add_betamat(beta,naxis,-weight, cioa2,cvoa2,cco2,ncoeff2);
#ifdef USE_THREADS
            if (*cioa2>=0)
              {
              QPTHREAD_MUTEX_UNLOCK(&matmutex2[imut]);
              }
#endif
            }
          cvoa += nicoeff;
          cioa += nicoeff;
          ccoa += nicoeff;
          cvoa2 += nicoeff2;
          cioa2 += nicoeff2;
          ccoa2 += nicoeff2;
	  }
        }
      }
    }

  if (nlsampmax)
    {
    free(coeffindex);
    free(coeffval);
    free(coeffconst);
    free(coeffindex2);
    free(coeffval2);
    free(coeffconst2);
    }

  return;
  }


/****** add_alphamat ******************************************************
PROTO	void add_alphamat(double *alpha, int naxis, int ncoefftot,
		double weight, int *ci, double *cv, int ncoeff,
		int *ci2, double *cv2, int ncoeff2)
PURPOSE	Add a product contribution to a normal equation matrix.
INPUT	ptr to the alpha matrix,
	number of dimensions,
	total size of the matrices
	weight of the contribution (can be negative for cross-products),
	index array for the first component,
	basis array for the first component,
	constraint array for the first component,
	size of the first component arrays,
	index array for the second component,
	basis array for the second component,
	size of the second component arrays,
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	20/12/2011
 ***/
static void add_alphamat(double *alpha, int naxis, int ncoefftot,
		double weight, int *ci, double *cv, int ncoeff,
		int *ci2, double *cv2, int ncoeff2)
  {
   double	*cvo,
		weight2;
   size_t	rowp;
   int		*cio,
		d, p, q, lim;

  if (*ci < 0 || *ci2 < 0)
    return;

  for (d=naxis; d--; ci2+=ncoeff2, cv2+=ncoeff2)
    for (p=ncoeff; p--; ci++, cv++)
      {
      cio = ci2;
      cvo = cv2;
      weight2 = weight**cv;
#ifdef MATSTORAGE_PACKED
      rowp = ((size_t)(2*ncoefftot-*ci-1)**ci)/2;
#else
      rowp = (size_t)ncoefftot**ci;
#endif
      for (q=ncoeff2; q--; cio++,cvo++)
        if (*cio>=*ci)
          alpha[rowp+(size_t)*cio] += weight2**cvo;
      }

  return;
  }


/****** add_betamat ******************************************************
PROTO	void add_betabat(double *alpha, int naxis,
		double weight, int *ci, double *cv, double *cc, int ncoeff)
PURPOSE	Add a product contribution to the right-hand matrix of normal
	equations.
INPUT	ptr to the right-hand matrix,
	number of dimensions,
	weight of the contribution (can be negative for cross-products),
	index array,
	basis array,
	constraint array,
	size of the first component arrays.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	30/12/2003
 ***/
static void add_betamat(double *beta, int naxis,
		double weight, int *ci, double *cv, double *cc, int ncoeff)
  {
   int		d, p;

  if (*ci < 0)
    return;

  for (d=naxis; d--;)
    for (p=ncoeff; p--;)
      beta[*(ci++)] -= weight**(cv++)**(cc++);

  return;
  }


#ifdef USE_THREADS

/****** pthread_fillastromatrix_thread ***************************************
PROTO   void *pthread_fillastromatrix_thread(void *arg)
PURPOSE thread that takes care of fill the astrometry matrix.
INPUT   Pointer to the thread number.
OUTPUT  -.
NOTES   Relies on global variables.
AUTHOR  E. Bertin (IAP)
VERSION 12/01/2005
 ***/
void    *pthread_fillastromatrix_thread(void *arg)
  {
   polystruct	*poly,*poly2;
   char		str[128];
   int		group2[NAXIS],
		findex, gindex, sindex, proc, naxis, groupdeg2,
		d;

  gindex = findex = -1;
  proc = *((int *)arg);
  naxis = pthread_fgroups[0]->naxis;
  for (d=0; d<naxis; d++)
    group2[d] = 1;
  groupdeg2 = 1;
  poly = poly_init(prefs.context_group, prefs.ncontext_name,
			prefs.group_deg, prefs.ngroup_deg);
  poly2 = poly_init(group2, naxis, &groupdeg2, 1);
  threads_gate_sync(pthread_startgate);
  while (!pthread_endflag)
    {
    QPTHREAD_MUTEX_LOCK(&fillastromutex);
    if (pthread_gindex<pthread_ngroup)
      {
      gindex = pthread_gindex;
      findex = pthread_findex;
      sindex = pthread_sindex++;
      if (pthread_sindex>=pthread_fgroups[gindex]->field[findex]->nset)
        {
        pthread_sindex = 0;
        pthread_findex++;
        }
      if (pthread_findex>=pthread_fgroups[gindex]->nfield)
        {
        pthread_findex = 0;
        pthread_gindex++;
        }
      if (!sindex)
        {
        sprintf(str, "Filling the global astrometry matrix: "
		"group %d/%d, field %d/%d",
		gindex+1, pthread_ngroup,
		findex+1, pthread_fgroups[gindex]->nfield);
        NFPRINTF(OUTPUT, str);
        }
      QPTHREAD_MUTEX_UNLOCK(&fillastromutex);

/*---- Match field */
      fill_astromatrix(pthread_fgroups[gindex]->field[findex]->set[sindex],
		pthread_alpha, pthread_beta, pthread_ncoefftot,
		poly, poly2, pthread_nconst);

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


/****** pthread_fill_astromatrix **********************************************
PROTO   void pthread_fill_astromatrix(fgroupstruct **fgroups, int ngroup,
		int ncoefftot, int *nconst, double *alpha, double *beta)
PURPOSE	Start multithreaded filling of the astrometry matrix
INPUT	ptr to the array of groups of fields,
	number of groups involved,
	matrix size,
	pointer to an array of constraint counts,
	pointer of pointer to the normal equation alpha matrix (output),
	pointer of pointer to the beta vector matrix (output).
OUTPUT	-.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	08/03/2007
 ***/
void	pthread_fill_astromatrix(fgroupstruct **fgroups, int ngroup,
		int ncoefftot, int *nconst, double *alpha, double *beta)
  {
   static pthread_attr_t	pthread_attr;
   int				*proc,
				p, i;

/* Number of active threads */
  nproc = prefs.nthreads;
  pthread_fgroups = fgroups;
  pthread_ngroup = ngroup;
  pthread_ncoefftot = ncoefftot;
  pthread_nconst = nconst;
  pthread_alpha = alpha;
  pthread_beta = beta;
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
	&pthread_fillastromatrix_thread, &proc[p]);
    }

  QPTHREAD_MUTEX_LOCK(&fillastromutex);
  pthread_gindex = pthread_findex = pthread_sindex = 0;
  pthread_endflag = 0;
  QPTHREAD_MUTEX_UNLOCK(&fillastromutex);
  NFPRINTF(OUTPUT, "Starting field matching...");
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
PROTO	void mat_to_wcs(polystruct *poly, double *mat, setstruct *set)
PURPOSE	Fill a WCS structure according to the matrix solution.
INPUT	ptr to the polynom structure,
	ptr to the solution vector,
	ptr to the WCS structure to be modified.
OUTPUT	-.
NOTES	Unfortunately the present WCS description of distortions is valid
	in 2D only.
AUTHOR	E. Bertin (IAP)
VERSION	15/12/2005
 ***/
void mat_to_wcs(polystruct *poly, double *mat, setstruct *set)
  {
  static const int pvindex[8][8] =	{{0,1,4,7,12,17,24,31},
					{2,5,8,13,18,25,32,0},
					{6,9,14,19,26,33,0,0},
					{10,15,20,27,34,0,0,0},
					{16,21,28,35,0,0,0,0},
					{22,29,36,0,0,0,0,0},
					{30,37,0,0,0,0,0,0},
					{38,0,0,0,0,0,0,0}};

   wcsstruct	*wcs;
   samplestruct	*sample;
   double	cd[NAXIS*NAXIS],
		*context, *cscale,*czero, *projp, *mat1, *mat2,
		dval, det;
   int		*powers, *powerst,
		c,cx,cy, d, e,ex,ey, p,s, lng,lat, naxis,ncoeff,ncontext, ival;

  wcs = set->wcs;
  naxis = wcs->naxis;
  lng = wcs->lng;
  lat = wcs->lat;
  ncoeff = set->ncoeff/naxis;
  ncontext = poly->ndim;
  cscale = set->contextscale;
  czero = set->contextoffset;
  projp = wcs->projp;
  mat1 = mat + set->index;	/* Start of the PV section in the solution */

/* Compute the array of exponents of polynom terms */
  powers = poly_powers(poly);

/* Precompute FITS or average parameters */
  QMALLOC(context, double, ncontext);
  cx = cy = -1;
  for (c=0; c<ncontext; c++)
    {
    if (!strcmp(set->contextname[c], "X_IMAGE")
	|| !strcmp(set->contextname[c], "XWIN_IMAGE")
	|| !strcmp(set->contextname[c], "XPSF_IMAGE")
	|| !strcmp(set->contextname[c], "XPEAK_IMAGE"))
      cx = c;
    else if (!strcmp(set->contextname[c], "Y_IMAGE")
	|| !strcmp(set->contextname[c], "YWIN_IMAGE")
	|| !strcmp(set->contextname[c], "YPSF_IMAGE")
	|| !strcmp(set->contextname[c], "YPEAK_IMAGE"))
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

  if (cx>=0 || cy>=0)
    {
    for (d=0; d<naxis; d++)
      if (!wcs->nprojp)
        projp[1+100*d] += 1.0;	/* Default to the linear term */
    powerst = powers;
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
      }

  free(context);
  free(powers);

/* Now take care of the field-dependent parameters */
  if (set->index2>=0)
    {
    mat2 = mat + set->index2;	/* Field-dependent section in the solution */
    projp[0+100*lng] += mat2[0+lng*(naxis+1)];
    projp[1+100*lng] += mat2[1+lat+lng*(naxis+1)];
    projp[2+100*lng] += mat2[1+lng+lat*(naxis+1)];
    projp[0+100*lat] += mat2[1+lng+lng*(naxis+1)];
    projp[1+100*lat] += mat2[1+lat+lat*(naxis+1)];
    projp[2+100*lat] += mat2[0+lat*(naxis+1)];
    }

  if (ncoeff <= naxis+1)
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


/****** compute jacobian ******************************************************
PROTO	void compute_jacobian(samplestruct *sample, double *dprojdred)
PURPOSE	Compute the Jacobian matrices required for solving the astrometric
	equations.
INPUT	ptr to the sample structure,
	ptr to the dproj/dred Jacobian (or NULL if no computation needed).
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	Memory must have been allocated (naxis*naxis*sizeof(double)) for the
	Jacobian array.
AUTHOR	E. Bertin (IAP)
VERSION	30/12/2003
 ***/
int	compute_jacobian(samplestruct *sample, double *dprojdred)
  {
   wcsstruct	*wcs;
   fgroupstruct	*fgroup;
   fieldstruct	*field;
   setstruct	*set;
   double	rawpos[NAXIS], rawpos2[NAXIS], redpos[NAXIS], wcspos[NAXIS],
		projpos1[NAXIS], projpos2[NAXIS];
   int		i,j, naxis;

  if (!(set=sample->set) || !(field=set->field) || !(fgroup=field->fgroup)
	|| !(wcs=fgroup->wcs))
    return RETURN_ERROR;
  naxis = wcs->naxis;
  for (i=0; i<naxis; i++)
    rawpos[i] = sample->rawpos[i];
  for (i=0; i<naxis; i++)
    {
/*-- Compute terms of the Jacobian dproj/dred matrix */
/*-- Things are a bit tortuous because it is computed with respect */
/*-- to reduced coordinates ("in projected degrees")*/
    raw_to_red(set->wcs, rawpos, redpos);
    redpos[i] -= 0.5*ARCSEC/DEG;
    red_to_raw(set->wcs, redpos, rawpos2);
    raw_to_wcs(set->wcs, rawpos2, wcspos);
    wcs_to_raw(wcs, wcspos, projpos1);
    redpos[i] += 1.0*ARCSEC/DEG;
    red_to_raw(set->wcs, redpos, rawpos2);
    raw_to_wcs(set->wcs, rawpos2, wcspos);
    wcs_to_raw(wcs, wcspos, projpos2);
    for (j=0; j<naxis; j++)
      *(dprojdred++) = 1.0*(projpos2[j] - projpos1[j])/(1*ARCSEC/DEG);
    }

  return RETURN_OK;
  }


/****** astrstats_fgroup *****************************************************
PROTO	void astrstat_fgroup(fgroupstruct *fgroup, fieldstruct *reffield,
			double hsn_thresh)
PURPOSE	Compute stats of a global astrometric solution among a group of fields.
INPUT	ptr to a group of fields pointers,
	ptr to a reference catalog,
	S/N threshold for the high-S/N sample.
OUTPUT	-.
NOTES	Input structures must have gone through crossid_fgroup() first.
AUTHOR	E. Bertin (IAP)
VERSION	31/01/2011
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
   int		i,f,n,s, naxis,nfield,nsamp, nsource,nsource_hsn,
		nmean,nmean_hsn, nmatch,nmatch_hsn;

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
            if ((samp2->sexflags & (OBJ_SATUR|OBJ_TRUNC)))
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
              if ((samp2->sexflags & (OBJ_SATUR|OBJ_TRUNC)))
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
                if ((samp3->sexflags & (OBJ_SATUR|OBJ_TRUNC)))
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
    if ((samp->sexflags & (OBJ_SATUR|OBJ_TRUNC)))
      continue;
    for (i=0; i<naxis; i++)
      mean[i] = mean_hsn[i] = 0.0;
    nmean = nmean_hsn = 0;
    for (samp2 = samp; (samp2=samp2->nextsamp);)
      {
      if ((samp2->sexflags & (OBJ_SATUR|OBJ_TRUNC)))
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
VERSION	27/03/2012
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
   int		i,f,n,s, naxis,nfield,nsamp, flag, nmeani,nmeanr, nclipi,nclipr;

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
            if ((samp2->sexflags & (OBJ_SATUR|OBJ_TRUNC)))
              continue;
            for (i=0; i<naxis; i++)
              meani[i] += samp2->projpos[i];
            nmeani++;
            }
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


/****** reproj_fgroup *******************************************************
PROTO	void reproj_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
PURPOSE	Perform re-projection of sources in a group of fields.
INPUT	ptr to the group of fields.
OUTPUT	-.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	09/09/2004
 ***/
void	reproj_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
  {
   fieldstruct	**field;
   wcsstruct	*wcs;
   setstruct	**pset,
		*set;
   samplestruct	*samp;
   double	*projposmin,*projposmax;
   int		i, f, nsamp,
		s, nfield, naxis;

  field = fgroup->field;
  nfield = fgroup->nfield;
  naxis = fgroup->naxis;

  wcs = fgroup->wcs;
  projposmin = fgroup->projposmin;
  projposmax = fgroup->projposmax;
  for (i=0; i<naxis; i++)
    {
    projposmin[i] = projposmin[i] = BIG;
    projposmax[i] = projposmax[i] = -BIG;
    }

  for (f=0; f<nfield; f++)
    {
    pset = field[f]->set;
    set = *(pset++);
    for (s=field[f]->nset; s--; set=*(pset++))
      {
      samp = set->sample;
      for (nsamp=set->nsample; nsamp--; samp++)
        {
/*------ Compute new projection */
        raw_to_wcs(set->wcs, samp->rawpos, samp->wcspos);
        wcs_to_raw(wcs, samp->wcspos, samp->projpos);
        }
/*---- Compute the boundaries of projected positions for the whole group */
      sort_samples(set);
      for (i=0; i<naxis; i++)
        {
        if (projposmin[i] > set->projposmin[i])
          projposmin[i] = set->projposmin[i];
        if (projposmax[i] < set->projposmax[i])
          projposmax[i] = set->projposmax[i];
	}
      }
    }

  if (reffield)
    {
    set = reffield->set[0];
    samp = set->sample;
    for (nsamp=set->nsample; nsamp--; samp++)
      wcs_to_raw(wcs, samp->wcspos, samp->projpos);
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
        for (n=set->nsample; n--; samp++)
          for (samp2=samp;(samp2=samp2->prevsamp);)
            if (samp2->set->field->astromlabel<0)
              nref++;
        set->weightfac = nref? (double)set->nsample / nref : 1.0;
        }
      }
    }

  return;
  }
