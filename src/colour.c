/*
*				colour.c
*
* Compute color indices and color dependencies.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2008-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "colour.h"
#include "merge.h"
#include "photsolve.h"
#include "prefs.h"
#include "samples.h"

#ifdef USE_THREADS
#include "threads.h"
#endif

/*------------------- global variables for multithreading -------------------*/
#ifdef USE_THREADS
#endif

/****** colour_fgroup ***************************************************
PROTO	void colour_fgroup(fgroupstruct **fgroups, int ngroup)
PURPOSE	Compute color index for sources in a group of fields.
INPUT	ptr to an array of groups,
	number of groups.
OUTPUT	-.
NOTES	Uses the global preferences. Input structures must have gone through
	reproj_fgroup() and crossid_fgroup() first, and preferably through
	astrsolve_fgroups and photsolve_fgroups() too.
AUTHOR	E. Bertin (IAP)
VERSION	12/11/2013
 ***/
void	colour_fgroup(fgroupstruct **fgroups, int ngroup)
  {
   fgroupstruct		*fgroup;
   fieldstruct		*field;
   setstruct		*set;
   msamplestruct	*msamp;
   samplestruct		*samp,*samp2;
   double		*colmat,*wcolmat, *mcol,*wmcol, *col,*wcol, *mag,*wmag,
			sum,wsum, weight, cweight, err2;
   float		colour;
   short		sexflagmask;
   unsigned int		imaflagmask;
   int			c,f,g,m,n,s, b1, b2, c1,c2, band, ncolour, ninstru;

  sexflagmask = (short)prefs.phot_sexflagsmask;
  imaflagmask = prefs.phot_imaflagsmask;
  ninstru = prefs.nphotinstrustr;
  ncolour = (ninstru * (ninstru-1)) / 2;

/* Allocate memory */
  QCALLOC(mag, double, ninstru);
  QCALLOC(wmag, double, ninstru);
  QCALLOC(mcol, double, ncolour);
  QCALLOC(wmcol, double, ncolour);
  QCALLOC(col, double, ncolour);
  QCALLOC(wcol, double, ncolour);

/* Compute average colours */
  for (g=0; g<ngroup; g++)
    {
    fgroup = fgroups[g];
    msamp = fgroup->msample;
    for (m=fgroup->nmsample; m--; msamp++)
      {
      samp = msamp->samp;
      if (!samp->prevsamp)
        continue;
      memset(mag, 0, ninstru*sizeof(double));
      memset(wmag, 0, ninstru*sizeof(double));
      for (samp2 = samp; samp2 && (band=samp2->set->field->photomlabel)>=0;
		samp2 = samp2->prevsamp)
        {
        if ((samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask)
		|| samp2->flux <= 0.0 
		|| (err2 = samp2->magerr*samp2->magerr)<=0.0)
          continue;
        weight = 1.0/(err2 + PHOTOM_MINMAGERR*PHOTOM_MINMAGERR);
        mag[band] += weight*samp2->mag;
        wmag[band] += weight;
        }
      for (band=0; band<ninstru; band++)
        if (wmag[band] > 0.0)
          mag[band] /= wmag[band];
/*---- Update the colour average vector */
      c = 0;
      for (b2=1; b2<ninstru; b2++)
        {
        if (wmag[b2] > 0.0)
          for (b1=0; b1<b2; b1++)
            if (wmag[b1] > 0.0)
              {
              cweight = 1.0;
              mcol[c+b1] += cweight*(mag[b1]-mag[b2]);
              wmcol[c+b1] += cweight;
              }
        c += b2;
        }
      }
    }

  for (c=0; c<ncolour; c++)
    if (wmcol[c]>0.0)
      mcol[c] /= wmcol[c]; 


/* Fill colour covariance matrix */
  QCALLOC(colmat, double, ncolour*ncolour);
  QCALLOC(wcolmat, double, ncolour*ncolour);

  for (g=0; g<ngroup; g++)
    {
    fgroup = fgroups[g];
    msamp = fgroup->msample;
    for (m=fgroup->nmsample; m--; msamp++)
      {
      samp = msamp->samp;
      if (!samp->prevsamp)
        {
        msamp->colour = 0.0;
        continue;
        }
      memset(mag, 0, ninstru*sizeof(double));
      memset(wmag, 0, ninstru*sizeof(double));
      memset(col, 0, ncolour*sizeof(double));
      memset(wcol, 0, ncolour*sizeof(double));
      for (samp2 = samp; samp2 && (band=samp2->set->field->photomlabel)>=0;
		samp2 = samp2->prevsamp)
        {
        if ((samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask)
		|| samp2->flux <= 0.0 
		|| (err2 = samp2->magerr*samp2->magerr)<=0.0)
          continue;
        weight = 1.0/( err2 + PHOTOM_MINMAGERR*PHOTOM_MINMAGERR);
        mag[band] += weight*samp2->mag;
        wmag[band] += weight;
        }
      for (band=0; band<ninstru; band++)
        if (wmag[band] > 0.0)
          mag[band] /= wmag[band];
/*---- Compute centred colour vector */
      c = 0;
      for (b2=1; b2<ninstru; b2++)
        {
        if (wmag[b2] > 0.0)
          for (b1=0; b1<b2; b1++)
            if (wmag[b1] > 0.0)
              {
              cweight = 1.0;
              col[c+b1] = cweight*(mag[b1]-mag[b2] - mcol[c+b1]);
              wcol[c+b1] = cweight;
              }
        c += b2;
        }
/*---- Update the covariance matrix */
      for (c2=0; c2<ncolour; c2++)
        {
        if (wcol[c2] > 0.0)
          for (c1=0; c1<=c2; c1++)
            if (wcol[c1] > 0.0)
              {
              cweight = 1.0;
              c = c2*ncolour+c1;
              colmat[c] += cweight*col[c1]*col[c2];
              wcolmat[c] += cweight;
              }
        }
      }
    }

/* Normalize and complete second half of covariance matrix */
  c = 0;
  for (c2=0; c2<ncolour; c2++)
    {
    for (c1=0; c1<=c2; c1++)
      {
      c = c2*ncolour+c1;
      if (wcolmat[c])
        colmat[c1*ncolour+c2] = (colmat[c] /= wcolmat[c]);
      }
    }

  free(wcolmat);

  colour_findpc(colmat, col, ncolour);
  free(colmat);

/* Compute the "colour index" of each source */
//printf("%g\t%g\t%g\t%g\t%g\t%g \n\n", mcol[0],mcol[1],mcol[2],mcol[3],mcol[4],mcol[5]);
//printf("%g\t%g\t%g\t%g\t%g\t%g \n\n", col[0],col[1],col[2],col[3],col[4],col[5]);
  for (g=0; g<ngroup; g++)
    {
    fgroup = fgroups[g];
    msamp = fgroup->msample;
    for (m=fgroup->nmsample; m--; msamp++)
      {
      samp = msamp->samp;
      if (!samp->prevsamp)
        continue;
      memset(mag, 0, ninstru*sizeof(double));
      memset(wmag, 0, ninstru*sizeof(double));
      for (samp2 = samp; samp2 && (band=samp2->set->field->photomlabel)>=0;
		samp2 = samp2->prevsamp)
        {
        if ((samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask)
		|| samp2->flux <= 0.0 
		|| (err2 = samp2->magerr*samp2->magerr)<=0.0)
          continue;
        weight = 1.0/(err2 + PHOTOM_MINMAGERR*PHOTOM_MINMAGERR);
        mag[band] += weight*samp2->mag;
        wmag[band] += weight;
        }
      for (band=0; band<ninstru; band++)
        if (wmag[band] > 0.0)
          mag[band] /= wmag[band];
/*---- Compute centred colour vector */
      c = 0;
      sum = wsum = 0.0;
      for (b2=1; b2<ninstru; b2++)
        {
        if (wmag[b2] > 0.0)
          for (b1=0; b1<b2; b1++)
            if (wmag[b1] > 0.0)
              {
              sum += col[c+b1]*(mag[b1]-mag[b2] - mcol[c+b1]);
              wsum += col[c+b1]*col[c+b1];
              }
        c += b2;
        }
/*---- Compute incomplete dot product */
      if (wsum>0.0)
        {
        colour = (float)(sum / wsum);
        msamp->colour = colour;
        }
      else
        {
        for (samp2 = samp; samp2 && samp2->set->field->photomlabel>=0;
		samp2 = samp2->prevsamp)
          samp2->scampflags |= SCAMP_PHOTNOCOLOR;
        msamp->colour = 0.0;
        msamp->scampflags |= SCAMP_PHOTNOCOLOR;
        }
      }
    }

/* Free memory */
  free(mag);
  free(wmag);
  free(mcol);
  free(wmcol);
  free(col);
  free(wcol);

  return;
  }


/****** colour_findpc *********************************************************
PROTO	double colour_findpc(double *covmat, float *vec, int nmat)
PURPOSE	Find the principal component (the one with the highest eigenvalue) and
	subtract its contribution from the covariance matrix, using the
	iterative "power" method.
INPUT	Covariance matrix,
	output vector,
	Number of principal components.
OUTPUT  Eigenvalue (variance) of the PC.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 07/02/2011
 ***/
double colour_findpc(double *covmat, double *vec, int nmat)
  {
   double	*tmat,*xmat, *c,*t,*x,
		dtval,dtnorm, xval,xnorm, lambda;
   int		i,j,n;

  QMALLOC(xmat, double, nmat);

/* First initialize the eigenvector to an arbitrary direction */
  tmat = vec;
  for (t=tmat,i=nmat; i--;)
    *(t++) = 1.0;

  dtnorm = 1.0;
  for (n=PCA_NITER; n-- && dtnorm>PCA_CONVEPS;)    
    {
/*-- Compute |x> = C|t> */
    xnorm = 0.0;
    for (c=covmat,x=xmat,j=nmat; j--;)
      {
      for (xval=0.0,t=tmat,i=nmat; i--;)
        xval += *(c++)**(t++);
      xnorm += xval*xval;
      *(x++) = xval;
      }
/*-- Compute |t> = |x>/||x|| and ||Delta t|| (for testing convergence) */
    xnorm = 1.0/sqrt(xnorm);
    dtnorm = 0.0;
    for (t=tmat,x=xmat,i=nmat; i--;)
      {
      dtval = *t;
      dtval -= (*(t++) = *(x++)*xnorm);
      dtnorm += dtval*dtval;
      }
    dtnorm = sqrt(dtnorm);
    }

  free(xmat);

/* Compute the eigenvalue lambda = <t|C|t> */
  lambda = 0.0;
  for (c=covmat,x=tmat,j=nmat; j--;)
    {
    for (xval=0.0,t=tmat,i=nmat; i--;)
      xval += *(c++)**(t++);
    lambda += xval**(x++);
    }

/* Finally subtract the contribution from the found PC: C -= lambda.|t><t| */
  for (c=covmat,x=tmat,j=nmat; j--;)
    for (xval=*(x++)*lambda,t=tmat,i=nmat; i--;)
      *(c++) -= xval**(t++);

  return lambda;
  }


