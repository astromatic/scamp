/*
                                  mosaic.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        SCAMP
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Routines related to the processing of mosaics
*
*       Last modify:    18/01/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
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
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "match.h"
#include "mosaic.h"
#include "misc.h"
#include "prefs.h"
#include "astrefcat.h"
#include "samples.h"
#include ATLAS_LAPACK_H
#ifdef USE_THREADS
#include "threads.h"
#endif

/*------------------- global variables for multithreading -------------------*/
#ifdef USE_THREADS
 pthread_t		*thread;
 pthread_mutex_t	adjmutex;
 threads_gate_t		*pthread_startgate, *pthread_stopgate;
 fieldstruct		**pthread_fields;
 int			*pthread_sviewflag,
			pthread_endflag, pthread_nset, pthread_nfield,
			pthread_lindex, pthread_sindex, pthread_sviewindex;
#endif

/****** adjust_mosaic ********************************************************
PROTO	void adjust_mosaic(fieldstruct **fields, int nfield)
PURPOSE	Adjust the relative positioning of mosaic elements (FITS extensions).
INPUT	ptr to the array of field pointers,
	number of fields.
OUTPUT	-.
NOTES	Uses the global preferences. All input fields must have a
	locate_field() up to date.
AUTHOR	E. Bertin (IAP)
VERSION	27/09/2004
 ***/
void	adjust_mosaic(fieldstruct **fields, int nfield)
  {
   fieldstruct	**matchfields,
		*field;
   setstruct	**sets;
   char		str[88];
   double	*wcsmean;
   int		d,f,l,s, nmatchfield, naxis, nlabel, flag;

  NFPRINTF(OUTPUT, "Making mosaic adjustments...");

  naxis = 0;		/* Avoid gcc -Wall warnings */
/* First step: convert CRVAL shifts to CRPIX variations */
  flag = 1;
  for (f=0; f<nfield; f++)
    {
    field = fields[f];
    if (prefs.mosaic_type[field->astromlabel] == MOSAIC_UNCHANGED)
      continue;
/*-- Use the average field center as a the new projection center */
    sprintf(str, "Computing set shifts for field %d/%d", f+1, nfield);
    NFPRINTF(OUTPUT, str);
    wcsmean = field->meanwcspos;
    naxis = field->naxis;
    sets = field->set;
    if (prefs.mosaic_type[field->astromlabel] == MOSAIC_SAMECRVAL)
      for (s=0; s<field->nset; s++)
        {
        for (d=0; d<naxis; d++)
          sets[s]->wcs->crval[d] = wcsmean[d];
        init_wcs(sets[s]->wcs);
        range_wcs(sets[s]->wcs);
        invert_wcs(sets[s]->wcs);
        }
    else
      for (s=0; s<field->nset; s++)
        crval_to_crpix(sets[s]->wcs, wcsmean);

    if (prefs.mosaic_type[field->astromlabel] == MOSAIC_FIXFOCALPLANE)
      flag = 0;

    }

/* Leave if every field has to be left unchanged */
  if (flag)
    return;

/* Do the job separately for each instrument */
  if ((nlabel=prefs.nastrinstrustr)>1)
    {
    QMALLOC(matchfields, fieldstruct *, nfield);
    nmatchfield = 0;
    } 
  else
    {
    nlabel = 1;
    matchfields = fields;
    nmatchfield = nfield;
    }

  for (l=0; l<nlabel; l++)
    {
    if (prefs.mosaic_type[l] != MOSAIC_FIXFOCALPLANE)
      continue;
    if (nlabel>1)
      {
      nmatchfield = 0;
      for (f=0; f<nfield; f++)
        if (fields[f]->astromlabel == l)
          matchfields[nmatchfield++] = fields[f];
      }

/*-- Compute median reduced coords between the optic center and ext. corners */
#ifdef USE_THREADS
    pthread_adjust_sets(matchfields, nmatchfield, l);
#else
    for (s=0; s<matchfields[0]->nset; s++)
      {
      sprintf(str, "Instrument A%-2d: Adjusting set %d/%d",
	l+1, s+1, matchfields[0]->nset);
      NFPRINTF(OUTPUT, str);
      adjust_set(matchfields, nmatchfield, s);
      }
#endif
    }

  if (nlabel>1)
    free(matchfields);

  return;
  }


/****** adjust_set ***********************************************************
PROTO	void adjust_set(fieldstruct **fields, int nfield, int s)
PURPOSE	Adjust the relative positioning of one mosaic element (FITS extension).
INPUT	ptr to the array of field pointers,
	number of fields,
	extension number.
OUTPUT	-.
NOTES	Uses the global preferences. All input fields must have a
	locate_field() up to date.
AUTHOR	E. Bertin (IAP)
VERSION	18/01/2006
 ***/
void	adjust_set(fieldstruct **fields, int nfield, int s)
  {
   double	x[(NAXIS+1)*NAXIS], cd[NAXIS*NAXIS],
		rawpos[NAXIS], redpos[NAXIS],
		*medstack;
   int		ipiv[NAXIS], *naxisn,
		d,e,f,i, naxis, noksets;


  naxisn = fields[0]->set[s]->wcs->naxisn;
  naxis = fields[0]->set[s]->wcs->naxis;
  QMALLOC(medstack, double, nfield*naxis);
  for (d=0; d<naxis; d++)
    rawpos[d] = 0.5;
  for (i=-1; i<naxis; i++)
    {
    if (i>=0)
      {
      rawpos[i] = naxisn[i]+0.5;	/* one corner*/
      if (i)
        rawpos[i-1] = 0.5;
      }
    noksets = 0;
    for (f=0; f<nfield; f++)
      {
      if (fields[f]->set[s]->nsample < prefs.fixfocalplane_nmin)
        continue;
      raw_to_red(fields[f]->set[s]->wcs, rawpos, redpos);
      for (d=0; d<naxis; d++)
        medstack[d*nfield+noksets] = redpos[d];
      noksets++;
      }
    if (noksets)
      for (d=0; d<naxis; d++)
/*------ Median of the corner coordinates */
        x[(i+1)*naxis+d] = fast_median(medstack+d*nfield, noksets);
    else
      {
/*---- No valid CCD, take the first one as a reference */
      raw_to_red(fields[0]->set[s]->wcs, rawpos, redpos);
      for (d=0; d<naxis; d++)
        x[(i+1)*naxis+d] = redpos[d];
      }
    }
/*-- Derive the CD's */
  for (d=0; d<naxis; d++)
    for (e=0; e<naxis; e++)
      cd[d*naxis+e] = (x[(e+1)*naxis+d] - x[d])/naxisn[e];
  for (f=0; f<nfield; f++)
    memcpy(fields[f]->set[s]->wcs->cd, cd,naxis*naxis*sizeof(double));
/* Derive the CRPIXs */
  clapack_dgesv(CblasRowMajor, naxis, 1, cd, naxis, ipiv, x, naxis);
  for (f=0; f<nfield; f++)
    for (d=0; d<naxis; d++)
      fields[f]->set[s]->wcs->crpix[d] = 0.5 - x[d];

/* Initialize other WCS structures */
  for (f=0; f<nfield; f++)
    {
    init_wcs(fields[f]->set[s]->wcs);
/*-- Find the range of coordinates */
    range_wcs(fields[f]->set[s]->wcs);
/*-- Invert projection corrections */
    invert_wcs(fields[f]->set[s]->wcs);
    }

  free(medstack);

  return;
  }


/****** crval_to_crpix *******************************************************
PROTO	void crval_to_crpix(wcsstruct *wcs, double *wcspos)
PURPOSE	Keep a WCS projection approx. the same but with a new CRVAL.
INPUT	ptr to the WCS structure,
	new CRVAL vector.
OUTPUT	-.
NOTES	The updated WCS is an approximation of the exact one.
AUTHOR	E. Bertin (IAP)
VERSION	19/02/2005
 ***/
void	crval_to_crpix(wcsstruct *wcs, double *wcspos)
  {
   double	rawpos[NAXIS], a[NAXIS*NAXIS],b[NAXIS*NAXIS],
		*c,*at,
		val, cas, sas, angle, dlng,dlat;
   int		i,j,k, lng,lat, naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;

  wcs_to_raw(wcs, wcspos, rawpos);

  if (lng != lat)
    {
/*-- Compute the angle difference towards the north pole induced by the shift*/
    dlng = wcspos[lng] - wcs->crval[lng];
    dlat = wcspos[lat] - wcs->crval[lat];
    angle = (dlng!=0.0 && dlat!= 0.0) ?
	(atan2(sin(dlng*DEG),
	cos(wcspos[lat]*DEG)*tan(wcs->crval[lat]*DEG)
	- sin(wcspos[lat]*DEG)*cos(dlng*DEG))
	+ atan2(sin(dlng*DEG),
	cos(wcs->crval[lat]*DEG)*tan(wcspos[lat]*DEG)
	- sin(wcs->crval[lat]*DEG)*cos(dlng*DEG)))/DEG - 180.0
	: 0.0;
/*-- A = B*C */
    c = wcs->cd;
/*-- The B matrix is made of 2 numbers */
    cas = cos(angle*DEG);
    sas = sin(angle*DEG);
    for (i=0; i<naxis; i++)
      b[i+i*naxis] = 1.0;
    b[lng+lng*naxis] = cas;
    b[lat+lng*naxis] = -sas;
    b[lng+lat*naxis] = sas;
    b[lat+lat*naxis] = cas;
    at = a;
    for (j=0; j<naxis; j++)
      for (i=0; i<naxis; i++)
        {
        val = 0.0;
        for (k=0; k<naxis; k++)
          val += b[k+j*naxis]*c[i+k*naxis];
        *(at++) = val;
        }

    at = a;

    for (i=0; i<naxis*naxis; i++)
      *(c++) = *(at++);
    }

  for (i=0; i<naxis; i++)
    {
    wcs->crval[i] = wcspos[i];
    wcs->crpix[i] = rawpos[i];
    }

/* Initialize other WCS structures */
  init_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);
  if (lng == lat)
    return;

  return;
  }


#ifdef USE_THREADS

/****** pthread_adjust_set ***************************************************
PROTO   void *pthread_adjust_set(void *arg)
PURPOSE thread that takes care of adjusting mosaic elements.
INPUT   Pointer to the thread number.
OUTPUT  -.
NOTES   Relies on global variables.
AUTHOR  E. Bertin (IAP)
VERSION 27/09/2004
 ***/
void    *pthread_adjust_set(void *arg)
  {
   char	str[80];
   int	sindex, proc;

  sindex = -1;
  proc = *((int *)arg);
  threads_gate_sync(pthread_startgate);
  while (!pthread_endflag)
    {
    QPTHREAD_MUTEX_LOCK(&adjmutex);
    if (sindex>-1)
/*---- Indicate that the field info is now suitable for viewing */
      pthread_sviewflag[sindex] = 1;
    while (pthread_sviewindex<pthread_nset
	&& pthread_sviewflag[pthread_sviewindex])
      {
      sprintf(str, "Instrument A%-2d: Adjusting set %d/%d",
	pthread_lindex+1, ++pthread_sviewindex, pthread_nset);
      NFPRINTF(OUTPUT, str);
      }
    if (pthread_sindex<pthread_nset)
      {
      sindex = pthread_sindex++;
      QPTHREAD_MUTEX_UNLOCK(&adjmutex);
/*---- Adjust set */
      adjust_set(pthread_fields, pthread_nfield, sindex);
      }
    else
      {
      QPTHREAD_MUTEX_UNLOCK(&adjmutex);
/*---- Wait for the input buffer to be updated */
      threads_gate_sync(pthread_stopgate);
/* ( Master thread process loads and saves new data here ) */
      threads_gate_sync(pthread_startgate);
      }
    }

  pthread_exit(NULL);

  return (void *)NULL;
  }


/****** pthread_adjust_sets ***************************************************
PROTO   void pthread_adjust_sets(fieldstruct **fields, int nfield, int l)
PURPOSE Adjust the relative positioning of mosaic elements in parallel using
	threads.
INPUT   Pointer to field structure pointers,
	number of set,
	number of fields,
	current astrometric instrument index (label).
OUTPUT  -.
NOTES   Relies on global variables.
AUTHOR  E. Bertin (IAP)
VERSION 27/09/2004
 ***/
void    pthread_adjust_sets(fieldstruct **fields, int nfield, int l)
  {
   static pthread_attr_t	pthread_attr;
   int				*proc,
				p;

/* Number of active threads */
  nproc = prefs.nthreads;
  pthread_fields = fields;
  pthread_nfield = nfield;
  pthread_nset = fields[0]->nset;
  QCALLOC(pthread_sviewflag, int, pthread_nset);
/* Set up multi-threading stuff */
  QMALLOC(proc, int, nproc);
  QMALLOC(thread, pthread_t, nproc);
  QPTHREAD_MUTEX_INIT(&adjmutex, NULL);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_startgate = threads_gate_init(nproc+1, NULL);
  pthread_stopgate = threads_gate_init(nproc+1, NULL);
/* Start the reading threads */
  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_adjust_set, &proc[p]);
    }
  QPTHREAD_MUTEX_LOCK(&adjmutex);
  pthread_sindex = pthread_sviewindex = 0;
  pthread_endflag = 0;
  pthread_lindex = l;
  QPTHREAD_MUTEX_UNLOCK(&adjmutex);
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
  QPTHREAD_MUTEX_DESTROY(&adjmutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  free(pthread_sviewflag);
  free(proc);
  free(thread);
  }

#endif

