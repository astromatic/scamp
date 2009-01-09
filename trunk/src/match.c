/*
                                  match.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        SCAMP
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Pattern matching routines
*
*       Last modify:    20/03/2007
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
#include "check.h"
#include "fft.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "match.h"
#include "mosaic.h"
#include "prefs.h"
#include "astrefcat.h"
#include "samples.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include "wcs/wcs.h"
#include ATLAS_LAPACK_H

/*------------------- global variables for multithreading -------------------*/
#ifdef USE_THREADS
 pthread_t		*thread;
 pthread_mutex_t	matchmutex,matchsortmutex;
 threads_gate_t		*pthread_startgate, *pthread_stopgate;
 fgroupstruct		**pthread_fgroups;
 fieldstruct		**pthread_reffields;
 int			**pthread_fviewflag,
			pthread_endflag, pthread_ngroup,
			pthread_gindex, pthread_findex,
			pthread_gviewindex,pthread_fviewindex;
#endif

 int	sort_coord;
 static int	compmag(const void *sample1, const void *sample2);
 static int	compraw(const void *sample1, const void *sample2);
 static void	putgauss(float *histo, int width,int height, double x,double y,
			double flux, double xsig, double ysig);

/****** match_field ***********************************************************
PROTO	void match_field(fieldstruct *field, fieldstruct *reffield)
PURPOSE	Perform (linear) pattern matching between a linked-set of catalogs and
	an astrometric reference catalog. The astrometric parameters are
	updated accordingly.
INPUT	ptr to the field to be matched,
	ptr to the reference field.
OUTPUT	-.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	28/09/2006
 ***/
void	match_field(fieldstruct *field, fieldstruct *reffield)
  {
   setstruct	*refset, *refset2, *set, *fieldset;
   samplestruct	*refsample;
   char		str[128];
   double	wcspos[NAXIS],
		angle,scale,dlng,dlat,sig, dlng2,dlat2,sig2, asig,
		dangle,dscale, ddlng,ddlat, shear,sangle, sheartot, sangletot,
		matchresol, refarea,area, refcrossec,crossec;
   int		i,j,k, nrefsample, lng,lat, reflng,reflat,
		naxis, nmax;

  sprintf(str, "Matching field %s...", field->rfilename);
  NFPRINTF(OUTPUT, str);

/* We are going to use FFT routines */
/* fft_init(); Moved to makeit() */

  naxis = field->naxis;
  lng = field->lng;
  lat = field->lat;

/* Prepare the WCS buffer by setting all elements to average values */
  for (k=0; k<naxis; k++)
    wcspos[k] = field->meanwcspos[k];
/* Convert the reference set to the current field coordinate labeling */
  reflng = reffield->lng;
  reflat = reffield->lat;
  refset = reffield->set[0];
  refsample = refset->sample;
  nrefsample = refset->nsample;
  for (j=0; j<nrefsample; j++)
    {
    wcspos[lng] = refsample[j].wcspos[reflng];
    wcspos[lat] = refsample[j].wcspos[reflat];
    for (k=0; k<naxis; k++)
      refsample[j].wcspos[k] = wcspos[k];
    }
  refset->naxis = reffield->naxis = naxis;
  refset->lng = reffield->lng = lng;
  refset->lat = reffield->lat = lat;

  if (field->mosaic_type == MOSAIC_LOOSE)
    {
/*-- Process each set independently */
    field->match_dscale = field->match_dangle = field->match_asig
			=  field->match_dlng = field->match_dlat
			= field->match_sig = 0.0;
    for (i=0; i<field->nset; i++)
      {
      set = field->set[i];
/*---- Only keep reference samples close enough to the current set */
      refset2 = frame_set(refset, set->wcs, set->wcspos,
		set->radius+prefs.radius_maxerr);
/*---- Compute the resolution used for MATCHing */
      if (prefs.match_resol)
        matchresol = prefs.match_resol;
      else
        {
/*------ Compute mean area per source and cross-section from (Gaussian) errors*/
        refarea = MATCH_CONFPROB * (PI*refset2->radius*refset2->radius)
		/ (refset2->nsample+1.0);
        area = MATCH_CONFPROB * set->wcs->naxisn[lng]*set->wcs->naxisn[lat]
		*set->wcsscale[lng]*set->wcsscale[lat] / (set->nsample+1.0);
        compute_rawpos(set->wcs, refset2->sample, refset2->nsample);
        refcrossec = -2.0*PI*mean_rawposvar(refset2)
		*set->wcsscale[lng]*set->wcsscale[lat]*log(MATCH_CONFPROB);
        crossec = -2.0*PI*mean_rawposvar(set)
		*set->wcsscale[lng]*set->wcsscale[lat]*log(MATCH_CONFPROB);
/*------ Choose the smallest of them all */
        if (refarea < area)
          area = refarea;
        if (refcrossec > crossec)
          crossec = refcrossec;
        matchresol = sqrt(crossec > area? crossec: area);
        }
/*-- Choose the number of stars used for MATCHing */
      if (prefs.nmatchmax>0)
        asig = match_setas(set, refset2, prefs.nmatchmax, matchresol,
			&angle, &scale);
      else
        for (nmax = AS_NSOURCESTART;
	nmax <= AS_NSOURCEMAX &&
	(asig = match_setas(set, refset2, nmax, matchresol, &angle,&scale))
		<MATCH_MINCONT && 
	(nmax < set->nsample || nmax < refset2->nsample);
	nmax*=2);
      locate_set(set);
      update_wcsas(set->wcs, -set->wcs->chirality*angle, scale);
      sig = match_setll(set, refset2, matchresol, &dlng, &dlat);
      if (prefs.posangle_maxerr > 90.0)
        {
/*------ The result might be rotated by more than 90 deg */
        update_wcsas(set->wcs, 180.0, 1.0);
        sig2 = match_setll(set, refset2, matchresol, &dlng2, &dlat2);
        if (sig2>sig)
          {
          sig = sig2;
          angle += 180.0;
          dlng = dlng2;
          dlat = dlat2;
          }
        else
          update_wcsas(set->wcs, 180.0, 1.0);
        }
      update_wcscc(set->wcs, dlng,dlat);
      end_set(refset2); 
      refset2 = frame_set(refset, set->wcs, set->wcspos,
		set->radius+prefs.radius_maxerr);
/*---- Now refine all parameters */
      match_refine(set, refset2, matchresol,
	&dangle, &dscale, &sangle, &shear, &ddlng, &ddlat);
/*---- Apply refinement */
      update_wcsas(set->wcs, set->wcs->chirality*dangle, dscale);
      update_wcsss(set->wcs, sangle, shear);
      update_wcscc(set->wcs, -ddlng, -ddlat);
      end_set(refset2);
      compute_wcsss(set->wcs, &sangletot, &sheartot);
      locate_set(set);

      field->match_dscale += (set->match_dscale = scale*dscale);
      field->match_dangle += (set->match_dangle = fmod(angle+dangle, 180.0));
      field->match_shear += (set->match_shear = sheartot);
      field->match_sangle += (set->match_sangle = sangletot);
      field->match_asig += (set->match_asig = asig);
      field->match_dlng += (set->match_dlng = dlng+ddlng);
      field->match_dlat += (set->match_dlat = dlat+ddlat);
      field->match_sig += (set->match_sig = sig);
      }
/*-- Update higher level localisation data */
    locate_field(field);

    field->match_dscale /= (float)field->nset;
    field->match_dangle /= (float)field->nset;
    field->match_shear /= (float)field->nset;
    field->match_sangle /= (float)field->nset;
    field->match_asig /= (float)field->nset;
    field->match_dlng /= (float)field->nset;
    field->match_dlat /= (float)field->nset;
    field->match_sig /= (float)field->nset;
    }
  else
    {
/*-- Reproject the whole thing to a single "virtual" image */
    fieldset = new_fieldset(field);

/*-- Keep only reference samples close enough to the current set */
    refset2 = frame_set(refset, fieldset->wcs, fieldset->wcspos,
		fieldset->radius+prefs.radius_maxerr);

/*-- Compute the resolution used for MATCHing */
    if (prefs.match_resol)
      matchresol = prefs.match_resol;
    else
      {
/*---- Compute mean area per source and cross-section from (Gaussian) errors*/
      refarea = MATCH_CONFPROB * (PI*refset2->radius*refset2->radius)
		/ (refset2->nsample+1.0);
      area = MATCH_CONFPROB * PI*fieldset->radius*fieldset->radius
		/ (fieldset->nsample+1.0);
      compute_rawpos(fieldset->wcs, refset2->sample, refset2->nsample);
      refcrossec = -2.0*PI*mean_rawposvar(refset2)
		*fieldset->wcsscale[lng]*fieldset->wcsscale[lat]
		*log(MATCH_CONFPROB);
      crossec = -2.0*PI*mean_rawposvar(fieldset)
		*fieldset->wcsscale[lng]*fieldset->wcsscale[lat]
		*log(MATCH_CONFPROB);
/*---- Choose the smallest of them all */
      if (refarea < area)
        area = refarea;
      if (refcrossec > crossec)
        crossec = refcrossec;
      matchresol = sqrt(crossec > area? crossec: area);
      }
/*-- Choose the number of stars used for MATCHing */
    if (prefs.nmatchmax>0)
      asig = match_setas(fieldset, refset2, prefs.nmatchmax, matchresol,
		&angle, &scale);
    else
      for (nmax = AS_NSOURCESTART;
	nmax <= AS_NSOURCEMAX &&
	(asig=match_setas(fieldset, refset2, nmax, matchresol, &angle,&scale))
		<MATCH_MINCONT
	&& (nmax < fieldset->nsample || nmax < refset2->nsample);
	nmax*=2);
    end_set(refset2); 
    update_wcsas(fieldset->wcs, -angle, scale);
/*-- Only keep reference samples close enough to the current set */
    refset2 = frame_set(refset, fieldset->wcs, fieldset->wcspos,
		fieldset->radius+prefs.radius_maxerr);
    sig = match_setll(fieldset, refset2, matchresol, &dlng, &dlat);
    if (prefs.posangle_maxerr > 90.0)
      {
      update_wcsas(fieldset->wcs, 180.0, 1.0);
/*---- Only keep reference samples close enough to the current set */
      end_set(refset2); 
      refset2 = frame_set(refset, fieldset->wcs, fieldset->wcspos,
		fieldset->radius+prefs.radius_maxerr);
      sig2 = match_setll(fieldset, refset2, matchresol, &dlng2, &dlat2);
      if (sig2>sig)
        {
        sig = sig2;
        angle = (angle>0.0)? angle - 180.0 : 180.0 + angle;
        dlng = dlng2;
        dlat = dlat2;
        }
      else
        update_wcsas(fieldset->wcs, 180.0, 1.0);
      }

    end_set(refset2);
    end_set(fieldset);

/*-- Apply basic corrections */
    for (i=0; i<field->nset; i++)
      {
      set = field->set[i];
      update_wcsas(set->wcs, -angle, scale);
      update_wcsll(set->wcs, dlng, dlat);
      locate_set(set);
      }
/* Update higher level localisation data */
    locate_field(field);

/*-- Now refine all parameters: */
/*-- Reproject once more the whole thing to a single "virtual" image */
    fieldset = new_fieldset(field);
    refset2 = frame_set(refset, fieldset->wcs, fieldset->wcspos,
		fieldset->radius+prefs.radius_maxerr);
    match_refine(fieldset, refset2, 2.0*matchresol,
	&dangle, &dscale, &sangle, &shear, &ddlng, &ddlat);
    update_wcsas(fieldset->wcs, -dangle, dscale);
    update_wcsss(fieldset->wcs, -sangle, shear);
    update_wcsll(fieldset->wcs, -ddlng, -ddlat);
    compute_wcsss(fieldset->wcs, &sangletot, &sheartot);
/*-- Apply refinements */
    for (i=0; i<field->nset; i++)
      {
      set = field->set[i];
      update_wcsas(set->wcs, dangle, dscale);
      update_wcsss(set->wcs, set->wcs->chirality*sangle, shear);
      update_wcsll(set->wcs, -ddlng, -ddlat);
      locate_set(set);
      }
    end_set(refset2); 
    end_set(fieldset);
/*-- Update higher level localisation data */
    locate_field(field);

    field->match_dscale = scale*dscale;
    field->match_dangle = fmod(angle+dangle, 180.0);
    field->match_shear = sheartot;
    field->match_sangle = sangletot;
    field->match_asig = asig;
    field->match_dlng = dlng+ddlng;
    field->match_dlat = dlat+ddlat;
    field->match_sig = sig;
    }

  return;
  }


/****** new_fieldset *********************************************************
PROTO	setstruct *new_fieldset(fieldstruct *field)
PURPOSE	Create a "virtual" set with default WCS parameters centered on the
	"field" content and containing all the sources belonging to "field".
INPUT	ptr to the original field.
OUTPUT	ptr to the virtual set.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/05/2006
 ***/
setstruct	*new_fieldset(fieldstruct *field)
  {
   setstruct	*set, *fieldset;
   samplestruct	*sample;
   wcsstruct	*wcs;
   char		*ctype[NAXIS];
   double	minpos[NAXIS],maxpos[NAXIS];
   int		i,d,s, naxis;

  naxis = field->naxis;
  fieldset = init_set();
/*
  wcs = fieldset->wcs = copy_wcs(field->set[0]->wcs);
*/
  wcs = field->set[0]->wcs;
  for (d=0; d<naxis; d++)
    ctype[d] = wcs->ctype[d];
  fieldset->wcs = create_wcs(ctype, wcs->crval, wcs->crpix, field->meanwcsscale,
	wcs->naxisn, wcs->naxis);
  wcs = fieldset->wcs;
  for (d=0; d<naxis; d++)
    maxpos[d] = -(minpos[d] = BIG);
  for (i=0; i<field->nset; i++)
    {
    set = field->set[i];
    copy_samples(set->sample, fieldset, set->nsample);
    fieldset->nsaturated += set->nsaturated;
    sample = fieldset->sample + (fieldset->nsample - set->nsample);
    for (s=set->nsample; s--; sample++)
      {
      raw_to_wcs(set->wcs, sample->rawpos, sample->wcspos);
      wcs_to_raw(wcs, sample->wcspos, sample->rawpos);
      for (d=0; d<naxis; d++)
	{
        if (sample->rawpos[d]<minpos[d])
          minpos[d] = sample->rawpos[d];
        if (sample->rawpos[d]>maxpos[d])
          maxpos[d] = sample->rawpos[d];
        }
      }
    }
  fieldset->field = field;
  for (d=0; d<naxis; d++)
    {
    wcs->crpix[d] -= minpos[d];
    wcs->naxisn[d] = (int)(maxpos[d]-minpos[d]+0.5);
    }
  init_wcs(wcs);
  range_wcs(wcs);
  invert_wcs(wcs);
  sample = fieldset->sample;
  for (s=fieldset->nsample; s--; sample++)
    for (d=0; d<naxis; d++)
      sample->rawpos[d] -= minpos[d];

  locate_set(fieldset);

  return fieldset;
  }

/****** match_setas **********************************************************
PROTO	double match_setas(setstruct *set, setstruct *refset, int nmax,
			double matchresol, double *angle, double *scale)
PURPOSE	Find position angle and scale of a set of catalogs relative to an
	astrometric reference catalog using pattern matching.
INPUT	ptr to the set to be matched,
	ptr to the reference set,
	max. number of non-saturated sources to be considered,
	resolution for MATCHing in deg,
	ptr to the position angle (output),
	ptr to the relative scale (output).
OUTPUT	Confidence level of the solution (in units of sigma).
NOTES	The output angle is expressed CCW in projected coordinates (not
	pixel coordinates) for homogeneity reasons.
	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	19/02/2007
 ***/
double	match_setas(setstruct *set, setstruct *refset, int nmax,
			double matchresol, double *angle, double *scale)
  {
   samplestruct	*refsample, *sample;
   wcsstruct	*wcs;
   char		str[MAXCHAR],
		*rfilename;
   double	lopass[NAXIS],hipass[NAXIS],
		thetafac,thetazero, rhofac,rhozero, dx,dy,dr2, theta,rho,
		theta2,rho2, rhosig,thetasig, sig,sig2, flux, dr2max,dr2mid,
		dr2mid2, flip;
   float	*histo, *histo2, *rhisto;
   int		csize[NAXIS],
		i,j,k,n, nrefsample,nsample,npair, lng,lat,
		naxis, cnaxis, cnpix, checknum, ival;

  histo2 = NULL;	/* to avoid gcc -Wall warnings */
  naxis = set->naxis;
  lng = set->lng;
  lat = set->lat;
  rfilename = set->field->rfilename;

/* Build a histogram of star pairs */
  cnaxis = 2;
/* Set the histogram dimensions (must be powers of 2) */
/* The relative position angle... */
  ival = (int)(1.4*PI*set->radius/matchresol);
  for (i=1; (ival>>=1); i<<=1);
  csize[0] = i;
  if (csize[0]>ANGLE_MAXSIZE)
    csize[0] = ANGLE_MAXSIZE;
/* ...and the relative scale (the 0.5 is a fudge factor for "typical" */
  ival = (int)(1.4*set->radius*log(SCALE_RANGE)/matchresol);
  for (i=1; (ival>>=1); i<<=1);  
  csize[1] = i;
  if (csize[1]>SCALE_MAXSIZE)
    csize[1] = SCALE_MAXSIZE;
  lopass[0] = lopass[1] = 1.0/ASLOPASS_LAMBDA;
  hipass[0] = hipass[1] = 1.0/(csize[1]*ASHIPASS_LAMBDA);
  cnpix = 1;
  for (i=0; i<cnaxis; i++)
    cnpix *= csize[i];
  QCALLOC(rhisto, float, cnpix);
  QCALLOC(histo, float, cnpix);
  thetafac = csize[0]/PI;
  thetazero = -PI/2.0;
  rhofac = csize[1]/log(SCALE_RANGE);
  rhozero = log(SCALE_ZERO);

  refsample = refset->sample;
/* Keep only the brightest sources in ref. set */
  if (refset->nsample>nmax)
    {
    qsort(refset->sample, refset->nsample, sizeof(samplestruct), compmag);
    nrefsample = nmax;
    }
  else
    nrefsample = refset->nsample;
  wcs = set->wcs;

/* Compute projected positions of the reference stars for the current set */
#ifndef USE_THREADS
  sprintf(str,"%.24s: projecting reference catalog...", rfilename);
  NFPRINTF(OUTPUT, str);
#endif
  compute_rawpos(wcs, refsample, nrefsample);

/* Compute maximum pair length (squared) */
  dr2max = wcs->naxisn[lng]<wcs->naxisn[lat]?
		wcs->naxisn[lng] : wcs->naxisn[lat]; 
  dr2max *= dr2max;
  dr2mid =dr2max/1.3;
  dr2mid2 = dr2max/(SCALE_RANGE/2.0);
/* First build histogram for the reference field */
  npair = nrefsample*(nrefsample-1)/2;
  n = 0;
  for (j=0; j<nrefsample; j++)
    for (k=j+1; k<nrefsample; k++)
      {
#ifdef USE_THREADS
      if (!(++n%100000) && prefs.nthreads<2)
#else
      if (!(++n%100000))
#endif
        {
        sprintf(str,"%.24s: reference pair #%d / %d processed",
		rfilename, n,npair);
        NFPRINTF(OUTPUT, str);
        }
      dx = refsample[j].rawpos[lng] - refsample[k].rawpos[lng];
      dy = refsample[j].rawpos[lat] - refsample[k].rawpos[lat];
      dr2 = dx*dx+dy*dy;
      if (dr2<dr2max && dr2>0.0)
	{
        theta = (dx!=0.0? (atan(dy/dx)-thetazero) : 0.0)*thetafac;
        rho=(0.5*log(dr2)-rhozero)*rhofac;
        sig = sqrt((refsample[j].rawposerr[lng]*refsample[j].rawposerr[lng]
		+refsample[k].rawposerr[lng]*refsample[k].rawposerr[lng]
		+refsample[j].rawposerr[lat]*refsample[j].rawposerr[lat]
		+refsample[k].rawposerr[lat]*refsample[k].rawposerr[lat])/dr2);
        thetasig = sig*thetafac;
        rhosig = sig*rhofac;
/*------ Take magnitudes into account */
        flux = -0.4*AS_FLUXEXP*(refsample[j].mag+refsample[k].mag);
        if (thetasig>GAUSS_MAXSIG)
          {
          flux *= (GAUSS_MAXSIG*GAUSS_MAXSIG) / (thetasig*thetasig);
          thetasig = GAUSS_MAXSIG;
          }
        if (rhosig>GAUSS_MAXSIG)
          {
          flux *= (GAUSS_MAXSIG*GAUSS_MAXSIG) / (rhosig*rhosig);
          rhosig = GAUSS_MAXSIG;
          }
        flux = (flux<70.0? (flux>-70.0? DEXP(flux) : 0.0) : BIG);
        if (dr2>dr2mid)
          flux *= (dr2max-dr2)/(dr2max-dr2mid);
        if (dr2<dr2mid2)
          flux *= dr2/dr2mid2;
        putgauss(rhisto, csize[0], csize[1], theta,rho, flux,
		thetasig,rhosig);
        }
      }

/* Check-image for angle-scale pair histogram of the reference catalog */
  if ((checknum=check_check(CHECK_ASREFPAIR))>=0)
    write_check(prefs.check_name[checknum], rhisto, csize[0], csize[1]);

/* Now the current field */
  sample = set->sample;
/* Keep only the brightest sources if too many samples in set */
  if (set->nsample>nmax)
    {
    qsort(set->sample, set->nsample, sizeof(samplestruct), compmag);
    nsample = nmax;
    }
  else
    nsample = set->nsample;

  npair = nsample*(nsample-1)/2;
  n = 0;
  for (j=0; j<nsample; j++)
    for (k=j+1; k<nsample; k++)
      {
#ifdef USE_THREADS
      if (!(++n%100000) && prefs.nthreads<2)
#else
      if (!(++n%100000))
#endif
        {
        sprintf(str,"%.24s: detection pair #%d / %d processed",
		rfilename, n,npair);
        NFPRINTF(OUTPUT, str);
        }
      dx = sample[j].rawpos[lng] - sample[k].rawpos[lng];
      dy = sample[j].rawpos[lat] - sample[k].rawpos[lat];
      dr2 = dx*dx+dy*dy;
      if (dr2<dr2max && dr2>0.0)
        {
        theta = (dx!=0.0? (atan(dy/dx)-thetazero) : 0.0)*thetafac;
        rho=(0.5*log(dr2)-rhozero)*rhofac;
        sig = sqrt((sample[j].rawposerr[lng]*sample[j].rawposerr[lng]
		+sample[k].rawposerr[lng]*sample[k].rawposerr[lng]
		+sample[j].rawposerr[lat]*sample[j].rawposerr[lat]
		+sample[k].rawposerr[lat]*sample[k].rawposerr[lat])/dr2);
        thetasig = sig*thetafac;
        rhosig = (sig+0.0002)*rhofac;
        flux = sample[j].flux*sample[k].flux;
        if (thetasig>GAUSS_MAXSIG)
          {
          flux *= (GAUSS_MAXSIG*GAUSS_MAXSIG) / (thetasig*thetasig);
          thetasig = GAUSS_MAXSIG;
          }
        if (rhosig>GAUSS_MAXSIG)
          {
          flux *= (GAUSS_MAXSIG*GAUSS_MAXSIG) / (rhosig*rhosig);
          rhosig = GAUSS_MAXSIG;
          }
        flux = (flux>=0.0? pow(flux,AS_FLUXEXP) : 0.0);
        if (dr2>dr2mid)
          flux *= (dr2max-dr2)/(dr2max-dr2mid);
        if (dr2<dr2mid2)
          flux *= dr2/dr2mid2;
        putgauss(histo, csize[0], csize[1], theta,rho, flux, thetasig,rhosig);
        }
      }

/* Check-image for angle-scale pair histogram of the reference catalog */
  if ((checknum=check_check(CHECK_ASPAIR))>=0)
    write_check(prefs.check_name[checknum], histo, csize[0], csize[1]);

/* Make the cross-correlation */
#ifndef USE_THREADS
  sprintf(str,"%.24s: finding field scale and position angle...", rfilename);
  NFPRINTF(OUTPUT, str);
#endif
/* Copy first histogram in case things go wrong */
  QMEMCPY(histo, histo2, float, cnpix);
  fastcorr(histo, rhisto, cnaxis, csize, lopass, hipass);
  sig = findcrosspeak(histo, csize[0], csize[1],
		prefs.posangle_maxerr*DEG*thetafac,
		log(prefs.pixscale_maxerr)*rhofac, 
		1.0/lopass[0], 1.0/lopass[1],
		&theta, &rho);
  flip = 1.0;
  if (prefs.flip_flag)
/* Flip the histogram and try again */
    {
    flipimage(histo2, csize[0], csize[1]);
    fastcorr(histo2, rhisto, cnaxis, csize, lopass, hipass);
    sig2 = findcrosspeak(histo2, csize[0], csize[1],
		prefs.posangle_maxerr*DEG*thetafac,
		log(prefs.pixscale_maxerr)*rhofac,
		1.0/lopass[0], 1.0/lopass[1],
		&theta2, &rho2);
    if (sig2>sig)
      {
      flip = -1.0;
      sig = sig2;
      free(histo);
      histo = histo2;
      rho = rho2;
      theta = theta2;
      }
    else
      free(histo2);
    }
  else
    free(histo2);

/* Check-image for angle-scale pair histogram of the reference catalog */
  if ((checknum=check_check(CHECK_ASXCORR))>=0)
    write_check(prefs.check_name[checknum], histo, csize[0], csize[1]);

  if (sig>1.0)
    {
    *angle = theta/thetafac/DEG;
    *scale = flip*exp(-rho/rhofac);
    }
  else
    {
    *angle = 0.0;
    *scale = 1.0;
    }

  free(rhisto);
  free(histo);

  return sig;
  }


/****** match_setll **********************************************************
PROTO	double match_setll(setstruct *set, setstruct *refset,
			double matchresol, double *dlng, double *dlat)
PURPOSE	Find the celestial coordinate shift of a set of catalogs relative to an
	astrometric reference catalog using pattern matching.
INPUT	ptr to the set to be matched,
	ptr to the reference set,
	resolution for matching in deg,
	ptr to the longitude shift (output),
	ptr to the latitude shift (output).
OUTPUT	Confidence level of the solution (in units of sigma).
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	19/02/2007
 ***/
double	match_setll(setstruct *set, setstruct *refset,
			double matchresol, double *dlng, double *dlat)
  {
   samplestruct	*refsample, *sample;
   wcsstruct	*wcs;
   char		str[MAXCHAR],
		*rfilename;
   double	lopass[NAXIS], hipass[NAXIS], rawpos[NAXIS], wcspos[NAXIS],
		xfac,xzero, yfac,yzero, x,y, xsig,ysig, sig, flux, poserr;
   float	*histo, *rhisto;
   int		csize[NAXIS],
		i,k,n, nrefsample, nsample, lng,lat,
		naxis, cnaxis, cnpix, checknum, ival;

  naxis = set->naxis;
  lng = set->lng;
  lat = set->lat;
  rfilename = set->field->rfilename;

/* Build a histogram of star positions */
  cnaxis = 2;
/* Set the histogram dimensions (must be powers of 2) */
/* The 0.5 factor is to allow for the 2x wrapping below */
  ival = (int)(1.4*2*set->radius/(matchresol*LL_WRAPPING));
  for (i=1; (ival>>=1); i<<=1);
  csize[0] = csize[1] = i;
  if (csize[0]>LL_MAXSIZE)
    csize[0] = csize[1] = LL_MAXSIZE;
  lopass[0] = lopass[1] = 1.0/LLLOPASS_LAMBDA;
  hipass[0] = hipass[1] = 1.0/(csize[1]*LLHIPASS_LAMBDA);
  cnpix = 1;
  for (i=0; i<cnaxis; i++)
    cnpix *= csize[i];
  QCALLOC(rhisto, float, cnpix);
  QCALLOC(histo, float, cnpix);
  poserr = refset->radius/set->wcsscale[lng] - set->wcs->naxisn[lng]/2.0;
  if (poserr<=0.0)
    poserr = 0.0;
/* We allow for a better precision by wrapping */
  xfac = LL_WRAPPING*(double)csize[0] / (set->wcs->naxisn[lng]+2*poserr);
  xzero = poserr*xfac;
  poserr = refset->radius/set->wcsscale[lat] - set->wcs->naxisn[lat]/2.0;
  if (poserr<=0.0)
    poserr = 0.0;
/* We allow for a better precision by wrapping */
  yfac = LL_WRAPPING*(double)csize[1] / (set->wcs->naxisn[lat]+2*poserr);
  yzero = poserr*yfac;
  refsample = refset->sample;
/* Keep only the brightest sources if too many samples in reference set */
  if (refset->nsample>LL_NSOURCEMAX)
    {
    qsort(refset->sample, refset->nsample, sizeof(samplestruct), compmag);
    nrefsample = LL_NSOURCEMAX;
    }
  else
    nrefsample = refset->nsample;

  wcs = set->wcs;
/* Compute projected positions of the reference stars for the current set */
#ifndef USE_THREADS
  sprintf(str,"%.24s: projecting reference catalog...", rfilename);
  NFPRINTF(OUTPUT, str);
#endif
  compute_rawpos(wcs, refsample, nrefsample);

/* First build histogram for the reference field */
  n = 0;
  for (i=0; i<nrefsample; i++)
    {
#ifdef USE_THREADS
    if (!(++n%100000) && prefs.nthreads<2)
#else
    if (!(++n%100000))
#endif
      {
      sprintf(str,"%.24s: reference source #%d / %d processed",
		rfilename, n,nrefsample);
      NFPRINTF(OUTPUT, str);
      }
    x = refsample[i].rawpos[lng]*xfac+xzero;
    y = refsample[i].rawpos[lat]*yfac+yzero;
    xsig = refsample[i].rawposerr[lng]*xfac;
    ysig = refsample[i].rawposerr[lat]*yfac;
/*-- Take magnitudes into account */
    flux = -0.4*LL_FLUXEXP*refsample[i].mag;
    if (xsig>GAUSS_MAXSIG)
      {
      flux *= GAUSS_MAXSIG / xsig;
      xsig = GAUSS_MAXSIG;
      }
    if (ysig>GAUSS_MAXSIG)
      {
      flux *= GAUSS_MAXSIG / ysig;
      ysig = GAUSS_MAXSIG;
      }
    flux = flux<70.0? (flux>-70.0? DEXP(flux) : 0.0) : BIG;
    putgauss(rhisto, csize[0], csize[1], x,y, flux, xsig,ysig);
    }

/* Check-image for longitude-latitude histogram of the reference catalog */
  if ((checknum=check_check(CHECK_LLREFPAIR))>=0)
    write_check(prefs.check_name[checknum], rhisto, csize[0], csize[1]);

/* Now the current field */
  sample = set->sample;
/* Keep only the brightest sources if too many samples in set */
  if (set->nsample>LL_NSOURCEMAX)
    {
    qsort(set->sample, set->nsample, sizeof(samplestruct), compmag);
    nsample = LL_NSOURCEMAX;
    }
  else
    nsample = set->nsample;

  n = 0;
  for (i=0; i<nsample; i++)
    {
#ifdef USE_THREADS
    if (!(++n%100000) && prefs.nthreads<2)
#else
    if (!(++n%100000))
#endif
      {
      sprintf(str,"%.24s: detection #%d / %d processed", rfilename, n,nsample);
      NFPRINTF(OUTPUT, str);
      }
    x = sample[i].rawpos[lng]*xfac+xzero;
    y = sample[i].rawpos[lat]*yfac+yzero;
    xsig = sample[i].rawposerr[lng]*xfac;
    ysig = sample[i].rawposerr[lat]*yfac;
/*-- Take magnitudes into account */
    flux = sample[i].flux>0.0? pow(sample[i].flux, LL_FLUXEXP): 0.0;
    putgauss(histo, csize[0], csize[1], x,y, flux, xsig,ysig);
    }

/* Check-image for longitude-latitude histogram of the reference catalog */
  if ((checknum=check_check(CHECK_LLPAIR))>=0)
    write_check(prefs.check_name[checknum], histo, csize[0], csize[1]);

/* Make the cross-correlation */
#ifndef USE_THREADS
  sprintf(str,"%.24s: finding field coordinate shift...", rfilename);
  NFPRINTF(OUTPUT, str);
#endif

  fastcorr(histo, rhisto, cnaxis, csize, lopass, hipass);

  sig = findcrosspeak(histo, csize[0], csize[1],
		prefs.position_maxerr[0]/set->wcsscale[lng]*xfac,
		prefs.position_maxerr[1]/set->wcsscale[lat]*yfac,
		1.0/lopass[0], 1.0/lopass[1],
		&x, &y);

/* Check-image for longitude-latitude correlation of the reference catalog */
  if ((checknum=check_check(CHECK_LLXCORR))>=0)
    write_check(prefs.check_name[checknum], histo, csize[0], csize[1]);

  if (sig>0.0)
    {
/*-- The new CRVAL corresponds approximately to a shifted raw coordinate */
    for (k=0; k<naxis; k++)
      rawpos[k] = wcs->crpix[k];
    rawpos[lng] -= x/xfac;
    rawpos[lat] -= y/yfac;
    raw_to_wcs(wcs, rawpos, wcspos);
    *dlng = wcspos[lng] - wcs->crval[lng];
    *dlat = wcspos[lat] - wcs->crval[lat];
    }
  else
    *dlng = *dlat = 0.0;

  free(rhisto);
  free(histo);

  return sig;
  }


/****** compmag *************************************************************
PROTO   int (*compmag(const void *sample1, const void *sample2))
PURPOSE Provide a sample magnitude comparison function for qsort()
INPUT   pointer to 1st sample,
	pointer to 2nd sample.
OUTPUT  <0 if mag1<mag2, >0 if mag1>mag2, 0 otherwise .
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 27/08/2002
*/
int	compmag(const void *sample1, const void *sample2)
  {
   float	dm;

  dm = ((samplestruct *)sample1)->mag - ((samplestruct *)sample2)->mag;

  return dm>0.0 ? 1: (dm<0.0? -1 : 0);
  }


/****** compraw **************************************************************
PROTO   int (*compraw_samples(const void *sample1, const void *sample2))
PURPOSE Provide a sample comparison function for qsort() based on
	raw coordinates.
INPUT   pointer to 1st sample,
	pointer to 2nd sample.
OUTPUT  <0 if y1<y2, >0 if y1>y2, 0 otherwise .
NOTES   Not fully reentrant because of the sort_coord global variable.
AUTHOR  E. Bertin (IAP)
VERSION 07/02/2005
*/
int	compraw(const void *sample1, const void *sample2)
  {
   double	dy;

  dy = ((samplestruct *)sample1)->rawpos[sort_coord]
	- ((samplestruct *)sample2)->rawpos[sort_coord];

  return dy>0.0 ? 1: (dy<0.0? -1 : 0);
  }


/****** putgauss ***********************************************************
PROTO	void putgauss(float *histo, int width, int height, double x, double y,
		double flux, double xsig, double ysig)
PURPOSE	Add Gaussian profiles to a 2D histogram.
INPUT	ptr to the histogram,
	histogram width,
	histogram height,
	x position of the Gaussian,
	y position of the Gaussian,
	sigma_x of the Gaussian,
	sigma_y of the Gaussian.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	22/09/2006
 ***/
void	putgauss(float *histo, int width, int height, double x, double y,
		double flux, double xsig, double ysig)
  {
   double	dx,dx0,dy, xalpha,yalpha;
   float	sb,a;
   int		i,j, vwidth,vheight, ix,ix0,iy;

/* Prepare the vignet */
  if (xsig<0.5)
    xsig = 0.5;
  vwidth = (int)(2*GAUSS_MAXNSIG*xsig);
  if (!(vwidth&1))
   vwidth++;
  if (ysig<0.5)
    ysig = 0.5;
  vheight = (int)(2*GAUSS_MAXNSIG*ysig);
  if (!(vheight&1))
   vheight++;

  if (vwidth<=0 || vwidth>width || vheight<=0 || vheight>height)
    return;

/* Generate the pseudo-Gaussian */
  xalpha = 1.0/(xsig*xsig);
  yalpha = 1.0/(ysig*ysig);
  ix = (int)(x-vwidth/2+0.5);
  iy = (int)(y-vheight/2+0.5);
  sb = flux/(xsig*ysig);
  dy = iy - y;
  dx0 = ix - x;
  ix0 = ix%width;
  if (ix0<0)
    ix0 += width;
  iy = iy%height;
  if (iy<0)
    iy += height;
  for (j=vheight; j--; dy+=1.0, iy = (iy+1)%height)
    {
    dx = dx0;
    ix = ix0;
    for (i=vwidth; i--; dx+=1.0, ix++)
      {
      if ((a = (GAUSS_MAXNSIG*GAUSS_MAXNSIG - xalpha*dx*dx-yalpha*dy*dy)) > 0.0)
        histo[iy*width+(ix%width)] += sb*a*a;
      }
    }


  return;
  }


/****** flipimage ***********************************************************
PROTO	void flipimage(float *histo, int width, int height)
PURPOSE	Flip an image over the X axis.
INPUT	ptr to the data,
	histogram width,
	histogram height.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	04/01/2003
 ***/
void	flipimage(float *histo, int width, int height)
  {
   unsigned int	*point0, *pointl, *pointh;
   int		i, j, w2, step;

  w2 = width/2;
  step = width-1;
  point0 = (unsigned int *)histo;
  for (j=height; j--; point0 += width)
    {
    pointl = point0;
    pointh = pointl+step;
    for (i=w2; i--; pointl++, pointh--)
      *pointl ^= (*pointh^=(*pointl^=*pointh));
    }

  return;
  }


/****** findcrosspeak ********************************************************
PROTO	double findcrosspeak(float *histo, int width, int height,
			double xrange, double yrange,
			double xresol, double yresol,
			double *xpeak, double *ypeak)
PURPOSE	Find the pixel with the peak value in a cross-correlation image
INPUT	ptr to the histogram,
	histogram width,
	histogram height,
	search window width,
	search window height,
	data resolution in x,
	data resolution in y,
	pointer to the x position of the maximum pixel (output),
	pointer to the y position of the maximum pixel (output),
OUTPUT	Significance of the peak.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	28/04/2006
 ***/
double	findcrosspeak(float *histo, int width, int height,
			double xrange, double yrange,
			double xresol, double yresol,
			double *xpeak, double *ypeak)

  {
   double	xstep,ystep, xmean,ymean, dx,dy, dxpeak,dypeak, dval, dsum;
   float	val, max,max1;
   int		i, ix,iy,iy2, ixmax,iymax, xmin,xmax,ymin,ymax, resol2,
		ixpeak,iypeak, idx,idy;

/* Apply a fudge factor to low-pass filter cut-off wavelength */
  resol2 = (int)(4.0*(xresol*xresol+yresol*yresol))+1;

/* Find the maximum */
  xmin = width/2 - (int)(xrange+0.499);
  if (xmin < 0)
    xmin = 0;
  xmax = width/2 + (int)(xrange+1.499);
  if (xmax > width)
    xmax = width;
  ymin = height/2 - (int)(yrange+0.499);
  if (ymin < 0)
    ymin = 0;
  ymax = height/2 + (int)(yrange+1.499);
  if (ymax > height)
    ymax = height;

  ix = ixmax = iymax = 0;	/* To avoid gcc -Wall warnings */
  max = -BIG;
  for (iy=ymin; iy<ymax; iy++)
    for (ix=xmin; ix<xmax; ix++)
      if ((val=histo[ix+iy*width]) > max)
        {
        max = val;
        ixmax = ix;
        iymax = iy;
        }
  max1 = max;
  ixpeak = ixmax;
  iypeak = iymax;

/* Fit the maximum */
  xstep = 1.0/xresol;	/* Includes a hidden fudge factor: sqrt(2) */
  ystep = 1.0/yresol;	/* Includes a hidden fudge factor: sqrt(2) */
  xmin = ixpeak - (idx=(int)(3.0/xstep+0.499));
  xmax = xmin + 2*idx + 1;
  ymin = iypeak - (idy=(int)(3.0/ystep+0.499));
  ymax = ymin + 2*idy + 1;
  dxpeak = dypeak = 0.0;
  for (i=0; i<PEAKFIND_NITER; i++)
    {
/*-- Find the barycenter */
    xmean = ymean = dsum = 0.0;
    dy = (ymin - iypeak - dypeak)*ystep;
    for (iy=ymin; iy<ymax; iy++)
      {
      dx = (xmin - ixpeak - dxpeak)*xstep;
      iy2 = ((iy+height)%height)*width;
      for (ix=xmin; ix<xmax; ix++)
        {
        dsum += (dval = histo[(ix+width)%width+iy2]*exp(-(dx*dx+dy*dy)));
        xmean += dval*dx;
        ymean += dval*dy;
        dx += xstep;
        }
      dy += ystep;
      }
    dxpeak = xmean/(xstep*dsum);
    dypeak = ymean/(ystep*dsum);
    }

  *xpeak = (double)(ixpeak - width/2) + dxpeak;
  *ypeak = (double)(iypeak - height/2) + dypeak;

  max = -BIG;
/* Find a second maximum */
  for (iy=0; iy<height; iy++)
    {
    for (ix=0; ix<width; ix++)
      if ((val=histo[ix+iy*width]) > max)
        {
        idx = fabs(ix-ixpeak);
        idy = fabs(iy-iypeak);
        if (idx*idx + idy*idy > resol2
		&& (width-idx)*(width-idx) + (height-idy)*(height-idy) > resol2)
          max = val;
        }
    }
  max = fabs(max);

  return(max>0.0 ? (double)max1/max : (double)max1);
  }


/****** match_refine *********************************************************
PROTO	void match_refine(setstruct *set, setstruct *refset, double matchresol,
			double *angle, double *scale, double *sangle,
			double *ratio, *dlng, double *dlat)
PURPOSE	Refine the celestial coordinate transformation of a set of catalogs 
	relative to an astrometric reference catalog using pattern matching.
INPUT	ptr to the set to be matched,
	ptr to the reference set,
	MATCHing resolution in arcsec,
	ptr to the rotation angle in deg (output),
	ptr to the scaling factor (output),
	ptr to shear angle in deg (output),
	ptr to the shear factor (output),
	ptr to the longitude shift (output),
	ptr to the latitude shift (output).
OUTPUT	Confidence level of the solution (in units of sigma).
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	20/03/2007
 ***/
void	match_refine(setstruct *set, setstruct *refset, double matchresol,
			double *angle, double *scale,
			double *sangle, double *shear,
			double *dlng,double *dlat)
  {
   wcsstruct	*wcs;
   samplestruct	*refsample, *samp1,*samp2,*samp2b,*samp2min;
#ifndef USE_THREADS
   char		str[MAXCHAR];
#endif
   double	alpha[9], blng[3], blat[3],
		rawpos[NAXIS], wcspos[NAXIS],crpix[NAXIS],
		sig, rlim, r2,r2min, reflng,reflat,
		lng1,lat1,lng2,lat2, latmin,latmax,
		dlng12,dlat12, wi, x,y, det,
		a11,a12,a21,a22, b11, c11,c12,c21,c22;
   int		d,i, nrefsample, nsample, nsamp2, nsamp2b, lng,lat,
		naxis;

  lng = set->lng;
  lat = set->lat;
  naxis = set->naxis;

#ifndef USE_THREADS
  sprintf(str,"%.24s: refining matching parameters...", set->field->rfilename);
  NFPRINTF(OUTPUT, str);
#endif

  wcs = set->wcs;
/* Compute local pixel scale and write down the crpix vector */
  for (d=0; d<naxis; d++)
    {
    rawpos[d] = wcs->naxisn[d]/2.0;
    crpix[d] = wcs->crpix[d];
    }
  matchresol /= sqrt(wcs_scale(wcs, rawpos));

/* Keep only the brightest sources if too many samples in reference set */
  if (refset->nsample>LL_NSOURCEMAX)
    {
    qsort(refset->sample, refset->nsample, sizeof(samplestruct), compmag);
    nrefsample = LL_NSOURCEMAX;
    }
  else
    nrefsample = refset->nsample;

/* Compute projected positions of the reference stars for the current set */
  refsample = refset->sample;
  for (i=nrefsample; i--; refsample++)
    {
    wcs_to_raw(wcs, refsample->wcspos, refsample->rawpos);
/*-- Compute position errors in pixel coordinates, too */
    refsample->rawposerr[lng] = refsample->rawposerr[lat]
	= sqrt(refsample->wcsposerr[lng]*refsample->wcsposerr[lat]
		/wcs_scale(wcs, refsample->rawpos));
    }

/* Sort samples to accelerate cross-matching */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&matchsortmutex);
#endif
  sort_coord = lat;
  qsort(refset->sample, nrefsample, sizeof(samplestruct), compraw);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&matchsortmutex);
#endif

/* Now the current field */
/* Keep only the brightest sources if too many samples in set */
  if (set->nsample>LL_NSOURCEMAX)
    {
    qsort(set->sample, set->nsample, sizeof(samplestruct), compmag);
    nsample = LL_NSOURCEMAX;
    }
  else
    nsample = set->nsample;

/* Sort samples to accelerate cross-matching */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&matchsortmutex);
#endif
  sort_coord = lat;
  qsort(set->sample, nsample, sizeof(samplestruct), compraw);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&matchsortmutex);
#endif

/* Zero the matrices */
  memset(alpha, 0, 9*sizeof(double));
  blng[0] = blng[1] = blng[2] = blat[0] = blat[1] = blat[2] = 0.0;

/* Cross-match! */
  rlim = matchresol;
  r2min = rlim*rlim;
  nsamp2b = nsample;
  samp1 = refset->sample;
  samp2b = set->sample;
  for (i=nrefsample; i--; samp1++)
    {
    reflng = samp1->rawpos[lng];
    reflat = samp1->rawpos[lat];
    lng1 = reflng - crpix[lng];
    lat1 = reflat - crpix[lat];
    latmin = reflat-rlim;
    latmax = reflat+rlim;
    samp2min = NULL;
    samp2 = samp2b;
/*-- Jump over sources that can't match in y */
    for (nsamp2=nsamp2b; nsamp2-- && samp2->rawpos[lat]<latmin; samp2++);
    samp2b = samp2;
    nsamp2b = ++nsamp2;
/*-- Find the brightest match */
    for (; nsamp2-- && samp2->rawpos[lat]<latmax; samp2++)
      {
      dlng12 = samp2->rawpos[lng] - reflng;
      dlat12 = samp2->rawpos[lat] - reflat;
      r2 = dlng12*dlng12 + dlat12*dlat12;
      if (r2<r2min && (!samp2min || samp2->flux>samp2min->flux))
        samp2min = samp2;
      }
/*-- Fill the normal equation matrix to fit trends in longitude and latitude */
    if (samp2min)
      {
      sig = samp1->rawposerr[lng]*samp1->rawposerr[lng]
	+ samp2min->rawposerr[lng]*samp2min->rawposerr[lng]
	+ 1e-2;	/* 0.1 pixel RMS error added to reduce the impact of glitches */
      if (sig>0.0)
        {
        wi = 1.0/sig;
        lng2 = samp2min->rawpos[lng] - crpix[lng];
        lat2 = samp2min->rawpos[lat] - crpix[lat];
        alpha[0] += wi*lng2*lng2;
        alpha[1] += wi*lng2*lat2;
        alpha[2] += wi*lng2;
        alpha[4] += wi*lat2*lat2;
        alpha[5] += wi*lat2;
        alpha[8] += wi;
        dlng12 = lng2 - lng1;
        dlat12 = lat2 - lat1;
        blng[0] += wi*lng2*dlng12;
        blng[1] += wi*lat2*dlng12;
        blng[2] += wi*dlng12;
        blat[0] += wi*lng2*dlat12;
        blat[1] += wi*lat2*dlat12;
        blat[2] += wi*dlat12;
        }
      }
    }

  alpha[3] = alpha[1];
  alpha[6] = alpha[2];
  alpha[7] = alpha[5];

/* Solve the system */
  if (clapack_dpotrf(CblasRowMajor, CblasUpper, 3, alpha, 3) == 0)
    {
    clapack_dpotrs(CblasRowMajor, CblasUpper, 3, 1, alpha, 3, blng, 3);
    clapack_dpotrs(CblasRowMajor, CblasUpper, 3, 1, alpha, 3, blat, 3);
    a11 = blng[0] + 1.0;
    a12 = blng[1];
    a21 = blat[0];
    a22 = blat[1] + 1.0;
    c11 = wcs->cd[lng*naxis+lng];
    c12 = wcs->cd[lng*naxis+lat];
    c21 = wcs->cd[lat*naxis+lng];
    c22 = wcs->cd[lat*naxis+lat];
    det = a11*a22 - a12*a21;
    if (fabs(det) < 1.0/BIG)
      {
      warning("Null CD determinant in ", "MATCHing refinement");
      det = 1.0;
      }

/*-- Replace A with the A' matrix such as A'.CD = CD.A, that is A' = CD.A.CD-1 */
    b11 = a22/det;
    a12 = -a12/det;
    a21 = -a21/det;
    a22 = a11/det;
    a11 = b11;

    x = sqrt((a11+a22)*(a11+a22) + (a12-a21)*(a12-a21));
    y = sqrt((a11-a22)*(a11-a22) + (a12+a21)*(a12+a21));
    *scale = 0.5*(x - y);
    *shear = 2.0 * y / (x - y);
    *angle = atan2(a21-a12, a11+a22) / DEG;
    *sangle = 0.5*atan2(2.0*(a11*a21+a12*a22),
		a11*a11-a12*a12-a22*a22+a21*a21)/ DEG;

/*-- The new CRVAL corresponds approximately to a shifted raw coordinate */
    for (d=0; d<naxis; d++)
      rawpos[d] = wcs->crpix[d];
    rawpos[lng] += blng[2];
    rawpos[lat] += blat[2];
    raw_to_wcs(wcs, rawpos, wcspos);

    *dlng = wcspos[lng] - wcs->crval[lng];
    *dlat = wcspos[lat] - wcs->crval[lat];
    }
  else
    {
/*-- No solution found: do not apply any change */
    *scale = 1.0;
    *shear = *angle = *sangle = *dlng = *dlat = 0.0;
    }

  return;
  }


/****** mean_rawposvar *******************************************************
PROTO	double mean_rawposvar(setstruct *set)
PURPOSE	Compute the mean variance in position of a set of sources.
INPUT	ptr to the input set.
OUTPUT	variance in position, or 0.0 if the input set is empty of sources.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	16/03/2006
 ***/
double	mean_rawposvar(setstruct *set)
  {
   samplestruct	*sample;
   double	poserr2;
   int		i, lng,lat;

  lng = set->lng;
  lat = set->lat;
  sample = set->sample;
  poserr2 = 0.0;
  for (i=set->nsample; i--; sample++)
    poserr2 += sample->rawposerr[lng]*sample->rawposerr[lat];

  return set->nsample? poserr2/set->nsample : 0.0;
  }


/****** compute_rawpos ******************************************************
PROTO	void compute_rawpos(wcsstruct *wcs, samplestruct *refsample,
			int nsample)
PURPOSE	Convert world positions and position errors into raw (pixel) units.
INPUT	ptr to the destination WCS,
	ptr to the input (reference) sample array,
	nb of samples.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	23/03/2006
 ***/
void	compute_rawpos(wcsstruct *wcs, samplestruct *refsample, int nrefsample)
  {
   int		i, lng,lat;

  lng = wcs->lng;
  lat = wcs->lat;
  for (i=nrefsample; i--; refsample++)
    {
    wcs_to_raw(wcs, refsample->wcspos, refsample->rawpos);
/*-- Compute position errors in pixel coordinates, too */
    refsample->rawposerr[lng] = refsample->rawposerr[lat]
	= sqrt(refsample->wcsposerr[lng]*refsample->wcsposerr[lat]
		/wcs_scale(wcs, refsample->rawpos));
   }

  return;
  }


/****** frame_set ***********************************************************
PROTO	setstruct *frame_set(setstruct *setin, wcsstruct *wcs,
		double *centerpos,double radius)
PURPOSE	Create a subset of the input set centered on given coordinates within
	a given radius.
INPUT	ptr to the input set,
	ptr to the array of coordinates,
	inclusion radius in degree.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	31/08/2002
 ***/
setstruct	*frame_set(setstruct *setin, wcsstruct *wcs,
			double *centerpos, double radius)
  {
   setstruct	*set;
   samplestruct	*sample;
   int		i, nsample;

/* Create a new set */
  set = init_set();
  sample = setin->sample;
  nsample = setin->nsample;
  for (i=nsample; i--; sample++)
   if (wcs_dist(wcs, centerpos, sample->wcspos) < radius)
     copy_samples(sample, set, 1);
  realloc_samples(set, set->nsample);
  set->radius = radius;

  return set;
  }


/****** update_wcsas *********************************************************
PROTO	void update_wcsas(wcsstruct *wcs, double angle, double scale)
PURPOSE	Update the scale, and position angle of a WCS projection.
INPUT	ptr to the WCS structure,
	position angle increment (CCW degrees in projected (not pixel)),
	scale factor.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/09/2006
 ***/
void	update_wcsas(wcsstruct *wcs, double angle, double scale)
  {
   double	a[NAXIS*NAXIS],b[NAXIS*NAXIS],
		*c,*at,
		val, cas, sas, fascale, sgn;
   int		i,j,k, lng,lat,naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;

  if (lng == lat)
    return;

  c = wcs->cd;

/* A = C*B */
/* The B matrix is made of 2 numbers */
  fascale = fabs(scale);
  sgn = (scale<0.0)? -1.0 : 1.0;
  if (wcs->chirality < 0)
    angle = -angle;
  cas = cos(angle*DEG)*fascale;
  sas = sin(angle*DEG)*fascale;
  for (i=0; i<naxis; i++)
    b[i+i*naxis] = 1.0;
  b[lng+lng*naxis] = cas;
  b[lat+lng*naxis] = -sas*sgn;
  b[lng+lat*naxis] = sas;
  b[lat+lat*naxis] = cas*sgn;
  at = a;
  for (j=0; j<naxis; j++)
    for (i=0; i<naxis; i++)
      {
      val = 0.0;
      for (k=0; k<naxis; k++)
        val += c[k+j*naxis]*b[i+k*naxis];
      *(at++) = val;
      }

  at = a;
  for (i=0; i<naxis*naxis; i++)
    *(c++) = *(at++);

/* Initialize other WCS structures */
  init_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);

  return;
  }


/****** update_wcsss *********************************************************
PROTO	void update_wcsss(wcsstruct *wcs, double sangle, double ratio)
PURPOSE	Update the shear angle, and aspect ratio of a WCS projection.
INPUT	ptr to the WCS structure,
	position angle of shear (CCW degrees in projected (not pixel)),
	expansion factor (<0.0 for contraction).
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	25/05/2005
 ***/
void	update_wcsss(wcsstruct *wcs, double sangle, double ratio)
  {
   double	a[NAXIS*NAXIS],b[NAXIS*NAXIS],
		*c,*at,
		val, cas, sas;
   int		i,j,k, lng,lat,naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;

  if (lng == lat)
    return;

  c = wcs->cd;
/* A = C*B */
/* The B matrix is made of 2 numbers */
  cas = cos(2*sangle*DEG);
  sas = sin(2*sangle*DEG);
  for (i=0; i<naxis; i++)
    b[i+i*naxis] = 1.0;
  b[lng+lng*naxis] = 1.0 + ratio*(1.0+cas)/2.0;
  b[lat+lng*naxis] = b[lng+lat*naxis] = ratio*sas/2.0;
  b[lat+lat*naxis] = 1.0 + ratio*(1.0-cas)/2.0;
  at = a;
  for (j=0; j<naxis; j++)
    for (i=0; i<naxis; i++)
      {
      val = 0.0;
      for (k=0; k<naxis; k++)
        val += c[k+j*naxis]*b[i+k*naxis];
      *(at++) = val;
      }

  at = a;
  for (i=0; i<naxis*naxis; i++)
    *(c++) = *(at++);

/* Initialize other WCS structures */
  init_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);

  return;
  }


/****** compute_wcsss *********************************************************
PROTO	void compute_wcsss(wcsstruct *wcs, double *sangle, double *shear)
PURPOSE	Compute the shear angle, and the amount of shear in a WCS projection.
INPUT	ptr to the WCS structure,
	position angle of shear (CCW degrees in projected (not pixel))(output),
	expansion factor (<0.0 for contraction) (output).
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	11/03/2005
 ***/
void	compute_wcsss(wcsstruct *wcs, double *sangle, double *shear)
  {
   double	a11,a12,a21,a22, x,y;
   int		lng,lat,naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;
  if (lng == lat)
    {
    lng = 0;
    lat = 1;
    }

  a11 = wcs->cd[lng*naxis+lng];
  a12 = wcs->cd[lng*naxis+lat];
  a21 = wcs->cd[lat*naxis+lng];
  a22 = wcs->cd[lat*naxis+lat];
/* Handle flipped case */
  if (a11*a22-a12*a21<0.0)
    {
    a11 = -a11;
    a12 = -a12;
    }
  x = sqrt((a11+a22)*(a11+a22) + (a12-a21)*(a12-a21));
  y = sqrt((a11-a22)*(a11-a22) + (a12+a21)*(a12+a21));
  *shear = 2.0 * y / (x - y);
  *sangle = fmod(90.0 + 0.5*atan2(2.0*(a11*a21+a12*a22),
		a11*a11+a12*a12-a22*a22-a21*a21) / DEG, 90.0);

  return;
  }


/****** update_wcscc *********************************************************
PROTO	void update_wcscc(wcsstruct *wcs, double drawlng, double drawlat)
PURPOSE	Update the projection center of a WCS projection.
INPUT	ptr to the WCS structure,
	longitude projected increment (pixel),
	latitude projected increment (pixel).
OUTPUT	-.
NOTES	The new celestial position is an approximation of the exact one.
AUTHOR	E. Bertin (IAP)
VERSION	23/09/2003
 ***/
void	update_wcscc(wcsstruct *wcs, double drawlng, double drawlat)
  {
   double	wcspos[NAXIS];

   int		i, lng,lat, naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;

  if (lng == lat)
    return;

  for (i=0; i<naxis; i++)
    wcspos[i] = wcs->crval[i];
  wcspos[lng] -= drawlng;
  wcspos[lat] -= drawlat;

  wcs_to_raw(wcs, wcspos, wcs->crpix);
/* Initialize other WCS structures */
  init_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);

  return;
  }


/****** update_wcsll *********************************************************
PROTO	void update_wcsll(wcsstruct *wcs, double dlng, double dlat)
PURPOSE	Update the celestial center of a WCS projection.
INPUT	ptr to the WCS structure,
	longitude raw increment (deg),
	latitude raw increment (deg).
OUTPUT	-.
NOTES	The new celestial position is an approximation of the exact one.
AUTHOR	E. Bertin (IAP)
VERSION	01/01/2004
 ***/
void	update_wcsll(wcsstruct *wcs, double dlng, double dlat)
  {
   double	wcspos[NAXIS], a[NAXIS*NAXIS],b[NAXIS*NAXIS],
		*c,*at,
		val, cas, sas, angle, dalpha;
   int		i,j,k, lng,lat, naxis;

  lng = wcs->lng;
  lat = wcs->lat;
  naxis = wcs->naxis;

  if (lng == lat)
    return;

  for (i=0; i<naxis; i++)
    wcspos[i] = wcs->crval[i];

  wcs->crval[lng] += dlng;
  wcs->crval[lat] += dlat;

  dalpha = (wcs->crval[lng] - wcspos[lng])*DEG;

/* Compute difference angle with the north axis between start and end */
  angle = (dlng!=0.0 && dlat!= 0.0) ?
	180.0 - (atan2(sin(dalpha),
	cos(wcs->crval[lat]*DEG)*tan(wcspos[lat]*DEG)
	- sin(wcs->crval[lat]*DEG)*cos(dalpha))
	+ atan2(sin(dalpha),
	cos(wcspos[lat]*DEG)*tan(wcs->crval[lat]*DEG)
	- sin(wcspos[lat]*DEG)*cos(dalpha)))/DEG
	: 0.0;

/* A = C*B */
  c = wcs->cd;
/* The B matrix is made of 2 numbers */
  cas = cos(-angle*DEG);
  sas = sin(-angle*DEG);

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
        val += c[k+j*naxis]*b[i+k*naxis];
      *(at++) = val;
      }

  at = a;
  for (i=0; i<naxis*naxis; i++)
    *(c++) = *(at++);

/* Initialize other WCS structures */
  init_wcs(wcs);
/* Find the range of coordinates */
  range_wcs(wcs);
/* Invert projection corrections */
  invert_wcs(wcs);

  return;
  }


/****** print_matchinfo *******************************************************
PROTO	void print_matchinfo(fieldstruct *field)
PURPOSE	Print field information after the pattern matching process.
INPUT	ptr to the field that has been matched.
OUTPUT	-.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	22/03/2006
 ***/
void	print_matchinfo(fieldstruct *field)
  {
   setstruct	*set;
   char		str[128];
   int		s;

  if (field->mosaic_type == MOSAIC_LOOSE)
    {
/*-- Print each set result independently */
    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT, "Catalog %s:\n", field->rfilename);
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      if (set->lat != set->lng)
        {
        sprintf(str, "%-11.11s[%2d/%-2d] A%-2d P%-2d%#+6.2f deg "
		"%#7.4g\"%#7.3g %#+8.2g\"%#+8.2g\"%#7.3g",
		field->rfilename, s+1, field->nset,
		field->astromlabel+1,field->photomlabel+1,
		set->match_dangle,
		set->wcsscale[0]/set->match_dscale*DEG/ARCSEC,
		set->match_asig,
		fabs(set->match_dlng)>1e-10?
		  set->match_dlng*cos(set->wcspos[set->lat]*DEG)*DEG/ARCSEC
		  : 0.0,
		fabs(set->match_dlat)>1e-10? set->match_dlat*DEG/ARCSEC : 0.0,
		set->match_sig);
        if (set->match_sig<MATCH_MINCONT)
          {
          QBPRINTF(OUTPUT, str);
          }
        else
          {
          QPRINTF(OUTPUT, "%s\n", str);
          }
        }
      else
        {
        sprintf(str, "%-11.11s[%2d/%-2d] A%-2d P%-2d%#+6.2f deg "
		"%#7.4gx%#7.3g %#+8.2g %#+8.2g %#7.3g",
		field->rfilename, s+1, field->nset,
		field->astromlabel+1,field->photomlabel+1,
		set->match_dangle,
		set->wcsscale[0]/set->match_dscale,
		set->match_asig,
		set->match_dlng,
		set->match_dlat,
		set->match_sig);
        if (set->match_sig<MATCH_MINCONT)
          {
          QBPRINTF(OUTPUT, str);
          }
        else
          {
          QPRINTF(OUTPUT, "%s\n", str);
          }
	}
      }
    }
  else
    {
    if (field->lat != field->lng)
      {
      sprintf(str, "%-18.18s A%-2d P%-2d%#+6.2f deg "
		"%#7.4g\"%#7.3g %#+8.2g\"%#+8.2g\"%#7.3g",
		field->rfilename,
		field->astromlabel+1,field->photomlabel+1,
		field->match_dangle,
		field->meanwcsscale[0]*DEG/ARCSEC,
		field->match_asig,
		fabs(field->match_dlng)>1e-10?
		 field->match_dlng*cos(field->meanwcspos[field->lat]*DEG)
		 *DEG/ARCSEC : 0.0,
		fabs(field->match_dlat)>1e-10?field->match_dlat*DEG/ARCSEC:0.0,
		field->match_sig);
      if (field->match_sig<MATCH_MINCONT)
        {
        QBPRINTF(OUTPUT, str);
        }
      else
        {
        QPRINTF(OUTPUT,"%s\n", str);
        }
      }
    else
      {
      sprintf(str, "%-18.18s A%-2d P%-2d%#+6.2f deg "
		"%#7.4gx%#7.3g %#+8.2g %#+8.2g %#7.3g",
		field->rfilename,
		field->astromlabel+1,field->photomlabel+1,
		field->match_dangle,
		field->meanwcsscale[0],
		field->match_asig,
		field->match_dlng,
		field->match_dlat,
		field->match_sig);
      if (field->match_sig<MATCH_MINCONT)
        {
        QBPRINTF(OUTPUT, str);
        }
      else
        {
        QPRINTF(OUTPUT, "%s\n", str);
        }
      }
    }

  return;
  }


#ifdef USE_THREADS

/****** pthread_match_field ***************************************************
PROTO   void *pthread_match_field(void *arg)
PURPOSE thread that takes care of matching fields.
INPUT   Pointer to the thread number.
OUTPUT  -.
NOTES   Relies on global variables.
AUTHOR  E. Bertin (IAP)
VERSION 27/12/2006
 ***/
void    *pthread_match_field(void *arg)
  {
   int	findex, gindex, proc;

  gindex = findex = -1;
  proc = *((int *)arg);
  threads_gate_sync(pthread_startgate);
  while (!pthread_endflag)
    {
    QPTHREAD_MUTEX_LOCK(&matchmutex);
    if (findex>-1)
/*---- Indicate that the field info is now suitable for viewing */
      pthread_fviewflag[gindex][findex] = 1;
    while (pthread_gviewindex<pthread_ngroup
	&& pthread_fviewflag[pthread_gviewindex][pthread_fviewindex])
      {
      if (!pthread_fviewindex)
        {
        NFPRINTF(OUTPUT, "");
        QPRINTF(OUTPUT, " Group %2d: %8d standard%s in %s (band %s)\n",
		pthread_gviewindex+1,
		pthread_reffields[pthread_gviewindex]->nsample,
		pthread_reffields[pthread_gviewindex]->nsample>1?"s":"",
		astrefcat[prefs.astrefcat].name,
		astrefcat[prefs.astrefcat].bandname);
        QIPRINTF(OUTPUT, "              instruments  pos.angle   scale    "
		"cont.        shift        cont.");
        }
      print_matchinfo(
	pthread_fgroups[pthread_gviewindex]->field[pthread_fviewindex++]);
      if (pthread_fviewindex >= pthread_fgroups[pthread_gviewindex]->nfield)
        {
        pthread_fviewindex = 0;
        pthread_gviewindex++;
        }
      }
    if (pthread_gindex<pthread_ngroup)
      {
      gindex = pthread_gindex;
      findex = pthread_findex++;
      if (pthread_findex>=pthread_fgroups[gindex]->nfield)
        {
        pthread_findex = 0;
        pthread_gindex++;
        }
      QPTHREAD_MUTEX_UNLOCK(&matchmutex);
/*---- Match field */
      match_field(pthread_fgroups[gindex]->field[findex],
		pthread_reffields[gindex]);
      }
    else
      {
      QPTHREAD_MUTEX_UNLOCK(&matchmutex);
/*---- Wait for the input buffer to be updated */
      threads_gate_sync(pthread_stopgate);
/* ( Master thread process loads and saves new data here ) */
      threads_gate_sync(pthread_startgate);
      }
    }

  pthread_exit(NULL);

  return (void *)NULL;
  }


/****** pthread_match_fields **************************************************
PROTO   void pthread_match_fields(fgroupstruct **fgroups,
				fieldstruct **reffields, int ngroup)
PURPOSE	Perform (linear) pattern matching between linked-sets of catalogs and
	an astrometric reference catalog in parallel. The astrometric
	parameters are updated accordingly.
INPUT	ptr to the array of groups of fields to be matched,
	ptr to the array of reference fields,
	number of groups involved.
OUTPUT	-.
NOTES	Uses the global preferences.
AUTHOR	E. Bertin (IAP)
VERSION	28/09/2004
 ***/
void	pthread_match_fields(fgroupstruct **fgroups,
				fieldstruct **reffields, int ngroup)
  {
   static pthread_attr_t	pthread_attr;
   int				*proc,
				g, p;

/* Number of active threads */
  nproc = prefs.nthreads;
  pthread_fgroups = fgroups;
  pthread_reffields = reffields;
  pthread_ngroup = ngroup;
/* Compute the total number of fields to consider */
  QMALLOC(pthread_fviewflag, int *, ngroup);
  for (g=0; g<ngroup; g++)
    {
    QCALLOC(pthread_fviewflag[g], int, fgroups[g]->nfield);
    }
/* Set up multi-threading stuff */
  QMALLOC(proc, int, nproc);
  QMALLOC(thread, pthread_t, nproc);
  QPTHREAD_MUTEX_INIT(&matchsortmutex, NULL);
  QPTHREAD_MUTEX_INIT(&matchmutex, NULL);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_startgate = threads_gate_init(nproc+1, NULL);
  pthread_stopgate = threads_gate_init(nproc+1, NULL);
/* Start the reading threads */
  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_match_field, &proc[p]);
    }
  QPTHREAD_MUTEX_LOCK(&matchmutex);
  pthread_findex=pthread_gindex = pthread_fviewindex=pthread_gviewindex = 0;
  pthread_endflag = 0;
  QPTHREAD_MUTEX_UNLOCK(&matchmutex);
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
  QPTHREAD_MUTEX_DESTROY(&matchsortmutex);
  QPTHREAD_MUTEX_DESTROY(&matchmutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  for (g=0; g<ngroup; g++)
    free(pthread_fviewflag[g]);
  free(pthread_fviewflag);
  free(proc);
  free(thread);
  }

#endif

