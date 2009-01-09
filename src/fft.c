/*
                                  fft.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        SCAMP
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Routines dealing with FFT.
*
*       Last modify:    29/09/2004
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

#include <fftw3.h>

#include "define.h"
#include "globals.h"
#include "fft.h"
#include "fits/fitscat.h"
#include "prefs.h"
#ifdef USE_THREADS
#include "threads.h"
#endif

 int	firsttimeflag;

#ifdef USE_THREADS
 pthread_mutex_t	fftmutex;
#endif

/****** fastcorr ************************************************************
PROTO	void fastcorr(float *data1, float *data2, int naxis, int *size,
		double *lambda_lopass, double *lambda_hipass)
PURPOSE	Optimized n-dimensional FFT cross-correlation using the FFTW library.
INPUT	ptr to the first data-cube,
	ptr to the second data-cube to correlate with,
	number of axes,
	ptr to the array of axis sizes,
	cut-wavelengths of the lopass-filter,
	cut-wavelengths of the highpass-filter.
OUTPUT	-.
NOTES	For data1 and data2, memory must be allocated for
	size[0]* ... * 2*(size[naxis-1]/2+1) floats (padding required).
	This function works in up to 10 dimensions.
AUTHOR	E. Bertin (IAP)
VERSION	29/09/2004
 ***/
void	fastcorr(float *data1, float *data2, int naxis, int *size,
		double *lambda_lopass, double *lambda_hipass)
  {
   fftwf_plan	plan;
   double	lambda_def[10],
		dx, r02l, r02h, r2l, r2h, passl,passh;
   float	*fdata1, *fdata2, *fdata1p,*fdata2p,
		real,imag, fac, pass;
   int		x[10], size2[10],
		d, i,j, npix,npix2, size0, x0;

  if (naxis>10)
    return;

/* Convert axis indexing to that of FFTW */
  npix = 1;
  for (j=0,i=naxis; i--;)
    npix *= (size2[j++] = size[i]);
  size0 = size[0];
  npix2 = npix/size0;
  npix2 *= 2 * (size0/2 + 1);

  /*  fft_init();*/	/* May have been called earlier (it does not matter) */

/* Forward FFT for data1 */
/* Unfortunately the whole process cannot be done totally in place... */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  QFFTWMALLOC(fdata1, float, npix2);
  plan = fftwf_plan_dft_r2c(naxis, (const int *)size2, data1,
	(fftwf_complex *)fdata1,
	FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftwf_destroy_plan(plan);

/*-- Forward FFT for data2 */
  QFFTWMALLOC(fdata2, float, npix2);
  plan = fftwf_plan_dft_r2c(naxis, (const int *)size2, data2,
	(fftwf_complex *)fdata2,
	FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftwf_destroy_plan(plan);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

/* Actual correlation (Fourier conjugate product) */
  fac = 1.0/(npix/2);  
  fdata1p = fdata1;
  fdata2p = fdata2;
  for (i=npix2/2; i--; fdata2p+=2)
    {
    real = *fdata1p**fdata2p + *(fdata1p+1)**(fdata2p+1);
    imag = *(fdata1p+1)**fdata2p - *fdata1p**(fdata2p+1);
    *(fdata1p++) = fac*real;
    *(fdata1p++) = fac*imag;
    }

/* Free one of the buffers */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  QFFTWFREE(fdata2);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

/* Filtering in (compressed) multidimensional Fourier space */
/* Set everything to 0.0 */
  for (d=0; d<naxis; d++)
    {
    x[d] = 0;
    lambda_def[d] = 0.0;
    }

  if (!lambda_lopass)
   lambda_lopass = lambda_def;
  if (!lambda_hipass)
    lambda_hipass = lambda_def;
  fdata1p = fdata1;
  for (j=npix/size0; j--;)
    {
    r02l = r02h = 0.0;
    for (d=1; d<naxis; d++)
      {
      if (x[d]<size[d]/2)
        {
        if (lambda_lopass[d]>0.0)
	  {
          dx = 2.0*x[d]/size[d] / lambda_lopass[d];
          r02l += dx*dx;
          }
        if (lambda_hipass[d]>0.0)
	  {
          dx = 2.0*x[d]/size[d] / lambda_hipass[d];
          r02h += dx*dx;
          }
        }
      else
        {
        if (lambda_lopass[d]>0.0)
	  {
          dx = 2.0*(1.0 - x[d]/(double)size[d]) / lambda_lopass[d];
          r02l += dx*dx;
          }
        if (lambda_hipass[d]>0.0)
	  {
          dx = 2.0*(1.0 - x[d]/(double)size[d]) / lambda_hipass[d];
          r02h += dx*dx;
          }
        }
      }

    x0 = 0;
    for (i=size0/2+1; i--; x0++)
      {
      r2l = r2h = 0.0;
      if (lambda_lopass[0]>0.0)
	{
        dx = 2.0*x0/size0 / lambda_lopass[0];
        r2l = r02l + dx*dx;
        }
      if (lambda_hipass[0]>0.0)
        {
        dx = 2.0*x0/size0 / lambda_hipass[0];
        r2h = r02h + dx*dx;
        }
      pass = 1.0;
      passl = LOPASS_SLOPE*(sqrt(r2l)-1.0);
      passl = (passl<70.0)? (passl>-70.0? (1/(1+exp(passl))) : 1.0) : 0.0;
      pass *= passl;
      passh = HIPASS_SLOPE*(1.0-sqrt(r2h));
      passh = (passh<70.0)? (passh>-70.0? (1/(1+exp(passh))) : 1.0) : 0.0;
      pass *= passh;
      *(fdata1p++) *= pass;
      *(fdata1p++) *= pass; 
      }
/*-- Update coordinate vector */
    for (d=1; d<naxis; d++)
      if ((++x[d])< size[d])
        break;
      else
        x[d] = 0;
    }

/* Reverse FFT */
#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  plan = fftwf_plan_dft_c2r(naxis, (const int *)size2,(fftwf_complex *)fdata1, 
	data1, 
	FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

  fftwf_execute(plan);

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&fftmutex);
#endif
  fftwf_destroy_plan(plan);

/* Free the fdata1 scratch array */
  QFFTWFREE(fdata1);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&fftmutex);
#endif

/* Put correlation 0-lag back at the cube center */
  shiftcube(data1, naxis, size);

  return;
  }


/****** shiftcube **********************************************************
PROTO	void shiftcube(float *data, int naxis, int *size)
PURPOSE	Shift an image on all axes over half their sizes.
INPUT	ptr to the data,
	number of axes,
	cu.
OUTPUT	-.
NOTES	Works only for width and height even.
AUTHOR	E. Bertin (IAP)
VERSION	09/11/2002
 ***/
void	shiftcube(float *data, int naxis, int *size)
  {
   unsigned int	*pointl, *pointh;
   int		d, i, j, nchunk, npix;

  for (d=1; d<=naxis; d++)
    {
    npix = 1;
    for (i=0; i<d;i++)
      npix *= size[i];
    nchunk = 1;
    for (i=d; i<naxis;i++)
      nchunk *= size[i];
    for (j=0;j<nchunk; j++)
      {
      pointl = (unsigned int *)data +j*npix;
      pointh = pointl+npix/2;
      for (i=npix/2; i--; pointh++, pointl++)
        *pointl ^= (*pointh^=(*pointl^=*pointh));
      }
    }


  return;
  }


/****** fft_init ************************************************************
PROTO	void fft_init(void)
PURPOSE	Initialize the FFT routines
INPUT	-.
OUTPUT	-.
NOTES	Global preferences are used for multhreading.
AUTHOR	E. Bertin (IAP)
VERSION	29/09/2004
 ***/
void	fft_init(void)
 {
  if (!firsttimeflag)
    {
#ifdef USE_THREADS
/*
    if (!fftwf_init_threads())
      error(EXIT_FAILURE, "*Error*: thread initialization failed in ", "FFTW");
    fftwf_plan_with_nthreads(prefs.nthreads);
*/
    QPTHREAD_MUTEX_INIT(&fftmutex, NULL);
#endif
    firsttimeflag = 1;
    }

  return;
  }


/****** fft_end ************************************************************
PROTO	void fft_init(void)
PURPOSE	Clear up stuff set by FFT routines
INPUT	-.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	29/09/2004
 ***/
void	fft_end(void)
 {

  if (firsttimeflag)
    {
    firsttimeflag = 0;
#ifdef USE_THREADS
/*
    fftwf_cleanup_threads();
*/
    QPTHREAD_MUTEX_DESTROY(&fftmutex);
#endif
    fftwf_cleanup();
    }

  return;
  }

