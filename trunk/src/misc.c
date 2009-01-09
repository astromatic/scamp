/*
                                  misc.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*       Part of:        SCAMP
*
*       Author:         E.BERTIN (IAP)
*
*       Contents:       Miscellaneous routines
*
*       Last modify:    05/08/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifdef HAVE_MATHIMF_H
#include <mathimf.h>
#else
#include <math.h>
#endif
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "types.h"
#include "globals.h"
#include "misc.h"

/******* fast_median **********************************************************
PROTO	double fast_median(double *arr, int n)
PURPOSE	Fast median from an input array, optimized version based on the
	select() routine (Numerical Recipes, 2nd ed. Section 8.5 and
	http://www.eso.org/~ndevilla/median/). If n is even, then the result
	is the average of the 2 "central" values.
INPUT	Input pixel array ptr,
	number of input pixels,
OUTPUT	Value of the median.
NOTES	n must be >0. Warning: changes the order of data (but does not sort
	them)!
AUTHOR	E. Bertin (IAP)
VERSION	04/08/2003
 ***/
#define MEDIAN_SWAP(a,b) \
	{val=a; a=b; b=val;}

double fast_median(double *arr, int n) 
  {
   double	*alow, *ahigh, *amedian, *amiddle, *all, *ahh,
		val, valmax, valmax2;
   int		i, nless;

  if (n==1)
    return *arr;
  else if (n==2)
      return 0.5*(arr[0]+arr[1]);
  else if (n==3)
    {
    if (arr[0]>arr[1])
      MEDIAN_SWAP(*arr, *(arr+1));
    if (arr[1]>arr[2])
      MEDIAN_SWAP(*(arr+1), *(arr+2));
    if (arr[0]>arr[1])
      MEDIAN_SWAP(*arr, *(arr+1));
    return arr[1];
    }
  alow = arr;
  ahigh = arr + n - 1;
  amedian = arr + n/2;
  while (ahigh > (all=alow + 1))
    {
/*-- Find median of low, middle and high items; swap into position low */
    amiddle = alow + (ahigh-alow)/2;
    if (*amiddle > *ahigh)
      MEDIAN_SWAP(*amiddle, *ahigh);
    if (*alow > *ahigh)
      MEDIAN_SWAP(*alow, *ahigh);
    if (*amiddle > *alow)
      MEDIAN_SWAP(*amiddle, *alow);

/*-- Swap low item (now in position middle) into position (low+1) */
    MEDIAN_SWAP(*amiddle, *all);

/*-- Nibble from each end towards middle, swapping items when stuck */
    ahh = ahigh;
    for (;;)
      {
      while (*alow > *(++all));
      while (*(--ahh) > *alow);

      if (ahh < all)
        break;

      MEDIAN_SWAP(*all, *ahh);
      }

/*-- Swap middle item (in position low) back into correct position */
    MEDIAN_SWAP(*alow, *ahh) ;

/*-- Re-set active partition */
    if (ahh <= amedian)
      alow = all;
    if (ahh >= amedian)
      ahigh = ahh - 1;
    }

/* One or two elements left */
  if (ahigh == all && *alow > *ahigh)
    MEDIAN_SWAP(*alow, *ahigh);

  if (n&1)
/*-- Odd case */
    return *amedian;
  else
    {
/*-- Even case */
    valmax2 = *amedian;
    valmax = -BIG;
    alow = arr;
    nless = 0;
    for (i=n/2;i--;)
      if ((val=*(alow++))<valmax2)
        {
        nless++;
        if (val > valmax)
          valmax = val;
        }
    return nless<n/2? *amedian : (*amedian+valmax)/2.0;
    }

  }

#undef MEDIAN_SWAP


