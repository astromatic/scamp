/*
*				fft.h
*
* Include file for fft.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2013 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		05/04/2013
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _FFT_H_
#define _FFT_H_

#include FFTW_H

/*---------------------------- Internal constants ---------------------------*/
#define		LOPASS_SLOPE	1000	/* Slope of the sigmoid filter */
#define		HIPASS_SLOPE	1000	/* Slope of the sigmoid filter */

/*------------------------------- Other Macros ------------------------------*/
#define	QFFTWF_MALLOC(ptr, typ, nel) \
		{ \
		if (!(ptr = (typ *)fftwf_malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
                 }

#define	QFFTWF_CALLOC(ptr, typ, nel) \
		{ \
		if (!(ptr = (typ *)fftwf_malloc((size_t)(nel)*sizeof(typ)))) \
		   { \
		   sprintf(gstr, #ptr " (" #nel "=%lu elements) " \
			"at line %d in module " __FILE__ " !", \
			(size_t)(nel)*sizeof(typ), __LINE__); \
		   error(EXIT_FAILURE, "Could not allocate memory for ", gstr);\
                   }; \
		memset(ptr, 0, (size_t)(nel)*sizeof(typ)); \
		}

#define	QFFTWF_FREE(ptr) \
		{fftwf_free(ptr); ptr=NULL;}

/*--------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern void	fastcorr(float *data1, float *data2, int naxis, int *size,
			double *lambda_lopass, double *lambda_hipass),
		fft_end(void),
		fft_init(int nthreads),
		shiftcube(float *data, int naxis, int *size);

#endif // _FFT_H_
