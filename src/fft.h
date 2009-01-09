 /*
 				fft.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for fft.c.
*
*	Last modify:	14/07/2003
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

/*---------------------------- Internal constants ---------------------------*/
#define		LOPASS_SLOPE	1000	/* Slope of the sigmoid filter */
#define		HIPASS_SLOPE	1000	/* Slope of the sigmoid filter */

/*------------------------------- Other Macros ------------------------------*/
#define	QFFTWMALLOC(ptr, typ, nel) \
		{if (!(ptr = (typ *)fftwf_malloc((size_t)(nel)*sizeof(typ)))) \
		  error(EXIT_FAILURE, "Not enough memory for ", \
			#ptr " (" #nel " elements) !");;}
#define	QFFTWFREE(ptr)	fftwf_free(ptr)

/*--------------------------- structure definitions -------------------------*/

/*---------------------------------- protos --------------------------------*/
extern void	fastcorr(float *data1, float *data2, int naxis, int *size,
			double *lambda_lopass, double *lambda_hipass),
		fft_end(void),
		fft_init(void),
		shiftcube(float *data, int naxis, int *size);
