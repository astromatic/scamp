 /*
 				field.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Handling of field structures.
*
*	Last modify:	26/10/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _SAMPLES_H_
#include "samples.h"
#endif

#ifndef _FIELD_H_
#define _FIELD_H_


/*----------------------------- Internal constants --------------------------*/

#define		MAXSET		1024		/* Max. number of sets/field */

/*--------------------------------- typedefs --------------------------------*/

typedef enum {MOSAIC_UNCHANGED, MOSAIC_SAMECRVAL, MOSAIC_SHAREPROJAXIS,
	      MOSAIC_FIXFOCALPLANE, MOSAIC_LOOSE} mosaicenum;
typedef enum {STABILITY_EXPOSURE, STABILITY_GROUP, STABILITY_INSTRUMENT,
		STABILITY_ALL}	stabilityenum;

typedef struct field
  {
/* ---- file parameters */
  char		path[MAXCHAR];		/* catalog path */
  char		filename[MAXCHAR];	/* catalog filename */
  char		*rfilename;		/* pointer to the reduced image name */
  char		hfilename[MAXCHAR];	/* header filename */
  int		headflag;		/* header found? */
  char		ident[MAXCHAR];		/* field identifier (read from FITS) */
/* ----  instrumental parameters */
  int		astromlabel;		/* Astrometry "instrument" label */
  mosaicenum	mosaic_type;		/* Mosaic type */
  stabilityenum	stability_type;		/* how stable the imager is */
/* ---- catalog parameters */
  struct set	**set;			/* Pointer to an array of sets */
  int		nset;			/* Number of sets */
  int		nsample;		/* Total number of samples */
/* ---- matching parameters */
  float		match_dscale;		/* Computed scale correction */
  float		match_dangle;		/* Computed angle correction */
  float		match_shear;		/* Computed shear amplitude */
  float		match_sangle;		/* Computed direction of the shear */
  float		match_dlng, match_dlat;	/* Positional corrections */
  float		match_asig, match_sig;	/* Contrasts for angle/scale and pos */
/* ---- astrometric parameters */
  double	epoch;			/* Epoch of observations */
  double	meanwcspos[NAXIS];	/* Mean pixel coordinate */
  double	meanwcsscale[NAXIS];	/* Mean pixel scale */
  double	maxradius;		/* Approximate radius of field (deg)*/
  int		naxis;			/* Number of axes */
  int		lng, lat;		/* Longitude and latitude indices */
  double	offset_ref[NAXIS];	/* Mean offset with reference cat. */
  double	offset_ref_hsn[NAXIS];	/* Mean offset with reference cat. */
  long double	chi2_int, chi2_int_hsn;	/* Internal astrometry chi2 */
  long long	nchi2_int,nchi2_int_hsn;/* Number of degrees of freedom */
  long double	sig_referr[NAXIS];	/* External astrometry std deviation */
  long double	sig_referr_hsn[NAXIS];	/* External astrometry std deviation */
  long double	sig_corr_ref;		/* X/Y correlation in ref. residuals */
  long double	sig_corr_ref_hsn;	/* X/Y correlation in ref. residuals*/
  long double	chi2_ref, chi2_ref_hsn;	/* External astrometry chi2 */
  long long	nchi2_ref,nchi2_ref_hsn;/* Number of degrees of freedom */
/* ---- photometric parameters */
  int		photomlabel;		/* Photometry "instrument" label */
  int		photomflag;		/* True if exposure is photometric */
  int		photomrank;		/* Rank as a photometric field */
  int		photomlink;		/* True if linked to a photom field*/
  double	dmagzero;		/* Magnitude zero-point correction */
  double	airmass;		/* Average field airmass */
  int		index;			/* CONTEXT index */  
  int		index2;			/* field-dependent index */  
/* ---- checkplot parameters */
  int		cplot_colour;		/* PLPLOT specific colour */
  struct fgroup	*fgroup;		/* Field group to which field belongs*/
  struct field	*prevfield;		/* Link to previous field */
  struct field	*nextfield;		/* Link to next field */
  }	fieldstruct;

/*------------------------------- functions ---------------------------------*/

extern fieldstruct	*inherit_field(char *filename, fieldstruct *reffield,
					int fflags),
			*init_field(fieldstruct **infield,
				fieldstruct **inwfield, int ninput,
				char *filename),
			*load_field(char *filename);

extern double		dhmedian(double *ra, int n);

extern void		end_field(fieldstruct *field),
			locate_field(fieldstruct *field),
			print_fieldinfo(fieldstruct *field),
			scale_field(fieldstruct *field, fieldstruct *reffield);
#ifdef USE_THREADS
void			pthread_end_fields(fieldstruct **fields, int nfield),
			*pthread_load_field(void *arg),
			pthread_load_fields(fieldstruct **fields, int nfield);
#endif
#endif
