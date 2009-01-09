 /*
 				fgroup.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Handling of field group structures.
*
*	Last modify:	31/10/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _FITSCAT_H_
#include "fits/fitscat.h"
#endif

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _FGROUP_H_
#define _FGROUP_H_


/*----------------------------- Internal constants --------------------------*/
#define		MAXASTRINSTRU	256	/* Max. number of astrom. instruments*/
#define		MAXPHOTINSTRU	256	/* Max. number of photom. instruments*/
#define		MAXNGROUP	65536	/* Max. number of groups (arbitrary) */

/*--------------------------------- typedefs --------------------------------*/
typedef struct fgroup
  {
  int		no;			/* Group ID (number) */
  struct field	**field;		/* Pointer to an array of fields */
  int		nfield;			/* Number of fields in the group */
  double	meanwcspos[NAXIS];	/* Mean field coordinates */
  double	meanwcsscale[NAXIS];	/* Mean pixel scale */
  double	projposmin[NAXIS];	/* Projected coordinate boundaries */
  double	projposmax[NAXIS];	/* Projected coordinate boundaries */
  int		naxis;			/* Number of dimensions */
  int		lng, lat;		/* Longitude and latitude indices */
  double	maxradius;		/* Maximum radius of group */
  struct wcs	*wcs;			/* Best WCS projection among fields */
  double	*intcolshiftscale[NAXIS];/* nfieldxnfield colour-shift scales */
  double	*intcolshiftzero[NAXIS];/* nfieldxnfield colour-shift zero-p.*/
  double	*colshiftscale[NAXIS];	/* nphotixnphoti colour-shift scales */
  double	*colshiftzero[NAXIS];	/* nphotixnphoti colour-shift zero-p.*/
  double	*refcolshiftscale[NAXIS];/* nfieldx1 colour-shift scales */
  double	*refcolshiftzero[NAXIS];/* nfieldx1 colour-shift zero-p.*/
  double	sig_interr[NAXIS];	/* Internal RMS relative error */
  double	chi2_int;		/* Chi2/d.o.f. for internal errors */
  int		nintmatch;		/* Number of internal matches */
  double	sig_referr[NAXIS];	/* RMS Error with respect to ref */
  double	chi2_ref;		/* Chi2/d.o.f. for reference errors */
  int		nrefmatch;		/* Number of matches with the ref cat*/
  double	sig_interr_hsn[NAXIS];	/* Internal RMS relative error */
  double	chi2_int_hsn;		/* Chi2/d.o.f. for internal errors */
  int		nintmatch_hsn;		/* Number of internal matches */
  double	sig_referr_hsn[NAXIS];	/* RMS Error with respect to ref */
  double	sig_corr_int;		/* proj. X/Y internal correl. coeff */ 
  double	sig_corr_int_hsn;	/* proj. X/Y internal correl. coeff */ 
  double	sig_corr_ref;		/* proj. X/Y ref correlation coeff */ 
  double	sig_corr_ref_hsn;	/* proj. X/Y ref correlation coeff */ 
  double	chi2_ref_hsn;		/* Chi2/d.o.f. for reference errors */
  int		nrefmatch_hsn;		/* Number of matches with the ref cat*/
  double	offset_ref[NAXIS];	/* Mean offset with reference cat */
  double	offset_ref_hsn[NAXIS];	/* Mean offset with reference cat */
  double	*sig_intmagerr;		/* Internal RMS relative mag errors */
  double	*chi2_intmag;		/* Chi2/d.o.f. for intern. mag errors*/
  int		*nintmagmatch;		/* Number of internal matches */
  double	*sig_intmagerr_hsn;	/* Internal RMS relative mag errors */
  double	*chi2_intmag_hsn;	/* Chi2/d.o.f. for intern. mag errors*/
  int		*nintmagmatch_hsn;	/* Number of internal matches */
  double	*sig_refmagerr;		/* RMS mag error with respect to refs*/
  double	*chi2_refmag;		/* Chi2/d.o.f. for ref. mag errors */
  int		*nrefmagmatch;		/* Number of matches with mag refs */
  double	*sig_refmagerr_hsn;	/* RMS mag error with respect to refs*/
  double	*chi2_refmag_hsn;	/* Chi2/d.o.f. for ref. mag errors */
  int		*nrefmagmatch_hsn;	/* Number of matches with mag refs */
  }	fgroupstruct;

/*------------------------------- functions ---------------------------------*/

extern fgroupstruct	**group_fields(fieldstruct **field, int nfield,
				int *nfgroup),
			*new_fgroup(void);

extern void		addfgroup_fgroup(fgroupstruct *fgroupin,
					fgroupstruct *fgroup),
			addfield_fgroup(fgroupstruct *fgroup,
				fieldstruct *field),
			end_fgroup(fgroupstruct *fgroup),
			locate_fgroup(fgroupstruct *fgroup),
			print_fgroupinfo(fgroupstruct **pfgroup, int nfgroup),
			print_instruinfo(void);

#endif
