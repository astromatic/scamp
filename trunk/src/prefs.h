 /*
 				prefs.h

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Include for prefs.c.
*
*	Last modify:	18/02/2007
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifndef _CATOUT_H_
#include "catout.h"
#endif

#ifndef _CHECK_H_
#include "check.h"
#endif

#ifndef _CPLOT_H_
#include "cplot.h"
#endif

#ifndef _FGROUP_H_
#include "fgroup.h"
#endif

#ifndef _FIELD_H_
#include "field.h"
#endif

#ifndef _FITSWCS_H_
#include "fitswcs.h"
#endif

#ifndef _MATCH_H_
#include "match.h"
#endif

#ifndef _ASTREFCAT_H_
#include "astrefcat.h"
#endif

#ifndef _SAMPLES_H_
#include "samples.h"
#endif

#ifndef _PREFS_H_
#define _PREFS_H_

/*----------------------------- Internal constants --------------------------*/

#define		MAXCHARL	16384	/* max. nb of chars in a string list */
#define		MAXLIST		256	/* max. nb of list members */

/*--------------------------------- typedefs --------------------------------*/

/*------------------------------- preferences -------------------------------*/
typedef struct
  {
  char		**command_line;		/* Command line */
  int		ncommand_line;		/* nb of params */
  char		prefs_name[MAXCHAR];	/* prefs filename*/
  char		*(file_name[MAXFILE]);	/* Filename(s) of input images */
  int		nfile;			/* Number of input images */
  char		ahead_global[MAXCHAR];	/* Global input FITS header filename */
  char		ahead_suffix[MAXCHAR];	/* Suffix for input FITS headers */
  char		head_suffix[MAXCHAR];	/* Suffix for output FITS headers */
  enum {NORMAL, FOCAL_PLANE}	header_type;	/* Output header type */
  double	sn_thresh[2];		/* S/N thresholds */
  int		nsn_thresh;		/* nb of params */
  double	fwhm_thresh[2];		/* Minimum and maximum FWHM allowed */
  int		nfwhm_thresh;		/* nb of params */
  int		flags_mask;		/* Rejection mask on SEx FLAGS */
  int		wflags_mask;		/* Rej. mask on SEx FLAGS_WEIGHT */
  int		imaflags_mask;		/* Rejection mask on SEx IMAFLAG_ISO */

/* Reference catalogs */
  astrefenum	astrefcat;		/* Reference catalog */
  char		*(ref_server[MAX_SERVER]);/* IP addresses of ref cat servers */
  int           nref_server;		/* nb of params */
  int		ref_port[MAX_SERVER];	/* Port of ref. cat. servers */
  int           nref_port;		/* nb of params */
  char		*(astref_name[MAXNGROUP]);/* Astrometric ref. cat. filenames */
  int           nastref_name;		/* nb of params */
  char		astref_bandname[MAXCHAR];/* Astrometric ref. band name */
  char		*(astrefcent_key[NAXIS]);/* Astrom ref. cat. centroid keywords*/
  int		nastrefcent_key;	/* nb of params */
  char		*(astreferr_key[(NAXIS*(NAXIS+1))/2]);/* err.ellipse keywords*/
  int		nastreferr_key;		/* nb of params */
  char		astrefmag_key[72];	/* Path for ref.catalog output files */
  int		outrefcat_flag;		/* Save a FITS-LDAC copy of ref.cats?*/
  char		outref_path[MAXCHAR];	/* Path for ref.catalog output files */

/* Merged output catalogs */
  char		mergedcat_name[MAXCHAR];/* Output filename */
  cattypenum	mergedcat_type;		/* Output catalog type */
  int		mergedcatpipe_flag;	/* Pipe output catalogs? */

/* Pattern-matching */
  int		match_flag;		/* Compute pattern-matching? */
  double	group_radius;	/* Max angular dist. between grouped fields */
  double	position_maxerr[NAXIS];	/* Maximum uncertainty along axes */
  int           nposition_maxerr;	/* nb of params */
  double	radius_maxerr;		/* Maximum radius uncertainty */
  double	posangle_maxerr;	/* Maximum uncertainty in pos. angle */
  double	pixscale_maxerr;	/* Maximum uncertainty in pix. scale */
  int		flip_flag;		/* Allow matching with flipped axes? */
  mosaicenum	mosaic_type[MAXASTRINSTRU];	/* Type of mosaic device */
  int           nmosaic_type;		/* nb of params */
  int		fixfocalplane_nmin;	/* Min.# of sources for FIX_FOCALPLANE*/
  double	match_resol;		/* Matching resolution */
  int		nmatchmax;		/* Max.# of sources for MATCHing */

/* Cross-identification */
  double	crossid_radius;		/* Cross-identification radius */

/* Photometric solution */
  int		solvphotom_flag;		/* Compute photometric sol.? */
  char		photflux_key[72];		/* Flux measurement cat.param.*/
  char		photfluxerr_key[72];		/* Flux error cat.parameter */
  char		*(photinstru_key[72]);		/* Photom instrument keywords*/
  int		nphotinstru_key;		/* nb of params */
  char		**photinstrustr;		/* Photom instrument labels */
  int		nphotinstrustr;			/* nb of params */
  char		airmass_key[72];	/* FITS keyword for air mass */
  char		expotime_key[72];	/* FITS keyword for expo. time */
  char		extcoeff_key[72];	/* FITS keyword for extinction coeff */
  char		magzero_key[72];	/* FITS keyword for zero-point */
  char		photomflag_key[72];	/* FITS keyword for photom. exposures*/
  double	magzero_out[MAXPHOTINSTRU];	/* Max nb of phot instruments*/
  int		nmagzero_out;			/* nb of params */
  double	photclip_nsig;			/* Photom. clipping threshold*/
  double	magzero_interr[MAXPHOTINSTRU];	/* Internal field ZP error RMS*/
  int		nmagzero_interr;		/* nb of params */
  double	magzero_referr[MAXPHOTINSTRU];	/* Photom.field ZP error RMS */
  int		nmagzero_referr;		/* nb of params */

/* Astrometric solution */
  int		solvastrom_flag;		/* Compute astrometric sol.? */
  char          *(centroid_key[NAXIS]);		/* Names of centroid measur. */
  int           ncentroid_key;			/* nb of params */
  char          *(centroiderr_key[(NAXIS*(NAXIS+1))/2]);/* err ellipse names */
  int           ncentroiderr_key;		/* nb of params */
  char		*(astrinstru_key[72]);		/* Astrom instrument keywords*/
  int		nastrinstru_key;		/* nb of params */
  char		**astrinstrustr;		/* Astrom instrument labels */
  int		nastrinstrustr;			/* nb of params */
  stabilityenum	stability_type[MAXASTRINSTRU];	/* How stable the devices are*/
  int           nstability_type;		/* nb of params */
  char          *(context_name[MAXCONTEXT]);	/* Names of context-keys */
  int           ncontext_name;			/* nb of params */
  int           context_group[MAXCONTEXT];	/* Context group */
  int           ncontext_group;			/* nb of params */
  int           group_deg[MAXCONTEXT];		/* Degree for each group */
  int           ngroup_deg;			/* nb of params */
  double	astrclip_nsig;			/* Astrom. clipping threshold*/
  double	astref_weight;			/* Ref.cat. relative weight */
/* Differential Chromatic Refraction (DRC) */
  int		colourshift_flag;		/* Correct for colour shifts? */

/* Check-plots */
  cplotenum	cplot_device[MAXCHECK];		/* check-plot format */
  int		ncplot_device;			/* nb of params */
  cplotenum	cplot_type[MAXCHECK];		/* check-plot types */
  int		ncplot_type;			/* nb of params */
  char		*(cplot_name[MAXCHECK]);	/* check-plot names */
  int		ncplot_name;			/* nb of params */
  int		cplot_flag;			/* = 0 if no check-plot */
  char		cplot_colourkey[72];		/* FITS keyword for colour */
  int		cplot_res[2];			/* X,Y check-plot resolution */
  int		ncplot_res;			/* nb of params */
  int		cplot_antialiasflag;		/* Anti-aliasing on/off */
/* Check-images */
  checkenum	check_type[MAXCHECK];		/* check-image types */
  int		ncheck_type;			/* nb of params */
  char		*(check_name[MAXCHECK]);	/* check-image names */
  int		ncheck_name;			/* nb of params */
/* Multithreading */
  int		nthreads;			/* Number of active threads */
/* Misc */
  char		cdsclient_path[MAXCHAR];	/* Path for CDSclient execs */
  enum {QUIET, NORM, LOG, FULL} verbose_type;	/* display type */
  int		xml_flag;			/* Write XML file? */
  char		xml_name[MAXCHAR];		/* XML file name */
  char		xsl_name[MAXCHAR];		/* XSL file name (or URL) */
  char		sdate_start[12];		/* SCAMP start date */
  char		stime_start[12];		/* SCAMP start time */
  char		sdate_end[12];			/* SCAMP end date */
  char		stime_end[12];			/* SCAMP end time */
  double	time_diff;			/* Execution time */
  }	prefstruct;

prefstruct	prefs;

/*-------------------------------- protos -----------------------------------*/
extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		endprefs(void),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);


#endif

