/*
*				prefs.h
*
* Include file for prefs.c
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		12/08/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _PREFS_H_
#define _PREFS_H_

#include "define.h"
#include "catout.h"
#include "check.h"
#include "cplot.h"
#include "fgroup.h"
#include "field.h"
#include "fitswcs.h"
#include "match.h"
#include "astrefcat.h"
#include "samples.h"

/*----------------------------- Internal constants --------------------------*/

#define		MAXCHARL	16384	/* max. nb of chars in a string list */
#define		MAXLIST		(MAXFILE)	/* max. nb of list members */
#define		MAXLISTSIZE	(100*MAXLIST)	/* max size of list */

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
  char		*(ahead_name[MAXFILE]);	/* Filename(s) of input FITS headers */
  int		nahead_name;		/* Number of input FITS headers */
  char		ahead_suffix[MAXCHAR];	/* Suffix for input FITS headers */
  char		*(head_name[MAXFILE]);	/* Filename(s) of output FITS headers */
  int		nhead_name;		/* Number of output FITS headers */
  char		head_suffix[MAXCHAR];	/* Suffix for output FITS headers */
  enum {NORMAL, FOCAL_PLANE}	header_type;	/* Output header type */
  double	sn_thresh[2];		/* S/N thresholds */
  int		nsn_thresh;		/* nb of params */
  double	fwhm_thresh[2];		/* Minimum and maximum FWHM allowed */
  int		nfwhm_thresh;		/* nb of params */
  double	maxellip;		/* Maximum ellipticity allowed */
  unsigned short sexflags_mask;		/* Rejection mask on SEx FLAGS */
  unsigned short wflags_mask;		/* Rej. mask on SEx FLAGS_WEIGHT */
  unsigned int	imaflags_mask;		/* Rejection mask on SEx IMAFLAG_ISO */

/* Reference catalogs */
  astrefenum	astrefcat;		/* Reference catalog */
  char		*(ref_server[MAX_SERVER]);/* IP addresses of ref cat servers */
  int           nref_server;		/* nb of params */
  int		ref_ntries[MAX_SERVER];	/* nb of tries per server */
  int           nref_ntries;		/* nb of params */
  double	ref_timeout[MAX_SERVER];/* Timeout of ref. cat. servers in s */
  int           nref_timeout;		/* nb of params */
  char		*(astref_name[MAXNGROUP]);/* Astrometric ref. cat. filenames */
  int           nastref_name;		/* nb of params */
  char		astref_bandname[MAXCHAR];/* Astrometric ref. band name */
  double	astref_maglim[2];	/* Selection magnitude range */
  int		nastref_maglim;		/* nb of params */
  char		*(astrefcent_key[NAXIS]);/* Astrom ref. cat. centroid keywords*/
  int		nastrefcent_key;	/* nb of params */
  char		*(astreferr_key[(NAXIS*(NAXIS+1))/2]);/* err.ellipse keywords*/
  int		nastreferr_key;		/* nb of params */
  char		*(astrefprop_key[NAXIS]);/* Astrom ref. cat. PM keywords*/
  int		nastrefprop_key;	/* nb of params */
  char		*(astrefproperr_key[NAXIS]);/* PM err keywords*/
  int		nastrefproperr_key;		/* nb of params */
  char		astrefmag_key[72];	/* Astrom ref. cat. mag. keyword */
  char		astrefmagerr_key[72];	/* Astrom ref. cat. mag. error keyword*/
  char		astrefobsdate_key[72];	/* Astrom ref. cat. obs. date keyword */
  int		outrefcat_flag;		/* Save a FITS-LDAC copy of ref.cats?*/
  char		outref_path[MAXCHAR];	/* Path for ref.catalog output files */

/* Merged output catalogs */
  char		mergedcat_name[MAXCHAR];/* Output filename */
  cattypenum	mergedcat_type;		/* Output catalog type */
  int		mergedcatpipe_flag;	/* Pipe output catalogs? */

/* Full output catalogs */
  char		fullcat_name[MAXCHAR];/* Output filename */
  cattypenum	fullcat_type;		/* Output catalog type */
  int		fullcatpipe_flag;	/* Pipe output catalogs? */
  int		spread_flag;		/* SPREAD_MODEL in input catalog(s)? */

/* Differential geometry maps */
  int		dgeomap_flag;		/* Compute and save Diff. Geom. maps? */
  char		dgeomap_name[MAXCHAR];	/* Output filename */
  int		dgeomap_nnearest;	/* Number of nearest neighbours */
  int		dgeomap_step;		/* Map sampling step */

/* Pattern-matching */
  int		match_flag;		/* Compute pattern-matching? */
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
  char		photflux_key[72];		/* Name of phot. flux key */
  char		photflux_rkey[72];		/* Reduced phot. flux key */
  int		photflux_num;			/* Phot. aperture # */
  char		photfluxerr_key[72];		/* Name of phot. flux err. key*/
  char		photfluxerr_rkey[72];		/* Reduced phot. flux err. key*/
  int		photfluxerr_num;		/* Phot.flux err. aperture # */
  char		*(photinstru_key[35]);		/* Photom instrument keywords */
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
  double	photaccuracy;			/* Photom. uncertainty floor */
  double	magzero_interr[MAXPHOTINSTRU];	/* Internal field ZP error RMS*/
  int		nmagzero_interr;		/* nb of params */
  double	magzero_referr[MAXPHOTINSTRU];	/* Photom.field ZP error RMS */
  int		nmagzero_referr;		/* nb of params */
  unsigned short phot_sexflagsmask;		/* Photom. mask on SEx FLAGS */
  unsigned int	phot_imaflagsmask;		/* Photom. mask on image flags*/

/* Astrometric solution */
  projenum	projection_type[MAXFILE];	/* Celestial projection type */
  int           nprojection_type;		/* nb of params */
  int		solvastrom_flag;		/* Compute astrometric sol.? */
  char          *(centroid_key[NAXIS]);		/* Names of centroid measur. */
  int           ncentroid_key;			/* nb of params */
  char          *(centroiderr_key[(NAXIS*(NAXIS+1))/2]);/* err ellipse names */
  int           ncentroiderr_key;		/* nb of params */
  char		*(astrinstru_key[35]);		/* Astrom instrument keywords */
  int		nastrinstru_key;		/* nb of params */
  char		**astrinstrustr;		/* Astrom instrument labels */
  int		nastrinstrustr;			/* nb of params */
  int		*nastrinstruext;		/* nb of extensions per instru*/
  stabilityenum	stability_type[MAXASTRINSTRU];	/* How stable the devices are*/
  int           nstability_type;		/* nb of params */
  char          *(context_name[MAXCONTEXT]);	/* Names of context-keys */
  int           ncontext_name;			/* nb of params */
  int           context_group[MAXCONTEXT];	/* Context group */
  int           ncontext_group;			/* nb of params */
  int           group_deg[MAXCONTEXT];		/* Degree for each group */
  int           ngroup_deg;			/* nb of params */
  int		focal_deg;			/* Focal coords polynom degree*/ 
  double	astrclip_nsig;			/* Astrom. clipping threshold*/
  double	astref_weight;			/* Ref.cat. relative weight */
  accuracyenum	astraccuracy_type;		/* Astrom. uncer. input type */
  char		astraccuracy_key[72];	/* Fits keyword for astrom. uncer. */
  double	astraccuracy;			/* Astrom. uncertainty floor */
  unsigned short astr_sexflagsmask;		/* Astrom. mask on SEx FLAGS */
  unsigned int	astr_imaflagsmask;		/* Astrom. mask on image flags*/
/* Parallaxes */
  int		parallax_flag;			/* Compute parallaxes? */
/* Proper motions */
  int		propmotion_flag;		/* Compute proper motions? */
  int		propmotioncorr_flag;		/* Correct proper motions? */
  int		astrefinprop_flag;		/* Use ref.catalog in pm? */
/* Differential Chromatic Refraction (DRC) */
  int		colourshiftcorr_flag;		/* Correct colour shifts?*/

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
  int		stats_maxmeshsize;		/* Max mesh size for stats */
/* Check-images */
  checkenum	check_type[MAXCHECK];		/* check-image types */
  int		ncheck_type;			/* nb of params */
  char		*(check_name[MAXCHECK]);	/* check-image names */
  int		ncheck_name;			/* nb of params */
/* Multithreading */
  int		nthreads;			/* Number of active threads */
/* Misc */
  enum {QUIET, NORM, LOG, FULL} verbose_type;	/* display type */
  int		xml_flag;			/* Write XML file? */
  char		xml_name[MAXCHAR];		/* XML file name */
  char		xsl_name[MAXCHAR];		/* XSL file name (or URL) */
  char		sdate_start[12];		/* SCAMP start date */
  char		stime_start[12];		/* SCAMP start time */
  char		sdate_end[12];			/* SCAMP end date */
  char		stime_end[12];			/* SCAMP end time */
  double	time_diff;			/* Execution time */
  int		ndets;				/* Final number of sources */
  }	prefstruct;

extern prefstruct	prefs;

/*-------------------------------- protos -----------------------------------*/
extern char	*list_to_str(char *listname);

extern int	cistrcmp(char *cs, char *ct, int mode);

extern void	dumpprefs(int state),
		endprefs(void),
		readprefs(char *filename,char **argkey,char **argval,int narg),
		useprefs(void);


#endif // _PREFS_H_
