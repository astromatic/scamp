/*
*				cplot.h
*
* Include file for cplot.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2012 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		09/07/2014
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef _CPLOT_H_
#define _CPLOT_H_

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include "fgroup.h"
#include "fitswcs.h"

/*------------------------------- constants ---------------------------------*/

#define	CPLOT_DEFRESX		800	/* Default X resol. for PNG and JPG */
#define	CPLOT_DEFRESY		600	/* Default X resol. for PNG and JPG */
#define	CPLOT_AAFAC		3	/* Anti-aliasing factor */
#define	CPLOT_NPOINTDEF		1024	/* Default number of points to plot */
#define	CPLOT_FGRIDLINES 	9	/* Number of grid lines per axis */
#define	CPLOT_NDISTGRID		32	/* # of distort steps in each CCD dim*/
#define	CPLOT_ASTNSUBPLOTS 	3	/* Number of subplot/dim/detector*/
#define	CPLOT_NTYPES	 	128 /* Number of CPLOT types (typedef below)*/
#define	CPLOT_NSHADES		32	/* Number of shading levels */
#define	CPLOT_NADERRHISTBIN	255	/* Number of int. 1D histo bins */
#define	CPLOT_NREFERRHISTBIN	127	/* Number of ref. histo bins */
#define	CPLOT_NPIXERRHISTBIN	127	/* Number of int. 1D histo bins */
#define	CPLOT_NSUBPIXERRHISTBIN	127	/* Number of int. 1D histo bins */
#define	CPLOT_PHOTERRNX		255	/* x resolution of phot 2D histogram */
#define	CPLOT_PHOTERRNY		31	/* y resolution of phot 2D histogram */
#define	CPLOT_PHOTERRNX_HSN	1023	/* x resolution of phot 2D histogram */
#define	CPLOT_PHOTERRNY_HSN	255	/* y resolution of phot 2D histogram */
#define	CPLOT_ADERR1DNX		255	/* x resolution of astrom 2D histo */
#define	CPLOT_ADERR1DNY		31	/* y resolution of astrom 2D histo */
#define	CPLOT_ADERR1DNX_HSN	1023	/* x resolution of astrom 2D histo */
#define	CPLOT_ADERR1DNY_HSN	255	/* y resolution of astrom 2D histo */
#define	CPLOT_ADERR2DN		63	/* resolution of astrom 2D histogram */
#define	CPLOT_ADERR2DN_HSN	1023	/* resolution of astrom 2D histogram */
#define	CPLOT_REFERR2DN		63	/* resolution of astrom 2D histogram */
#define	CPLOT_REFERR2DN_HSN	511	/* resolution of astrom 2D histogram */
#define	CPLOT_PIXERR1DNX	127	/* x resolution of astrom 2D histo*/
#define	CPLOT_PIXERR1DNY	31	/* y resolution of astrom 2D histo */
#define	CPLOT_PIXERR1DNX_HSN	1023	/* x resolution of astrom 2D histo */
#define	CPLOT_PIXERR1DNY_HSN	255	/* y resolution of astrom 2D histo */
#define	CPLOT_SUBPIXERR1DNX	63	/* x resolution of astrom 2D histo*/
#define	CPLOT_SUBPIXERR1DNY	31	/* y resolution of astrom 2D histo */
#define	CPLOT_SUBPIXERR1DNX_HSN	1023	/* x resolution of astrom 2D histo*/
#define	CPLOT_SUBPIXERR1DNY_HSN	255	/* y resolution of astrom 2D histo*/
#define	CPLOT_ASTRCOLSHIFT1DNX	31	/* x resolution of colshift 1D histo*/
#define	CPLOT_ASTRCOLSHIFT1DNY	31	/* x resolution of colshift 1D histo*/
#define	CPLOT_ASTRCOLSHIFT1DNX_HSN 511	/* x resol of colshift 1D histo*/
#define	CPLOT_ASTRCOLSHIFT1DNY_HSN 511	/* x resol of colshift 1D histo*/
#define	CPLOT_REFPROPN		63	/* resolution of prop.motion histo */
#define	CPLOT_REFPROPN_HSN	511	/* resolution of prop.motion histo */
#define	CPLOT_ASTREFPROPMINSN	20.0	/* minimum S/N on ref. prop motion */
#define	CPLOT_NPIXERRGRID	32	/* resolution of pix.error 2D histo */

/*---------------------------- return messages ------------------------------*/
/*-------------------------------- macros -----------------------------------*/
// Work-around to emulate the plwid() function replaced in later versions of
// the PLPlot library.
#ifdef HAVE_PLPLOT
 #ifndef __PLPLOT_H__
  #include	PLPLOT_H
 #endif

 #ifndef __PLPLOTP_H__
  #include	PLPLOTP_H
  #endif

 #ifdef plwidth
  #define	CPLOT_PLWID(wid)	plwidth((PLFLT)(wid))
 #else
  #define	CPLOT_PLWID		plwid
 #endif
#endif

/*--------------------------------- typedefs --------------------------------*/
typedef enum {CPLOT_NONE, CPLOT_ALLSKY, CPLOT_FGROUPS, CPLOT_PHOTOM,
	CPLOT_ADERROR1D, CPLOT_REFERROR1D, CPLOT_PIXERROR1D,
	CPLOT_SUBPIXERROR1D,
	CPLOT_DISTORT, CPLOT_SHEAR, CPLOT_PHOTERROR, CPLOT_PHOTERRORVSMAG,
	CPLOT_PHOTZP, CPLOT_CHI2, CPLOT_ADERROR2D, CPLOT_REFERROR2D,
	CPLOT_PHOTZP3D, CPLOT_ASTRCOLSHIFT1D, CPLOT_REFPROP,
	CPLOT_ADSYSMAP2D, CPLOT_REFSYSMAP2D, CPLOT_ASTREPOCH3D,
	CPLOT_ADPROP2D, CPLOT_XPIXERROR2D, CPLOT_YPIXERROR2D}
		cplotenum;

typedef enum {CPLOT_NULL, CPLOT_XWIN, CPLOT_TK, CPLOT_PLMETA, CPLOT_PS,
	CPLOT_PSC, CPLOT_XFIG, CPLOT_PNG, CPLOT_JPEG, CPLOT_PSTEX, CPLOT_AQT,
	CPLOT_PDF, CPLOT_SVG} cplotdevenum;

typedef struct {wcsstruct	*wcsin,*wcsout;
		int		ngridx,ngridy;}
		distortstruct;

typedef struct {cplotdevenum device; char *devname; char *extension;}
		devicestruct;

/*---------------------------------- svgp -----------------------------------*/
/*------------------------------- functions ---------------------------------*/

extern int		cplot_aderrhisto1d(fgroupstruct *fgroup,
					double hsn_thresh),
			cplot_aderrhisto2d(fgroupstruct *fgroup,
					double hsn_thresh),
			cplot_allsky(fgroupstruct **fgroups, int ngroup),
			cplot_astintsysmap(fgroupstruct **fgroups, int ngroup,
					int instru, double hsn_thresh),
			cplot_adprophisto2d(fgroupstruct *fgroup,
					double hsn_thresh),
			cplot_astrefsysmap(fgroupstruct **fgroups, int ngroup,
					int instru, double hsn_thresh),
			cplot_astrcolshift1d(fgroupstruct *fgroup,
					double hsn_thresh),
			cplot_astrepoch3d(fgroupstruct *fgroup),
			cplot_astrefprop(fgroupstruct *fgroup,
					fieldstruct *reffield,
					double hsn_thresh),
			cplot_check(cplotenum cplottype),
			cplot_chi2(fgroupstruct *fgroup),
			cplot_distort(fieldstruct *field),
			cplot_drawbounds(wcsstruct *wcsin, wcsstruct *wcsout),
			cplot_drawcoordgrid(wcsstruct *wcs, double xmin,
					double xmax, double ymin, double ymax),
			cplot_drawloccoordgrid(wcsstruct *wcs, double xmin,
					double xmax, double ymin, double ymax),
			cplot_drawfgrid(wcsstruct *wcsin, wcsstruct *wcsout),
			cplot_end(cplotenum cplottype),
			cplot_fgroup(fgroupstruct *fgroup,
				fieldstruct *reffield),
			cplot_init(int nx, int ny, cplotenum cplottype),
			cplot_photerrhisto(fgroupstruct *fgroup,
				fieldstruct *reffield, double hsn_thresh),
			cplot_photerrhistomag(fgroupstruct *fgroup,
				fieldstruct *reffield, double hsn_thresh),
			cplot_photom(fgroupstruct **fgroups, int ngroup,
				fieldstruct **reffields),
			cplot_photzp(fgroupstruct *fgroup),
			cplot_photzp3d(fgroupstruct *fgroup),
			cplot_pixerrhisto1d(fgroupstruct **fgroups,
				int ngroup, int instru, double hsn_thresh),
			cplot_shear(fgroupstruct **fgroups, int ngroup,
					int instru),
			cplot_subpixerrhisto1d(fgroupstruct **fgroups,
				int ngroup, int instru, double hsn_thresh),
			cplot_referrhisto1d(fgroupstruct *fgroup,
					fieldstruct *reffield,
					double hsn_thresh),
			cplot_referrhisto2d(fgroupstruct *fgroup,
					fieldstruct *reffield,
					double hsn_thresh),
			cplot_xpixerrhisto2d(fgroupstruct **fgroups, int ngroup,
				int instru),
			cplot_ypixerrhisto2d(fgroupstruct **fgroups, int ngroup,
				int instru);


char			*cplot_degtosexal(char *str, double alpha,double step),
			*cplot_degtosexde(char *str, double delta,double step);
#endif // _CPLOT_H_
