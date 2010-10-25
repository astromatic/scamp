/*
*				astrcplot.c
*
* Produce astrometric check plots.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include	<math.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>

#include	PLPLOT_H
#include	PLPLOTP_H

#include	"define.h"
#include	"globals.h"
#include	"cplot.h"
#include	"fgroup.h"
#include	"field.h"
#include	"fitswcs.h"
#include	"prefs.h"
#include	"samples.h"

extern devicestruct	cplot_device[];
extern int		plotaaflag;

static void	distort_map(PLFLT x,PLFLT y, PLFLT *tx,PLFLT *ty,
			    void *pltr_data); 


/****** cplot_allsky *******************************************************
PROTO	int cplot_allsky(fgroupstruct **fgroups, int ngroup)
PURPOSE	Plot the position of field groups on the celestial sphere.
INPUT	Pointer to the array of field groups,
	number of groups.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	19/07/2004
 ***/
int	cplot_allsky(fgroupstruct **fgroups, int ngroup)
  {
   wcsstruct	*wcsin,*wcsout;
   fieldstruct	**field;
   setstruct	**set;
   PLFLT	*x, *y;
   PLINT	mark[2],space[2],
		lwid;
   char		str[32];
   double	wcspos[NAXIS],wcspos2[NAXIS], rawpos[NAXIS],rawpos2[NAXIS],
		rawposc[NAXIS],
		rad, scale;
   int		nfield,nset,lng,lat,
		i, g,f,s;

  if (cplot_init(1,1, CPLOT_ALLSKY) == RETURN_ERROR)
    {
    cplot_end(CPLOT_ALLSKY);
    return RETURN_OK;
    }
  pladv(0);
  plvpas(0.0,1.0,0.0,1.0,0.5);
  plwind(-180.0,180.0, 90.0,-90.0);

  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;

/* Prepare an Aitoff projection */
  wcsout = create_wcs(NULL,NULL,NULL,NULL,NULL,2);

  QMALLOC(x, PLFLT, CPLOT_NPOINTDEF);
  QMALLOC(y, PLFLT, CPLOT_NPOINTDEF);
/* Draw meridians and parallels */
  plfont(2);
  plschr(0.0, 0.3);
  plwid(0);
  plcol(7);
  mark[0] = 500;
  space[0] = 500;
  for (wcspos[0] = -180.0; wcspos[0]<=180.0; wcspos[0]+=14.999)
    {
    i=0;
    for (wcspos[1]=-90.0; wcspos[1]<90.01 && i<CPLOT_NPOINTDEF;
	wcspos[1]+=2.5, i++)
      {
      wcs_to_raw(wcsout, wcspos, rawpos);
      x[i] = rawpos[0];
      y[i] = rawpos[1];
      if (wcspos[1]>-0.1 && wcspos[1]<0.1 && wcspos[0] >-179.0)
        {
        sprintf(str, "%d#uh", ((int)(wcspos[0]/15.0+0.1+24.0))%24);
        plptex(x[i], y[i]-2.0, 1.0, 0.0, -0.2, str);
        }
      }
    if (wcspos[0]>-179.0 && wcspos[0]<179.0)
      plstyl(1, mark,space);
    else
      pllsty(1);
    plline(i, x, y);
    }

  plstyl(1, mark,space);
  for (wcspos[1] =-90.0; wcspos[1]<=90.01; wcspos[1]+=10.0)
    {
    i=0;
    for (wcspos[0]=-180.0; wcspos[0]<180.0 && i<CPLOT_NPOINTDEF;
	wcspos[0]+=3.999999, i++)
      {
      wcs_to_raw(wcsout, wcspos, rawpos);
      x[i] = rawpos[0];
      y[i] = rawpos[1];
      if (wcspos[0] ==-180.0)
        {
        sprintf(str, "%d#uo", (int)(wcspos[1]+(wcspos[1]>=0.0?0.1:-0.1)));
        plptex(x[i]*1.03, y[i]*1.04, 1.0, 0.0, 0.5, str);
        }
      }
    plline(i, x, y);
    }

/* Plot fields */
  pllsty(1);
  plwid(lwid);
  for (g=0; g<ngroup; g++)
    {
    field = fgroups[g]->field;
    nfield = fgroups[g]->nfield;
    plcol(15);
    for (f=0; f<nfield; f++)
      {
      set = field[f]->set;
      nset = field[f]->nset;
/*---- Draw image boundaries */
      for (s=0; s<nset; s++)
        cplot_drawbounds(set[s]->wcs, wcsout);
      }

/*-- Plot a circle around the position */
    wcsin = fgroups[g]->wcs;
    lng = wcsin->lng;
    lat = wcsin->lat;
    wcs_to_raw(wcsin, fgroups[g]->meanwcspos, rawposc);
    rad = fgroups[g]->maxradius+2.0;	/* Add a 2deg margin */
    if (rad>45.0)			/* Put some limitation! */
      rad = 45.0;
    scale = sqrt(wcs_scale(wcsin, rawposc));
    if (scale>0.0)
      rad /= scale;
    for (i=0; i<wcsin->naxis; i++)
      rawpos[i] = rawposc[i];
    for (i=0; i<37; i++)
      {
      rawpos[lng] = rawposc[lng] + rad*cos(i*10.0*DEG);
      rawpos[lat] = rawposc[lat] + rad*sin(i*10.0*DEG);
      raw_to_wcs(wcsin, rawpos, wcspos);
      wcspos2[0] = wcspos[lng];
      wcspos2[1] = wcspos[lat];
      wcs_to_raw(wcsout, wcspos2, rawpos2);
      x[i] = rawpos2[0];
      y[i] = rawpos2[1];
      }
    plcol(3);
    plline(37, x,y);
    }

  plend();

  free(x);
  free(y);

  cplot_allsky(fgroups, ngroup);

  return RETURN_OK;
  }


/****** cplot_drawbounds *****************************************************
PROTO	int cplot_drawbounds(wcsstruct *wcsin, wcsstruct *wcsout)
PURPOSE	Draw the projected image boundaries in a given projection.
INPUT	Pointer to the image WCS structure,
	pointer to the plot WCS structure.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	18/10/2002
 ***/
int	cplot_drawbounds(wcsstruct *wcsin, wcsstruct *wcsout)
  {
   PLFLT	x[5],y[5];
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS];
   int		i;

/* Initialize the input coordinates to an "average" value */
  for (i=0; i<wcsin->naxis; i++)
    rawpos2[i] = wcsin->naxisn[i]/2.0;

/* 1st corner */
  rawpos2[wcsin->lng] = 0.0;
  rawpos2[wcsin->lat] = 0.0;
  raw_to_wcs(wcsin, rawpos2, wcspos2);
  wcspos[wcsout->lng] = wcspos2[wcsin->lng];
  wcspos[wcsout->lat] = wcspos2[wcsin->lat];
  wcs_to_raw(wcsout, wcspos, rawpos);
  x[4] = x[0] = rawpos[wcsout->lng];
  y[4] = y[0] = rawpos[wcsout->lat];
/* 2nd corner */
  rawpos2[wcsin->lng] = wcsin->naxisn[wcsin->lng]-1.0;
  raw_to_wcs(wcsin, rawpos2, wcspos2);
  wcspos[wcsout->lng] = wcspos2[wcsin->lng];
  wcspos[wcsout->lat] = wcspos2[wcsin->lat];
  wcs_to_raw(wcsout, wcspos, rawpos);
  x[1] = rawpos[wcsout->lng];
  y[1] = rawpos[wcsout->lat];
/* 3rd corner */
  rawpos2[wcsin->lat] = wcsin->naxisn[wcsin->lat]-1.0;
  raw_to_wcs(wcsin, rawpos2, wcspos2);
  wcspos[wcsout->lng] = wcspos2[wcsin->lng];
  wcspos[wcsout->lat] = wcspos2[wcsin->lat];
  wcs_to_raw(wcsout, wcspos, rawpos);
  x[2] = rawpos[wcsout->lng];
  y[2] = rawpos[wcsout->lat];
/* Last corner */
  rawpos2[wcsin->lng] = 0.0;
  raw_to_wcs(wcsin, rawpos2, wcspos2);
  wcspos[wcsout->lng] = wcspos2[wcsin->lng];
  wcspos[wcsout->lat] = wcspos2[wcsin->lat];
  wcs_to_raw(wcsout, wcspos, rawpos);
  x[3] = rawpos[wcsout->lng];
  y[3] = rawpos[wcsout->lat];
/* Draw */
  plline(5, x,y);

  return RETURN_OK;
  }


/****** cplot_drawfgrids *****************************************************
PROTO	int cplot_drawfgrids(wcsstruct *wcsin, wcsstruct *wcsout)
PURPOSE	Draw a projected image grid in a given projection.
INPUT	Pointer to the image WCS structure,
	pointer to the plot WCS structure.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	15/12/2003
 ***/
int	cplot_drawfgrids(wcsstruct *wcsin, wcsstruct *wcsout)
  {
   PLFLT	x[CPLOT_FGRIDLINES],y[CPLOT_FGRIDLINES];
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		xpos, ypos, xstep,ystep;
   int		i,j, lngin,latin, lngout,latout, fgridlines;

  lngin = wcsin->lng;
  latin = wcsin->lat;
  lngout = wcsout->lng;
  latout = wcsout->lat;
  fgridlines = CPLOT_FGRIDLINES - 1;
  xstep = (wcsin->naxisn[lngin]-1)/(double)fgridlines;
  ystep = (wcsin->naxisn[latin]-1)/(double)fgridlines;

/* Initialize the input coordinates to an "average" value */
  for (i=0; i<wcsin->naxis; i++)
    rawpos2[i] = wcsin->naxisn[i]/2.0;

/* Horizontal lines */
  for (j=1; j<fgridlines; j++)
    {
    ypos = j*ystep+1.0;
    for (i=0; i<CPLOT_FGRIDLINES; i++)
      {
      rawpos2[lngin] = i*xstep+1.0;
      rawpos2[latin] = ypos;
      raw_to_wcs(wcsin, rawpos2, wcspos2);
      wcspos[lngout] = wcspos2[lngin];
      wcspos[latout] = wcspos2[latin];
      wcs_to_raw(wcsout, wcspos, rawpos);
      x[i] = rawpos[lngout];
      y[i] = rawpos[latout];
      }
/*-- Draw */
    plline(CPLOT_FGRIDLINES, x,y);
    }

/* Vertical lines */
  for (i=1; i<fgridlines; i++)
    {
    xpos = i*xstep+1.0;
    for (j=0; j<CPLOT_FGRIDLINES; j++)
      {
      rawpos2[lngin] = xpos;
      rawpos2[latin] = j*ystep+1.0;
      raw_to_wcs(wcsin, rawpos2, wcspos2);
      wcspos[lngout] = wcspos2[lngin];
      wcspos[latout] = wcspos2[latin];
      wcs_to_raw(wcsout, wcspos, rawpos);
      x[j] = rawpos[lngout];
      y[j] = rawpos[latout];
      }
/*-- Draw */
    plline(CPLOT_FGRIDLINES, x,y);
    }

  return RETURN_OK;
  }


/****** cplot_drawcoordgrid **************************************************
PROTO	int cplot_drawcoordgrid(wcsstruct *wcs, double xmin, double xmax,
				double ymin, double ymax)
PURPOSE	Draw an atlas-like grid of angular coordinates.
INPUT	Pointer to the WCS projection structure,
	left projected coordinate limit,
	right projected coordinate limit,
	bottom projected coordinate limit,
	top projected coordinate limit.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	26/10/2005
 ***/
int	cplot_drawcoordgrid(wcsstruct *wcs, double xmin, double xmax,
				double ymin, double ymax)
  {
  wcsstruct	*wcs2;
  char		str[32];
  PLFLT		*x,*y;
  double	rawpos[NAXIS], wcspos[NAXIS],
		dx,dy, alphabeg,deltabeg, alphaend,deltaend,
		alphastep,deltastep, xd,yd,xdo,ydo, xm,ym, xmd,ymd, xmu,ymu;
  int		i, lng,lat;

  lng = wcs->lng;
  lat = wcs->lat;
/* Exit if we are not dealing with angular coordinates */
  if (lng==lat)
    return RETURN_ERROR;

/* A small trick to compute the min and max longitudes and latitudes */
  wcs2 = copy_wcs(wcs);
  dx = xmax - xmin;
  dy = ymax - ymin;
  wcs2->naxisn[lng] = (int)(dx+1.0);
  wcs2->naxisn[lat] = (int)(dy+1.0);
  wcs2->crpix[lng] -= xmin;
  wcs2->crpix[lat] -= ymin;
  init_wcs(wcs2);
  range_wcs(wcs2);
  alphabeg = fmod(wcs2->wcsmin[lng]+360.0, 360.0);
  alphaend = fmod(wcs2->wcsmax[lng]+360.0, 360.0);
  while (alphaend<alphabeg)
    alphaend += 360.0;
  deltabeg = wcs2->wcsmin[lat];
  deltaend = wcs2->wcsmax[lat];
  end_wcs(wcs2);

  alphastep = alphaend - alphabeg;
/* Quantize at the hour level */
  if (alphastep*DEG>60.0*DEG)
     alphastep = 15.0*DEG/DEG;
/* Quantize at the 15 minutes level */
  else if (alphastep*DEG>15.0*DEG)
    alphastep = 3.75*DEG/DEG;
/* Quantize at the 5 minutes level */
  else if (alphastep*DEG>3.75*DEG)
    alphastep = 75.0*ARCMIN/DEG;
/* Quantize at the 1 minute level */
  else if (alphastep*DEG>75.0*ARCMIN)
    alphastep = 30.0*ARCMIN/DEG;
/* Quantize at the 20 seconds level */
  else if (alphastep*DEG>30.0*ARCMIN)
    alphastep = 10.0*ARCMIN/DEG;
/* Quantize at the 10 seconds level */
  else if (alphastep*DEG>10.0*ARCMIN)
    alphastep = 150.0*ARCSEC/DEG;
/* Quantize at the 2 seconds level */
  else if (alphastep*DEG>150.0*ARCSEC)
    alphastep = 30.0*ARCSEC/DEG;
/* Quantize at the 0.5 seconds level */
  else if (alphastep*DEG>37.5*ARCSEC)
    alphastep = 7.5*ARCSEC/DEG;
/* Quantize at the 0.1 seconds level */
  else if (alphastep*DEG>7.5*ARCSEC)
    alphastep = 1.5*ARCSEC/DEG;
/* Quantize at the 0.01 seconds level */
  else
    alphastep = 0.15*ARCSEC/DEG;
  alphabeg -= fmod(alphabeg, alphastep) + alphastep;
  alphabeg = fmod (alphabeg+360.0, 360.0);
  alphaend += alphastep;
  while (alphaend<alphabeg)
    alphaend += 360.0;

  deltastep = deltaend - deltabeg;
/* Quantize at the 15 degrees level */
  if (deltastep*DEG>45.0*DEG)
    deltastep = 15.0*DEG/DEG;
/* Quantize at the 5 degrees level */
  else if (deltastep*DEG>15.0*DEG)
    deltastep = 5.0*DEG/DEG;
/* Quantize at the 1 degree level */
  else if (deltastep*DEG>5.0*DEG)
    deltastep = 1.0*DEG/DEG;
/* Quantize at the 20 arcmin level */
  else if (deltastep*DEG>1.0*DEG)
    deltastep = 20.0*ARCMIN/DEG;
/* Quantize at the 5 arcmin level */
  else if (deltastep*DEG>20.0*ARCMIN)
    deltastep = 5.0*ARCMIN/DEG;
/* Quantize at the 1 arcmin level */
  else if (deltastep*DEG>5.0*ARCMIN)
    deltastep = 1.0*ARCMIN/DEG;
/* Quantize at the 20 arcsec level */
  else if (deltastep*DEG>1.0*ARCMIN)
    deltastep = 20.0*ARCSEC/DEG;
/* Quantize at the 5 arcsec level */
  else if (deltastep*DEG>20.0*ARCSEC)
    deltastep = 5.0*ARCSEC/DEG;
/* Quantize at the 1 arcsec level */
  else if (deltastep*DEG>5.0*ARCSEC)
    deltastep = 1.0*ARCSEC/DEG;
/* Quantize at the 0.25 arcsec level */
  else if (deltastep*DEG>1.0*ARCSEC)
    deltastep = 0.25*ARCSEC/DEG;
  else
    deltastep = 0.1*ARCSEC/DEG;
  deltabeg -= fmod(deltabeg, deltastep) + deltastep;
  deltaend += deltastep;
  if (deltabeg<-90.0)
    deltabeg = -90.0;
  if (deltaend>90.0)
    deltabeg = 90.0;

  QMALLOC(x, PLFLT, CPLOT_NPOINTDEF);
  QMALLOC(y, PLFLT, CPLOT_NPOINTDEF);

/* Draw meridians */
  plschr(0.0, 0.33);
  plwid(0);
  pllsty(2);
  xmd = xmu = xdo = -0.5;
  ymd = ymu = ydo = -0.5;
  for (wcspos[0] = alphabeg; wcspos[0]<=alphaend; wcspos[0] += alphastep)
    {
    i=0;
    for (wcspos[1]=deltabeg; wcspos[1]<deltaend && i<CPLOT_NPOINTDEF;
	wcspos[1]+=deltastep, i++)
      {
      wcs_to_raw(wcs, wcspos, rawpos);
      x[i] = (PLFLT)(xd = rawpos[0]);
      y[i] = (PLFLT)(yd = rawpos[1]);
      if (i>0)
        {
        if ((yd-ymin)*(ydo-ymin) < 0.0 
		&& (xm = (xd-xmin)/dx) > 0.0 && xm < 1.0
		&& fabs(xm-xmd) > 0.1)
          {
          plmtex("b", 2.0, (PLFLT)xm, 0.5, cplot_degtosexal(str,wcspos[0],
		alphastep));
          xmd = xm;
          }
        if ((yd-ymax)*(ydo-ymax) < 0.0
		&& (xm = (xd-xmin)/dx) > 0.0 && xm < 1.0
		&& fabs(xm-xmu) > 0.1)
          {
          plmtex("t", 1.5, (PLFLT)xm, 0.5, cplot_degtosexal(str,wcspos[0],
		alphastep));
          xmu = xm;
          }
        }
      xdo = xd;
      ydo = yd;
      }
    plline(i,x,y);
    }

/* Draw parallels */
  for (wcspos[1] = deltabeg; wcspos[1]<=deltaend; wcspos[1] += deltastep)
    {
    i=0;
    for (wcspos[0]=alphabeg; wcspos[0]<alphaend && i<CPLOT_NPOINTDEF;
	wcspos[0]+=alphastep, i++)
      {
      wcs_to_raw(wcs, wcspos, rawpos);
      x[i] = (PLFLT)(xd = rawpos[0]);
      y[i] = (PLFLT)(yd = rawpos[1]);
      if (i>0)
        {
        if ((xd-xmin)*(xdo-xmin) < 0.0
		&& (ym = (yd-ymin)/dy) > 0.0 && ym < 1.0
		&& fabs(ym-ymd) > 0.1)
          {
          plmtex("lv", 1.0, (PLFLT)ym, 1.0, cplot_degtosexde(str,wcspos[1],
		deltastep));
          ymd = ym;
          }
        if ((xd-xmax)*(xdo-xmax) < 0.0
		&& (ym = (yd-ymin)/dy) > 0.0 && ym < 1.0
		&& fabs(ym-ymu) > 0.1)
          {
          plmtex("rv", 1.0, (PLFLT)ym, 0.0, cplot_degtosexde(str,wcspos[1],
		deltastep));
          ymu = ym;
          }
	}
      xdo = xd;
      ydo = yd;
      }
    plline(i,x,y);
    }

  free(x);
  free(y);

  return RETURN_OK;

  }

/****** cplot_drawloccoordgrid ************************************************
PROTO	int cplot_drawloccoordgrid(wcsstruct *wcs, double xmin, double xmax,
				double ymin, double ymax)
PURPOSE	Draw an atlas-like grid of angular coordinates with respect to
	projection center.
INPUT	Pointer to the WCS projection structure,
	left projected coordinate limit,
	right projected coordinate limit,
	bottom projected coordinate limit,
	top projected coordinate limit.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	26/10/2005
 ***/
int	cplot_drawloccoordgrid(wcsstruct *wcs, double xmin, double xmax,
				double ymin, double ymax)
  {
  wcsstruct	*wcs2;
  char		str[32];
  PLFLT		*x,*y;
  double	rawpos[NAXIS], wcspos[NAXIS],
		dx,dy, alphabeg,deltabeg, alphaend,deltaend,
		alphastep,deltastep, xd,yd,xdo,ydo, xm,ym, xmd,ymd, xmu,ymu;
  int		i, lng,lat;

  lng = wcs->lng;
  lat = wcs->lat;
/* Exit if we are not dealing with angular coordinates */
  if (lng==lat)
    return RETURN_ERROR;

/* A small trick to compute the min and max longitudes and latitudes */
  wcs2 = copy_wcs(wcs);
  dx = xmax - xmin;
  dy = ymax - ymin;
  wcs2->naxisn[lng] = (int)(dx+1.0);
  wcs2->naxisn[lat] = (int)(dy+1.0);
  wcs2->crpix[lng] -= xmin;
  wcs2->crpix[lat] -= ymin;
  wcs2->crval[lng] = 0.0;
  wcs2->crval[lat] = 0.0;
  init_wcs(wcs2);
  range_wcs(wcs2);
  alphabeg = fmod(wcs2->wcsmin[lng], 360.0);
  alphaend = fmod(wcs2->wcsmax[lng], 360.0);
  while (alphaend<alphabeg)
    alphaend += 360.0;
  deltabeg = wcs2->wcsmin[lat];
  deltaend = wcs2->wcsmax[lat];

  alphastep = alphaend - alphabeg;
/* Quantize at the 15 degrees level */
  if (alphastep*DEG>45.0*DEG)
     alphastep = 15.0*DEG/DEG;
/* Quantize at the 5 degrees level */
  else if (alphastep*DEG>15.0*DEG)
    alphastep = 5.0*DEG/DEG;
/* Quantize at the 1 degree level */
  else if (alphastep*DEG>5.0*DEG)
    alphastep = 1.0*DEG/DEG;
/* Quantize at the 20 arcmin level */
  else if (alphastep*DEG>1.0*DEG)
    alphastep = 20.0*ARCMIN/DEG;
/* Quantize at the 5 arcmin level */
  else if (alphastep*DEG>20.0*ARCMIN)
    alphastep = 5.0*ARCMIN/DEG;
/* Quantize at the 1 arcmin level */
  else if (alphastep*DEG>5.0*ARCMIN)
    alphastep = 1.0*ARCMIN/DEG;
/* Quantize at the 20 arcsec level */
  else if (alphastep*DEG>1.0*ARCMIN)
    alphastep = 20.0*ARCSEC/DEG;
/* Quantize at the 5 arcsec level */
  else if (alphastep*DEG>20.0*ARCSEC)
    alphastep = 5.0*ARCSEC/DEG;
/* Quantize at the 1 arcsec level */
  else if (alphastep*DEG>5.0*ARCSEC)
    alphastep = 1.0*ARCSEC/DEG;
/* Quantize at the 0.25 arcsec level */
  else if (alphastep*DEG>1.0*ARCSEC)
    alphastep = 0.25*ARCSEC/DEG;
  else
    alphastep = 0.1*ARCSEC/DEG;
  alphabeg -= fmod(alphabeg, alphastep) + alphastep;
  alphabeg = fmod (alphabeg, 360.0);
  alphaend += alphastep;
  while (alphaend<alphabeg)
    alphaend += 360.0;

  deltastep = deltaend - deltabeg;
/* Quantize at the 15 degrees level */
  if (deltastep*DEG>45.0*DEG)
    deltastep = 15.0*DEG/DEG;
/* Quantize at the 5 degrees level */
  else if (deltastep*DEG>15.0*DEG)
    deltastep = 5.0*DEG/DEG;
/* Quantize at the 1 degree level */
  else if (deltastep*DEG>5.0*DEG)
    deltastep = 1.0*DEG/DEG;
/* Quantize at the 20 arcmin level */
  else if (deltastep*DEG>1.0*DEG)
    deltastep = 20.0*ARCMIN/DEG;
/* Quantize at the 5 arcmin level */
  else if (deltastep*DEG>20.0*ARCMIN)
    deltastep = 5.0*ARCMIN/DEG;
/* Quantize at the 1 arcmin level */
  else if (deltastep*DEG>5.0*ARCMIN)
    deltastep = 1.0*ARCMIN/DEG;
/* Quantize at the 20 arcsec level */
  else if (deltastep*DEG>1.0*ARCMIN)
    deltastep = 20.0*ARCSEC/DEG;
/* Quantize at the 5 arcsec level */
  else if (deltastep*DEG>20.0*ARCSEC)
    deltastep = 5.0*ARCSEC/DEG;
/* Quantize at the 1 arcsec level */
  else if (deltastep*DEG>5.0*ARCSEC)
    deltastep = 1.0*ARCSEC/DEG;
/* Quantize at the 0.25 arcsec level */
  else if (deltastep*DEG>1.0*ARCSEC)
    deltastep = 0.25*ARCSEC/DEG;
  else
    deltastep = 0.1*ARCSEC/DEG;
  deltabeg -= fmod(deltabeg, deltastep) + deltastep;
  deltaend += deltastep;
  if (deltabeg<-90.0)
    deltabeg = -90.0;
  if (deltaend>90.0)
    deltabeg = 90.0;

  QMALLOC(x, PLFLT, CPLOT_NPOINTDEF);
  QMALLOC(y, PLFLT, CPLOT_NPOINTDEF);

/* Draw meridians */
  plschr(0.0, 0.33);
  plwid(0);
  pllsty(2);
  xmd = xmu = xdo = -0.5;
  ymd = ymu = ydo = -0.5;
  for (wcspos[0] = alphabeg; wcspos[0]<=alphaend; wcspos[0] += alphastep)
    {
    i=0;
    for (wcspos[1]=deltabeg; wcspos[1]<deltaend && i<CPLOT_NPOINTDEF;
	wcspos[1]+=deltastep, i++)
      {
      wcs_to_raw(wcs2, wcspos, rawpos);
      x[i] = (PLFLT)(xd = rawpos[0]);
      y[i] = (PLFLT)(yd = rawpos[1]);
      if (i>0)
        {
        if ((yd-ymin)*(ydo-ymin) < 0.0 
		&& (xm = (xd-xmin)/dx) > 0.0 && xm < 1.0
		&& fabs(xm-xmd) > 0.1)
          {
          plmtex("b", 2.0, (PLFLT)xm, 0.5, cplot_degtosexde(str,wcspos[0],
		alphastep));
          xmd = xm;
          }
        if ((yd-ymax)*(ydo-ymax) < 0.0
		&& (xm = (xd-xmin)/dx) > 0.0 && xm < 1.0
		&& fabs(xm-xmu) > 0.1)
          {
          plmtex("t", 1.5, (PLFLT)xm, 0.5, cplot_degtosexde(str,wcspos[0],
		alphastep));
          xmu = xm;
          }
        }
      xdo = xd;
      ydo = yd;
      }
    plline(i,x,y);
    }

/* Draw parallels */
  for (wcspos[1] = deltabeg; wcspos[1]<=deltaend; wcspos[1] += deltastep)
    {
    i=0;
    for (wcspos[0]=alphabeg; wcspos[0]<alphaend && i<CPLOT_NPOINTDEF;
	wcspos[0]+=alphastep, i++)
      {
      wcs_to_raw(wcs2, wcspos, rawpos);
      x[i] = (PLFLT)(xd = rawpos[0]);
      y[i] = (PLFLT)(yd = rawpos[1]);
      if (i>0)
        {
        if ((xd-xmin)*(xdo-xmin) < 0.0
		&& (ym = (yd-ymin)/dy) > 0.0 && ym < 1.0
		&& fabs(ym-ymd) > 0.1)
          {
          plmtex("lv", 1.0, (PLFLT)ym, 1.0, cplot_degtosexde(str,wcspos[1],
		deltastep));
          ymd = ym;
          }
        if ((xd-xmax)*(xdo-xmax) < 0.0
		&& (ym = (yd-ymin)/dy) > 0.0 && ym < 1.0
		&& fabs(ym-ymu) > 0.1)
          {
          plmtex("rv", 1.0, (PLFLT)ym, 0.0, cplot_degtosexde(str,wcspos[1],
		deltastep));
          ymu = ym;
          }
	}
      xdo = xd;
      ydo = yd;
      }
    plline(i,x,y);
    }

  free(x);
  free(y);
  end_wcs(wcs2);

  return RETURN_OK;

  }


/****** cplot_fgroup *******************************************************
PROTO	int cplot_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
PURPOSE	Plot photometric relations between fgroups and the reference field
	(if present).
INPUT	Pointer to the field group,
	Reference field to the group.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	19/03/2007
 ***/
int	cplot_fgroup(fgroupstruct *fgroup, fieldstruct *reffield)
  {
   wcsstruct	*wcs;
   fieldstruct	*field;
   setstruct	*set, *refset;
   samplestruct	*samp;
   char		str[64];
   PLFLT	*x,*y, *x2,*y2, *xt,*yt, *x2t,*y2t, psize;
   PLINT	lwid, lsize;
   double	dx,dy, xmin,xmax,ymin,ymax;
   int		f,s,n, lng,lat, nsamp, npoint;

  x = y = x2 = y2 = NULL;	/* to avoid gcc -Wall warnings */
  psize = 0;			/* to avoid gcc -Wall warnings */
  if (cplot_init(1,1, CPLOT_FGROUPS) == RETURN_ERROR)
    {
    cplot_end(CPLOT_FGROUPS);
    return RETURN_OK;
    }

  wcs = fgroup->wcs;
  lng = fgroup->lng;
  lat = fgroup->lat;
  refset = reffield->set[0];
  dx = fgroup->projposmax[lng] - fgroup->projposmin[lng];
  dy = fgroup->projposmax[lat] - fgroup->projposmin[lat];
  if (dy>dx)
    dx = dy;
  xmin = 0.5*(fgroup->projposmin[lng]+fgroup->projposmax[lng]) - 0.55*dx;
  xmax = 0.5*(fgroup->projposmin[lng]+fgroup->projposmax[lng]) + 0.55*dx;
  ymin = 0.5*(fgroup->projposmin[lat]+fgroup->projposmax[lat]) - 0.55*dx;
  ymax = 0.5*(fgroup->projposmin[lat]+fgroup->projposmax[lat]) + 0.55*dx;

/* Use the total number of points to scale symbols */
  nsamp = 0;
  for (f=0; f<fgroup->nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      nsamp += field->set[s]->nsample;
    }
  lsize = 0;
  if (nsamp>0)
    {
    psize = 30.0/sqrt(nsamp);
    if (psize<0.2)
      psize = 0.2;
    if (psize>1.0)
      {
      lsize = (int)sqrt(psize);
      psize = 1.0;
      }
    }
  else
    psize = 1.0;

  plfont(2);
  plcol(15);
  plschr(0.0, 0.67);
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uGroup ##%d: detections", fgroup->no);
  pllab("","", str);
  nsamp = refset->nsample;
  if (nsamp)
    {
    QMALLOC(x, PLFLT, nsamp);
    QMALLOC(y, PLFLT, nsamp);
    QMALLOC(x2, PLFLT, nsamp);
    QMALLOC(y2, PLFLT, nsamp);
    }
  xt = x;
  yt = y;
  x2t = x2;
  y2t = y2;
  samp = refset->sample;
  npoint = 0;
  for (n=nsamp; n--; samp++)
    {
    if (samp->nextsamp)
      {
      *(xt++) = (PLFLT)samp->projpos[lng];
      *(yt++) = (PLFLT)samp->projpos[lat];
      npoint++;
      }
    else
      {
      *(x2t++) = (PLFLT)samp->projpos[lng];
      *(y2t++) = (PLFLT)samp->projpos[lat];
      }
    }
  plssym(0.0, psize);
  plwid(lsize);
  plcol(3);
  plpoin((PLINT)npoint, x,y, 11);
  plcol(1);
  plpoin((PLINT)(nsamp-npoint), x2,y2, 0);
  free(x);
  free(y);
  free(x2);
  free(y2);

  for (f=0; f<fgroup->nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      if (nsamp)
        {
        QMALLOC(x, PLFLT, nsamp);
        QMALLOC(y, PLFLT, nsamp);
        xt = x;
        x2t = x+nsamp;
        yt = y;
        y2t = y+nsamp;
        samp = set->sample;
        for (n=nsamp; n--; samp++)
          if (samp->nextsamp)
            {
            *(--x2t) = (PLFLT)samp->projpos[lng];
            *(--y2t) = (PLFLT)samp->projpos[lat];
            }
          else if (!samp->prevsamp)
            {
            *(xt++) = (PLFLT)samp->projpos[lng];
            *(yt++) = (PLFLT)samp->projpos[lat];
            }
        plwid(lsize);
        plcol(8);
        if (lsize)
          plpoin((PLINT)(xt-x), x,y, 17);
        else
          plpoin((PLINT)(xt-x), x,y, 1);
        plcol(4);
        plpoin((PLINT)(x+nsamp-x2t), x2t,y2t, 2);
        free(x);
        free(y);
        }
      plwid(2*lwid);
      if (field->cplot_colour==15)
        {
        plcol(15);
        cplot_drawbounds(set->wcs, wcs);
        }
      }
    }

/* End with fields coloured differently */
  for (f=0; f<fgroup->nfield; f++)
    {
    field = fgroup->field[f];
    if (field->cplot_colour!=15)
      {
      plwid(3*lwid);
      plcol(field->cplot_colour);
      for (s=0; s<field->nset; s++)
        {
        set = field->set[s];
        cplot_drawbounds(set->wcs, wcs);
        }
      }
    }

  plcol(7);
  plwid(lwid);
  cplot_drawcoordgrid(wcs, xmin, xmax, ymin, ymax);
  plflush();

  plend();

  cplot_fgroup(fgroup, reffield);	/* Recursive stuff */

  return RETURN_OK;
  }

/****** distort_map *******************************************************   
PROTO   void distort_map(PLFLT x,PLFLT y,PLFLT *tx,PLFLT *ty, void *pltr_data)
PURPOSE Astrometric mapping function for shade plots in cplot_distort().
INPUT   TBW
OUTPUT  -.
NOTES   see http://plplot.sourceforge.net/docbook-manual/
            plplot-html-5.5.3/contour-plots.html#contour-plots-c
AUTHOR  E. Bertin (IAP)
VERSION 30/11/2005
***/
static void distort_map(PLFLT x,PLFLT y,PLFLT *tx,PLFLT *ty, void *pltr_data)
{
 wcsstruct	**wcs;
 double	 	rawpos[NAXIS], wcspos[NAXIS], wcspos2[NAXIS];
 int		i, lng, lat, naxis;


  wcs = (wcsstruct **)pltr_data;
  if ((naxis=wcs[0]->naxis) > 2)
    for (i=0; i<naxis; i++)
      rawpos[i]= wcs[0]->naxisn[i]/2.0 + 0.5;
  lng = wcs[0]->lng;
  lat = wcs[0]->lat;
  rawpos[lng] = (x+0.5)*wcs[0]->naxisn[lng]/CPLOT_NDISTGRID + 0.5;
  rawpos[lat] = (y+0.5)*wcs[0]->naxisn[lat]/CPLOT_NDISTGRID + 0.5;
  raw_to_wcs(wcs[0], rawpos, wcspos);
  wcspos2[wcs[1]->lng] = wcspos[lng];
  wcspos2[wcs[1]->lat] = wcspos[lat];
  wcs_to_raw(wcs[1], wcspos2, rawpos);
  *tx = rawpos[wcs[1]->lng];
  *ty = rawpos[wcs[1]->lat];

  return;
  }

/****** cplot_distort *******************************************************
PROTO	int cplot_distort(fieldstruct *field)
PURPOSE	Plot astrometric mappings for a given field (instrument).
INPUT	Pointer to the field.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	02/12/2005
 ***/
int	cplot_distort(fieldstruct *field)
  {
   wcsstruct	*wcsptr[2],
		*wcs, *wcsin;
   setstruct	*set;
   PLFLT	**scale,
		clevel[CPLOT_NSHADES], cpoint[3], r[3],g[3],b[3],
		cscale, scalemin,scalemax, mscale,dscale;
   PLINT	lwid;
   char		*ctype[NAXIS],
		str[64];
   double	crpix[NAXIS], cdelt[NAXIS], raw[NAXIS],
		xmin,ymin,xmax,ymax, xstep,ystep;
   int		naxisn[NAXIS],
		i,j, s, lng,lat, naxis;

  if (cplot_init(1,1, CPLOT_DISTORT) == RETURN_ERROR)
    {
    cplot_end(CPLOT_DISTORT);
    return RETURN_OK;
    }

  set = field->set[0];
  wcsin = set->wcs;
  if (!wcsin || wcsin->naxis<2)
    return RETURN_ERROR;
  naxis = wcsin->naxis;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(ctype[i], char, 16); 
    strncpy(ctype[i],wcsin->ctype[i], 16);
    crpix[i] = 50.0;
    cdelt[i] = field->maxradius/50.0;
    if (i==set->lng)
      cdelt[i] = -cdelt[i];	/* Put East to the left */
    naxisn[i] = 100;
    }
  wcs = create_wcs(ctype,field->meanwcspos,crpix,cdelt,naxisn, naxis);

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uInstrument A%d: distortion map", field->astromlabel+1);
  pllab("","", str);
  plwid(0);
  plcol(7);
  cplot_drawloccoordgrid(wcs, xmin, xmax, ymin, ymax);

  pllsty(1);
  plcol(15);
  plscmap1n(256);

  scalemin = BIG;
  scalemax = -BIG;

/* First pass through the data to find min and max pixel scales */
  for (s=0; s<field->nset; s++)
    {
    set = field->set[s];
    for (i=0; i<naxis; i++)
      raw[i] = set->wcs->naxisn[i]/2.0 + 0.5;
    lng = set->lng;
    lat = set->lat;
    xstep = set->wcs->naxisn[lng] / CPLOT_NDISTGRID;
    ystep = set->wcs->naxisn[lat] / CPLOT_NDISTGRID;
    raw[lat] = ystep / 2.0 + 0.5;
    for (j=0; j<CPLOT_NDISTGRID; j++)
      {
      raw[lng] = xstep / 2.0 + 0.5;
      for (i=0; i<CPLOT_NDISTGRID; i++)
        {
        cscale = sqrt(wcs_scale(set->wcs, raw));
        if (cscale<scalemin)
          scalemin = cscale;
        if (cscale>scalemax)
          scalemax = cscale;
        raw[lng] += xstep;
        }
      raw[lat] += ystep;
      }
    }

/* Lower bound to variability in pixel scale is 1e-6 */
  if ((mscale=(scalemin+scalemax)/2.0) < 1.0e-10*ARCSEC/DEG
       || (dscale=(scalemax-scalemin))/mscale < 1.0e-6)
    dscale = 1.0e-6;
  for (i=0; i<CPLOT_NSHADES; i++)
    clevel[i] = scalemin + (i-0.5) * dscale / (CPLOT_NSHADES-2);
  cpoint[0] = 0.0; r[0] = 0.5; g[0] = 0.5; b[0] = 1.0;
  cpoint[1] = 0.5; r[1] = 0.5; g[1] = 1.0; b[1] = 0.5;
  cpoint[2] = 1.0; r[2] = 1.0; g[2] = 0.5; b[2] = 0.5;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);

  plAlloc2dGrid(&scale, CPLOT_NDISTGRID, CPLOT_NDISTGRID);

/* Now the real 2D pixel-scale mapping */
  for (s=0; s<field->nset; s++)
    {
    set = field->set[s];
    for (i=0; i<naxis; i++)
      raw[i] = set->wcs->naxisn[i]/2.0 + 0.5;
    lng = set->lng;
    lat = set->lat;
    xstep = set->wcs->naxisn[lng] / CPLOT_NDISTGRID;
    ystep = set->wcs->naxisn[lat] / CPLOT_NDISTGRID;
    raw[lat] = ystep / 2.0 + 0.5;
    for (j=0; j<CPLOT_NDISTGRID; j++)
      {
      raw[lng] = xstep / 2.0 + 0.5;
      for (i=0; i<CPLOT_NDISTGRID; i++)
        {
        scale[i][j] = sqrt(wcs_scale(set->wcs, raw));
        raw[lng] += xstep;
        }
      raw[lat] += ystep;
      }

    wcsptr[0] = set->wcs;
    wcsptr[1] = wcs;
    plshades(scale, CPLOT_NDISTGRID, CPLOT_NDISTGRID, NULL,
	     xstep/2.0+0.5, set->wcs->naxisn[lng]-xstep/2.0+0.5,
             ystep/2.0+0.5, set->wcs->naxisn[lat]-ystep/2.0+0.5,
	     clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, distort_map, wcsptr);
/*
    plwid(0);
    cplot_drawfgrids(set->wcs, wcs);
*/
    plcol(7);
    plwid(lwid);
    cplot_drawbounds(set->wcs, wcs);
    }

  plFree2dGrid(scale, CPLOT_NDISTGRID, CPLOT_NDISTGRID);

/* Draw left colour scale */
  plAlloc2dGrid(&scale, 2, CPLOT_NSHADES);
  for (j=0; j<CPLOT_NSHADES; j++)
    scale[0][j] = scale[1][j] = scalemin + j * dscale/(CPLOT_NSHADES-1);

  plvpor(0.91,0.935,0.115,0.885);
  plwind(0.0,1.0,scalemin*DEG/ARCSEC,scalemax*DEG/ARCSEC);
  plshades(scale, 2, CPLOT_NSHADES, NULL, 0.0, 1.0,
	   scalemin*DEG/ARCSEC,scalemax*DEG/ARCSEC, clevel,
	   CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
  plcol(15);
  plschr(0.0, 0.5);
  plbox("bc", 0.0, 0, "bnstv", 0.0, 0);
  sprintf(str, "%s", mscale < 0.09*ARCSEC/DEG?
	  "mas" : (mscale < ARCMIN/DEG?
		   "arcsec" : (mscale < 1.0? "arcmin": "deg")));
  plschr(0.0, 0.6);
  plmtex("l", 5.0, 0.5, 0.0, str);
  plmtex("b", 2.0, 0.5, 0.5, "pixel scale");

/* Draw right colour scale */
  mscale /= 100.0;			/* convert to percentage */
  scalemin = scalemin/mscale - 100.0;
  scalemax = scalemax/mscale - 100.0;
  dscale /= mscale;
  plwind(0.0,1.0,scalemin,scalemax);
  plschr(0.0, 0.5);
  plbox("", 0.0, 0, "cmstv", 0.0, 0);
  plschr(0.0, 0.5);
  plmtex("r", 5.0, 0.5, 0.0, "%");

  plFree2dGrid(scale, 2, CPLOT_NSHADES);
  plend();
  end_wcs(wcs);
  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_distort(field);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_astintsysmap ****************************************************
PROTO	int cplot_astintsysmap(fgroupstruct **fgroups, int ngroup, int instru,
		double hsn_thresh)
PURPOSE	Plot a map of internal astrometric systematics for an instrument.
INPUT	Pointer to the field group,
	Reference field to the group.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	24/04/2008
 ***/
int	cplot_astintsysmap(fgroupstruct **fgroups, int ngroup, int instru,
		double hsn_thresh)
  {
   fgroupstruct	*fgroup, *fgroup0;
   fieldstruct	*field, *field0;
   setstruct	*set, *set0;
   samplestruct	*samp,*samp2;
   wcsstruct	*wcs, *wcs0;
   PLINT	lwid;
   PLFLT	x[2], y[2],
		scalel;
   char		*ctype[NAXIS],
		str[64], dispunit[32];
   double	vecshift[CPLOT_ASTNSUBPLOTS*CPLOT_ASTNSUBPLOTS][NAXIS],
		vecvar[CPLOT_ASTNSUBPLOTS*CPLOT_ASTNSUBPLOTS][NAXIS][NAXIS],
		rawpos[NAXIS],vecpos[NAXIS],vecpos2[NAXIS],wcspos[NAXIS],
		crpix[NAXIS], cdelt[NAXIS],
		mean[NAXIS], covar[NAXIS][NAXIS], dx[NAXIS],
		xmin,ymin,xmax,ymax, xoff,yoff,xstep,ystep, minscale,scale,
		dispscale;
   int		naxisn[NAXIS], nvec[CPLOT_ASTNSUBPLOTS*CPLOT_ASTNSUBPLOTS],
		d,d2,d3,f,g,i,n,s,  fx,fy, nx,ny, lng,lat,
		count, naxis,nsamp, nset;

  if (cplot_init(1,1, CPLOT_ADSYSMAP2D) == RETURN_ERROR)
    {
    cplot_end(CPLOT_ADSYSMAP2D);
    return RETURN_OK;
    }

/* Compute instrument projection (for display only) */
  wcs0 = NULL;
  fgroup0 = NULL;	/* to avoid gcc -Wall warnings */
  field0 = NULL;	/* to avoid gcc -Wall warnings */
  set0 = NULL;		/* to avoid gcc -Wall warnings */
  naxis = 0;		/* to avoid gcc -Wall warnings */
  nset = 0;

  for (g=0; g<ngroup; g++)
    {
    fgroup = fgroups[g];
    for (f=0; f<fgroup->nfield; f++)
      {
      field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
      if (field->astromlabel == instru && !wcs0)
        {
        field0 = field;
        fgroup0 = fgroup;
        nset = field0->nset;
        set0 = field0->set[0];
        wcs = set0->wcs;
        if (!wcs || wcs->naxis<2)
          return RETURN_ERROR;
        naxis = wcs->naxis;
        for (i=0; i<naxis; i++)
          {
          QMALLOC(ctype[i], char, 16); 
          strncpy(ctype[i],wcs->ctype[i], 16);
          crpix[i] = 50.0;
          cdelt[i] = field0->maxradius/50.0;
          if (i==set0->lng)
            cdelt[i] = -cdelt[i];	/* Put East to the left */
          naxisn[i] = 100;
          }
        wcs0 = create_wcs(ctype,field0->meanwcspos,crpix,cdelt,naxisn, naxis);
        f = fgroup->nfield;
        g = ngroup;
        }
      }
    }

  if (!nset)
    return RETURN_ERROR;

/* The scaling factor is computed in a way to avoid overlaps between vectors */
  minscale = BIG;
  for (s=0; s<nset; s++)
    {
    set = field0->set[s];
    lng = set->lng;
    lat = set->lat;
    scale = set->wcs->naxisn[lng]/(double)CPLOT_ASTNSUBPLOTS/2.0
		* field0->meanwcsscale[lng] / (0.1*fgroup0->sig_interr[lng]);
    if (minscale > scale)
      minscale = scale;
    scale = set->wcs->naxisn[lat] / (double)CPLOT_ASTNSUBPLOTS / 2.0
		* field0->meanwcsscale[lat] / (0.1*fgroup0->sig_interr[lat]);
    if (minscale > scale)
      minscale = scale;
    }

/* Compute the displayed scale */
  for (d=0; d<naxis; d++)
    rawpos[d] = wcs0->naxisn[d]/2.0;
  dispscale = sqrt(wcs_scale(wcs0, rawpos));			/* in deg */
  if ((scalel = 10.0*minscale/dispscale)     < 30.0)
    strcpy(dispunit, "10 #(2218)");
  else if ((scalel = 1.0*minscale/dispscale)  < 30.0)	/* 1 deg */
    strcpy(dispunit, "1 #(2218)");
  else if ((scalel = minscale/dispscale/6.0)  < 30.0)	/* 10 ' */
    strcpy(dispunit, "10 #(2216)");
  else if ((scalel = minscale/dispscale/60.0) < 30.0)	/* 1 ' */
    strcpy(dispunit, "1 #(2216)");
  else if ((scalel = minscale/dispscale/360.0) < 30.0)	/* 10 '' */
    strcpy(dispunit, "10 #(2216)#(2216)");
  else if ((scalel = minscale/dispscale/3600.0) < 30.0)	/* 1 '' */
    strcpy(dispunit, "1 #(2216)#(2216)");
  else if ((scalel = minscale/dispscale/36000.0) < 30.0)/* 100 mas */
    strcpy(dispunit, "100 mas");
  else if ((scalel = minscale/dispscale/3.6e5) < 30.0)	/* 10 mas */
    strcpy(dispunit, "10 mas");
  else if ((scalel = minscale/dispscale/3.6e6) < 30.0)	/* 1 mas */
    strcpy(dispunit, "1 mas");
  else if ((scalel = minscale/dispscale/3.6e7) < 30.0)	/* 100 muas */
    strcpy(dispunit, "100 #gmas");
  else if ((scalel = minscale/dispscale/3.6e8) < 30.0)	/* 10 muas */
    strcpy(dispunit, "10 #gmas");
  else
    {
    scalel = minscale/dispscale/3.6e9;			/* 1 muas */
    strcpy(dispunit, "1 #gmas");
    }

  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  sprintf(str, "#uInstrument A%d: map of astrometric systematics (internal)", 
	field0->astromlabel+1);
  pllab("","", str);
  plwid(0);
  plcol(7);
  cplot_drawloccoordgrid(wcs0, xmin, xmax, ymin, ymax);
  y[0] = y[1] = 3.0;
  x[0] = 7.0;
  x[1] = x[0] + scalel;
  pllsty(1);
  plcol(15);
  plwid(lwid*3);
  plline(2,x,y);
  plwid(lwid);
  plschr(0.0, 0.5);
  plptex((x[0] + x[1]) / 2.0, y[0] + 2.5, 1.0, 0.0, 0.5, dispunit);
  x[0] = x[1] = 7.0;
  y[0] = 2.5;
  y[1] = 3.5;
  plline(2,x,y);
  x[0] = x[1] = 7.0 + scalel;
  y[0] = 2.5;
  y[1] = 3.5;
  plline(2,x,y);
  for (s=0; s<nset; s++)
    {
    set0 = field0->set[s];
    plcol(7);
    cplot_drawbounds(set0->wcs, wcs0);
    lng = set0->lng;
    lat = set0->lat;
    nx = set0->wcs->naxisn[lng];
    ny = set0->wcs->naxisn[lat];
    xstep = nx / (double)CPLOT_ASTNSUBPLOTS;
    xoff = xstep/2.0 + 0.5;
    ystep = ny / (double)CPLOT_ASTNSUBPLOTS;
    yoff = ystep/2.0 + 0.5;
    for (g=0; g<ngroup; g++)
      {
      fgroup = fgroups[g];
      for (d3=0; d3<CPLOT_ASTNSUBPLOTS*CPLOT_ASTNSUBPLOTS; d3++)
        {
        for (d=0; d<naxis; d++)
          {
          vecshift[d3][d] = 0.0;
          for (d2=0; d2<naxis; d2++)
            vecvar[d3][d][d2] = 0.0;
          }
        nvec[d3] = 0;
        }
      for (f=0; f<fgroup->nfield; f++)
        {
        field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
        if (field->astromlabel != instru)
          continue;
        set = field->set[s];
        samp = set->sample;
        nsamp = set->nsample;
        for (n=nsamp; n--; samp++)
          {
          if (samp->flags & (OBJ_SATUR|OBJ_TRUNC))
            continue;
/*-------- Reset mean and count */
          for (d=0; d<naxis; d++)
            {
            mean[d]  = 0.0;
            for (d2=0; d2<naxis; d2++)
              covar[d][d2] = 0.0;
            }
          count = 0;
/*-------- Explore forward */
          samp2 = samp;
          while ((samp2=samp2->nextsamp))
            {
            for (d=0; d<naxis; d++)
              {
              dx[d] = (samp2->projpos[d] - samp->projpos[d]);
              mean[d] += dx[d];
              for (d2=0; d2<=d; d2++)
                covar[d][d2] += dx[d]*dx[d2];
              }
            count++;
            }
/*-------- Explore backward */
          samp2 = samp;
          while ((samp2=samp2->prevsamp) && samp2->set->field->astromlabel>=0)
            {
            if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
              continue;
            for (d=0; d<naxis; d++)
              {
              dx[d] = (samp2->projpos[d] - samp->projpos[d]);
              mean[d] += dx[d];
              for (d2=0; d2<=d; d2++)
                covar[d][d2] += dx[d]*dx[d2];
              }
            count++;
            }
          if (count)
            {
            fx = (int)((samp->rawpos[lng]-0.5)/xstep);
            if (fx<0)
              fx = 0;
            else if (fx >= CPLOT_ASTNSUBPLOTS)
              fx = CPLOT_ASTNSUBPLOTS - 1;
            fy = (int)((samp->rawpos[lat]-0.5)/ystep);
            if (fy<0)
              fy = 0;
            else if (fy >= CPLOT_ASTNSUBPLOTS)
              fy = CPLOT_ASTNSUBPLOTS - 1;
            d3 = fy*CPLOT_ASTNSUBPLOTS + fx;
            for (d=0; d<naxis; d++)
              {
              vecshift[d3][d] += mean[d]/count;
              for (d2=0; d2<=d; d2++)
                vecvar[d3][d][d2] += (covar[d][d2]-mean[d]*mean[d2])
				/(count>1?(count-1):1);
              }
            nvec[d3]++;
            }
          }
        }
/*---- Compute coordinates of vectors */
      for (d=0; d<naxis; d++)
        rawpos[d] = set0->wcs->naxisn[d]/2.0;
      for (fy=0; fy<CPLOT_ASTNSUBPLOTS; fy++)
        for (fx=0; fx<CPLOT_ASTNSUBPLOTS; fx++)
          {
          d3 = fy*CPLOT_ASTNSUBPLOTS + fx;
          rawpos[lng] = xoff + xstep*fx;
          rawpos[lat] = yoff + ystep*fy;
/*-------- Compute application point coordinates in the display system */
          raw_to_wcs(set0->wcs, rawpos, wcspos);
          wcs_to_raw(wcs0, wcspos, vecpos);
/*-------- Compute application point coordinates in the group system */
          wcs_to_raw(fgroup->wcs, wcspos, rawpos);
/*-------- Add shift */
          if (nvec[d3])
            for (d=0; d<naxis; d++)
              rawpos[d] += vecshift[d3][d]/nvec[d3];
/*-------- Compute vector extremity coordinates in the display system */
          raw_to_wcs(fgroup->wcs, rawpos, wcspos);
          wcs_to_raw(wcs0, wcspos, vecpos2);
          x[0] = vecpos[lng];
          x[1] = vecpos2[lng] + minscale*(vecpos2[lng] - vecpos[lng]);
          y[0] = vecpos[lat];
          y[1] = vecpos2[lat] + minscale*(vecpos2[lat] - vecpos[lat]);
          plwid(lwid*2);
          plcol(3);
          plline(2, x,y);
          plcol(15);
          plpoin(1,x,y,1);
          plwid(lwid);
          }
      }
    }

  plend();

  end_wcs(wcs0);

  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_astintsysmap(fgroups, ngroup, instru, hsn_thresh);/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_astrefsysmap ****************************************************
PROTO	int cplot_astrefsysmap(fgroupstruct **fgroups, int ngroup, int instru,
				double hsn_thresh)
PURPOSE	Plot a map of referenc astrometric systematics for an instrument.
INPUT	Pointer to the field group,
	Reference field to the group.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	24/04/2008
 ***/
int	cplot_astrefsysmap(fgroupstruct **fgroups, int ngroup, int instru,
			double hsn_thresh)
  {
   fgroupstruct	*fgroup, *fgroup0;
   fieldstruct	*field, *field0;
   setstruct	*set, *set0;
   samplestruct	*samp,*samp2;
   wcsstruct	*wcs, *wcs0;
   PLINT	lwid;
   PLFLT	x[2], y[2],
		scalel;
   char		*ctype[NAXIS],
		str[64], dispunit[32];
   double	vecshift[CPLOT_ASTNSUBPLOTS*CPLOT_ASTNSUBPLOTS][NAXIS],
		vecvar[CPLOT_ASTNSUBPLOTS*CPLOT_ASTNSUBPLOTS][NAXIS][NAXIS],
		rawpos[NAXIS],vecpos[NAXIS],vecpos2[NAXIS],wcspos[NAXIS],
		crpix[NAXIS], cdelt[NAXIS],
		mean[NAXIS], covar[NAXIS][NAXIS], dx[NAXIS],
		xmin,ymin,xmax,ymax, xoff,yoff,xstep,ystep, minscale,scale,
		dispscale;
   int		naxisn[NAXIS], nvec[CPLOT_ASTNSUBPLOTS*CPLOT_ASTNSUBPLOTS],
		d,d2,d3,f,g,i,n,s,  fx,fy, nx,ny, lng,lat,
		count, naxis,nsamp, nset;

  if (cplot_init(1,1, CPLOT_REFSYSMAP2D) == RETURN_ERROR)
    {
    cplot_end(CPLOT_REFSYSMAP2D);
    return RETURN_OK;
    }

/* Compute instrument projection (for display only) */
  wcs0 = NULL;
  fgroup0 = NULL;	/* to avoid gcc -Wall warnings */
  field0 = NULL;	/* to avoid gcc -Wall warnings */
  set0 = NULL;		/* to avoid gcc -Wall warnings */
  naxis = 0;		/* to avoid gcc -Wall warnings */
  nset = 0;

  for (g=0; g<ngroup; g++)
    {
    fgroup = fgroups[g];
    for (f=0; f<fgroup->nfield; f++)
      {
      field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
      if (field->astromlabel == instru && !wcs0)
        {
        field0 = field;
        fgroup0 = fgroup;
        nset = field0->nset;
        set0 = field0->set[0];
        wcs = set0->wcs;
        if (!wcs || wcs->naxis<2)
          return RETURN_ERROR;
        naxis = wcs->naxis;
        for (i=0; i<naxis; i++)
          {
          QMALLOC(ctype[i], char, 16); 
          strncpy(ctype[i],wcs->ctype[i], 16);
          crpix[i] = 50.0;
          cdelt[i] = field0->maxradius/50.0;
          if (i==set0->lng)
            cdelt[i] = -cdelt[i];	/* Put East to the left */
          naxisn[i] = 100;
          }
        wcs0 = create_wcs(ctype,field0->meanwcspos,crpix,cdelt,naxisn, naxis);
        f = fgroup->nfield;
        g = ngroup;
        }
      }
    }

  if (!nset)
    return RETURN_ERROR;

/* The scaling factor is computed in a way to avoid overlaps between vectors */
  minscale = BIG;
  for (s=0; s<nset; s++)
    {
    set = field0->set[s];
    lng = set->lng;
    lat = set->lat;
    scale = set->wcs->naxisn[lng]/(double)CPLOT_ASTNSUBPLOTS/2.0
		* field0->meanwcsscale[lng] / (0.333*fgroup0->sig_referr[lng]);
    if (minscale > scale)
      minscale = scale;
    scale = set->wcs->naxisn[lat] / (double)CPLOT_ASTNSUBPLOTS / 2.0
		* field0->meanwcsscale[lat] / (0.333*fgroup0->sig_referr[lat]);
    if (minscale > scale)
      minscale = scale;
    }

/* Compute the displayed scale */
  for (d=0; d<naxis; d++)
    rawpos[d] = wcs0->naxisn[d]/2.0;
  dispscale = sqrt(wcs_scale(wcs0, rawpos));	/* in deg */
  if ((scalel = 10.0*minscale/dispscale)     < 30.0)
    strcpy(dispunit, "10 #(2218)");
  else if ((scalel = 1.0*minscale/dispscale)  < 30.0)	/* 1 deg */
    strcpy(dispunit, "1 #(2218)");
  else if ((scalel = minscale/dispscale/6.0)  < 30.0)	/* 10 ' */
    strcpy(dispunit, "10 #(2216)");
  else if ((scalel = minscale/dispscale/60.0) < 30.0)	/* 1 ' */
    strcpy(dispunit, "1 #(2216)");
  else if ((scalel = minscale/dispscale/360.0) < 30.0)	/* 10 '' */
    strcpy(dispunit, "10 #(2216)#(2216)");
  else if ((scalel = minscale/dispscale/3600.0) < 30.0)	/* 1 '' */
    strcpy(dispunit, "1 #(2216)#(2216)");
  else if ((scalel = minscale/dispscale/36000.0) < 30.0)/* 100 mas */
    strcpy(dispunit, "100 mas");
  else if ((scalel = minscale/dispscale/3.6e5) < 30.0)	/* 10 mas */
    strcpy(dispunit, "10 mas");
  else if ((scalel = minscale/dispscale/3.6e6) < 30.0)	/* 1 mas */
    strcpy(dispunit, "1 mas");
  else if ((scalel = minscale/dispscale/3.6e7) < 30.0)	/* 100 muas */
    strcpy(dispunit, "100 #gmas");
  else if ((scalel = minscale/dispscale/3.6e8) < 30.0)	/* 10 muas */
    strcpy(dispunit, "10 #gmas");
  else
    {
    scalel = minscale/dispscale/3.6e9;			/* 1 muas */
    strcpy(dispunit, "1 #gmas");
    }
  xmin = 0.5;
  xmax = 100.5;
  ymin = 0.5;
  ymax = 100.5;
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plfont(2);
  plcol(15);
  plenv((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 1, -1);
  plschr(0.0, 0.8);
  sprintf(str, "#uInstrument A%d: map of astrometric systematics (reference)",
	field0->astromlabel+1);
  pllab("","", str);
  plwid(0);
  plcol(7);
  cplot_drawloccoordgrid(wcs0, xmin, xmax, ymin, ymax);
  y[0] = y[1] = 3.0;
  x[0] = 7.0;
  x[1] = x[0] + scalel;
  pllsty(1);
  plcol(15);
  plwid(lwid*3);
  plline(2,x,y);
  plwid(lwid);
  plschr(0.0, 0.5);
  plptex((x[0] + x[1]) / 2.0, y[0] + 2.5, 1.0, 0.0, 0.5, dispunit);
  x[0] = x[1] = 7.0;
  y[0] = 2.5;
  y[1] = 3.5;
  plline(2,x,y);
  x[0] = x[1] = 7.0 + scalel;
  y[0] = 2.5;
  y[1] = 3.5;
  plline(2,x,y);
  for (s=0; s<nset; s++)
    {
    set0 = field0->set[s];
    plcol(7);
    cplot_drawbounds(set0->wcs, wcs0);
    lng = set0->lng;
    lat = set0->lat;
    nx = set0->wcs->naxisn[lng];
    ny = set0->wcs->naxisn[lat];
    xstep = nx / (double)CPLOT_ASTNSUBPLOTS;
    xoff = xstep/2.0 + 0.5;
    ystep = ny / (double)CPLOT_ASTNSUBPLOTS;
    yoff = ystep/2.0 + 0.5;
    for (g=0; g<ngroup; g++)
      {
      fgroup = fgroups[g];
      for (d3=0; d3<CPLOT_ASTNSUBPLOTS*CPLOT_ASTNSUBPLOTS; d3++)
        {
        for (d=0; d<naxis; d++)
          {
          vecshift[d3][d] = 0.0;
          for (d2=0; d2<naxis; d2++)
            vecvar[d3][d][d2] = 0.0;
          }
        nvec[d3] = 0;
        }
      for (f=0; f<fgroup->nfield; f++)
        {
        field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
        if (field->astromlabel != instru)
          continue;
        set = field->set[s];
        samp = set->sample;
        nsamp = set->nsample;
        for (n=nsamp; n--; samp++)
          {
          if (samp->flags & (OBJ_SATUR|OBJ_TRUNC))
            continue;
/*-------- Explore backward */
          samp2 = samp;
          while ((samp2=samp2->prevsamp) && samp2->set->field->astromlabel>=0);
          if (!samp2 || (samp2->flags & (OBJ_SATUR|OBJ_TRUNC)))
            continue;
/*-------- Reset mean and count */
          for (d=0; d<naxis; d++)
            {
            mean[d]  = 0.0;
            for (d2=0; d2<naxis; d2++)
              covar[d][d2] = 0.0;
            }
          count = 0;
          for (; samp2 && samp2->set->field->astromlabel < 0;
            samp2=samp2->prevsamp)
            {
            for (d=0; d<naxis; d++)
              {
              dx[d] = (samp2->projpos[d] - samp->projpos[d]);
              mean[d] += dx[d];
              for (d2=0; d2<=d; d2++)
                covar[d][d2] += dx[d]*dx[d2];
              }
            count++;
            }
          if (count)
            {
            fx = (int)((samp->rawpos[lng]-0.5)/xstep);
            if (fx<0)
              fx = 0;
            else if (fx >= CPLOT_ASTNSUBPLOTS)
              fx = CPLOT_ASTNSUBPLOTS - 1;
            fy = (int)((samp->rawpos[lat]-0.5)/ystep);
            if (fy<0)
              fy = 0;
            else if (fy >= CPLOT_ASTNSUBPLOTS)
              fy = CPLOT_ASTNSUBPLOTS - 1;
            d3 = fy*CPLOT_ASTNSUBPLOTS + fx;
            for (d=0; d<naxis; d++)
              {
              vecshift[d3][d] += mean[d]/count;
              for (d2=0; d2<=d; d2++)
                vecvar[d3][d][d2] += (covar[d][d2]-mean[d]*mean[d2])
				/(count>1?(count-1):1);
              }
            nvec[d3]++;
            }
          }
        }
/*---- Compute coordinates of vectors */
      for (d=0; d<naxis; d++)
        rawpos[d] = set0->wcs->naxisn[d]/2.0;
      for (fy=0; fy<CPLOT_ASTNSUBPLOTS; fy++)
        for (fx=0; fx<CPLOT_ASTNSUBPLOTS; fx++)
          {
          d3 = fy*CPLOT_ASTNSUBPLOTS + fx;
          rawpos[lng] = xoff + xstep*fx;
          rawpos[lat] = yoff + ystep*fy;
/*-------- Compute application point coordinates in the display system */
          raw_to_wcs(set0->wcs, rawpos, wcspos);
          wcs_to_raw(wcs0, wcspos, vecpos);
/*-------- Compute application point coordinates in the group system */
          wcs_to_raw(fgroup->wcs, wcspos, rawpos);
/*-------- Add shift */
          if (nvec[d3])
            for (d=0; d<naxis; d++)
              rawpos[d] += vecshift[d3][d]/nvec[d3];
/*-------- Compute vector extremity coordinates in the display system */
          raw_to_wcs(fgroup->wcs, rawpos, wcspos);
          wcs_to_raw(wcs0, wcspos, vecpos2);
          x[0] = vecpos[lng];
          x[1] = vecpos2[lng] + minscale*(vecpos2[lng] - vecpos[lng]);
          y[0] = vecpos[lat];
          y[1] = vecpos2[lat] + minscale*(vecpos2[lat] - vecpos[lat]);
          plwid(lwid*2);
          plcol(1);
          plline(2, x,y);
          plcol(15);
          plpoin(1,x,y,1);
          plwid(lwid);
          }
      }
    }

  plend();

  end_wcs(wcs0);

  for (i=0; i<naxis; i++)
    free(ctype[i]);

  cplot_astrefsysmap(fgroups, ngroup, instru, hsn_thresh);/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_chi2 *********************************************************
PROTO	int cplot_chi2(fgroupstruct *fgroup)
PURPOSE	Plot internal and reference chi2/d.o.f for each field.
INPUT	Pointer to the field group.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	astrstats_fgroup() must have been run on group first.
AUTHOR	E. Bertin (IAP)
VERSION	12/06/2005
 ***/
int	cplot_chi2(fgroupstruct *fgroup)
  {
   fieldstruct	**fields;
   PLFLT	xl[2], yl[2],
		x,y, lim, ymin,ymin2, ymax,ymax2;
   PLINT	lwid;
   char		xlabel[80], ylabel[80], str[80];
   int		f,n, nx,ny, nfield;

  nx = 1;
  ny = 2;

  if (cplot_init(nx, ny , CPLOT_CHI2) == RETURN_ERROR)
    {
    cplot_end(CPLOT_CHI2);
    return RETURN_OK;
    }

  fields = fgroup->field;
  nfield = fgroup->nfield;

/* Compute chi2 range for plot limits */
  ymax = ymax2 = 1.0;
  for (f=0; f<nfield; f++)
    {
    lim = fields[f]->chi2_int;
    if (lim>ymax)
      ymax = lim;
    lim = fields[f]->chi2_ref;
    if (lim>ymax2)
      ymax2 = lim;
    }
  ymin = -ymax*0.05;
  ymax *= 1.05;
  ymin2 = -ymax2*0.05;
  ymax2 *= 1.05;

  yl[0] = yl[1] = 0.0;
  xl[0] = 0.0;
  xl[1] = nfield+1.0;

/*-- Now plot! */
/* Internal chi2/d.o.f. first */
  plcol(15);
  plschr(0.0,0.5);
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plenv(0.0, nfield + 1.0, ymin, ymax, 0, 0);
  sprintf(xlabel, "Field ##");
  sprintf(ylabel, "#gx#u2#d / d.o.f.");
  sprintf(str, "Group ##%d: chi2 / d.o.f. (astrometry, internal)", fgroup->no);
  pllab(xlabel, ylabel, str);
  plssym(0.0,1.0);
  plschr(0.0,0.5);
  n = 0;
  for (f=0; f<nfield; f++)
    {
    x = (PLFLT)++n;
    y = fields[f]->chi2_int;
    plpoin((PLINT)1, &x,&y, 5);
    plptex(x,y, 0.0, -1.0, -0.1, fields[f]->rfilename);
    }
  pllsty(2);
  plline(2, xl, yl);
  pllsty(1);
/* Reference chi2/d.o.f. */
  plschr(0.0,0.5);
  plwid(lwid);
  plenv(0.0, nfield + 1.0, ymin2, ymax2, 0, 0);
  sprintf(xlabel, "Field ##");
  sprintf(ylabel, "#gx#u2#d / d.o.f.");
  sprintf(str, "Group ##%d: chi2 / d.o.f. (astrometry, reference)", fgroup->no);
  pllab(xlabel, ylabel, str);
  plssym(0.0,1.0);
  plschr(0.0,0.5);
  n = 0;
  for (f=0; f<nfield; f++)
    {
    x = (PLFLT)++n;
    y = fields[f]->chi2_ref;
    plpoin((PLINT)1, &x,&y, 5);
    plptex(x,y, 0.0, -1.0, -0.1, fields[f]->rfilename);
    }
  pllsty(2);
  plline(2, xl, yl);
  pllsty(1);

  plend();

  cplot_chi2(fgroup);		/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_shear *********************************************************
PROTO	int cplot_shear(fgroupstruct **fgroups, int ngroup, int instru)
PURPOSE	Plot internal and reference chi2/d.o.f for each field.
INPUT	Pointer to the array of field groups,
	Number of groups,
	Astrometric instrument index.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	astrstats_fgroup() must have been run on group first.
AUTHOR	E. Bertin (IAP)
VERSION	11/03/2006
 ***/
int	cplot_shear(fgroupstruct **fgroups, int ngroup, int instru)
  {
   fgroupstruct	*fgroup;
   fieldstruct	*field;
   PLFLT	xl[200], yl[200],
		x,y, shear,shearmin,shearmax,airmassmin,airmassmax;
   PLINT	lwid;
   double	cste, temp;
   char		xlabel[80], ylabel[80], str[80];
   int		f,g,i,j, nx,ny, nfield;

  nx = 1;
  ny = 1;

  if (cplot_init(nx, ny , CPLOT_SHEAR) == RETURN_ERROR)
    {
    cplot_end(CPLOT_SHEAR);
    return RETURN_OK;
    }

  nfield = 0;

/* First, find limits in shear and airmass */
  shearmax = airmassmax = 0.0;
  shearmin = airmassmin = 100.0;
  for (g=0; g<ngroup; g++)
    {
    fgroup = fgroups[g];
    for (f=0; f<fgroup->nfield; f++)
      {
      field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
      if (field->astromlabel == instru)
        {
        nfield++;
        if (field->match_shear < shearmin)
          shearmin = field->match_shear;
        if (field->match_shear > shearmax)
          shearmax = field->match_shear;
        if (field->airmass < airmassmin)
          airmassmin = field->airmass;
        if (field->airmass > airmassmax)
          airmassmax = field->airmass;
        }
      }
    }
  if (airmassmin>=airmassmax)
    {
    airmassmax += 0.1;
    airmassmin = airmassmax - 0.2;
    }

  if (shearmin>=shearmax)
    {
    shearmax += 0.001;
    shearmin = shearmax - 0.002;
    if (shearmin<0.0)
      shearmin = 0.0;
    }

  airmassmin -= (airmassmax-airmassmin)*(0.05+0.5/nfield);
  airmassmax += (airmassmax-airmassmin)*(0.15+0.5/nfield)/(1.05+0.5/nfield);
  shearmin -= (shearmax-shearmin)*0.1;
  shearmax += (shearmax-shearmin)*0.05/1.1;
  shearmin *= 100.0;	/* in % */
  shearmax *= 100.0;	/* in % */
  shear = (airmassmin<0.0) ? 0.0 : airmassmin*airmassmin*2.790e-2/2.0;
  if (shearmin>shear)
    shearmin = shear;
  shear = airmassmax*airmassmax*2.790e-2;
  if (shearmax<shear)
    shearmax = shear;

/*-- Now plot! */
/* Internal chi2/d.o.f. first */
  plcol(15);
  plschr(0.0,0.5);
  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plenv(airmassmin, airmassmax, shearmin, shearmax, 0, 0);
  sprintf(xlabel, "Airmass");
  sprintf(ylabel, "Contraction (%%)");
  sprintf(str, "Instrument A%d: image flattening as a fonction of airmass",
	instru+1);
  pllab(xlabel, ylabel, str);
  plschr(0.0,0.33);
  for (j=0; j<9; j++)
    {
    temp = 293.0-5.0*j;
    cste = 2.790e-2*(288.0/temp)*exp(-68.0*j/(temp+288.0)); /* in % */
    for (i=0; i<200; i++)
      {
      xl[i] = airmassmin+i/199.0*(airmassmax-airmassmin);
      yl[i] = (xl[i]>1.0? xl[i]*xl[i]*cste : cste);
      }
    plline(200, xl,yl);
    sprintf(str, "altitude %d km", j);
    plptex(xl[195],yl[195]-(shearmax-shearmin)/80.0,
	xl[195]-xl[194], yl[195]-yl[194], 1.0, str);
    }
  plssym(0.0,1.0);
  for (g=0; g<ngroup; g++)
    {
    fgroup = fgroups[g];
    for (f=0; f<fgroup->nfield; f++)
      {
      field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
      if (field->astromlabel == instru)
        {
        x = (PLFLT)field->airmass;
        y = (PLFLT)field->match_shear*100.0;
        plpoin((PLINT)1, &x,&y, 5);
        plptex(x,y, 0.0, -1.0, -0.1, field->rfilename);
        }
      }
    }
  pllsty(1);

  plend();

  cplot_shear(fgroups, ngroup, instru);		/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_aderrhisto1d ****************************************************
PROTO	int cplot_aderrhisto1d(fgroupstruct *fgroup, double hsn_thresh)
PURPOSE	Plot an astrometric difference histogram between star pairs along 1
	dimension.
INPUT	Pointer to the field group,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_aderrhisto1d(fgroupstruct *fgroup, double hsn_thresh)
  {
   fieldstruct	*field;
   wcsstruct	*wcs;
   setstruct	*set;
   samplestruct	*samp, *samp2;
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		pixscale[NAXIS], xscale[NAXIS],xscale_hsn[NAXIS],
		xoffset[NAXIS],xoffset_hsn[NAXIS],
		yscale,yscale_hsn, yoffset,yoffset_hsn, boffset,bscale, dy,
		maxwidth,margin;
   PLFLT	**histo[NAXIS*NAXIS],**histo_hsn[NAXIS*NAXIS],
		*cuty[NAXIS*NAXIS],*cuty_hsn[NAXIS*NAXIS],
		*clevel,*cutbin,
		xl[2], yl[2],zmax[NAXIS*NAXIS],zmax_hsn[NAXIS*NAXIS],
		cutymax[NAXIS*NAXIS], cutymax_hsn[NAXIS*NAXIS],
		r[2],g[2],b[2],cpoint[2],
		lim,maxlim, cy, z;
   PLINT	lwid;
   char		xlabel[80], ylabel[80], str[80];
   int		d,d2,d3, f,i,s,n, naxis, nsamp, firstflag, ix,iy;

  if (cplot_init(1,fgroup->naxis*fgroup->naxis, CPLOT_ADERROR1D)==RETURN_ERROR)
    {
    cplot_end(CPLOT_ADERROR1D);
    return RETURN_OK;
    }

  wcs = fgroup->wcs;
  if (!wcs)
    return RETURN_ERROR;
  naxis = fgroup->naxis;

  for (d=0; d<naxis; d++)
    rawpos2[d] = rawpos[d] = wcs->naxisn[d]/2.0;
  raw_to_wcs(wcs, rawpos, wcspos);

  QMALLOC(cutbin, PLFLT, CPLOT_NADERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);
  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  for (d=0; d<naxis; d++)
    {
    rawpos2[d] += 1.0;
    raw_to_wcs(wcs, rawpos2, wcspos2);
    pixscale[d] = wcs_dist(wcs, wcspos, wcspos2)*DEG/ARCSEC;	/* in arcsec */
    rawpos2[d] -= 1.0;
    for (d2=0; d2<naxis; d2++)
      {
      d3 = d2*naxis+d;
      plAlloc2dGrid(&histo[d3], CPLOT_ADERR1DNX, CPLOT_ADERR1DNY);
      plAlloc2dGrid(&histo_hsn[d3], CPLOT_ADERR1DNX_HSN, CPLOT_ADERR1DNY_HSN);
      QCALLOC(cuty[d3], PLFLT, CPLOT_NADERRHISTBIN);
      QCALLOC(cuty_hsn[d3], PLFLT, CPLOT_NADERRHISTBIN);
      cutymax[d3] = cutymax_hsn[d3] = zmax[d3] = zmax_hsn[d3] = 0.0;
      }
    xoffset[d] = xoffset_hsn[d] = fgroup->projposmin[d];
    xscale[d] = CPLOT_ADERR1DNX/(fgroup->projposmax[d]-fgroup->projposmin[d]);
    xscale_hsn[d] = CPLOT_ADERR1DNX_HSN
	/ (fgroup->projposmax[d] - fgroup->projposmin[d]);
    }

  maxlim = 0.0;
  for (d2=0; d2<naxis; d2++)
    if ((lim=fgroup->sig_interr[d2]*DEG/ARCSEC*3.0) > maxlim)
      maxlim = lim;
  if (maxlim<=0.0)
    maxlim = 1.0;
  boffset = -maxlim;
  bscale = CPLOT_NADERRHISTBIN / (2.0*maxlim);
  yoffset = yoffset_hsn = -maxlim;
  yscale = CPLOT_ADERR1DNY/(2.0*maxlim);
  yscale_hsn = CPLOT_ADERR1DNY_HSN/(2.0*maxlim);

  for (f=0; f<fgroup->nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      samp = set->sample;
      for (n=nsamp; n--; samp++)
        if (!samp->nextsamp && samp->prevsamp
		&& !(samp->flags & (OBJ_SATUR|OBJ_TRUNC)))
	  {
          samp2 = samp;
          while ((samp2=samp2->prevsamp)
		&& samp2->set->field->astromlabel>=0)
            {
            if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
              continue;
            for (d2=0; d2<naxis; d2++)
              {
              ix = (int)((samp2->projpos[d2]-xoffset[d2])*xscale[d2]);
              for (d=0; d<naxis; d++)
                {
                d3 = d2*naxis+d;
                dy = (samp2->projpos[d]-samp->projpos[d])*pixscale[d];
                iy = (int)((dy -yoffset)*yscale);
                if (ix>=0 && ix<CPLOT_ADERR1DNX
			&& iy>=0 && iy<CPLOT_ADERR1DNY)
                  {
                  z = (histo[d3][ix][iy] += 1.0);
                  if (z>zmax[d3])
                    zmax[d3] = z;
                  }
                iy = (int)((dy - boffset)*bscale);
                if (iy>=0 && iy<CPLOT_NADERRHISTBIN)
                  {
                  cy = (cuty[d3][iy] += 1.0);
                  if (cy>cutymax[d3])
                  cutymax[d3] = cy;
                  }
                }
	      }
            if (samp2->flux/samp2->fluxerr>=hsn_thresh)
              {
              for (d2=0; d2<naxis; d2++)
                {
                ix = (int)((samp2->projpos[d2]-xoffset_hsn[d2])*xscale_hsn[d2]);
                for (d=0; d<naxis; d++)
                  {
                  d3 = d2*naxis+d;
                  dy = (samp2->projpos[d]-samp->projpos[d])*pixscale[d];
                  iy = (int)((dy -yoffset_hsn)*yscale_hsn);
                  if (ix>=0 && ix<CPLOT_ADERR1DNX_HSN
			&& iy>=0 && iy<CPLOT_ADERR1DNY_HSN)
                    {
                    z = (histo_hsn[d3][ix][iy] += 1.0);
                    if (z>zmax_hsn[d3])
                      zmax_hsn[d3] = z;
                    }
                  iy = (int)((dy - boffset)*bscale);
                  if (iy>=0 && iy<CPLOT_NADERRHISTBIN)
                    {
                    cy = (cuty_hsn[d3][iy] += 1.0);
                    if (cy>cutymax_hsn[d3])
                    cutymax_hsn[d3] = cy;
                    }
                  }
                }
	      }
            }
          }
      }
    }

/* Now plot! */

  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  for (i=0; i<CPLOT_NADERRHISTBIN; i++)
    cutbin[i] = boffset+(i+0.5)/bscale;

  firstflag = 1;
  plschr(0.0,0.67);
  for (d2=0; d2<fgroup->naxis; d2++)
    for (d=0; d<fgroup->naxis; d++)
      {
      d3 = d2*naxis+d;
      maxwidth = fgroup->projposmax[d2]-fgroup->projposmin[d2];
      margin = 0.1*maxwidth;
/* Adjust histogram to fit in the displayed box */
      for (i=0; i<CPLOT_NADERRHISTBIN; i++)
        {
        cuty[d3][i] = fgroup->projposmin[d2] - margin
			+ cuty[d3][i]/cutymax[d3] * 0.9*margin;
        cuty_hsn[d3][i] = fgroup->projposmin[d2] - margin
			+ cuty_hsn[d3][i]/cutymax_hsn[d3] * 0.9*margin;
        }
      plwid(lwid);
      plenv(fgroup->projposmin[d2]-margin, fgroup->projposmax[d2],
		-maxlim, maxlim, 0, -2);
/* Use a non-linear shade level distribution */
      if (zmax[d3]>=1.0)
        {
        for (i=0; i<CPLOT_NSHADES; i++)
          clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax[d3]+0.5;
        r[0] = 0.96; g[0] = 1.0; b[0] = 0.96;
        r[1] = 0.2; g[1] = 0.3; b[1] = 0.2;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plshades(histo[d3], CPLOT_ADERR1DNX, CPLOT_ADERR1DNY, NULL,
		fgroup->projposmin[d2], fgroup->projposmax[d2], -maxlim, maxlim,
		clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
        }
      else
        {
        plcol(1);
        plptex((fgroup->projposmin[d2] - margin + fgroup->projposmax[d2])/2.0,
		maxlim/2.0, 1.0, 0.0, 0.5, "No overlapping detections!");
        }
      if (zmax_hsn[d3]>=1.0)
        {
        r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
        r[1] = 0.7; g[1] = 0.7; b[1] = 0.7;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plimage(histo_hsn[d3], CPLOT_ADERR1DNX_HSN, CPLOT_ADERR1DNY_HSN,
		fgroup->projposmin[d2], fgroup->projposmax[d2], -maxlim, maxlim,
		0.5, zmax_hsn[d3],
		fgroup->projposmin[d2], fgroup->projposmax[d2], -maxlim, maxlim);
        }
      sprintf(xlabel, "AXIS%d [pixels]", d2+1);
      sprintf(ylabel, "#gDAXIS%d [\"]", d+1);
      plscolbg(255,255,255);	/* Force the background colour to white */
      plscol0(15, 0,0,0);	/* Force the foreground colour to black */
/* 1D histograms */
      plcol(3);
      plwid(2*lwid);
      plline(CPLOT_NADERRHISTBIN, cuty[d3], cutbin);
      plcol(7);
      plline(CPLOT_NADERRHISTBIN, cuty_hsn[d3], cutbin);
      plwid(lwid);
      plcol(15);
      plwid(lwid);
      xl[0] = fgroup->projposmin[d2] - margin;
      xl[1] = fgroup->projposmax[d2];
      yl[0] = yl[1] = 0.0;
      pllsty(2);
      plline(2, xl, yl);
      pllsty(1);
      plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
      if (firstflag)
        {
        sprintf(str, "Group ##%d: 1D internal astrometric errors", fgroup->no);
        pllab(xlabel, ylabel, str);
        }
      else
        pllab(xlabel, ylabel, "");
      firstflag = 0;
      }

/* Free array of points */
  free(clevel);
  free(cutbin);
  for (d3=0; d3<naxis*naxis; d3++)
    {
    plFree2dGrid(histo[d3], CPLOT_ADERR1DNX, CPLOT_ADERR1DNY); 
    plFree2dGrid(histo_hsn[d3], CPLOT_ADERR1DNX_HSN, CPLOT_ADERR1DNY_HSN);
    free(cuty[d3]);
    free(cuty_hsn[d3]);
    }

  plend();

  cplot_aderrhisto1d(fgroup, hsn_thresh);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_aderrhisto2d ****************************************************
PROTO	int cplot_aderrhisto2d(fgroupstruct *fgroup, double hsn_thresh)
PURPOSE	Plot astrometric difference between star pairs as a 2D histogram.
INPUT	Pointer to the field group,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_aderrhisto2d(fgroupstruct *fgroup, double hsn_thresh)
  {
   fieldstruct	*field;
   wcsstruct	*wcs;
   setstruct	*set;
   samplestruct	*samp, *samp2;
   char		str[64];
   double	offset,offset_hsn, scale,scale_hsn, boffset,bscale,
		cutxmax,cutxmax_hsn, cutymax,cutymax_hsn, cx,cy, dx,dy;
   PLFLT	**histo,**histo_hsn,
		rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		pixscale[NAXIS],
		xl[5], yl[5],r[2],g[2],b[2],cpoint[2],
		*clevel,*cutbin,*cutx,*cutx_hsn,*cuty,*cuty_hsn,
		lim,maxlim, z,zmax,zmax_hsn;
   PLINT	lwid;
   int		d,d2, f,i,s,n, ix,iy,
		nsamp;

  if (cplot_init(1,1, CPLOT_ADERROR2D) == RETURN_ERROR)
    {
    cplot_end(CPLOT_ADERROR2D);
    return RETURN_OK;
    }

  wcs = fgroup->wcs;
  if (!wcs)
    return RETURN_ERROR;
  for (d=0; d<fgroup->naxis; d++)
    rawpos2[d] = rawpos[d] = wcs->naxisn[d]/2.0;
  raw_to_wcs(wcs, rawpos, wcspos);

  for (d=0; d<fgroup->naxis; d++)
    {
    rawpos2[d] += 1.0;
    raw_to_wcs(wcs, rawpos2, wcspos2);
    pixscale[d] = wcs_dist(wcs, wcspos, wcspos2)*DEG/ARCSEC;	/* in arcsec */
    rawpos2[d] -= 1.0;
    }

  maxlim = 0.0;
  for (d2=0; d2<fgroup->naxis; d2++)
    if ((lim=fgroup->sig_interr[d2]*DEG/ARCSEC*3.0) > maxlim)
      maxlim = lim;
  if (maxlim<=0.0)
    maxlim = 1.0;
  plAlloc2dGrid(&histo, CPLOT_ADERR2DN, CPLOT_ADERR2DN);
  plAlloc2dGrid(&histo_hsn, CPLOT_ADERR2DN_HSN, CPLOT_ADERR2DN_HSN);
  offset = offset_hsn = -maxlim;
  scale = CPLOT_ADERR2DN / (2.0*maxlim);
  scale_hsn = CPLOT_ADERR2DN_HSN / (2.0*maxlim);

  boffset = -maxlim;
  bscale = CPLOT_NADERRHISTBIN / (2.0*maxlim);

  QMALLOC(cutbin, PLFLT, CPLOT_NADERRHISTBIN);
  QCALLOC(cutx, PLFLT, CPLOT_NADERRHISTBIN);
  QCALLOC(cuty, PLFLT, CPLOT_NADERRHISTBIN);
  QCALLOC(cutx_hsn, PLFLT, CPLOT_NADERRHISTBIN);
  QCALLOC(cuty_hsn, PLFLT, CPLOT_NADERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);

  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  cutxmax = cutxmax_hsn = cutymax = cutymax_hsn = 0.0;;
  zmax = zmax_hsn = 0.0;
  for (f=0; f<fgroup->nfield; f++)
    {
    field = fgroup->field[f];
    for (s=0; s<field->nset; s++)
      {
      set = field->set[s];
      nsamp = set->nsample;
      samp = set->sample;
      for (n=nsamp; n--; samp++)
        if (!samp->nextsamp && samp->prevsamp
		&& !(samp->flags & (OBJ_SATUR|OBJ_TRUNC)))
          {
          samp2 = samp;
          while ((samp2=samp2->prevsamp)
		&& samp2->set->field->astromlabel>=0)
              {
              if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
                continue;
              dx = (samp2->projpos[0]-samp->projpos[0])*pixscale[0];
              dy = (samp2->projpos[1]-samp->projpos[1])*pixscale[1];
              ix = (int)((dx - offset)*scale);
              iy = (int)((dy - offset)*scale);
              if (ix>=0 && ix<CPLOT_ADERR2DN && iy>=0 && iy<CPLOT_ADERR2DN)
                {
                z = (histo[ix][iy] += 1.0);
                if (z>zmax)
                  zmax = z;
                }
              ix = (int)((dx - boffset)*bscale);
              if (ix>=0 && ix<CPLOT_NADERRHISTBIN)
                {
                cx = (cutx[ix] += 1.0);
                if (cx>cutxmax)
                  cutxmax = cx;
                }
              iy = (int)((dy - boffset)*bscale);
              if (iy>=0 && iy<CPLOT_NADERRHISTBIN)
                {
                cy = (cuty[iy] += 1.0);
                if (cy>cutymax)
                  cutymax = cy;
                }
              if (samp2->flux/samp2->fluxerr >= hsn_thresh)
                {
                ix = (int)((dx - offset_hsn)*scale_hsn);
                iy = (int)((dy - offset_hsn)*scale_hsn);
                if (ix>=0 && ix<CPLOT_ADERR2DN_HSN
			&& iy>=0 && iy<CPLOT_ADERR2DN_HSN)
                  {
                  z = (histo_hsn[ix][iy] += 1.0);
                  if (z>zmax_hsn)
                    zmax_hsn = z;
                  }
                ix = (int)((dx - boffset)*bscale);
                if (ix>=0 && ix<CPLOT_NADERRHISTBIN)
                  {
                  cx = (cutx_hsn[ix] += 1.0);
                  if (cx>cutxmax_hsn)
                    cutxmax_hsn = cx;
                  }
                iy = (int)((dy - boffset)*bscale);
                if (iy>=0 && iy<CPLOT_NADERRHISTBIN)
                  {
                  cy = (cuty_hsn[iy] += 1.0);
                  if (cy>cutymax_hsn)
                    cutymax_hsn = cy;
                  }
                }
	      }
          }
        }
    }

/* Adjust histogram to fit in the displayed box */
  for (i=0; i<CPLOT_NADERRHISTBIN; i++)
    {
    cutbin[i] = boffset+(i+0.5)/bscale;
    cutx[i] = -maxlim + cutx[i]/cutxmax*maxlim/2.0;
    cutx_hsn[i] = -maxlim + cutx_hsn[i]/cutxmax_hsn*maxlim/2.0;
    cuty[i] = -maxlim + cuty[i]/cutymax*maxlim/2.0;
    cuty_hsn[i] = -maxlim + cuty_hsn[i]/cutymax_hsn*maxlim/2.0;
    }

  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plenv(-maxlim, maxlim, -maxlim, maxlim, 1, -2);

/* Use a non-linear shade level distribution */
  if (zmax>=1.0)
    {
    for (i=0; i<CPLOT_NSHADES; i++)
      clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax+0.5;
    r[0] = 0.96; g[0] = 1.0; b[0] = 0.96;
    r[1] = 0.2; g[1] = 0.3; b[1] = 0.2;
    plscmap1l(1, 2, cpoint, r, g, b, NULL);
    plshades(histo, CPLOT_ADERR2DN, CPLOT_ADERR2DN, NULL,
	-maxlim,maxlim, -maxlim,maxlim,
	clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
    }
  else
    {
    plcol(1);
    plptex(0.0, maxlim/2.0, 1.0, 0.0, 0.5, "No overlapping detections!");
    }

  if (zmax_hsn>=1.0)
    {
    r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
    r[1] = 0.8; g[1] = 0.8; b[1] = 0.8;
    plscmap1l(1, 2, cpoint, r, g, b, NULL);
    plimage(histo_hsn, CPLOT_ADERR2DN_HSN, CPLOT_ADERR2DN_HSN,
	-maxlim,maxlim, -maxlim, maxlim,
	0.5, zmax_hsn,
	-maxlim, maxlim, -maxlim, maxlim);
    }

  plscolbg(255,255,255);	/* Force the background colour to white */
  plscol0(15, 0,0,0);		/* Force the foreground colour to black */
  plschr(0.0,0.5);
/* Pixel footprint */
  plcol(15);
  pllsty(3);
  xl[0] = xl[1] = xl[4] = pixscale[0]/2.0;
  yl[0] = yl[3] = yl[4] = pixscale[1]/2.0;
  xl[2] = xl[3] = -xl[0];
  yl[1] = yl[2] = -yl[0];
  plline(5, xl, yl);
  pllsty(1);
/* 1D histograms */
  plcol(3);
  plwid(2*lwid);
  plline(CPLOT_NADERRHISTBIN, cutbin, cutx);
  plline(CPLOT_NADERRHISTBIN, cuty, cutbin);
  plcol(7);
  plline(CPLOT_NADERRHISTBIN, cutbin, cutx_hsn);
  plline(CPLOT_NADERRHISTBIN, cuty_hsn, cutbin);
  plwid(lwid);
  plcol(15);
  plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
  sprintf(str, "Group ##%d: 2D internal astrometric errors", fgroup->no);
  pllab( "#gDAXIS1 [\"]", "#gDAXIS2 [\"]", str);
/* reticulus */
  pllsty(2);
  xl[0] = -maxlim;
  xl[1] = maxlim;
  yl[0] = yl[1] = 0.0;
  plline(2, xl, yl);
  xl[0] = xl[1] = 0.0;
  yl[0] = -maxlim;
  yl[1] = maxlim;
  plline(2, xl, yl);
  pllsty(1);

/*-- Free array of points */
  plFree2dGrid(histo, CPLOT_ADERR2DN, CPLOT_ADERR2DN); 
  plFree2dGrid(histo_hsn, CPLOT_ADERR2DN_HSN, CPLOT_ADERR2DN_HSN); 
  free(clevel);
  free(cutbin);
  free(cutx);
  free(cuty);
  free(cutx_hsn);
  free(cuty_hsn);

  plend();

  cplot_aderrhisto2d(fgroup, hsn_thresh);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_referrhisto1d ****************************************************
PROTO	int cplot_referrhisto1d(fgroupstruct *fgroup, fieldstruct *reffield,
		double hsn_thresh)
PURPOSE	Plot an astrometric difference histogram between star pairs along 1
	dimension.
INPUT	Pointer to the field group,
	pointer to the reference field,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_referrhisto1d(fgroupstruct *fgroup, fieldstruct *reffield,
		double hsn_thresh)
  {
   wcsstruct	*wcs;
   setstruct	*set;
   samplestruct	*samp, *samp2;
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		pixscale[NAXIS], xscale[NAXIS],xscale_hsn[NAXIS],
		xoffset[NAXIS],xoffset_hsn[NAXIS],
		yscale,yscale_hsn, yoffset,yoffset_hsn, boffset,bscale, dy,
		maxwidth,margin;
   PLFLT	**histo[NAXIS*NAXIS],**histo_hsn[NAXIS*NAXIS],
		*cuty[NAXIS*NAXIS],*cuty_hsn[NAXIS*NAXIS],
		*clevel,*cutbin,
		xl[2], yl[2],zmax[NAXIS*NAXIS],zmax_hsn[NAXIS*NAXIS],
		cutymax[NAXIS*NAXIS], cutymax_hsn[NAXIS*NAXIS],
		r[2],g[2],b[2],cpoint[2],
		lim,maxlim, cy, z;
   PLINT	lwid;
   char		xlabel[80], ylabel[80], str[80];
   int		d,d2,d3, i,s,n, naxis, nsamp, firstflag, ix,iy;

  if (cplot_init(1,fgroup->naxis*fgroup->naxis, CPLOT_REFERROR1D)==RETURN_ERROR)
    {
    cplot_end(CPLOT_REFERROR1D);
    return RETURN_OK;
    }

  wcs = fgroup->wcs;
  if (!wcs)
    return RETURN_ERROR;
  naxis = fgroup->naxis;

  for (d=0; d<naxis; d++)
    rawpos2[d] = rawpos[d] = wcs->naxisn[d]/2.0;
  raw_to_wcs(wcs, rawpos, wcspos);

  QMALLOC(cutbin, PLFLT, CPLOT_NREFERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);
  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  for (d=0; d<naxis; d++)
    {
    rawpos2[d] += 1.0;
    raw_to_wcs(wcs, rawpos2, wcspos2);
    pixscale[d] = wcs_dist(wcs, wcspos, wcspos2)*DEG/ARCSEC;	/* in arcsec */
    rawpos2[d] -= 1.0;
    for (d2=0; d2<naxis; d2++)
      {
      d3 = d2*naxis+d;
      plAlloc2dGrid(&histo[d3], CPLOT_ADERR1DNX, CPLOT_ADERR1DNY);
      plAlloc2dGrid(&histo_hsn[d3], CPLOT_ADERR1DNX_HSN, CPLOT_ADERR1DNY_HSN);
      QCALLOC(cuty[d3], PLFLT, CPLOT_NREFERRHISTBIN);
      QCALLOC(cuty_hsn[d3], PLFLT, CPLOT_NREFERRHISTBIN);
      cutymax[d3] = cutymax_hsn[d3] = zmax[d3] = zmax_hsn[d3] = 0.0;
      }
    xoffset[d] = xoffset_hsn[d] = fgroup->projposmin[d];
    xscale[d] = CPLOT_ADERR1DNX/(fgroup->projposmax[d]-fgroup->projposmin[d]);
    xscale_hsn[d] = CPLOT_ADERR1DNX_HSN
	/ (fgroup->projposmax[d] - fgroup->projposmin[d]);
    }

  maxlim = 0.0;
  for (d2=0; d2<naxis; d2++)
    if ((lim=fgroup->sig_referr[d2]*DEG/ARCSEC*3.0) > maxlim)
      maxlim = lim;
  if (maxlim<=0.0)
    maxlim = 1.0;
  boffset = -maxlim;
  bscale = CPLOT_NREFERRHISTBIN / (2.0*maxlim);
  yoffset = yoffset_hsn = -maxlim;
  yscale = CPLOT_ADERR1DNY/(2.0*maxlim);
  yscale_hsn = CPLOT_ADERR1DNY_HSN/(2.0*maxlim);

  for (s=0; s<reffield->nset; s++)
    {
    set = reffield->set[s];
    nsamp = set->nsample;
    samp = set->sample;
    for (n=nsamp; n--; samp++)
      if (samp->nextsamp && !(samp->flags & (OBJ_SATUR|OBJ_TRUNC)))
        {
        samp2 = samp;
        while ((samp2=samp2->nextsamp))
          {
          if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
            continue;
          for (d2=0; d2<naxis; d2++)
            {
            ix = (int)((samp2->projpos[d2]-xoffset[d2])*xscale[d2]);
            for (d=0; d<naxis; d++)
              {
              d3 = d2*naxis+d;
              dy = (samp2->projpos[d]-samp->projpos[d])*pixscale[d];
              iy = (int)((dy -yoffset)*yscale);
              if (ix>=0 && ix<CPLOT_ADERR1DNX
			&& iy>=0 && iy<CPLOT_ADERR1DNY)
                {
                z = (histo[d3][ix][iy] += 1.0);
                if (z>zmax[d3])
                  zmax[d3] = z;
                }
              iy = (int)((dy - boffset)*bscale);
              if (iy>=0 && iy<CPLOT_NREFERRHISTBIN)
                {
                cy = (cuty[d3][iy] += 1.0);
                if (cy>cutymax[d3])
                cutymax[d3] = cy;
                }
              }
	    }
          if (samp2->flux/samp2->fluxerr>=hsn_thresh)
            {
            for (d2=0; d2<naxis; d2++)
              {
              ix = (int)((samp2->projpos[d2]-xoffset_hsn[d2])*xscale_hsn[d2]);
              for (d=0; d<naxis; d++)
                {
                d3 = d2*naxis+d;
                dy = (samp2->projpos[d]-samp->projpos[d])*pixscale[d];
                iy = (int)((dy -yoffset_hsn)*yscale_hsn);
                if (ix>=0 && ix<CPLOT_ADERR1DNX_HSN
			&& iy>=0 && iy<CPLOT_ADERR1DNY_HSN)
                  {
                  z = (histo_hsn[d3][ix][iy] += 1.0);
                  if (z>zmax_hsn[d3])
                    zmax_hsn[d3] = z;
                  }
                iy = (int)((dy - boffset)*bscale);
                if (iy>=0 && iy<CPLOT_NREFERRHISTBIN)
                  {
                  cy = (cuty_hsn[d3][iy] += 1.0);
                  if (cy>cutymax_hsn[d3])
                    cutymax_hsn[d3] = cy;
                  }
                }
              }
            }
          }
        }
    }

/* Now plot! */

  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  for (i=0; i<CPLOT_NREFERRHISTBIN; i++)
    cutbin[i] = boffset+(i+0.5)/bscale;

  firstflag = 1;
  plschr(0.0,0.67);
  for (d2=0; d2<fgroup->naxis; d2++)
    for (d=0; d<fgroup->naxis; d++)
      {
      d3 = d2*naxis+d;
      maxwidth = fgroup->projposmax[d2]-fgroup->projposmin[d2];
      margin = 0.1*maxwidth;
/*---- Adjust histogram to fit in the displayed box */
      for (i=0; i<CPLOT_NREFERRHISTBIN; i++)
        {
        cuty[d3][i] = fgroup->projposmin[d2] - margin
			+ cuty[d3][i]/cutymax[d3] * 0.9*margin;
        cuty_hsn[d3][i] = fgroup->projposmin[d2] - margin
			+ cuty_hsn[d3][i]/cutymax_hsn[d3] * 0.9*margin;
        }
      plwid(lwid);
      plenv(fgroup->projposmin[d2] - margin, fgroup->projposmax[d2],
		-maxlim, maxlim, 0, -2);
/*---- Use a non-linear shade level distribution */
      if (zmax[d3]>=1.0)
        {
        for (i=0; i<CPLOT_NSHADES; i++)
          clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax[d3]+0.5;
        r[0] = 1.0; g[0] = 0.98; b[0] = 0.98;
        r[1] = 0.6; g[1] = 0.1; b[1] = 0.1;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plshades(histo[d3], CPLOT_ADERR1DNX, CPLOT_ADERR1DNY, NULL,
		fgroup->projposmin[d2], fgroup->projposmax[d2], -maxlim, maxlim,
		clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
        }
      else
        {
        plcol(1);
        plptex((fgroup->projposmin[d2] - margin + fgroup->projposmax[d2])/2.0,
		maxlim/2.0, 1.0, 0.0, 0.5, "No match with a reference!");
        }
      if (zmax_hsn[d3]>=1.0)
        {
        r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
        r[1] = 0.8; g[1] = 0.8; b[1] = 0.8;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plimage(histo_hsn[d3], CPLOT_ADERR1DNX_HSN, CPLOT_ADERR1DNY_HSN,
		fgroup->projposmin[d2], fgroup->projposmax[d2], -maxlim, maxlim,
		0.5, zmax_hsn[d3],
		fgroup->projposmin[d2], fgroup->projposmax[d2], -maxlim, maxlim);
        }
      plscolbg(255,255,255);	/* Force the background colour to white */
      plscol0(15, 0,0,0);	/* Force the foreground colour to black */
      sprintf(xlabel, "AXIS%d [pixels]", d2+1);
      sprintf(ylabel, "#gDAXIS%d [\"]", d+1);
/*---- 1D histograms */
      plcol(1);
      plwid(2*lwid);
      plline(CPLOT_NREFERRHISTBIN, cuty[d3], cutbin);
      plcol(7);
      plline(CPLOT_NREFERRHISTBIN, cuty_hsn[d3], cutbin);
      plwid(lwid);
      plcol(15);
      plwid(lwid);
      xl[0] = fgroup->projposmin[d2] - margin;
      xl[1] = fgroup->projposmax[d2];
      yl[0] = yl[1] = 0.0;
      pllsty(2);
      plline(2, xl, yl);
      pllsty(1);
      plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
      if (firstflag)
        {
        sprintf(str, "Group ##%d: 1D reference astrometric errors", fgroup->no);
        pllab(xlabel, ylabel, str);
        }
      else
        pllab(xlabel, ylabel, "");
      firstflag = 0;
      }

/* Free array of points */
  free(clevel);
  free(cutbin);
  for (d3=0; d3<naxis*naxis; d3++)
    {
    plFree2dGrid(histo[d3], CPLOT_ADERR1DNX, CPLOT_ADERR1DNY); 
    plFree2dGrid(histo_hsn[d3], CPLOT_ADERR1DNX_HSN, CPLOT_ADERR1DNY_HSN);
    free(cuty[d3]);
    free(cuty_hsn[d3]);
    }

  plend();

  cplot_referrhisto1d(fgroup, reffield, hsn_thresh);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_referrhisto2d ***************************************************
PROTO	int cplot_referrhisto2d(fgroupstruct *fgroup, fieldstruct *reffield,
			double hsn_thresh)
PURPOSE	Plot astrometric difference between star pairs as a 2D histogram.
INPUT	Pointer to the field group,
	pointer to the reference field,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_referrhisto2d(fgroupstruct *fgroup, fieldstruct *reffield,
			double hsn_thresh)
  {
   wcsstruct	*wcs;
   setstruct	*set;
   samplestruct	*samp, *samp2;
   char		str[64];
   double	offset,offset_hsn, scale,scale_hsn, boffset,bscale,
		cutxmax,cutxmax_hsn, cutymax,cutymax_hsn, cx,cy, dx,dy;
   PLFLT	**histo,**histo_hsn,
		rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		pixscale[NAXIS],
		xl[5], yl[5],r[2],g[2],b[2],cpoint[2],
		*clevel,*cutbin,*cutx,*cutx_hsn,*cuty,*cuty_hsn,
		lim,maxlim, z,zmax,zmax_hsn;
   PLINT	lwid;
   int		d,d2, i,s,n, ix,iy,
		nsamp;

  if (cplot_init(1,1, CPLOT_REFERROR2D) == RETURN_ERROR)
    {
    cplot_end(CPLOT_REFERROR2D);
    return RETURN_OK;
    }

  wcs = fgroup->wcs;
  if (!wcs)
    return RETURN_ERROR;
  for (d=0; d<fgroup->naxis; d++)
    rawpos2[d] = rawpos[d] = wcs->naxisn[d]/2.0;
  raw_to_wcs(wcs, rawpos, wcspos);

  for (d=0; d<fgroup->naxis; d++)
    {
    rawpos2[d] += 1.0;
    raw_to_wcs(wcs, rawpos2, wcspos2);
    pixscale[d] = wcs_dist(wcs, wcspos, wcspos2)*DEG/ARCSEC;	/* in arcsec */
    rawpos2[d] -= 1.0;
    }

  maxlim = 0.0;
  for (d2=0; d2<fgroup->naxis; d2++)
    if ((lim=fgroup->sig_referr[d2]*DEG/ARCSEC*3.0) > maxlim)
      maxlim = lim;
  if (maxlim<=0.0)
    maxlim = 1.0;
  plAlloc2dGrid(&histo, CPLOT_REFERR2DN, CPLOT_REFERR2DN);
  plAlloc2dGrid(&histo_hsn, CPLOT_REFERR2DN_HSN, CPLOT_REFERR2DN_HSN);
  offset = offset_hsn = -maxlim;
  scale = CPLOT_REFERR2DN / (2.0*maxlim);
  scale_hsn = CPLOT_REFERR2DN_HSN / (2.0*maxlim);

  boffset = -maxlim;
  bscale = CPLOT_NREFERRHISTBIN / (2.0*maxlim);

  QMALLOC(cutbin, PLFLT, CPLOT_NREFERRHISTBIN);
  QCALLOC(cutx, PLFLT, CPLOT_NREFERRHISTBIN);
  QCALLOC(cuty, PLFLT, CPLOT_NREFERRHISTBIN);
  QCALLOC(cutx_hsn, PLFLT, CPLOT_NREFERRHISTBIN);
  QCALLOC(cuty_hsn, PLFLT, CPLOT_NREFERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);

  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  cutxmax = cutxmax_hsn = cutymax = cutymax_hsn = 0.0;;
  zmax = zmax_hsn = 0.0;
  for (s=0; s<reffield->nset; s++)
    {
    set = reffield->set[s];
    nsamp = set->nsample;
    samp = set->sample;
    for (n=nsamp; n--; samp++)
      if (samp->nextsamp && !samp->prevsamp
		&& !(samp->flags & (OBJ_SATUR|OBJ_TRUNC)))
        {
        samp2 = samp;
        while ((samp2=samp2->nextsamp))
          {
          if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
            continue;
          dx = (samp2->projpos[0]-samp->projpos[0])*pixscale[0];
          dy = (samp2->projpos[1]-samp->projpos[1])*pixscale[1];
          ix = (int)((dx - offset)*scale);
          iy = (int)((dy - offset)*scale);
          if (ix>=0 && ix<CPLOT_REFERR2DN && iy>=0 && iy<CPLOT_REFERR2DN)
            {
            z = (histo[ix][iy] += 1.0);
            if (z>zmax)
              zmax = z;
            }
          ix = (int)((dx - boffset)*bscale);
          if (ix>=0 && ix<CPLOT_NREFERRHISTBIN)
            {
            cx = (cutx[ix] += 1.0);
            if (cx>cutxmax)
              cutxmax = cx;
            }
          iy = (int)((dy - boffset)*bscale);
          if (iy>=0 && iy<CPLOT_NREFERRHISTBIN)
            {
            cy = (cuty[iy] += 1.0);
            if (cy>cutymax)
              cutymax = cy;
            }
          if (samp2->flux/samp2->fluxerr>=hsn_thresh)
            {
            ix = (int)((dx - offset_hsn)*scale_hsn);
            iy = (int)((dy - offset_hsn)*scale_hsn);
            if (ix>=0 && ix<CPLOT_REFERR2DN_HSN
		&& iy>=0 && iy<CPLOT_REFERR2DN_HSN)
              {
               z = (histo_hsn[ix][iy] += 1.0);
              if (z>zmax_hsn)
                zmax_hsn = z;
              }
            ix = (int)((dx - boffset)*bscale);
            if (ix>=0 && ix<CPLOT_NREFERRHISTBIN)
              {
              cx = (cutx_hsn[ix] += 1.0);
              if (cx>cutxmax_hsn)
                cutxmax_hsn = cx;
              }
            iy = (int)((dy - boffset)*bscale);
            if (iy>=0 && iy<CPLOT_NREFERRHISTBIN)
              {
              cy = (cuty_hsn[iy] += 1.0);
              if (cy>cutymax_hsn)
                cutymax_hsn = cy;
              }
	    }
          }
        }
    }

/* Adjust histogram to fit in the displayed box */
  for (i=0; i<CPLOT_NREFERRHISTBIN; i++)
    {
    cutbin[i] = boffset+(i+0.5)/bscale;
    cutx[i] = -maxlim + cutx[i]/cutxmax*maxlim/2.0;
    cutx_hsn[i] = -maxlim + cutx_hsn[i]/cutxmax_hsn*maxlim/2.0;
    cuty[i] = -maxlim + cuty[i]/cutymax*maxlim/2.0;
    cuty_hsn[i] = -maxlim + cuty_hsn[i]/cutymax_hsn*maxlim/2.0;
    }

  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  plwid(lwid);
  plenv(-maxlim, maxlim, -maxlim, maxlim, 1, -2);

/* Use a non-linear shade level distribution */
  if (zmax>=1.0)
    {
    for (i=0; i<CPLOT_NSHADES; i++)
      clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax+0.5;
    r[0] = 1.0; g[0] = 0.98; b[0] = 0.98;
    r[1] = 0.6; g[1] = 0.1; b[1] = 0.1;
    plscmap1l(1, 2, cpoint, r, g, b, NULL);
    plshades(histo, CPLOT_REFERR2DN, CPLOT_REFERR2DN, NULL,
	-maxlim,maxlim, -maxlim,maxlim,
	clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
    }
  else
    {
    plcol(1);
    plptex(0.0, maxlim/2.0, 1.0, 0.0, 0.5, "No match with a reference!");
    }
  if (zmax_hsn>=1.0)
    {
    r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
    r[1] = 0.8; g[1] = 0.8; b[1] = 0.8;
    plscmap1l(1, 2, cpoint, r, g, b, NULL);
    plimage(histo_hsn, CPLOT_REFERR2DN_HSN, CPLOT_REFERR2DN_HSN,
	-maxlim,maxlim, -maxlim, maxlim,
	0.5, zmax_hsn,
	-maxlim, maxlim, -maxlim, maxlim);
    }
  plscolbg(255,255,255);	/* Force the background colour to white */
  plscol0(15, 0,0,0);		/* Force the foreground colour to black */
  plschr(0.0,0.5);
/* Pixel footprint */
  plcol(15);
  pllsty(3);
  xl[0] = xl[1] = xl[4] = pixscale[0]/2.0;
  yl[0] = yl[3] = yl[4] = pixscale[1]/2.0;
  xl[2] = xl[3] = -xl[0];
  yl[1] = yl[2] = -yl[0];
  plline(5, xl, yl);
  pllsty(1);
/* 1D histograms */
  plcol(1);
  plwid(2*lwid);
  plline(CPLOT_NREFERRHISTBIN, cutbin, cutx);
  plline(CPLOT_NREFERRHISTBIN, cuty, cutbin);
  plcol(7);
  plline(CPLOT_NREFERRHISTBIN, cutbin, cutx_hsn);
  plline(CPLOT_NREFERRHISTBIN, cuty_hsn, cutbin);
  plwid(lwid);
  plcol(15);
  plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
  sprintf(str, "Group ##%d: 2D reference astrometric errors", fgroup->no);
  pllab( "#gDAXIS1 [\"]", "#gDAXIS2 [\"]", str);
/* reticulus */
  pllsty(2);
  xl[0] = -maxlim;
  xl[1] = maxlim;
  yl[0] = yl[1] = 0.0;
  plline(2, xl, yl);
  xl[0] = xl[1] = 0.0;
  yl[0] = -maxlim;
  yl[1] = maxlim;
  plline(2, xl, yl);
  pllsty(1);

/*-- Free array of points */
  plFree2dGrid(histo, CPLOT_REFERR2DN, CPLOT_REFERR2DN); 
  plFree2dGrid(histo_hsn, CPLOT_REFERR2DN_HSN, CPLOT_REFERR2DN_HSN); 
  free(clevel);
  free(cutbin);
  free(cutx);
  free(cuty);
  free(cutx_hsn);
  free(cuty_hsn);

  plend();

  cplot_referrhisto2d(fgroup, reffield, hsn_thresh);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_pixerrhisto1d ***************************************************
PROTO	int cplot_pixerrhisto1d(fgroupstruct **fgroups, int ngroup,
		int instru, double hsn_thresh)
PURPOSE	Plot an astrometric difference histogram between star pairs along 1
	dimension as a function of image pixel coordinates for a given
	astrometric intrument.
INPUT	Pointer to an array of field group pointers,
	Number of field groups,
	Astrometric intrument index,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_pixerrhisto1d(fgroupstruct **fgroups, int ngroup, int instru,
		double hsn_thresh)
  {
   fgroupstruct	*fgroup;
   fieldstruct	*field;
   wcsstruct	*wcs;
   setstruct	*set;
   samplestruct	*samp, *samp2;
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		pixscale[NAXIS], xscale[NAXIS],xscale_hsn[NAXIS],
		xoffset[NAXIS],xoffset_hsn[NAXIS], mean[NAXIS],
		yscale,yscale_hsn, yoffset,yoffset_hsn, boffset,bscale, dy,
		margin, sig2;
   PLFLT	**histo[NAXIS*NAXIS],**histo_hsn[NAXIS*NAXIS],
		*cuty[NAXIS*NAXIS],*cuty_hsn[NAXIS*NAXIS],
		*line[NAXIS*NAXIS],*weight[NAXIS*NAXIS],
		*clevel,*cutbin,*cutx,
		xl[2], yl[2],zmax[NAXIS*NAXIS],zmax_hsn[NAXIS*NAXIS],
		cutymax[NAXIS*NAXIS], cutymax_hsn[NAXIS*NAXIS],
		r[2],g[2],b[2],cpoint[2],
		lim,maxlim,maxwidth, cy, z;
   PLINT	lwid;
   char		xlabel[80], ylabel[80], str[80];
   int		d,d2,d3, f,fg,i,s,n, naxis, nsamp, firstflag, ix,iy, count;

  if (cplot_init(1,fgroups[0]->naxis*fgroups[0]->naxis,
	CPLOT_PIXERROR1D)==RETURN_ERROR)
    {
    cplot_end(CPLOT_PIXERROR1D);
    return RETURN_OK;
    }

  naxis = fgroups[0]->naxis;

  QMALLOC(cutbin, PLFLT, CPLOT_NPIXERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);
  QMALLOC(cutx, PLFLT, CPLOT_PIXERR1DNX+1);
  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  maxlim = maxwidth = 0.0;
  for (fg=0; fg<ngroup; fg++)
    {
    for (d2=0; d2<naxis; d2++)
      if ((lim=fgroups[fg]->sig_interr_hsn[d2]/fgroups[fg]->meanwcsscale[d2])
		> maxlim)
        maxlim = lim;
    if (maxlim<=0.0)
      maxlim = 1.0;
    fgroup = fgroups[fg];
    for (f=0; f<fgroup->nfield; f++)
      {
      field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
      if (field->astromlabel != instru)
        continue;
      for (s=0; s<field->nset; s++)
        {
        set = field->set[s];
        wcs = set->wcs;
        if (!wcs)
          continue;
        for (d2=0; d2<set->naxis; d2++)
          if (wcs->naxisn[d2]>maxwidth)
            maxwidth = (PLFLT)wcs->naxisn[d2];
        }
      }
    }

  maxwidth += 0.5;	/* This is the extreme limit of the frame */

  boffset = -maxlim;
  bscale = CPLOT_NPIXERRHISTBIN / (2.0*maxlim);
  yoffset = yoffset_hsn = -maxlim;
  yscale = CPLOT_PIXERR1DNY/(2.0*maxlim);
  yscale_hsn = CPLOT_PIXERR1DNY_HSN/(2.0*maxlim);

  for (d=0; d<naxis; d++)
    {
    for (d2=0; d2<naxis; d2++)
      {
      d3 = d2*naxis+d;
      plAlloc2dGrid(&histo[d3], CPLOT_PIXERR1DNX, CPLOT_PIXERR1DNY);
      plAlloc2dGrid(&histo_hsn[d3], CPLOT_PIXERR1DNX_HSN, CPLOT_PIXERR1DNY_HSN);
      QCALLOC(cuty[d3], PLFLT, CPLOT_NPIXERRHISTBIN);
      QCALLOC(cuty_hsn[d3], PLFLT, CPLOT_NPIXERRHISTBIN);
      QCALLOC(line[d3], PLFLT, CPLOT_PIXERR1DNX+1);
      QCALLOC(weight[d3], PLFLT, CPLOT_PIXERR1DNX);
      cutymax[d3] = cutymax_hsn[d3] = zmax[d3] = zmax_hsn[d3] = 0.0;
      }
    xoffset[d] = xoffset_hsn[d] = 0.5;
    xscale[d] = CPLOT_PIXERR1DNX/maxwidth;
    xscale_hsn[d] = CPLOT_PIXERR1DNX_HSN / maxwidth;
    }

  for (fg=0; fg<ngroup; fg++)
    {
    fgroup = fgroups[fg];
    wcs = fgroup->wcs;
    if (!wcs)
      return RETURN_ERROR;
    for (d=0; d<naxis; d++)
      rawpos2[d] = rawpos[d] = wcs->naxisn[d]/2.0;
    raw_to_wcs(wcs, rawpos, wcspos);
    for (d=0; d<naxis; d++)
      {
      rawpos2[d] += 1.0;
      raw_to_wcs(wcs, rawpos2, wcspos2);
      pixscale[d] = wcs_dist(wcs, wcspos, wcspos2)*DEG/ARCSEC;	/* in arcsec */
      rawpos2[d] -= 1.0;
      }
    for (f=0; f<fgroup->nfield; f++)
      {
      field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
      if (field->astromlabel != instru)
        continue;
      for (s=0; s<field->nset; s++)
        {
        set = field->set[s];
        nsamp = set->nsample;
        samp = set->sample;
        for (n=nsamp; n--; samp++)
          {
          if (samp->flags & (OBJ_SATUR|OBJ_TRUNC))
            continue;
/*-------- Reset mean and count */
          for (d2=0; d2<naxis; d2++)
            mean[d2] = 0.0;
          count = 0;
/*-------- Explore forward */
          samp2 = samp;
          while ((samp2=samp2->nextsamp))
            {
            if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
              continue;
            for (d2=0; d2<naxis; d2++)
              mean[d2] += samp2->projpos[d2] - samp->projpos[d2];
            count++;
            }
/*-------- Explore backward */
          samp2 = samp;
          while ((samp2=samp2->prevsamp) && samp2->set->field->astromlabel>=0)
            {
            if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
              continue;
            for (d2=0; d2<naxis; d2++)
              mean[d2] += samp2->projpos[d2] - samp->projpos[d2];
            count++;
            }
          if (count)
            {
            for (d2=0; d2<naxis; d2++)
              mean[d2] = mean[d2]/count + samp->projpos[d2];
/*---------- Convert to raw local coordinates */
            raw_to_wcs(wcs, mean, wcspos);
            wcs_to_raw(samp->set->wcs, wcspos, rawpos);
            for (d2=0; d2<naxis; d2++)
              {
              ix = (int)((rawpos[d2]-xoffset[d2])*xscale[d2]);
              for (d=0; d<naxis; d++)
                {
                d3 = d2*naxis+d;
                dy = rawpos[d]-samp->rawpos[d];
                iy = (int)((dy -yoffset)*yscale);
                if (ix>=0 && ix<CPLOT_PIXERR1DNX+1)
                  {
                  sig2 = 1.0/*fabs(samp->wcsposerr[d])*/;
                  if (sig2>0.0 && samp->flux/samp->fluxerr>=hsn_thresh)
                    {
                    line[d3][ix] += dy/sig2;
                    weight[d3][ix] += 1.0/sig2;
                    }
                  if (ix<CPLOT_PIXERR1DNX && iy>=0 && iy<CPLOT_PIXERR1DNY)
                    {
                    z = (histo[d3][ix][iy] += 1.0);
                    if (z>zmax[d3])
                      zmax[d3] = z;
                    }
                  }
                iy = (int)((dy - boffset)*bscale);
                if (iy>=0 && iy<CPLOT_NPIXERRHISTBIN)
                  {
                  cy = (cuty[d3][iy] += 1.0);
                  if (cy>cutymax[d3])
                  cutymax[d3] = cy;
                  }
                }
	      }
            if (samp->flux/samp->fluxerr>=hsn_thresh)
              {
              for (d2=0; d2<naxis; d2++)
                {
                ix =(int)((rawpos[d2]-xoffset_hsn[d2])*xscale_hsn[d2]);
                for (d=0; d<naxis; d++)
                  {
                  d3 = d2*naxis+d;
                  dy = rawpos[d]-samp->rawpos[d];
                  iy = (int)((dy -yoffset_hsn)*yscale_hsn);
                  if (ix>=0 && ix<CPLOT_PIXERR1DNX_HSN
			&& iy>=0 && iy<CPLOT_PIXERR1DNY_HSN)
                    {
                    z = (histo_hsn[d3][ix][iy] += 1.0);
                    if (z>zmax_hsn[d3])
                      zmax_hsn[d3] = z;
                    }
                  iy = (int)((dy - boffset)*bscale);
                  if (iy>=0 && iy<CPLOT_NPIXERRHISTBIN)
                    {
                    cy = (cuty_hsn[d3][iy] += 1.0);
                    if (cy>cutymax_hsn[d3])
                    cutymax_hsn[d3] = cy;
                    }
                  }
                }
	      }
            }
          }
        }
      }
    }

/* Now plot! */

  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  for (i=0; i<CPLOT_NPIXERRHISTBIN; i++)
    cutbin[i] = boffset+(i+0.5)/bscale;
  for (i=0; i<=CPLOT_PIXERR1DNX+1; i++)
    cutx[i] = 0.5 + (PLFLT)i/CPLOT_PIXERR1DNX*(maxwidth-0.5);

  firstflag = 1;
  plschr(0.0,0.67);
  for (d2=0; d2<naxis; d2++)
    for (d=0; d<naxis; d++)
      {
      d3 = d2*naxis+d;
      margin = 0.1*(maxwidth-0.5);
/* Adjust histogram to fit in the displayed box */
      for (i=0; i<CPLOT_NPIXERRHISTBIN; i++)
        {
        cuty[d3][i] = 0.5 - margin
			+ cuty[d3][i]/cutymax[d3] * 0.9*margin;
        cuty_hsn[d3][i] = 0.5 - margin
			+ cuty_hsn[d3][i]/cutymax_hsn[d3] * 0.9*margin;
        }
      plwid(lwid);
      plenv(-0.5-margin, maxwidth, -maxlim, maxlim, 0, -2);
/* Use a non-linear shade level distribution */
      if (zmax[d3]>=1.0)
        {
/*
        for (i=0; i<CPLOT_NSHADES; i++)
          clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax[d3]+0.5;
        r[0] = 0.96; g[0] = 1.0; b[0] = 0.96;
        r[1] = 0.3; g[1] = 0.4; b[1] = 0.3;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plshades(histo[d3], CPLOT_PIXERR1DNX, CPLOT_PIXERR1DNY, NULL,
		0.5, maxwidth, -maxlim, maxlim,
		clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
*/
        }
      else
        {
        plcol(1);
        plptex(-margin/2.0, maxlim/2.0, 1.0, 0.0, 0.5,
		"No overlapping detections!");
        }
      if (zmax_hsn[d3]>=1.0)
        {
        r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
        r[1] = 0.7; g[1] = 0.7; b[1] = 0.7;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plimage(histo_hsn[d3], CPLOT_PIXERR1DNX_HSN, CPLOT_PIXERR1DNY_HSN,
		0.5, maxwidth, -maxlim, maxlim,
		0.5, zmax_hsn[d3],
		0.5, maxwidth, -maxlim, maxlim);
        }
      sprintf(xlabel, "AXIS%d [pixels]", d2+1);
      sprintf(ylabel, "#gDAXIS%d [pixels]", d+1);
      plscolbg(255,255,255);	/* Force the background colour to white */
      plscol0(15, 0,0,0);	/* Force the foreground colour to black */
/* 1D histograms */
      plcol(3);
      plwid(2*lwid);
      plline(CPLOT_NPIXERRHISTBIN, cuty[d3], cutbin);
      plcol(7);
      plline(CPLOT_NPIXERRHISTBIN, cuty_hsn[d3], cutbin);
      if (zmax[d3]>=1.0)
        {
        plcol(15);
        for (i=0; i<CPLOT_PIXERR1DNX+1; i++)
          if (weight[d3][i]>0.0)
            line[d3][i] /= weight[d3][i];
        plwid(6*lwid);
        plcol(15);
        plline(CPLOT_PIXERR1DNX+1, cutx, line[d3]);
        plwid(3*lwid);
        plcol(3);
        plline(CPLOT_PIXERR1DNX+1, cutx, line[d3]);
        }
      plwid(lwid);
      plcol(15);
      plwid(lwid);
      xl[0] = 0.5 - margin;
      xl[1] = maxwidth;
      yl[0] = yl[1] = 0.0;
      pllsty(2);
      plline(2, xl, yl);
      pllsty(1);
      plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
      if (firstflag)
        {
        sprintf(str, "Instrument A%d: mean 1D astrometric residuals vs "
		"pixel coordinates", instru+1);
        pllab(xlabel, ylabel, str);
        }
      else
        pllab(xlabel, ylabel, "");
      firstflag = 0;
      }

/* Free array of points */
  free(clevel);
  free(cutbin);
  free(cutx);
  for (d3=0; d3<naxis*naxis; d3++)
    {
    plFree2dGrid(histo[d3], CPLOT_PIXERR1DNX, CPLOT_PIXERR1DNY); 
    plFree2dGrid(histo_hsn[d3], CPLOT_PIXERR1DNX_HSN, CPLOT_PIXERR1DNY_HSN);
    free(cuty[d3]);
    free(cuty_hsn[d3]);
    free(line[d3]);
    }

  plend();

/* Recursive stuff*/
  cplot_pixerrhisto1d(fgroups, ngroup, instru, hsn_thresh);

  return RETURN_OK;
  }


/****** cplot_subpixerrhisto1d ************************************************
PROTO	int cplot_subpixerrhisto1d(fgroupstruct **fgroups, int ngroup,
		int instru, double hsn_thresh)
PURPOSE	Plot an astrometric difference histogram between star pairs along 1
	dimension as a function of sub-pixel coordinates for a given
	astrometric intrument.
INPUT	Pointer to an array of field group pointers,
	Number of field groups,
	Astrometric intrument index,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_subpixerrhisto1d(fgroupstruct **fgroups, int ngroup, int instru,
		double hsn_thresh)
  {
   fgroupstruct	*fgroup;
   fieldstruct	*field;
   wcsstruct	*wcs;
   setstruct	*set;
   samplestruct	*samp, *samp2;
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		pixscale[NAXIS], xscale[NAXIS],xscale_hsn[NAXIS],
		xoffset[NAXIS],xoffset_hsn[NAXIS], mean[NAXIS],
		yscale,yscale_hsn, yoffset,yoffset_hsn, boffset,bscale, dy,
		maxwidth,margin, sig2;
   PLFLT	**histo[NAXIS*NAXIS],**histo_hsn[NAXIS*NAXIS],
		*cuty[NAXIS*NAXIS],*cuty_hsn[NAXIS*NAXIS],
		*line[NAXIS*NAXIS],*weight[NAXIS*NAXIS],
		*clevel,*cutbin,*cutx,
		xl[2], yl[2],zmax[NAXIS*NAXIS],zmax_hsn[NAXIS*NAXIS],
		cutymax[NAXIS*NAXIS], cutymax_hsn[NAXIS*NAXIS],
		r[2],g[2],b[2],cpoint[2],
		lim,maxlim, cy, z;
   PLINT	lwid;
   char		xlabel[80], ylabel[80], str[80];
   int		d,d2,d3, f,fg,i,s,n, naxis, nsamp, firstflag, ix,iy, count;

  if (cplot_init(1,fgroups[0]->naxis*fgroups[0]->naxis,
	CPLOT_SUBPIXERROR1D)==RETURN_ERROR)
    {
    cplot_end(CPLOT_SUBPIXERROR1D);
    return RETURN_OK;
    }

  naxis = fgroups[0]->naxis;

  QMALLOC(cutbin, PLFLT, CPLOT_NSUBPIXERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);
  QMALLOC(cutx, PLFLT, CPLOT_SUBPIXERR1DNX+1);
  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  for (d=0; d<naxis; d++)
    {
    for (d2=0; d2<naxis; d2++)
      {
      d3 = d2*naxis+d;
      plAlloc2dGrid(&histo[d3], CPLOT_SUBPIXERR1DNX, CPLOT_SUBPIXERR1DNY);
      plAlloc2dGrid(&histo_hsn[d3], CPLOT_SUBPIXERR1DNX_HSN,
				CPLOT_SUBPIXERR1DNY_HSN);
      QCALLOC(cuty[d3], PLFLT, CPLOT_NSUBPIXERRHISTBIN);
      QCALLOC(cuty_hsn[d3], PLFLT, CPLOT_NSUBPIXERRHISTBIN);
      QCALLOC(line[d3], PLFLT, CPLOT_SUBPIXERR1DNX+1);
      QCALLOC(weight[d3], PLFLT, CPLOT_SUBPIXERR1DNX);
      cutymax[d3] = cutymax_hsn[d3] = zmax[d3] = zmax_hsn[d3] = 0.0;
      }
    xoffset[d] = xoffset_hsn[d] = -0.5;
    xscale[d] = CPLOT_SUBPIXERR1DNX/1.0;
    xscale_hsn[d] = CPLOT_SUBPIXERR1DNX_HSN / 1.0;
    }

  maxlim = 0.0;
  for (fg=0; fg<ngroup; fg++)
    {
    for (d2=0; d2<naxis; d2++)
      if ((lim=fgroups[fg]->sig_interr_hsn[d2]/fgroups[fg]->meanwcsscale[d2])
		> maxlim)
        maxlim = lim;
    if (maxlim<=0.0)
      maxlim = 1.0;
    }
  boffset = -maxlim;
  bscale = CPLOT_NSUBPIXERRHISTBIN / (2.0*maxlim);
  yoffset = yoffset_hsn = -maxlim;
  yscale = CPLOT_SUBPIXERR1DNY/(2.0*maxlim);
  yscale_hsn = CPLOT_SUBPIXERR1DNY_HSN/(2.0*maxlim);

  for (fg=0; fg<ngroup; fg++)
    {
    fgroup = fgroups[fg];
    wcs = fgroup->wcs;
    if (!wcs)
      return RETURN_ERROR;
    for (d=0; d<naxis; d++)
      rawpos2[d] = rawpos[d] = wcs->naxisn[d]/2.0;
    raw_to_wcs(wcs, rawpos, wcspos);
    for (d=0; d<naxis; d++)
      {
      rawpos2[d] += 1.0;
      raw_to_wcs(wcs, rawpos2, wcspos2);
      pixscale[d] = wcs_dist(wcs, wcspos, wcspos2)*DEG/ARCSEC;	/* in arcsec */
      rawpos2[d] -= 1.0;
      }
    for (f=0; f<fgroup->nfield; f++)
      {
      field = fgroup->field[f];
/*---- Skip field if not observed with the right astrometric instrument */
      if (field->astromlabel != instru)
        continue;
      for (s=0; s<field->nset; s++)
        {
        set = field->set[s];
        nsamp = set->nsample;
        samp = set->sample;
        for (n=nsamp; n--; samp++)
          {
          if (samp->flags & (OBJ_SATUR|OBJ_TRUNC))
            continue;
/*-------- Reset mean and count */
          for (d2=0; d2<naxis; d2++)
            mean[d2] = 0.0;
          count = 0;
/*-------- Explore forward */
          samp2 = samp;
          while ((samp2=samp2->nextsamp))
            {
            if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
              continue;
            for (d2=0; d2<naxis; d2++)
              mean[d2] += samp2->projpos[d2] - samp->projpos[d2];
            count++;
            }
/*-------- Explore backward */
          samp2 = samp;
          while ((samp2=samp2->prevsamp) && samp2->set->field->astromlabel>=0)
            {
            if (samp2->flags & (OBJ_SATUR|OBJ_TRUNC))
              continue;
            for (d2=0; d2<naxis; d2++)
              mean[d2] += samp2->projpos[d2] - samp->projpos[d2];
            count++;
            }
          if (count)
            {
            for (d2=0; d2<naxis; d2++)
              mean[d2] = mean[d2]/count + samp->projpos[d2];
/*---------- Convert to raw local coordinates */
            raw_to_wcs(wcs, mean, wcspos);
            wcs_to_raw(samp->set->wcs, wcspos, rawpos);
            for (d2=0; d2<naxis; d2++)
              {
              ix = (int)((rawpos[d2]-(int)(rawpos[d2]+0.4999999)
			-xoffset[d2])*xscale[d2]);
              for (d=0; d<naxis; d++)
                {
                d3 = d2*naxis+d;
                dy = rawpos[d]-samp->rawpos[d];
                iy = (int)((dy -yoffset)*yscale);
                if (ix>=0 && ix<CPLOT_SUBPIXERR1DNX)
                  {
                  sig2 = 1.0/*fabs(samp->wcsposerr[d])*/;
                  if (sig2>0.0 && samp->flux/samp->fluxerr>=hsn_thresh)
                    {
                    line[d3][ix] += dy/sig2;
                    weight[d3][ix] += 1.0/sig2;
                    }
                  if (iy>=0 && iy<CPLOT_SUBPIXERR1DNY)
                    {
                    z = (histo[d3][ix][iy] += 1.0);
                    if (z>zmax[d3])
                      zmax[d3] = z;
                    }
                  }
                iy = (int)((dy - boffset)*bscale);
                if (iy>=0 && iy<CPLOT_NSUBPIXERRHISTBIN)
                  {
                  cy = (cuty[d3][iy] += 1.0);
                  if (cy>cutymax[d3])
                  cutymax[d3] = cy;
                  }
                }
	      }
            if (samp->flux/samp->fluxerr>=hsn_thresh)
              {
              for (d2=0; d2<naxis; d2++)
                {
                ix =(int)((rawpos[d2]-(int)(rawpos[d2]+0.4999999)
			-xoffset_hsn[d2])*xscale_hsn[d2]);
                for (d=0; d<naxis; d++)
                  {
                  d3 = d2*naxis+d;
                  dy = rawpos[d]-samp->rawpos[d];
                  iy = (int)((dy -yoffset_hsn)*yscale_hsn);
                  if (ix>=0 && ix<CPLOT_SUBPIXERR1DNX_HSN
			&& iy>=0 && iy<CPLOT_SUBPIXERR1DNY_HSN)
                    {
                    z = (histo_hsn[d3][ix][iy] += 1.0);
                    if (z>zmax_hsn[d3])
                      zmax_hsn[d3] = z;
                    }
                  iy = (int)((dy - boffset)*bscale);
                  if (iy>=0 && iy<CPLOT_NSUBPIXERRHISTBIN)
                    {
                    cy = (cuty_hsn[d3][iy] += 1.0);
                    if (cy>cutymax_hsn[d3])
                    cutymax_hsn[d3] = cy;
                    }
                  }
                }
	      }
            }
          }
        }
      }
    }

/* Now plot! */

  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
  for (i=0; i<CPLOT_NSUBPIXERRHISTBIN; i++)
    cutbin[i] = boffset+(i+0.5)/bscale;
  for (i=0; i<=CPLOT_SUBPIXERR1DNX; i++)
    cutx[i] = -0.5 + (PLFLT)i/CPLOT_SUBPIXERR1DNX;

  firstflag = 1;
  plschr(0.0,0.67);
  for (d2=0; d2<naxis; d2++)
    for (d=0; d<naxis; d++)
      {
      d3 = d2*naxis+d;
      maxwidth = 1.0;
      margin = 0.1*maxwidth;
/* Adjust histogram to fit in the displayed box */
      for (i=0; i<CPLOT_NSUBPIXERRHISTBIN; i++)
        {
        cuty[d3][i] = -0.5 - margin
			+ cuty[d3][i]/cutymax[d3] * 0.9*margin;
        cuty_hsn[d3][i] = -0.5 - margin
			+ cuty_hsn[d3][i]/cutymax_hsn[d3] * 0.9*margin;
        }
      plwid(lwid);
      plenv(-0.5-margin, 0.5, -maxlim, maxlim, 0, -2);
/* Use a non-linear shade level distribution */
      if (zmax[d3]>=1.0)
        {
/*
        for (i=0; i<CPLOT_NSHADES; i++)
          clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax[d3]+0.5;
        r[0] = 0.96; g[0] = 1.0; b[0] = 0.96;
        r[1] = 0.3; g[1] = 0.4; b[1] = 0.3;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plshades(histo[d3], CPLOT_SUBPIXERR1DNX, CPLOT_SUBPIXERR1DNY, NULL,
		-0.5, 0.5, -maxlim, maxlim,
		clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
*/
        }
      else
        {
        plcol(1);
        plptex(-margin/2.0, maxlim/2.0, 1.0, 0.0, 0.5,
		"No overlapping detections!");
        }
      if (zmax_hsn[d3]>=1.0)
        {
        r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
        r[1] = 0.7; g[1] = 0.7; b[1] = 0.7;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plimage(histo_hsn[d3], CPLOT_SUBPIXERR1DNX_HSN,
		CPLOT_SUBPIXERR1DNY_HSN,
		-0.5, 0.5, -maxlim, maxlim,
		0.5, zmax_hsn[d3],
		-0.5, 0.5, -maxlim, maxlim);
        }
      sprintf(xlabel, "AXIS%d [sub-pixel]", d2+1);
      sprintf(ylabel, "#gDAXIS%d [pixels]", d+1);
      plscolbg(255,255,255);	/* Force the background colour to white */
      plscol0(15, 0,0,0);	/* Force the foreground colour to black */
/* 1D histograms */
      plcol(3);
      plwid(2*lwid);
      plline(CPLOT_NSUBPIXERRHISTBIN, cuty[d3], cutbin);
      plcol(7);
      plline(CPLOT_NSUBPIXERRHISTBIN, cuty_hsn[d3], cutbin);
      if (zmax[d3]>=1.0)
        {
        plcol(15);
        for (i=0; i<CPLOT_SUBPIXERR1DNX; i++)
          if (weight[d3][i]>0.0)
            line[d3][i] /= weight[d3][i];
        line[d3][CPLOT_SUBPIXERR1DNX] = line[d3][0];
        plwid(6*lwid);
        plcol(15);
        plline(CPLOT_SUBPIXERR1DNX+1, cutx, line[d3]);
        plwid(3*lwid);
        plcol(3);
        plline(CPLOT_SUBPIXERR1DNX+1, cutx, line[d3]);
        }
      plwid(lwid);
      plcol(15);
      plwid(lwid);
      xl[0] = -0.5 - margin;
      xl[1] = 0.5;
      yl[0] = yl[1] = 0.0;
      pllsty(2);
      plline(2, xl, yl);
      pllsty(1);
      plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
      if (firstflag)
        {
        sprintf(str, "Instrument A%d: mean 1D astrometric residuals vs "
		"sub-pixel coordinates", instru+1);
        pllab(xlabel, ylabel, str);
        }
      else
        pllab(xlabel, ylabel, "");
      firstflag = 0;
      }

/* Free array of points */
  free(clevel);
  free(cutbin);
  free(cutx);
  for (d3=0; d3<naxis*naxis; d3++)
    {
    plFree2dGrid(histo[d3], CPLOT_SUBPIXERR1DNX, CPLOT_SUBPIXERR1DNY); 
    plFree2dGrid(histo_hsn[d3],CPLOT_SUBPIXERR1DNX_HSN,CPLOT_SUBPIXERR1DNY_HSN);
    free(cuty[d3]);
    free(cuty_hsn[d3]);
    free(line[d3]);
    }

  plend();

/* Recursive stuff*/
  cplot_subpixerrhisto1d(fgroups, ngroup, instru, hsn_thresh);

  return RETURN_OK;
  }


/****** cplot_astrcolshift1d **************************************************
PROTO	int cplot_astrcolshift1d(fgroupstruct *fgroup, double hsn_thresh)
PURPOSE	Plot an astrometric difference histogram between star pairs along 1
	dimension as a function of magnitude-difference (colour).
INPUT	Pointer to the field group,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	astrcolshift_fgroup() must have been run on group first.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_astrcolshift1d(fgroupstruct *fgroup, double hsn_thresh)
  {
   fieldstruct	*field;
   wcsstruct	*wcs;
   setstruct	*set;
   samplestruct	*samp, *samp1, *samp2;
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		pixscale[NAXIS], xscale,xscale_hsn, xoffset,xoffset_hsn,
		yscale,yscale_hsn, yoffset,yoffset_hsn, boffset,bscale, dy,
		maxwidth,margin, dmag, mdmag, mdmag2, sdmag, ndmag,
		ymin, ymax;
   PLFLT	**histo[NAXIS],**histo_hsn[NAXIS],
		**histot,
		*cuty[NAXIS],*cuty_hsn[NAXIS],
		*clevel,*cutbin,
		cutymax[NAXIS], cutymax_hsn[NAXIS],
		xl[2], yl[2],zmax[NAXIS],zmax_hsn[NAXIS],
		r[2],g[2],b[2],cpoint[2],
		lim,maxlim, cy, z, dmagmin, dmagmax;
   PLINT	lwid;
   char		xlabel[80], ylabel[80];
   int		d, f,f1,f2,ff, i,s,n, naxis, nsamp, ninstru, nfield, ix,iy, gra,
		instru1,instru2;

  ninstru = prefs.nphotinstrustr;
  naxis = fgroup->naxis;
  nfield = fgroup->nfield;

  if (cplot_init(ninstru,ninstru*naxis,CPLOT_ASTRCOLSHIFT1D) == RETURN_ERROR)
    {
    cplot_end(CPLOT_ASTRCOLSHIFT1D);
    return RETURN_OK;
    }

  wcs = fgroup->wcs;
  if (!wcs)
    return RETURN_ERROR;

  for (d=0; d<naxis; d++)
    rawpos2[d] = rawpos[d] = wcs->naxisn[d]/2.0;
  raw_to_wcs(wcs, rawpos, wcspos);

  QMALLOC(cutbin, PLFLT, CPLOT_NADERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);
  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  maxlim = 0.0;
  for (d=0; d<naxis; d++)
    if ((lim=fgroup->sig_interr_hsn[d]*DEG/ARCSEC*2.0) > maxlim)
      maxlim = lim;
  if (maxlim<=0.0)
    maxlim = 1.0;
  boffset = -maxlim;
  bscale = CPLOT_NADERRHISTBIN / (2.0*maxlim);
  yoffset = yoffset_hsn = -maxlim;
  yscale = CPLOT_ASTRCOLSHIFT1DNY/(2.0*maxlim);
  yscale_hsn = CPLOT_ASTRCOLSHIFT1DNY_HSN/(2.0*maxlim);

  for (d=0; d<naxis; d++)
    {
    rawpos2[d] += 1.0;
    raw_to_wcs(wcs, rawpos2, wcspos2);
    pixscale[d] = wcs_dist(wcs, wcspos, wcspos2)*DEG/ARCSEC;	/* in arcsec */
    rawpos2[d] -= 1.0;
    plAlloc2dGrid(&histo[d], CPLOT_ASTRCOLSHIFT1DNX, CPLOT_ASTRCOLSHIFT1DNY);
    plAlloc2dGrid(&histo_hsn[d], CPLOT_ASTRCOLSHIFT1DNX_HSN,
				CPLOT_ASTRCOLSHIFT1DNY_HSN);
    QMALLOC(cuty[d], PLFLT, CPLOT_NADERRHISTBIN);
    QMALLOC(cuty_hsn[d], PLFLT, CPLOT_NADERRHISTBIN);
    }

  for (i=0; i<CPLOT_NADERRHISTBIN; i++)
    cutbin[i] = boffset+(i+0.5)/bscale;

  for (instru1=0; instru1<ninstru; instru1++)
    {
    for (instru2=0; instru2<ninstru; instru2++)
      {
/*----- Initialize histograms */
      for (d=0; d<fgroup->naxis; d++)
        {
        histot = histo[d];
        for (ix=0; ix<CPLOT_ASTRCOLSHIFT1DNX; ix++)
          memset(histot[ix], 0, CPLOT_ASTRCOLSHIFT1DNY*sizeof(PLFLT));
        histot = histo_hsn[d];
        for (ix=0; ix<CPLOT_ASTRCOLSHIFT1DNX_HSN; ix++)
          memset(histot[ix], 0, CPLOT_ASTRCOLSHIFT1DNY_HSN*sizeof(PLFLT));
        memset(cuty[d], 0, CPLOT_NADERRHISTBIN*sizeof(PLFLT));
        memset(cuty_hsn[d], 0, CPLOT_NADERRHISTBIN*sizeof(PLFLT));
        cutymax[d] = cutymax_hsn[d] = zmax[d] = zmax_hsn[d] = 0.0;
        }

/*---- Find the range in Delta-mag. */
      mdmag = mdmag2 = ndmag = 0.0;
      for (f=0; f<nfield; f++)
        {
        field = fgroup->field[f];
        for (s=0; s<field->nset; s++)
          {
          set = field->set[s];
          nsamp = set->nsample;
          samp1 = set->sample;
          for (n=nsamp; n--; samp1++)
            if (!samp1->nextsamp && samp1->prevsamp
		&& !(samp1->flags & (OBJ_SATUR|OBJ_TRUNC)))
              {
/*------------ Look for a counterpart from the right photometric instrument */
              for (samp = samp1; samp && samp->set->field->photomlabel>=0;
			samp=samp->prevsamp)
                {
/*-------------- Don't bother if field is a different instru or photom. ref */
/*-------------- or the flux is negative */
                if (samp->set->field->photomlabel != instru1
			|| samp->flux <= 0.0
			|| (samp->flags & (OBJ_SATUR|OBJ_TRUNC)))
                  continue;
                for (samp2=samp1; samp2 && samp2->set->field->photomlabel>=0;
			samp2=samp2->prevsamp)
                  {
/*---------------- Don't bother if field is a different instru or photom. ref*/
/*---------------- or if the flux is negative */
                  if (samp2==samp || samp2->set->field->photomlabel != instru2
			|| samp2->flux <= 0.0
			|| (samp2->flags & (OBJ_SATUR|OBJ_TRUNC)))
                    continue;
                  dmag = samp2->mag - samp->mag;
                  mdmag += dmag;
                  mdmag2 += dmag*dmag;
                  ndmag += 1.0;
                  }
                }
              }
          }
        }

      if (ndmag)
        {
        mdmag /= ndmag;
        sdmag = sqrt(fabs(mdmag2-mdmag*mdmag)/ndmag);
        }
      else
        sdmag = mdmag = 0.0;
      dmagmin = mdmag - 2.0*sdmag;
      dmagmax = mdmag + 2.0*sdmag;
      xoffset = xoffset_hsn = dmagmin;
      xscale = CPLOT_ASTRCOLSHIFT1DNX/(dmagmax - dmagmin);
      xscale_hsn = CPLOT_ASTRCOLSHIFT1DNX_HSN/(dmagmax - dmagmin);

      for (f=0; f<nfield; f++)
        {
        field = fgroup->field[f];
        for (s=0; s<field->nset; s++)
          {
          set = field->set[s];
          nsamp = set->nsample;
          samp1 = set->sample;
          for (n=nsamp; n--; samp1++)
            if (!samp1->nextsamp && samp1->prevsamp)
              {
/*------------ Look for a counterpart from the right photometric instrument */
              for (samp = samp1; samp && samp->set->field->photomlabel>=0;
			samp=samp->prevsamp)
                {
/*-------------- Don't bother if field is a different instru or photometric ref */
/*-------------- or the flux is negative */
                if (samp->set->field->photomlabel != instru1
			|| samp->flux <= 0.0
			|| (samp->flags & (OBJ_SATUR|OBJ_TRUNC)))
                  continue;
                for (samp2=samp1; samp2 && samp2->set->field->photomlabel>=0;
			samp2=samp2->prevsamp)
                  {
/*---------------- Don't bother if field is a different instru or photom. ref*/
/*---------------- or if the flux is negative */
                  if (samp2==samp || samp2->set->field->photomlabel != instru2
			|| samp2->flux <= 0.0
			|| (samp2->flags & (OBJ_SATUR|OBJ_TRUNC)))
                    continue;
                  ix = (int)((samp2->mag - samp->mag - xoffset)*xscale);
                  for (d=0; d<naxis; d++)
                    {
                    dy = (samp2->projpos[d]-samp->projpos[d])*pixscale[d];
                    iy = (int)((dy -yoffset)*yscale);
                    if (ix>=0 && ix<CPLOT_ASTRCOLSHIFT1DNX
			&& iy>=0 && iy<CPLOT_ASTRCOLSHIFT1DNY)
                      {
                      z = (histo[d][ix][iy] += 1.0);
                      if (z>zmax[d])
                        zmax[d] = z;
                      }
                    iy = (int)((dy - boffset)*bscale);
                    if (iy>=0 && iy<CPLOT_NADERRHISTBIN)
                      {
                      cy = (cuty[d][iy] += 1.0);
                      if (cy>cutymax[d])
                        cutymax[d] = cy;
                      }
                    }
                  if (samp->flux/samp->fluxerr < hsn_thresh
			|| samp2->flux/samp2->fluxerr < hsn_thresh)
                    continue;
                  ix = (int)((samp2->mag - samp->mag - xoffset_hsn)*xscale_hsn);
                  for (d=0; d<naxis; d++)
                    {
                    dy = (samp2->projpos[d]-samp->projpos[d])*pixscale[d];
                    iy = (int)((dy -yoffset_hsn)*yscale_hsn);
                    if (ix>=0 && ix<CPLOT_ASTRCOLSHIFT1DNX_HSN
			&& iy>=0 && iy<CPLOT_ASTRCOLSHIFT1DNY_HSN)
                      {
                      z = (histo_hsn[d][ix][iy] += 1.0);
                      if (z>zmax_hsn[d])
                        zmax_hsn[d] = z;
                      }
                    iy = (int)((dy - boffset)*bscale);
                    if (iy>=0 && iy<CPLOT_NADERRHISTBIN)
                      {
                      cy = (cuty_hsn[d][iy] += 1.0);
                      if (cy>cutymax_hsn[d])
                        cutymax_hsn[d] = cy;
                      }
                    }
                  }
                }
              }
          }
        }      

/*---- Now plot! */

      lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
      plschr(0.0,0.5);
      for (d=0; d<naxis; d++)
        {
        if (dmagmin == dmagmax)
          dmagmin = -(dmagmax = 1.0);
        maxwidth = dmagmax - dmagmin;
        margin = 0.1*maxwidth;
/*------ Adjust histogram to fit in the displayed box */
        for (i=0; i<CPLOT_NADERRHISTBIN; i++)
          {
          cuty[d][i] = dmagmin - margin
			+ cuty[d][i]/cutymax[d] * 0.9*margin;
          cuty_hsn[d][i] = dmagmin - margin
			+ cuty_hsn[d][i]/cutymax_hsn[d] * 0.9*margin;
          }
        plwid(lwid);
        gra = instru2 + (instru1*naxis + d)*ninstru;
        if (gra)
          pladv(gra);
        plenv(dmagmin - margin, dmagmax, -maxlim, maxlim, 0, -2);
/*------ Use a non-linear shade level distribution */
        if (zmax[d]>=1.0)
          {
          for (i=0; i<CPLOT_NSHADES; i++)
            clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax[d]+0.5;
          r[0] = 1.0; g[0] = 0.96; b[0] = 1.0;
          r[1] = 0.3; g[1] = 0.2; b[1] = 0.3;
          plscmap1l(1, 2, cpoint, r, g, b, NULL);
          plshades(histo[d], CPLOT_ASTRCOLSHIFT1DNX, CPLOT_ASTRCOLSHIFT1DNY,
		NULL, dmagmin, dmagmax, -maxlim, maxlim,
		clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
          }
        else
          {
          plcol(1);
          plptex((PLFLT)mdmag, maxlim/2.0, 1.0, 0.0, 0.5,
		"No overlapping detections!");
          }
        if (zmax_hsn[d]>=1.0)
          {
          r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
          r[1] = 0.7; g[1] = 0.7; b[1] = 0.7;
          plscmap1l(1, 2, cpoint, r, g, b, NULL);
          plimage(histo_hsn[d], CPLOT_ASTRCOLSHIFT1DNX_HSN,
		CPLOT_ASTRCOLSHIFT1DNY_HSN,
		dmagmin, dmagmax, -maxlim, maxlim, 0.5, zmax_hsn[d],
		dmagmin, dmagmax, -maxlim, maxlim);
          }
        sprintf(xlabel, "P%d - P%d [mag]", instru2+1, instru1+1);
        sprintf(ylabel, "#gDAXIS%d [\"]", d+1);
        plscolbg(255,255,255);	/* Force the background colour to white */
        plscol0(15, 0,0,0);	/* Force the foreground colour to black */
/*------ 1D histograms */
        plcol(13);
        plwid(2*lwid);
        plline(CPLOT_NADERRHISTBIN, cuty[d], cutbin);
        plcol(7);
        plline(CPLOT_NADERRHISTBIN, cuty_hsn[d], cutbin);
        plwid(lwid);
        plcol(15);
        plwid(lwid);
        xl[0] = dmagmin;
        xl[1] = dmagmax;
        ymin = ymax = 0.0;
        n = 0;
        for (f1 = 0; f1<nfield; f1++)
          if (fgroup->field[f1]->photomlabel == instru1)
            for (f2 = 0; f2<nfield; f2++)
              if (fgroup->field[f2]->photomlabel == instru2)
                {
                ff = f1*nfield+f2;
                ymin += fgroup->intcolshiftzero[d][ff]
			+ dmagmin*fgroup->intcolshiftscale[d][ff];
                ymax += fgroup->intcolshiftzero[d][ff]
			+ dmagmax*fgroup->intcolshiftscale[d][ff];
                n++;
                }
        if (n)
          {
          yl[0] = ymin*pixscale[d] / n;
          yl[1] = ymax*pixscale[d] / n;
          plline(2, xl, yl);
          }
        pllsty(2);
        xl[0] = dmagmin - margin;
        xl[1] = dmagmax;
        yl[0] = yl[1] = 0.0;
        plline(2, xl, yl);
        xl[0] = xl[1] = 0.0;
        yl[0] = -(yl[1] = maxlim);
        plline(2, xl, yl);
        pllsty(1);
        plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
        pllab(xlabel, ylabel, "");
        }
      }
    }

/* Free array of points */
  free(clevel);
  free(cutbin);
  for (d=0; d<naxis; d++)
    {
    plFree2dGrid(histo[d], CPLOT_ASTRCOLSHIFT1DNX, CPLOT_ASTRCOLSHIFT1DNY); 
    plFree2dGrid(histo_hsn[d], CPLOT_ASTRCOLSHIFT1DNX_HSN,
		CPLOT_ASTRCOLSHIFT1DNY_HSN); 
    free(cuty[d]);
    free(cuty_hsn[d]);
    }

  plend();

  cplot_astrcolshift1d(fgroup, hsn_thresh);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_astrefprop ******************************************************
PROTO	int cplot_astrefprop(fgroupstruct *fgroup, fieldstruct *reffield,
			double hsn_thresh)
PURPOSE	Plot recovered vs ref catalog proper motions as a 2D histogram.
INPUT	Pointer to the field group,
	pointer to the reference field,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() and astrprop_fgroup() must have been run on all groups
	first.
AUTHOR	E. Bertin (IAP)
VERSION	10/09/2009
 ***/
int	cplot_astrefprop(fgroupstruct *fgroup, fieldstruct *reffield,
			double hsn_thresh)
  {
   wcsstruct	*wcs;
   setstruct	*set;
   samplestruct	*samp, *samp2;
   char		str[64];
   double	offset,offset_hsn, scale,scale_hsn, boffset,bscale, err;
   PLFLT	**histo[2],**histo_hsn[2],
		xl[5], yl[5],r[2],g[2],b[2],cpoint[2],zmax[2],zmax_hsn[2],
		*clevel,
		lim,maxlim, z;
   PLINT	lwid;
   int		d,d2, i,s,n, ix,iy, nsamp, lng,lat;

  if (cplot_init(1,1, CPLOT_REFPROP) == RETURN_ERROR)
    {
    cplot_end(CPLOT_REFPROP);
    return RETURN_OK;
    }

  wcs = fgroup->wcs;
  lng = wcs->lng;
  lat = wcs->lat;
  if (!wcs || lng<0 || lat<0)
    return RETURN_ERROR;

  maxlim = 0.0;
  for (d2=0; d2<fgroup->naxis; d2++)
    if ((lim=fgroup->sig_interr[d2]*DEG/MAS*2.0) > maxlim)
      maxlim = lim;
  if (maxlim<=0.0)
    maxlim = 1.0;

  for (d=0; d<2; d++)
    {
    plAlloc2dGrid(&histo[d], CPLOT_REFPROPN, CPLOT_REFPROPN);
    plAlloc2dGrid(&histo_hsn[d], CPLOT_REFPROPN_HSN, CPLOT_REFPROPN_HSN);
    zmax[d] = zmax_hsn[d] = 0.0;
    }
  offset = offset_hsn = -maxlim;
  scale = CPLOT_REFERR2DN / (2.0*maxlim);
  scale_hsn = CPLOT_REFERR2DN_HSN / (2.0*maxlim);

  boffset = -maxlim;
  bscale = CPLOT_NREFERRHISTBIN / (2.0*maxlim);

  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);

  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  for (s=0; s<reffield->nset; s++)
    {
    set = reffield->set[s];
    nsamp = set->nsample;
    samp = set->sample;
    for (n=nsamp; n--; samp++)
      if ((samp2=samp->nextsamp) && !samp->prevsamp
	&& !(samp->flags & (OBJ_SATUR|OBJ_TRUNC))
	&& !(samp2->flags & (OBJ_SATUR|OBJ_TRUNC)))
        {
/*------ Do not plot objects with bad S/N on proper motions */
        err = samp->wcsproperr[lng]*samp->wcsproperr[lng]
		+ samp->wcsproperr[lat]*samp->wcsproperr[lat];
/*
        if (err == 0.0 || (samp->wcsprop[0]*samp->wcsprop[0]
		+ samp->wcsprop[1]*samp->wcsprop[1])/err
		< CPLOT_ASTREFPROPMINSN*CPLOT_ASTREFPROPMINSN)
          continue;
*/
        ix = (int)((samp->wcsprop[lng]*DEG/MAS - offset)*scale);
        iy = (int)((samp2->wcsprop[lng]*DEG/MAS - offset)*scale);
        if (ix>=0 && ix<CPLOT_REFPROPN && iy>=0 && iy<CPLOT_REFPROPN)
          {
          z = (histo[0][ix][iy] += 1.0);
          if (z>zmax[0])
            zmax[0] = z;
          }
        if (samp2->flux/samp2->fluxerr>=hsn_thresh)
          {
          ix = (int)((samp->wcsprop[lng]*DEG/MAS - offset_hsn)*scale_hsn);
          iy = (int)((samp2->wcsprop[lng]*DEG/MAS - offset_hsn)*scale_hsn);
          if (ix>=0 && ix<CPLOT_REFPROPN_HSN
		&& iy>=0 && iy<CPLOT_REFPROPN_HSN)
            {
            z = (histo_hsn[0][ix][iy] += 1.0);
            if (z>zmax_hsn[0])
              zmax_hsn[0] = z;
            }
          }
        ix = (int)((samp->wcsprop[lat]*DEG/MAS - offset)*scale);
        iy = (int)((samp2->wcsprop[lat]*DEG/MAS - offset)*scale);
        if (ix>=0 && ix<CPLOT_REFPROPN && iy>=0 && iy<CPLOT_REFPROPN)
          {
          z = (histo[1][ix][iy] += 1.0);
          if (z>zmax[1])
            zmax[1] = z;
          }
        if (samp2->flux/samp2->fluxerr>=hsn_thresh)
          {
          ix = (int)((samp->wcsprop[lat]*DEG/MAS - offset_hsn)*scale_hsn);
          iy = (int)((samp2->wcsprop[lat]*DEG/MAS - offset_hsn)*scale_hsn);
          if (ix>=0 && ix<CPLOT_REFPROPN_HSN
		&& iy>=0 && iy<CPLOT_REFPROPN_HSN)
            {
            z = (histo_hsn[1][ix][iy] += 1.0);
            if (z>zmax_hsn[1])
              zmax_hsn[1] = z;
            }
          }
        }
    }

  lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;

  plschr(0.0,0.5);
  for (d=0; d<2; d++)
    {
    plwid(lwid);
    pladv(1);
    plvpas(0.07+d*0.5,0.48+d*0.5,0.1,0.9,1.0);
    plwind(-maxlim,maxlim,-maxlim,maxlim);
/*
    plenv(-maxlim, maxlim, -maxlim, maxlim, 1, -2);
*/
/*-- Use a non-linear shade level distribution */
    if (zmax[d]>=1.0)
      {
      for (i=0; i<CPLOT_NSHADES; i++)
        clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax[d]+0.5;
      r[0] = 1.0; g[0] = 0.96; b[0] = 1.0;
      r[1] = 0.3; g[1] = 0.2; b[1] = 0.3;
      plscmap1l(1, 2, cpoint, r, g, b, NULL);
      plshades(histo[d], CPLOT_REFPROPN, CPLOT_REFPROPN, NULL,
	-maxlim,maxlim, -maxlim,maxlim,
	clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
      }
    else
      {
      plcol(1);
      plptex(0.0, maxlim/2.0, 1.0, 0.0, 0.5, "No match with a reference!");
      }
    if (zmax_hsn[d]>=1.0)
      {
      r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
      r[1] = 0.8; g[1] = 0.8; b[1] = 0.8;
      plscmap1l(1, 2, cpoint, r, g, b, NULL);
      plimage(histo_hsn[d], CPLOT_REFERR2DN_HSN, CPLOT_REFERR2DN_HSN,
	-maxlim,maxlim, -maxlim, maxlim,
	0.5, zmax_hsn[d],
	-maxlim, maxlim, -maxlim, maxlim);
      }
    plscolbg(255,255,255);	/* Force the background colour to white */
    plscol0(15, 0,0,0);		/* Force the foreground colour to black */
    plcol(15);
    plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
    pllab(d?"#gm#d#gd#u(ref) [mas/yr]":"#gm#d#ga#u(ref) [mas/yr]",
	d?"#gm#d#gd#u(SCAMP) [mas/yr]":"#gm#d#ga#u(SCAMP) [mas/yr]", "");
/* reticulus */
    pllsty(2);
    xl[0] = -maxlim;
    xl[1] = maxlim;
    yl[0] = yl[1] = 0.0;
    plline(2, xl, yl);
    xl[0] = xl[1] = 0.0;
    yl[0] = -maxlim;
    yl[1] = maxlim;
    plline(2, xl, yl);
    pllsty(1);
    }
  pladv(1);
  plvpor(0.0,1.0,0.0,1.0);
  sprintf(str,"Group ##%d: "
	"Proper motions - measured #fivs#fn ref. astrometric catalog",
	fgroup->no);
  plmtex("t", -4.0, 0.5,0.5, str);

/*-- Free array of points */
  for (d=0; d<2; d++)
    {
    plFree2dGrid(histo[d], CPLOT_REFERR2DN, CPLOT_REFERR2DN); 
    plFree2dGrid(histo_hsn[d], CPLOT_REFERR2DN_HSN, CPLOT_REFERR2DN_HSN); 
    }
  free(clevel);

  plend();

  cplot_astrefprop(fgroup, reffield, hsn_thresh);	/* Recursive stuff */

  return RETURN_OK;
  }

