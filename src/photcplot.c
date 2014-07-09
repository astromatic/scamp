/*
*				photcplot.c
*
* Produce photometric check plots.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2014 Emmanuel Bertin -- IAP/CNRS/UPMC
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
struct	focplanestruct {PLFLT x[5], y[5], z[5]; PLINT colour; char *str;};
extern int		plotaaflag;

int		comp_focz(const void *focplane1, const void *focplane2);

/****** cplot_photom *******************************************************
PROTO	int cplot_photom(fgroupstruct **fgroups, int ngroup,
		fieldstruct **reffields)
PURPOSE	Plot photometric relations between fgroups and the reference field
	(if present).
INPUT	Pointer to the array of field groups,
	number of groups,
	pointer to the array of reference field to the groups.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2004
 ***/
int	cplot_photom(fgroupstruct **fgroups, int ngroup,
		fieldstruct **reffields)
  {
   setstruct	*refset;
   samplestruct	*refsamp, *samp;
   PLFLT	*x,*y;
   int		g,n, npoint, nrefsamp, ngw,ngh, lng,lat;

  if (cplot_init(1,1, CPLOT_PHOTOM) == RETURN_ERROR)
    {
    cplot_end(CPLOT_PHOTOM);
    return RETURN_OK;
    }
  x = y = NULL;		/* to avoid gcc -Wall warnings */
  ngw = (int)(sqrt(ngroup)+0.49);
  ngh = (int)(ngroup/ngw+0.49);
  for (g=0; g<ngroup; g++)
    {
    refset = reffields[g]->set[0];
    lng = refset->lng;
    lat = refset->lat;
    refsamp = refset->sample;
    nrefsamp = refset->nsample;
    if (nrefsamp)
      {
      QMALLOC(x, PLFLT, nrefsamp);
      QMALLOC(y, PLFLT, nrefsamp);
      }
    pladv(0);
    npoint = 0;
    plenv(10.0,22.0,10.0,22.0, 1, 0);
    for (n=nrefsamp; n--; refsamp++)
      {
      if ((samp=refsamp->nextsamp))
        {
        x[npoint] = refsamp->mag;
        y[npoint] = samp->mag+30.0;
        npoint++;
        }
      }
    if (nrefsamp)
      {
      plpoin((PLINT)npoint, x,y, 20);
      free(x);
      free(y);
      }
    plflush();
    }

  plend();

  cplot_photom(fgroups, ngroup, reffields);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_photzp *******************************************************
PROTO	int cplot_photzp(fgroupstruct *fgroup)
PURPOSE	Plot photometric zero-point corrections.
INPUT	Pointer to the field group.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	09/07/2014
 ***/
int	cplot_photzp(fgroupstruct *fgroup)
  {
   fieldstruct	**fields;
   PLFLT	xl[2], yl[2],
		*x,*y,*xt,*yt,dy, lim, ymin,ymax;
   PLINT	lwid;
   char		xlabel[80], ylabel[80], str[80];
   int		f,n, npointmax,
		firstflag, instru, ninstru, npinstru, nx,ny;

  ninstru = prefs.nphotinstrustr;
  npinstru = 0;
  for (instru=0; instru<ninstru; instru++)
    if (fgroup->chi2_intmag[instru]!=0.0)
      npinstru++;
  if (!npinstru)
    npinstru = 1;

  nx = 1;
  ny = npinstru;

  if (cplot_init(nx, ny , CPLOT_PHOTZP) == RETURN_ERROR)
    {
    cplot_end(CPLOT_PHOTZP);
    return RETURN_OK;
    }

  fields = fgroup->field;
  for (instru=0; instru<ninstru; instru++)
    {
    if (fgroup->chi2_intmag[instru]==0.0)
      continue;
    npointmax = 0;
    for (f=0; f<fgroup->nfield; f++)
      if (fields[f]->photomlabel == instru)
        npointmax++;
    QMALLOC(x, PLFLT, npointmax);
    QMALLOC(y, PLFLT, npointmax);
    xt = x;
    yt = y;
    ymin = BIG;
    ymax = -BIG;
    n = 0;
    for (f=0; f<fgroup->nfield; f++)
      {
      if (fields[f]->photomlabel == instru)
        {
        *(xt++) = (PLFLT)++n;
        *(yt++) = lim = fields[f]->dmagzero;
        if (lim<ymin)
  	  ymin = lim;
        if (lim>ymax)
  	  ymax = lim;
        }
      }
    dy = ymax - ymin;
    ymin -= dy*0.2;
    ymax += dy*0.05;

/*-- Now plot! */
    firstflag = 1;
    yl[0] = yl[1] = 0.0;
    plcol0(15);
    plschr(0.0,0.5);
    lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
    CPLOT_PLWID(lwid);
    plenv(0.0, npointmax + 1.0, ymin, ymax, 0, 0);
    sprintf(xlabel, "Field ##");
    sprintf(ylabel, "#gDZP [mag]");
    if (firstflag)
      {
      sprintf(str, "Group ##%d / Instrument P%d: Zero-point corrections",
		fgroup->no, instru+1);
      pllab(xlabel, ylabel, str);
      }
    else
      pllab(xlabel, ylabel, "");
    plssym(0.0,1.0);
    plschr(0.0,0.5);
    n = 0;
    for (f=0; f<fgroup->nfield; f++)
      if (fields[f]->photomlabel == instru)
        {
        if (fields[f]->photomflag==1)
 	  plcol0(9);
        else
 	  plcol0(8);
        plpoin((PLINT)1, x+n,y+n, 5);
        plptex(x[n],y[n], 0.0, -1.0, -0.1, fields[f]->rfilename);
        n++;
        }
    xl[0] = 0.0;
    xl[1] = npointmax+1.0;
    pllsty(2);
    plcol0(15);
    plline(2, xl, yl);
    pllsty(1);
    firstflag = 0;

/*-- Free array of points */
    free(x);
    free(y);
    }

  plend();

  cplot_photzp(fgroup);		/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_photerrhisto ****************************************************
PROTO	int cplot_photerrhisto(fgroupstruct *fgroup, fieldstruct *reffield,
				double hsn_thresh)
PURPOSE	Plot a 2D histogram of photometric difference between star pairs.
INPUT	Pointer to the field group,
	reference field to the group,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	09/07/2014
 ***/
int	cplot_photerrhisto(fgroupstruct *fgroup, fieldstruct *reffield,
				double hsn_thresh)
  {
   fieldstruct	*field;
   setstruct	*set;
   samplestruct	*samp,*samp1,*samp2;
   double	xscale[NAXIS],xscale_hsn[NAXIS],
		xoffset[NAXIS],xoffset_hsn[NAXIS],
		yscale,yscale_hsn, yoffset,yoffset_hsn, boffset,bscale, dy,
		maxwidth,margin;

   PLFLT	**histo[NAXIS], **histo_hsn[NAXIS], **histot,
		*cuty[NAXIS],*cuty_hsn[NAXIS],
		*clevel,*cutbin,
		xl[2], yl[2],zmax[NAXIS],zmax_hsn[NAXIS],
		cutymax[NAXIS], cutymax_hsn[NAXIS],
		r[2],g[2],b[2],cpoint[2],
		maxlim, cy, z;
   PLINT	lwid;

   char		xlabel[80], ylabel[80], str[80];
   short	sexflagmask;
   unsigned int	imaflagmask;
   int		d,f,i, s,n,
		nsamp, firstflag, instru, ninstru, npinstru,
		ix,iy, nx,ny;

  sexflagmask = (short)prefs.phot_sexflagsmask;
  imaflagmask = prefs.phot_imaflagsmask;
  ninstru = prefs.nphotinstrustr;
  npinstru = 0;
  for (instru=0; instru<ninstru; instru++)
    if (fgroup->chi2_intmag[instru]!=0.0)
      npinstru++;
  if (!npinstru)
    npinstru = 1;

  nx = 1;
  ny = npinstru*fgroup->naxis;

  if (cplot_init(nx, ny , CPLOT_PHOTERROR) == RETURN_ERROR)
    {
    cplot_end(CPLOT_PHOTERROR);
    return RETURN_OK;
    }

  for (d=0; d<fgroup->naxis; d++)
    {
    plAlloc2dGrid(&histo[d], CPLOT_PHOTERRNX, CPLOT_PHOTERRNY);
    plAlloc2dGrid(&histo_hsn[d], CPLOT_PHOTERRNX_HSN, CPLOT_PHOTERRNY_HSN);
    QMALLOC(cuty[d], PLFLT, CPLOT_NADERRHISTBIN);
    QMALLOC(cuty_hsn[d], PLFLT, CPLOT_NADERRHISTBIN);
    xoffset[d] = xoffset_hsn[d] = fgroup->projposmin[d];
    xscale[d] = CPLOT_PHOTERRNX / (fgroup->projposmax[d]-fgroup->projposmin[d]);
    xscale_hsn[d] = CPLOT_PHOTERRNX_HSN
	/ (fgroup->projposmax[d] - fgroup->projposmin[d]);
    }

  QMALLOC(cutbin, PLFLT, CPLOT_NADERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);
  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  for (instru=0; instru<ninstru; instru++)
    {
    maxlim = fgroup->sig_intmagerr[instru]*4;
    if (maxlim<=0.0)
      maxlim = 1.0;
    boffset = -maxlim;
    bscale = CPLOT_NADERRHISTBIN / (2.0*maxlim);
    yoffset = yoffset_hsn = -maxlim;
    yscale = CPLOT_PHOTERRNY/(2.0*maxlim);
    yscale_hsn = CPLOT_PHOTERRNY_HSN/(2.0*maxlim);
/*-- Initialize histograms */
    for (d=0; d<fgroup->naxis; d++)
      {
      histot = histo[d];
      for (ix=0; ix<CPLOT_PHOTERRNX; ix++)
        {
        memset(histot[ix], 0, CPLOT_PHOTERRNY*sizeof(PLFLT));
        }
      histot = histo_hsn[d];
      for (ix=0; ix<CPLOT_PHOTERRNX_HSN; ix++)
        {
        memset(histot[ix], 0, CPLOT_PHOTERRNY_HSN*sizeof(PLFLT));
        }
      memset(cuty[d], 0, CPLOT_NADERRHISTBIN*sizeof(PLFLT));
      memset(cuty_hsn[d], 0, CPLOT_NADERRHISTBIN*sizeof(PLFLT));
      cutymax[d] = cutymax_hsn[d] = zmax[d] = zmax_hsn[d] = 0.0;
      }
    for (f=0; f<fgroup->nfield; f++)
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
            samp2 = NULL;
/*---------- Look for a counterpart from the right photometric instrument */
            for (samp = samp1; samp && samp->set->field->photomlabel>=0;
		samp=samp->prevsamp)
              if (samp->set->field->photomlabel == instru
		&& samp->flux > 0.0
		&& !(samp->sexflags & sexflagmask)
		&& !(samp->imaflags & imaflagmask))
                {
                samp2 = samp;
                break;
                }
/*---------- No right instrument found: skip */
            if (!samp2)
              continue;
            for (samp2=samp2->prevsamp;
		samp2 && samp2->set->field->photomlabel>=0;
		samp2=samp2->prevsamp)
	      {
/*------------ Don't bother if field is a different instru or photometric ref */
/*------------ or the flux is negative or the source is cropped/saturated */
              if (samp2->set->field->photomlabel != instru
		|| samp2->flux <= 0.0
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
                continue;
              dy = samp2->mag - samp->mag;
              for (d=0; d<fgroup->naxis; d++)
                {
                ix = (int)((samp2->projpos[d]-xoffset[d])*xscale[d]);
                iy = (int)((dy - yoffset)*yscale);
                if (ix>=0 && ix<CPLOT_PHOTERRNX
			&& iy>=0 && iy<CPLOT_PHOTERRNY)
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
              if (samp->flux/samp->fluxerr >= hsn_thresh)
                {
                for (d=0; d<fgroup->naxis; d++)
                  {
                  ix = (int)((samp2->projpos[d]-xoffset_hsn[d])
			*xscale_hsn[d]);
                  iy = (int)((dy - yoffset_hsn)*yscale_hsn);
                  if (ix>=0 && ix<CPLOT_PHOTERRNX_HSN
			&& iy>=0 && iy<CPLOT_PHOTERRNY_HSN)
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


/*-- Now plot! */
    lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
    for (i=0; i<CPLOT_NADERRHISTBIN; i++)
      cutbin[i] = boffset+(i+0.5)/bscale;
    firstflag = 1;
    yl[0] = yl[1] = 0.0;
    plschr(0.0,0.5);
    for (d=0; d<fgroup->naxis; d++)
      {
      maxwidth = fgroup->projposmax[d]-fgroup->projposmin[d];
      margin = 0.1*maxwidth;
/*---- Adjust histogram to fit in the displayed box */
      for (i=0; i<CPLOT_NADERRHISTBIN; i++)
        {
        cuty[d][i] = fgroup->projposmin[d] - margin
			+ cuty[d][i]/cutymax[d] * 0.9*margin;
        cuty_hsn[d][i] = fgroup->projposmin[d] - margin
			+ cuty_hsn[d][i]/cutymax_hsn[d] * 0.9*margin;
        }
      CPLOT_PLWID(lwid);
      plenv(fgroup->projposmin[d] - margin, fgroup->projposmax[d],
		-maxlim, maxlim, 0, -2); 
/* Use a non-linear shade level distribution */
      if (zmax[d]>=1.0)
        {
        for (i=0; i<CPLOT_NSHADES; i++)
          clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax[d]+0.5;
        r[0] = 0.98; g[0] = 0.98; b[0] = 1.0;
        r[1] = 0.3; g[1] = 0.3; b[1] = 0.4;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
        plshades((const PLFLT **)histo[d],
		CPLOT_PHOTERRNX, CPLOT_PHOTERRNY, NULL,
		fgroup->projposmin[d], fgroup->projposmax[d], -maxlim, maxlim,
		clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
        }
      else
        {
        plcol0(1);
        plptex((fgroup->projposmin[d] - margin + fgroup->projposmax[d])/2.0,
		maxlim/2.0, 1.0, 0.0, 0.5, "No overlapping detections!");
        }
      if (zmax_hsn[d]>=1.0)
        {
        r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
        r[1] = 0.8; g[1] = 0.8; b[1] = 0.8;
        plscmap1l(1, 2, cpoint, r, g, b, NULL);
       plimage((const PLFLT **)histo_hsn[d],
		CPLOT_PHOTERRNX_HSN, CPLOT_PHOTERRNY_HSN,
		fgroup->projposmin[d], fgroup->projposmax[d], -maxlim, maxlim,
		0.5, zmax_hsn[d],
		fgroup->projposmin[d], fgroup->projposmax[d], -maxlim, maxlim);
        }
      plscolbg(255,255,255);	/* Force the background colour to white */
      plscol0(15, 0,0,0);	/* Force the foreground colour to black */
      plcol0(9);
      CPLOT_PLWID(2*lwid);
      plline(CPLOT_NADERRHISTBIN, cuty[d], cutbin);
      plcol0(7);
      plline(CPLOT_NADERRHISTBIN, cuty_hsn[d], cutbin);
      CPLOT_PLWID(lwid);
      plcol0(15);
      xl[0] = fgroup->projposmin[d] - margin;
      xl[1] = fgroup->projposmax[d];
      pllsty(2);
      plline(2, xl, yl);
      pllsty(1);
      plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
      sprintf(xlabel, "AXIS%d [pixels]", d+1);
      sprintf(ylabel, "#gDmag");
      if (firstflag)
        {
        sprintf(str, "Group ##%d / Instrument P%d: Internal photometric error",
		fgroup->no, instru+1);
        pllab(xlabel, ylabel, str);
        }
      else
        pllab(xlabel, ylabel, "");
      firstflag = 0;
      }
    }

/*-- Free array of points */
  free(clevel);
  free(cutbin);
  for (d=0; d<fgroup->naxis; d++)
    {
    plFree2dGrid(histo[d], CPLOT_PHOTERRNX, CPLOT_PHOTERRNY); 
    plFree2dGrid(histo_hsn[d], CPLOT_PHOTERRNX_HSN, CPLOT_PHOTERRNY_HSN); 
    free(cuty[d]);
    free(cuty_hsn[d]);
    }

  plend();

  cplot_photerrhisto(fgroup, reffield, hsn_thresh);	/* Recursive stuff */

  return RETURN_OK;
  }

/****** cplot_photerrhistomag *************************************************
PROTO	int cplot_photerrhistomag(fgroupstruct *fgroup, fieldstruct *reffield,
				double hsn_thresh)
PURPOSE	Plot a 2D histogram of photometric difference between star pairs versus
	magnitude.
INPUT	Pointer to the field group,
	reference field to the group,
	S/N threshold for the high-S/N sample.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP), C. MARMO (IAP)
VERSION	09/07/2014
 ***/
int	cplot_photerrhistomag(fgroupstruct *fgroup, fieldstruct *reffield,
				double hsn_thresh)
  {
   fieldstruct	*field;
   setstruct	*set;
   samplestruct	*samp,*samp1,*samp2;
   double	xscale,xscale_hsn,
		xoffset,xoffset_hsn,
		yscale,yscale_hsn, yoffset,yoffset_hsn, boffset,bscale, dy,
		maxwidth,margin, magmax;
   PLFLT	**histo, **histo_hsn, **histot,
		*cuty,*cuty_hsn,
		*clevel,*cutbin,
		xl[2], yl[2],zmax,zmax_hsn,
		cutymax, cutymax_hsn,
		r[2],g[2],b[2],cpoint[2],
		maxlim, cy, z;
   PLINT	lwid;
   char		xlabel[80], ylabel[80], str[80];
   short	sexflagmask;
   unsigned int	imaflagmask;
   int		f,i, s,n,
		nsamp, firstflag, instru, ninstru, npinstru,
		ix,iy, nx,ny;

  sexflagmask = (short)prefs.phot_sexflagsmask;
  imaflagmask = prefs.phot_imaflagsmask;
  ninstru = prefs.nphotinstrustr;
  npinstru = 0;
  for (instru=0; instru<ninstru; instru++)
    if (fgroup->chi2_intmag[instru]!=0.0)
      npinstru++;
  if (!npinstru)
    npinstru = 1;

  nx = 1;
  ny = npinstru;
  if (cplot_init(nx, ny , CPLOT_PHOTERRORVSMAG) == RETURN_ERROR)
    {
    cplot_end(CPLOT_PHOTERRORVSMAG);
    return RETURN_OK;
    }

  plAlloc2dGrid(&histo, CPLOT_PHOTERRNX, CPLOT_PHOTERRNY);
  plAlloc2dGrid(&histo_hsn, CPLOT_PHOTERRNX_HSN, CPLOT_PHOTERRNY_HSN);
  QMALLOC(cuty, PLFLT, CPLOT_NADERRHISTBIN);
  QMALLOC(cuty_hsn, PLFLT, CPLOT_NADERRHISTBIN);

  QMALLOC(cutbin, PLFLT, CPLOT_NADERRHISTBIN);
  QMALLOC(clevel, PLFLT, CPLOT_NSHADES);
  plscmap1n(256);
  cpoint[0] = 0.0;
  cpoint[1] = 1.0;

  for (instru=0; instru<ninstru; instru++)
    {
/* First, find limits and bin width in magnitude */
    xoffset = BIG;
    magmax = -BIG;
    for (f=0; f<fgroup->nfield; f++)
      {
      field = fgroup->field[f];
/*---- Skip field if not observed with the right photometric instrument */
      if (field->photomlabel == instru)
        {
        for (s=0; s<field->nset; s++)
          {
          set = field->set[s];
          nsamp = set->nsample;
          samp = set->sample;
          for (n=nsamp; n--; samp++)
            {
            if (xoffset > samp->mag)
              xoffset = samp->mag;
            if (magmax< samp->mag)
              magmax = samp->mag;
            }
          }
        }
      }

    xoffset_hsn = xoffset;
    xscale = CPLOT_PHOTERRNX / (magmax-xoffset);
    xscale_hsn = CPLOT_PHOTERRNX_HSN / (magmax-xoffset);

    maxlim = fgroup->sig_intmagerr[instru]*4;
    if (maxlim<=0.0)
      maxlim = 1.0;
    boffset = -maxlim;
    bscale = CPLOT_NADERRHISTBIN / (2.0*maxlim);
    yoffset = yoffset_hsn = -maxlim;
    yscale = CPLOT_PHOTERRNY/(2.0*maxlim);
    yscale_hsn = CPLOT_PHOTERRNY_HSN/(2.0*maxlim);
/*-- Initialize histograms */
    histot = histo;
    for (ix=0; ix<CPLOT_PHOTERRNX; ix++)
      {
      memset(histot[ix], 0, CPLOT_PHOTERRNY*sizeof(PLFLT));
      }
    histot = histo_hsn;
    for (ix=0; ix<CPLOT_PHOTERRNX_HSN; ix++)
      {
      memset(histot[ix], 0, CPLOT_PHOTERRNY_HSN*sizeof(PLFLT));
      }
    memset(cuty, 0, CPLOT_NADERRHISTBIN*sizeof(PLFLT));
    memset(cuty_hsn, 0, CPLOT_NADERRHISTBIN*sizeof(PLFLT));
    cutymax = cutymax_hsn = zmax = zmax_hsn = 0.0;

    for (f=0; f<fgroup->nfield; f++)
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
            samp2 = NULL;
/*---------- Look for a counterpart from the right photometric instrument */
            for (samp = samp1; samp && samp->set->field->photomlabel>=0;
		samp=samp->prevsamp)
              if (samp->set->field->photomlabel == instru
		&& samp->flux > 0.0
		&& !(samp->sexflags & sexflagmask)
		&& !(samp->imaflags & imaflagmask))
                {
                samp2 = samp;
                break;
                }
/*---------- No right instrument found: skip */
            if (!samp2)
              continue;
            for (samp2=samp2->prevsamp;
		samp2 && samp2->set->field->photomlabel>=0;
		samp2=samp2->prevsamp)
	      {
/*------------ Don't bother if field is a different instru or photometric ref */
/*------------ or the flux is negative */
              if (samp2->set->field->photomlabel != instru
		|| samp2->flux <= 0.0
		|| (samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
               continue;
              dy = samp2->mag - samp->mag;
              ix = (int)((samp->mag-xoffset)*xscale);
              iy = (int)((dy - yoffset)*yscale);
              if (ix>=0 && ix<CPLOT_PHOTERRNX
			&& iy>=0 && iy<CPLOT_PHOTERRNY)
                {
                z = (histo[ix][iy] += 1.0);
                if (z>zmax)
                  zmax = z;
                }
              iy = (int)((dy - boffset)*bscale);
              if (iy>=0 && iy<CPLOT_NADERRHISTBIN)
                {
                cy = (cuty[iy] += 1.0);
                if (cy>cutymax)
                  cutymax = cy;
                }
              if (samp->flux/samp->fluxerr >= hsn_thresh)
                {
                ix = (int)((samp->mag-xoffset_hsn)
			*xscale_hsn);
                iy = (int)((dy - yoffset_hsn)*yscale_hsn);
                if (ix>=0 && ix<CPLOT_PHOTERRNX_HSN
			&& iy>=0 && iy<CPLOT_PHOTERRNY_HSN)
                  {
                  z = (histo_hsn[ix][iy] += 1.0);
                  if (z>zmax_hsn)
                    zmax_hsn = z;
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

/*-- Now plot! */
    lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
    for (i=0; i<CPLOT_NADERRHISTBIN; i++)
      cutbin[i] = boffset+(i+0.5)/bscale;
    firstflag = 1;
    yl[0] = yl[1] = 0.0;
    plschr(0.0,0.5);
    maxwidth = magmax-xoffset;
    margin = 0.1*maxwidth;
/*---- Adjust histogram to fit in the displayed box */
    for (i=0; i<CPLOT_NADERRHISTBIN; i++)
      {
      cuty[i] = xoffset - margin + cuty[i]/cutymax * 0.9*margin;
      cuty_hsn[i] = xoffset - margin + cuty_hsn[i]/cutymax_hsn * 0.9*margin;
      }
    CPLOT_PLWID(lwid);
    plenv(xoffset - margin, magmax, -maxlim, maxlim, 0, -2);
/* Use a non-linear shade level distribution */
    if (zmax>=1.0)
      {
      for (i=0; i<CPLOT_NSHADES; i++)
        clevel[i] = pow(i/(CPLOT_NSHADES-1.0),1.8)*zmax+0.5;
      r[0] = 0.98; g[0] = 0.98; b[0] = 1.0;
      r[1] = 0.3; g[1] = 0.3; b[1] = 0.4;
      plscmap1l(1, 2, cpoint, r, g, b, NULL);
      plshades((const PLFLT **)histo, CPLOT_PHOTERRNX, CPLOT_PHOTERRNY, NULL,
		xoffset, magmax, -maxlim, maxlim,
		clevel, CPLOT_NSHADES, 1, 0, 0, plfill, 1, NULL, NULL);
      }
    else
      {
      plcol0(1);
      plptex((xoffset - margin + magmax)/2.0,
		maxlim/2.0, 1.0, 0.0, 0.5, "No overlapping detections!");
      }
    if (zmax_hsn>=1.0)
      {
      r[0] = 0.0; g[0] = 0.0; b[0] = 0.0;
      r[1] = 0.8; g[1] = 0.8; b[1] = 0.8;
      plscmap1l(1, 2, cpoint, r, g, b, NULL);
     plimage((const PLFLT **)histo_hsn,
		CPLOT_PHOTERRNX_HSN, CPLOT_PHOTERRNY_HSN,
		xoffset, magmax, -maxlim, maxlim,
		0.5, zmax_hsn,
		xoffset, magmax, -maxlim, maxlim);
      }
    plscolbg(255,255,255);	/* Force the background colour to white */
    plscol0(15, 0,0,0);	/* Force the foreground colour to black */
    plcol0(9);
    CPLOT_PLWID(2*lwid);
    plline(CPLOT_NADERRHISTBIN, cuty, cutbin);
    plcol0(7);
    plline(CPLOT_NADERRHISTBIN, cuty_hsn, cutbin);
    CPLOT_PLWID(lwid);
    plcol0(15);
    xl[0] = xoffset - margin;
    xl[1] = magmax;
    pllsty(2);
    plline(2, xl, yl);
    pllsty(1);
    plbox("bcnst", 0.0, 0.0, "bcnst", 0.0, 0.0);
    sprintf(xlabel, "mag");
    sprintf(ylabel, "#gDmag");
    sprintf(str, "Group ##%d / Instrument P%d: Internal photometric error",
	fgroup->no, instru+1);
    pllab(xlabel, ylabel, str);
    }

/*-- Free array of points */
  free(clevel);
  free(cutbin);
  plFree2dGrid(histo, CPLOT_PHOTERRNX, CPLOT_PHOTERRNY); 
  plFree2dGrid(histo_hsn, CPLOT_PHOTERRNX_HSN, CPLOT_PHOTERRNY_HSN); 
  free(cuty);
  free(cuty_hsn);

  plend();

  cplot_photerrhistomag(fgroup, reffield, hsn_thresh);	/* Recursive stuff */

  return RETURN_OK;
  }


/****** cplot_photzp3d *******************************************************
PROTO	int cplot_photzp3d(fgroupstruct *fgroup)
PURPOSE	Plot photometric zero-point corrections in 3D.
INPUT	Pointer to the field group.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	crossid_fgroup() must have been run on all groups first.
AUTHOR	E. Bertin (IAP)
VERSION	09/07/2014
 ***/
int	cplot_photzp3d(fgroupstruct *fgroup)
  {
   struct focplanestruct *focplane, *focplanet;
   fieldstruct	**fields,
		*field;
   wcsstruct	*wcsin,*wcsout;
   PLFLT	cpoint[3], r[3],g[3],b[3], xl[3], yl[3], zl[3],
		dx,dy,dz, lim, xmin,xmax, ymin,ymax, zmin,zmax, x1,y1, col;
   PLINT	lwid;
   double	rawpos[NAXIS],rawpos2[NAXIS], wcspos[NAXIS],wcspos2[NAXIS],
		raw,rawmax, dd,ddmax;
   char		str[80];
   int		f,i,imax,s,stxt, npointmax,
		instru, ninstru, npinstru, nx,ny, nset, lng,lat;

  ninstru = prefs.nphotinstrustr;
  npinstru = 0;
  for (instru=0; instru<ninstru; instru++)
    if (fgroup->chi2_intmag[instru]!=0.0)
      npinstru++;
  if (!npinstru)
    npinstru = 1;

  ny = (int)(sqrt(npinstru+1.0));
  nx = 1+(npinstru-1)/ny;

  if (cplot_init(nx, ny , CPLOT_PHOTZP3D) == RETURN_ERROR)
    {
    cplot_end(CPLOT_PHOTZP3D);
    return RETURN_OK;
    }

  plscmap1n(256);
  cpoint[0] = 0.0; r[0] = 1.0; g[0] = 0.0; b[0] = 0.0;
  cpoint[1] = 0.5; r[1] = 1.0; g[1] = 1.0; b[1] = 0.0;
  cpoint[2] = 1.0; r[2] = 0.0; g[2] = 1.0; b[2] = 0.0;
  plscmap1l(1, 3, cpoint, r, g, b, NULL);
  plschr(0.0,0.5);

  wcsout = fgroup->wcs;
  lng = wcsout->lng;
  lat = wcsout->lat;
  fields = fgroup->field;
  for (instru=0; instru<ninstru; instru++)
    {
    if (fgroup->chi2_intmag[instru]==0.0)
      continue;
    npointmax = nset = 0;
    for (f=0; f<fgroup->nfield; f++)
      if (fields[f]->photomlabel == instru)
        {
        nset = fields[f]->nset;
        npointmax++;
        }
    if (!npointmax)
      continue;
    QMALLOC(focplane, struct focplanestruct, npointmax*nset);
    focplanet = focplane;
    zmin = BIG;
    zmax = -BIG;
    for (f=0; f<fgroup->nfield; f++)
      {
      field= fields[f];
      if (field->photomlabel == instru)
        {
        lim = field->dmagzero;
        if (lim<zmin)
          zmin = lim;
        if (lim>zmax)
          zmax = lim;
/*------ Find the set with highest projpos[lng] and lowest projpos[lat] */
        rawmax = -BIG;
        stxt = 0;
        for (s=0; s<field->nset; s++)
          {
          wcsin = field->set[s]->wcs;
          for (i=0; i<wcsin->naxis; i++)
            rawpos2[i] = wcsin->naxisn[i]/2.0;
          raw_to_wcs(wcsin, rawpos2, wcspos2);
          wcspos[lng] = wcspos2[wcsin->lng];
          wcspos[lat] = wcspos2[wcsin->lat];
          wcs_to_raw(wcsout, wcspos, rawpos);
          if ((raw=rawpos[0]-rawpos[1]) > rawmax)
            {
            rawmax = raw;
            stxt = s;
            }
          }
        for (s=0; s<field->nset; s++)
          {
          wcsin = field->set[s]->wcs;
/*-------- Initialize the input coordinates to an "average" value */
          for (i=0; i<wcsin->naxis; i++)
            rawpos2[i] = wcsin->naxisn[i]/2.0;

/*-------- 1st corner */
          rawpos2[wcsin->lng] = 0.0;
          rawpos2[wcsin->lat] = 0.0;
          raw_to_wcs(wcsin, rawpos2, wcspos2);
          wcspos[lng] = wcspos2[wcsin->lng];
          wcspos[lat] = wcspos2[wcsin->lat];
          wcs_to_raw(wcsout, wcspos, rawpos);
          focplanet->x[4] = focplanet->x[0] = rawpos[lng];
          focplanet->y[4] = focplanet->y[0] = rawpos[lat];
/*-------- 2nd corner */
          rawpos2[wcsin->lng] = wcsin->naxisn[wcsin->lng]-1.0;
          raw_to_wcs(wcsin, rawpos2, wcspos2);
          wcspos[lng] = wcspos2[wcsin->lng];
          wcspos[lat] = wcspos2[wcsin->lat];
          wcs_to_raw(wcsout, wcspos, rawpos);
          focplanet->x[1] = rawpos[lng];
          focplanet->y[1] = rawpos[lat];
/*-------- 3rd corner */
          rawpos2[wcsin->lat] = wcsin->naxisn[wcsin->lat]-1.0;
          raw_to_wcs(wcsin, rawpos2, wcspos2);
          wcspos[lng] = wcspos2[wcsin->lng];
          wcspos[lat] = wcspos2[wcsin->lat];
          wcs_to_raw(wcsout, wcspos, rawpos);
          focplanet->x[2] = rawpos[lng];
          focplanet->y[2] = rawpos[lat];
/*-------- Last corner */
          rawpos2[wcsin->lng] = 0.0;
          raw_to_wcs(wcsin, rawpos2, wcspos2);
          wcspos[lng] = wcspos2[wcsin->lng];
          wcspos[lat] = wcspos2[wcsin->lat];
          wcs_to_raw(wcsout, wcspos, rawpos);
          focplanet->x[3] = rawpos[lng];
          focplanet->y[3] = rawpos[lat];
          focplanet->z[0] = focplanet->z[1] = focplanet->z[2]
		= focplanet->z[3] = focplanet->z[4] = lim;
          focplanet->colour = field->cplot_colour;
          if (field->photomflag==1 && field->cplot_colour == 15)
            focplanet->colour = 9;
          focplanet->str = (s==stxt? field->rfilename : NULL);
          focplanet++;
          }
        }
      }
    dz = zmax - zmin;
    zmin -= dz*0.2;
    zmax += dz*0.05;
/*-- Sort fields by increasing z */
    qsort(focplane, npointmax*nset, sizeof(struct focplanestruct), comp_focz);
/*-- Now plot! */
    yl[0] = yl[1] = 0.0;
    plcol0(15);
    lwid = plotaaflag? ((CPLOT_AAFAC+1)/2) : 1;
    CPLOT_PLWID(lwid);
    pladv(0);
    plvpor(0.0, 1.0, 0.0, 0.95);
    plwind(-0.75, 0.75, -0.6, 1.1);
    dx = fgroup->projposmax[lng] - fgroup->projposmin[lng];
    dy = fgroup->projposmax[lat] - fgroup->projposmin[lat];
    if (dy>dx)
    dx = dy;
    xmin = 0.5*(fgroup->projposmin[lng]+fgroup->projposmax[lng]) - 0.55*dx;
    xmax = 0.5*(fgroup->projposmin[lng]+fgroup->projposmax[lng]) + 0.55*dx;
    ymin = 0.5*(fgroup->projposmin[lat]+fgroup->projposmax[lat]) - 0.55*dx;
    ymax = 0.5*(fgroup->projposmin[lat]+fgroup->projposmax[lat]) + 0.55*dx;
    plw3d(1.0, 1.0, 1.0, xmin, xmax, ymin, ymax, zmin, zmax, 50.0, -20.0);
    plbox3("bfnstu", "AXIS1", 0.0, 0,
	"bfnstu", "AXIS2", 0.0, 0,
	"bcdmnstuv", "#gDZP [mag]", 0.0, 0);
    xl[0] = xl[1] = xmin;
    xl[2] = xmax;
    yl[0] = ymin;
    yl[1] = yl[2] = ymax;
    zl[0] = zl[1] = zl[2];
    pllsty(2);
    plcol0(15);
    plline3(3, xl, yl, zl);
    pllsty(1);
    plpsty(0);
    focplanet = focplane;
    for (s=npointmax*nset; s--; focplanet++)
      {
      col = (focplanet->z[0] - zmin)/(zmax - zmin);
      plcol1(col);
      plfill3(5, focplanet->x, focplanet->y, focplanet->z);
      if (focplanet->colour!=15)
        {
        CPLOT_PLWID(2*lwid);
        plcol0(focplanet->colour);
        }
      else
        plcol0(15);
      plline3(5, focplanet->x, focplanet->y, focplanet->z); 
      CPLOT_PLWID(lwid);
/*---- Find the point with highest projpos[lng] and lowest projpos[lat] */
      if (focplanet->str)
        {
        ddmax = -BIG;
        imax = 0;
        for (i=0; i<5; i++)
          if ((dd = focplanet->x[i]-focplanet->y[i])> ddmax)
            {
            ddmax = dd;
            imax = i;
            }
/* 
       x1=plP_w3wcx(focplanet->x[imax],focplanet->y[imax],focplanet->z[imax]);
        y1=plP_w3wcy(focplanet->x[imax],focplanet->y[imax],focplanet->z[imax]);
        plptex(x1, y1-0.02, 1.0, -0.26, 1.1, focplanet->str);
*/
        }
      plcol0(15);
      }
    xl[0] = xmin;
    xl[1] = xl[2] = xmax;
    yl[0] = yl[1] = ymin;
    yl[2] = ymax;
    zl[0] = zl[1] = zl[2];
    pllsty(2);
    plline3(3, xl, yl, zl);
    pllsty(1);
    sprintf(str, "Group ##%d / Instrument P%d: Zero-point corrections",
		fgroup->no, instru+1);
    pllab("", "", str);

/*-- Free array of points */
    free(focplane);
    }

  plend();

  cplot_photzp3d(fgroup);		/* Recursive stuff */

  return RETURN_OK;
  }


/****** comp_focz ************************************************************
PROTO   int (*comp_focz(const void *focplane1, const void *focplane2))
PURPOSE Provide a focplane z comparison function for qsort()
INPUT   pointer to 1st focplane,
	pointer to 2nd focplane.
OUTPUT  <0 if z1<z2, >0 if z1>z2, 0 otherwise .
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 25/01/2005
*/
int	comp_focz(const void *focplane1, const void *focplane2)
  {
   PLFLT	dz;

  dz = ((struct focplanestruct *)focplane1)->z[0]
	- ((struct focplanestruct *)focplane2)->z[0];

  return dz>0.0 ? 1: (dz<0.0? -1 : 0);
  }

