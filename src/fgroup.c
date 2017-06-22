/*
*				fgroup.c
*
* Manage group of fields.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2017 IAP/CNRS/UPMC
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
*	Last modified:		22/06/2017
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "fgroup.h"
#include "field.h"
#include "prefs.h"
#include "samples.h"


/****** group_fields ********************************************************
PROTO   fgroupstruct	**group_fields(fieldstruct **fields, int nfield,
			int *nfgroup)
PURPOSE	Group fields that are likely to overlap.
INPUT   Pointer to field structure pointers,
	number of fields,
	Pointer to the total number of groups found (filled by group_field()).
OUTPUT  Pointer to the array of groups found.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 22/06/2017
*/
fgroupstruct	**group_fields(fieldstruct **fields, int nfield, int *nfgroup) {

   fgroupstruct	**fgroups;
   wcsstruct	*wcs1;
   fieldstruct	*field1, *field2;
   int		*gflag,
		f1,f2, g,g2,g3, s1,s2, nset1, nset2, testflag, ngroup;

  if (!nfield)
    return 0;

  ngroup = 0;
/* Allocate memory */
  QMALLOC(gflag, int, nfield);
  QMALLOC(fgroups, fgroupstruct *, nfield);
  for (f1 = 0; f1 < nfield; f1++) {
    if (ngroup)
      memset(gflag, 0, ngroup*sizeof(int));
    testflag = 1;
    field1 = fields[f1];
    nset1 = field1->nset;
    for (s1 = 0; s1 < nset1; s1++) {
      wcs1 = field1->set[s1]->wcs;
      for (g = 0; g < ngroup; g++)
        for (f2 = 0; f2 < fgroups[g]->nfield && !gflag[g]; f2++) {
          field2 = fields[f2];
          nset2 = field2->nset;
          for (s2 = 0; s2 < nset2; s2++) {
            if (frame_wcs(wcs1, field2->set[s2]->wcs)) {
              testflag = 0;
              gflag[g] = 1;
              break;
  	    }
          }
        }
    }

    if (testflag) {
/*---- field too far: Create a new group */
      fgroups[ngroup] = new_fgroup();
      addfield_fgroup(fgroups[ngroup], fields[f1]);
      ngroup++;
    } else {
/*---- Add to an existing group */
      for (g=0; !gflag[g] && g<ngroup; g++);        
      g2 = g;
      addfield_fgroup(fgroups[g], fields[f1]);
/*---- Check that the newcomer doesn't link groups together */
      for (g=g2+1; g<ngroup; g++)
        if (gflag[g]) {
/*-------- Fusion this group with the first one */          
          addfgroup_fgroup(fgroups[g], fgroups[g2]);
          end_fgroup(fgroups[g]);
          for (g3=g+1; g3<ngroup; g3++)
            fgroups[g3-1] = fgroups[g3];
          ngroup--;
        }
    }
  }
  free(gflag);
  QREALLOC(fgroups, fgroupstruct *, ngroup);

/* Number groups */
  for (g=0; g<ngroup; g++)
    fgroups[g]->no = g+1;

/* Update astrometric stuff */
  for (g=0; g<ngroup; g++)
    locate_fgroup(fgroups[g]);

  *nfgroup = ngroup;
  return fgroups;
}


/****** locate_fgroup ********************************************************
PROTO   void locate_fgroup(fgroupstruct *fgroup)
PURPOSE Compute statiscal informations gathered for a field group.
INPUT   Field group pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 26/07/2012
*/
void	locate_fgroup(fgroupstruct *fgroup)
  {
   fieldstruct	**fields,
		*field;
   wcsstruct	*wcs;
   char		*ctype[NAXIS];
   double	*scale[NAXIS],*scalet[NAXIS],
		crpix[NAXIS], cdelt[NAXIS],
		*wcsmean, *wcspos,
		cosalpha,sinalpha,sindelta, dist, maxradius, dscale,
		epoch,epochmin,epochmax;
   int		naxisn[NAXIS],
		f,i, lat,lng, naxis,nfield, nepoch;

  if (!fgroup->nfield)
    return;

/* Set angular coordinates (we assume there are the same in all fields) */
  fields = fgroup->field;
  field = fields[0];
  wcs = field->set[0]->wcs;
  lng = fgroup->lng = field->lng;
  lat = fgroup->lat = field->lat;
  naxis = fgroup->naxis = field->naxis;
  nfield = fgroup->nfield;
  wcsmean = fgroup->meanwcspos;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(scale[i], double, nfield);
    QMALLOC(ctype[i], char, 16);
    scalet[i] = scale[i];
    wcsmean[i] = 0.0;
    }

/* Compute min/max/mean epoch and average coordinates for the whole group */
  epoch = 0.0;
  epochmax = -(epochmin = BIG);
  nepoch = 0;
  cosalpha = sinalpha = sindelta = 0.0;
  for (f=0; f<nfield; f++)
    {
    field = fields[f];
    wcspos = field->meanwcspos;
    if (lng != lat)
      {
      cosalpha += cos(wcspos[lng]*DEG);
      sinalpha += sin(wcspos[lng]*DEG);
      sindelta += sin(wcspos[lat]*DEG);
      }
    for (i=0; i<naxis; i++)
      {
      if (lat==lng || (i!=lng && i!=lat))
        wcsmean[i] += wcspos[i];
      *(scalet[i]++) = field->meanwcsscale[i];
      }
    if (field->epochmin != 0.0 && field->epochmin < epochmin)
      epochmin = field->epochmin;
    if (field->epochmax != 0.0 && field->epochmax > epochmax)
      epochmax = field->epochmax;
    if (field->epoch != 0.0)
      {
      epoch += field->epoch;
      nepoch++;
      }
    }

  fgroup->epoch = (nepoch)? epoch / nepoch : 0.0;
  fgroup->epochmin = epochmin < BIG/2? epochmin : 0.0;
  fgroup->epochmax = epochmax > -BIG/2? epochmax : 0.0;

/* Now make the stats on each axis */
  for (i=0; i<naxis; i++)
    {
    if (lat!=lng && (i==lng))
      {
      wcsmean[i] = atan2(sinalpha/nfield, cosalpha/nfield)/DEG;
      wcsmean[i] = fmod(wcsmean[i]+360.0, 360.0);
      }
    else if (lat!=lng && (i==lat))
      wcsmean[i] = asin(sindelta/nfield)/DEG;
    else
      wcsmean[i] /= (double)nfield;
    fgroup->meanwcsscale[i] = dhmedian(scale[i], nfield);
    }

/* Compute (average) angular pixel scale */
  dscale = 0.0;
  for (i=0; i<naxis; i++)
    if (lat==lng || i==lng || i==lat)
      dscale += fgroup->meanwcsscale[i];
  dscale /= (double)naxis;
  for (i=0; i<naxis; i++)
    if (lat==lng || i==lng || i==lat)
      fgroup->meanwcsscale[i] = dscale;

/* Find the projection center closest to the group center */
/* ...and compute the radius of the surveyed region */
  maxradius = 0.0;
  for (f=0; f<nfield; f++)
    {
    field = fields[f];
    dist = wcs_dist(wcs,field->meanwcspos,fgroup->meanwcspos)+ field->maxradius;
    if (dist>maxradius)
      maxradius = dist;
    }

  fgroup->maxradius = maxradius;
  field = fields[0];
  for (i=0; i<naxis; i++)
    {
    cdelt[i] = fgroup->meanwcsscale[i];
    if (lat!=lng && i==lng)
      cdelt[i] = -cdelt[i];	/* To look like a sky image */
    if (lat!=lng && (i==lng || i==lat))
      {
      naxisn[i] = (int)(2.0*fgroup->maxradius / fgroup->meanwcsscale[i]+1);
      sprintf(ctype[i], "%5.5sSTG", field->set[0]->wcs->ctype[i]);
      }
    else
      {
      naxisn[i] = (lat==lng)?
		(int)(2.0*fgroup->maxradius / fgroup->meanwcsscale[i] + 1.0)
		: field->set[0]->wcs->naxisn[i];
      strcpy(ctype[i], field->set[0]->wcs->ctype[i]);
      }
    crpix[i] = (naxisn[i]+1.0) / 2.0;
    }

/* Create new fgroup STG projection */
  fgroup->wcs = create_wcs((char **)ctype, wcsmean, crpix, cdelt,naxisn,naxis);

/* Free memory */
  for (i=0; i<naxis; i++)
    {
    free(ctype[i]);
    free(scale[i]);
    }

  return;
  }


/****** addfield_fgroup ******************************************************
PROTO   void addfield_fgroup(fgroupstruct *fgroup, fieldstruct *field)
PURPOSE Add a field to a field group.
INPUT   Field group pointer,
	Field pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 30/12/2003
*/
void	addfield_fgroup(fgroupstruct *fgroup, fieldstruct *field)
  {
  if (fgroup->nfield)
    {
    QREALLOC(fgroup->field, fieldstruct *, fgroup->nfield+1);
    }
  else
    {
    QMALLOC(fgroup->field, fieldstruct *, 1);
    }

  fgroup->field[fgroup->nfield++] = field;
  field->fgroup = fgroup;

  return;
  }


/****** addfgroup_fgroup ******************************************************
PROTO   void addfgroup_fgroup(fgroupstruct *fgroup, fgroupstruct *fgroupin)
PURPOSE Add a group to a field group.
INPUT   Field group pointer (dest),
	Field group pointer (source).
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 10/12/2003
*/
void	addfgroup_fgroup(fgroupstruct *fgroupin, fgroupstruct *fgroup)
  {
   int	f, nfield;

  if (!fgroupin->nfield)
    return;

  nfield = fgroup->nfield;
  fgroup->nfield += fgroupin->nfield;

  if (nfield)
    {
    QREALLOC(fgroup->field, fieldstruct *, fgroup->nfield);
    }
  else
    {
    QMALLOC(fgroup->field, fieldstruct *, fgroup->nfield);
    }

  for (f=0; f<fgroupin->nfield; f++)
    fgroup->field[nfield+f] = fgroupin->field[f];

  for (f=0; f<fgroupin->nfield; f++)
    fgroup->field[nfield+f]->fgroup = fgroup;

  return;
  }


/****** new_fgroup ***********************************************************
PROTO   fgroupstruct	*new_fgroup(void)
PURPOSE Create a new field group.
INPUT   -.
OUTPUT  Pointer to an allocated field group.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 25/05/2012
*/
fgroupstruct	*new_fgroup(void)
  {
   fgroupstruct	*fgroup;

  QCALLOC(fgroup, fgroupstruct, 1);
  QCALLOC(fgroup->sig_intmagerr, double, prefs.nphotinstrustr);
  QCALLOC(fgroup->chi2_intmag, double, prefs.nphotinstrustr);
  QCALLOC(fgroup->nintmagmatch, int, prefs.nphotinstrustr);
  QCALLOC(fgroup->sig_intmagerr_hsn, double, prefs.nphotinstrustr);
  QCALLOC(fgroup->chi2_intmag_hsn, double, prefs.nphotinstrustr);
  QCALLOC(fgroup->nintmagmatch_hsn, int, prefs.nphotinstrustr);
  QCALLOC(fgroup->sig_refmagerr, double, prefs.nphotinstrustr);
  QCALLOC(fgroup->chi2_refmag, double, prefs.nphotinstrustr);
  QCALLOC(fgroup->nrefmagmatch, int, prefs.nphotinstrustr);
  QCALLOC(fgroup->sig_refmagerr_hsn, double, prefs.nphotinstrustr);
  QCALLOC(fgroup->chi2_refmag_hsn, double, prefs.nphotinstrustr);
  QCALLOC(fgroup->nrefmagmatch_hsn, int, prefs.nphotinstrustr);

  return fgroup;
  }


/****** end_fgroup ***********************************************************
PROTO   void end_fgroup(fgroupstruct *fgroup)
PURPOSE Deallocate field group data.
INPUT   Field group pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 28/01/2013
*/
void	end_fgroup(fgroupstruct *fgroup)
  {
   int	i;

  free(fgroup->field);
  free(fgroup->msample);
  free(fgroup->sig_intmagerr);
  free(fgroup->chi2_intmag);
  free(fgroup->nintmagmatch);
  free(fgroup->sig_intmagerr_hsn);
  free(fgroup->chi2_intmag_hsn);
  free(fgroup->nintmagmatch_hsn);
  free(fgroup->sig_refmagerr);
  free(fgroup->chi2_refmag);
  free(fgroup->nrefmagmatch);
  free(fgroup->sig_refmagerr_hsn);
  free(fgroup->chi2_refmag_hsn);
  free(fgroup->nrefmagmatch_hsn);
  for (i=0; i<NAXIS; i++)
    {
    free(fgroup->intcolshiftscale[i]);
    free(fgroup->intcolshiftzero[i]);
    free(fgroup->colshiftscale[i]);
    free(fgroup->colshiftzero[i]);
    free(fgroup->refcolshiftscale[i]);
    free(fgroup->refcolshiftzero[i]);
/*-- Maps */
    free(fgroup->sig_interr_map[i]);
    }
  if (fgroup->wcs)
    end_wcs(fgroup->wcs);
  free(fgroup);

  return;
  }


/****** print_fgroupinfo *****************************************************
PROTO	void print_fgroupinfo(fgroupstruct **pfgroup, int nfgroup)
PURPOSE	Print info about a series of fgroups.
INPUT	Pointer to an array of fgroup pointers,
	Number of fgroups.
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	03/08/2010
 ***/
void	print_fgroupinfo(fgroupstruct **pfgroup, int nfgroup)
  {
   char		str1[32], str2[32];
   fgroupstruct	*fgroup;
   fieldstruct	**fields,
		*field;
   int		i,f, lng,lat;

  if (!pfgroup || !nfgroup)
    return;


  QPRINTF(OUTPUT, "\n----- %d field %s found:\n",
	nfgroup, nfgroup>1? "groups":"group");

  for (i=0; i<nfgroup; i++)
    {
    fgroup = pfgroup[i];
    if (fgroup->lat != fgroup->lng)
      {
      QPRINTF(OUTPUT,
	"\n Group %2d: %d field%s at %s %s with radius %.4g'\n",
        i+1,
        fgroup->nfield, fgroup->nfield>1 ? "s":"",
	      degtosexal(fgroup->meanwcspos[fgroup->lng], str1),
	      degtosexde(fgroup->meanwcspos[fgroup->lat], str2),
        fgroup->maxradius*DEG/ARCMIN);
      }
    else
      {
      QPRINTF(OUTPUT,
	"\n Group %2d: %d field%s at %.3g %.3g with radius %.4g\n",
        i+1,
        fgroup->nfield, fgroup->nfield>1 ? "s":"",
        fgroup->meanwcspos[0],
        fgroup->naxis>1? fgroup->meanwcspos[1] : 0.0,
        fgroup->maxradius);
      }
    QIPRINTF(OUTPUT,
	"                  instruments  epoch      center coordinates "
	"    radius   scale ");
    fields = fgroup->field;
    for (f=fgroup->nfield; f--;)
      {
      field = *(fields++);
      lng = field->lng;
      lat = field->lat;
      if (lat != lng)
        {
        QPRINTF(OUTPUT,
		"%-20.20s A%-2d P%-2d %c %#6.1f  %s %s %#7.4g' %#7.4g\"\n",
		field->rfilename,
		field->astromlabel+1,field->photomlabel+1,
		field->photomflag? '*': ' ',
		field->epoch,
		degtosexal(field->meanwcspos[lng], str1),
		degtosexde(field->meanwcspos[lat], str2),
		field->maxradius*DEG/ARCMIN,
		field->meanwcsscale[lng]*DEG/ARCSEC);
        }
      else
        {
        QPRINTF(OUTPUT,
	"%-20.20s A%-2d P%-2d           %#+11.4e %#+11.4e %#7.4g %#7.4g\n",
		field->rfilename,
		field->astromlabel+1,field->photomlabel+1,
		field->meanwcspos[0],
		field->naxis>1? field->meanwcspos[1] : 0.0,
		field->maxradius,
		field->meanwcsscale[0]);
        }
      }
    }

  QPRINTF(OUTPUT, "\n");

  return;
  }


/****** print_instruinfo *****************************************************
PROTO	void print_instruinfo(void)
PURPOSE	Print info about astrometric and photometric instruments found
INPUT	-.
OUTPUT	-.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	27/04/2010
 ***/
void	print_instruinfo(void)
  {
   int		i,l,len;

  QPRINTF(OUTPUT, "\n----- %d %s found for astrometry:\n",
	prefs.nastrinstrustr, prefs.nastrinstrustr>1? "instruments":"instrument");
  for (i=0; i<prefs.nastrinstrustr; i++)
    {
    QPRINTF(OUTPUT, "\nInstrument A%-2d:\n", i+1);
    QPRINTF(OUTPUT, "%d extensions\n", prefs.nastrinstruext[i]);
    len = fitsfind(prefs.astrinstrustr[i], "END     ");
    for (l=0; l<len; l++)
      {
      QPRINTF(OUTPUT, "%.80s\n", prefs.astrinstrustr[i]+l*80);
      }
    }


  QPRINTF(OUTPUT, "\n----- %d %s found for photometry:\n",
	prefs.nphotinstrustr, prefs.nphotinstrustr>1? "instruments":"instrument");
  for (i=0; i<prefs.nphotinstrustr; i++)
    {
    QPRINTF(OUTPUT, "\nInstrument P%-2d:\n", i+1);
    len = fitsfind(prefs.photinstrustr[i], "END     ");
    for (l=0; l<len; l++)
      {
      QPRINTF(OUTPUT, "%.80s\n", prefs.photinstrustr[i]+l*80);
      }
    }

  return;
  }

