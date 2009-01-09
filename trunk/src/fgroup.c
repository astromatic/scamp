 /*
				fgroup.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Handle groups of fields.
*
*	Last modify:	30/01/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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
PROTO   fgroupstruct	**group_fields(fieldstruct **field, int nfield,
			int *nfgroup)
PURPOSE	Group fields which are sufficiently close on the sky to share a common
	projection.
INPUT   Pointer to field structure pointers,
	number of fields,
	Pointer to the total number of groups found (filled by group_field()).
OUTPUT  Pointer to the array of groups found.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 05/07/2004
*/
fgroupstruct	**group_fields(fieldstruct **field, int nfield, int *nfgroup)
  {
   fgroupstruct	**fgroup;
   int		*gflag,
		lng,lat, naxis, i, f1,f2, g,g2,g3,testflag, ngroup;
   double	*pos1, *pos2,
		dist,maxdist,dx;

  if (!nfield)
    return 0;

  maxdist = prefs.group_radius;
/* Set angular coordinates (we assume there are the same in all fields */
  lng = field[0]->set[0]->wcs->lng;
  lat = field[0]->set[0]->wcs->lat;
  naxis = field[0]->set[0]->wcs->naxis;
  ngroup = 0;
/* Allocate memory */
  QMALLOC(gflag, int, nfield);
  QMALLOC(fgroup, fgroupstruct *, nfield);
  for (f1=0; f1<nfield; f1++)
    {
    pos1 = field[f1]->meanwcspos;
    testflag = 1;
    if (ngroup)
      memset(gflag, 0, ngroup*sizeof(int));
    for (g=0; g<ngroup; g++)
      for (f2=0; f2<fgroup[g]->nfield && !gflag[g]; f2++)
        {
        pos2 = fgroup[g]->field[f2]->meanwcspos;
        if (lng != lat)
	  {
          dist = sin(pos1[lat]*DEG)*sin(pos2[lat]*DEG)
		+cos(pos1[lat]*DEG)*cos(pos2[lat]*DEG)
			*cos((pos2[lng]-pos1[lng])*DEG);
          dist = dist>-1.0? (dist<1.0 ? acos(dist)/DEG : 0.0) : 180.0;
          }
        else
	  {
          dist = 0.0;
          for (i=0; i<naxis; i++)
	    {
            dx = pos2[i] - pos1[i];
            dist += dx*dx;
            }
          dist = sqrt(dist);
          }
/*------ Check whether it is close enough */
        if (dist<maxdist)
	  {
          testflag = 0;
          gflag[g] = 1;
	  }
	}

    if (testflag)
/*---- field too far: Create a new group */
      {
      fgroup[ngroup] = new_fgroup();
      addfield_fgroup(fgroup[ngroup], field[f1]);
      ngroup++;
      }
    else
      {
/*---- Add to an existing group */
      for (g=0; !gflag[g] && g<ngroup; g++);        
      g2 = g;
      addfield_fgroup(fgroup[g], field[f1]);
/*---- Check that the newcomer doesn't link groups together */
      for (g=g2+1; g<ngroup; g++)
        if (gflag[g])
	  {
/*-------- Fusion this group with the first one */          
          addfgroup_fgroup(fgroup[g], fgroup[g2]);
          end_fgroup(fgroup[g]);
          for (g3=g+1; g3<ngroup; g3++)
            fgroup[g3-1] = fgroup[g3];
          ngroup--;
          }
      }
    }
  free(gflag);
  QREALLOC(fgroup, fgroupstruct *, ngroup);

/* Number groups */
  for (g=0; g<ngroup; g++)
    fgroup[g]->no = g+1;

/* Update astrometric stuff */
  for (g=0; g<ngroup; g++)
    locate_fgroup(fgroup[g]);

  *nfgroup = ngroup;
  return fgroup;
  }


/****** locate_fgroup ********************************************************
PROTO   void locate_fgroup(fgroupstruct *fgroup)
PURPOSE Compute statiscal informations gathered for a field group.
INPUT   Field group pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/09/2004
*/
void	locate_fgroup(fgroupstruct *fgroup)
  {
   fieldstruct	**field;
   wcsstruct	*wcs;
   char		*ctype[NAXIS];
   double	*scale[NAXIS],*scalet[NAXIS],
		crpix[NAXIS], cdelt[NAXIS],
		*wcsmean, *wcspos,
		cosalpha,sinalpha,sindelta, dist, maxradius, dscale;
   int		naxisn[NAXIS],
		lat,lng, naxis,nfield, f,i;

  if (!fgroup->nfield)
    return;

/* Set angular coordinates (we assume there are the same in all fields) */
  field = fgroup->field;
  wcs = field[0]->set[0]->wcs;
  lng = fgroup->lng = field[0]->lng;
  lat = fgroup->lat = field[0]->lat;
  naxis = fgroup->naxis = field[0]->naxis;
  nfield = fgroup->nfield;
  wcsmean = fgroup->meanwcspos;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(scale[i], double, nfield);
    QMALLOC(ctype[i], char, 16);
    scalet[i] = scale[i];
    wcsmean[i] = 0.0;
    }

/* Compute average angular coordinates */
  cosalpha = sinalpha = sindelta = 0.0;
  for (f=0; f<nfield; f++)
    {
    wcspos = field[f]->meanwcspos;
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
      *(scalet[i]++) = field[f]->meanwcsscale[i];
      }
    }

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
    dist = wcs_dist(wcs, field[f]->meanwcspos, fgroup->meanwcspos)
		+ field[f]->maxradius;
    if (dist>maxradius)
      maxradius = dist;
    }

  fgroup->maxradius = maxradius;
  for (i=0; i<naxis; i++)
    {
    cdelt[i] = fgroup->meanwcsscale[i];
    if (lat!=lng && i==lng)
      cdelt[i] = -cdelt[i];	/* To look like a sky image */
    if (lat!=lng && (i==lng || i==lat))
      {
      naxisn[i] = (int)(2.0*fgroup->maxradius / fgroup->meanwcsscale[i]+1);
      sprintf(ctype[i], "%5.5sSTG", field[0]->set[0]->wcs->ctype[i]);
      }
    else
      {
      naxisn[i] = (lat==lng)?
		(int)(2.0*fgroup->maxradius / fgroup->meanwcsscale[i] + 1.0)
		: field[0]->set[0]->wcs->naxisn[i];
      strcpy(ctype[i], field[0]->set[0]->wcs->ctype[i]);
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
VERSION 10/09/2004
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
VERSION 30/01/2006
*/
void	end_fgroup(fgroupstruct *fgroup)
  {
   int	i;

  if (fgroup->field)
    free(fgroup->field);
  if (fgroup->sig_intmagerr)
    free(fgroup->sig_intmagerr);
  if (fgroup->chi2_intmag)
    free(fgroup->chi2_intmag);
  if (fgroup->nintmagmatch)
    free(fgroup->nintmagmatch);
  if (fgroup->sig_intmagerr_hsn)
    free(fgroup->sig_intmagerr_hsn);
  if (fgroup->chi2_intmag_hsn)
    free(fgroup->chi2_intmag_hsn);
  if (fgroup->nintmagmatch_hsn)
    free(fgroup->nintmagmatch_hsn);
  if (fgroup->sig_refmagerr)
    free(fgroup->sig_refmagerr);
  if (fgroup->chi2_refmag)
    free(fgroup->chi2_refmag);
  if (fgroup->nrefmagmatch)
    free(fgroup->nrefmagmatch);
  if (fgroup->sig_refmagerr_hsn)
    free(fgroup->sig_refmagerr_hsn);
  if (fgroup->chi2_refmag_hsn)
    free(fgroup->chi2_refmag_hsn);
  if (fgroup->nrefmagmatch_hsn)
    free(fgroup->nrefmagmatch_hsn);
  for (i=0; i<NAXIS; i++)
    {
    if (fgroup->intcolshiftscale[i])
      free(fgroup->intcolshiftscale[i]);
    if (fgroup->intcolshiftzero[i])
      free(fgroup->intcolshiftzero[i]);
    if (fgroup->colshiftscale[i])
      free(fgroup->colshiftscale[i]);
    if (fgroup->colshiftzero[i])
      free(fgroup->colshiftzero[i]);
    if (fgroup->refcolshiftscale[i])
      free(fgroup->refcolshiftscale[i]);
    if (fgroup->refcolshiftzero[i])
      free(fgroup->refcolshiftzero[i]);
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
VERSION	25/09/2004
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
	"\nGroup %2d: %d field%s at %s %s with radius %.4g'\n",
        i+1,
        fgroup->nfield, fgroup->nfield>1 ? "s":"",
	      degtosexal(fgroup->meanwcspos[fgroup->lng], str1),
	      degtosexde(fgroup->meanwcspos[fgroup->lat], str2),
        fgroup->maxradius*DEG/ARCMIN);
      }
    else
      {
      QPRINTF(OUTPUT,
	"\nGroup %2d: %d field%s at %.3g %.3g with radius %.4g\n",
        i+1,
        fgroup->nfield, fgroup->nfield>1 ? "s":"",
        fgroup->meanwcspos[0],
        fgroup->naxis>1? fgroup->meanwcspos[1] : 0.0,
        fgroup->maxradius);
      }
    QIPRINTF(OUTPUT,
	"                  instruments     center coordinates "
	"    radius    scale");
    fields = fgroup->field;
    for (f=fgroup->nfield; f--;)
      {
      field = *(fields++);
      lng = field->lng;
      lat = field->lat;
      if (lat != lng)
        {
        QPRINTF(OUTPUT,"%-20.20s A%-2d P%-2d %c %s %s %#7.4g' %#7.4g\"\n",
		field->rfilename,
		field->astromlabel+1,field->photomlabel+1,
		field->photomflag? '*': ' ',
		degtosexal(field->meanwcspos[lng], str1),
		degtosexde(field->meanwcspos[lat], str2),
		field->maxradius*DEG/ARCMIN,
		field->meanwcsscale[lng]*DEG/ARCSEC);
        }
      else
        {
        QPRINTF(OUTPUT,
		"%-20.20s A%-2d P%-2d %#+11.3g %#+11.3g %#7.4g  %#7.4g\n",
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
VERSION	11/07/2004
 ***/
void	print_instruinfo(void)
  {
   int		i,l,len;

  QPRINTF(OUTPUT, "\n----- %d %s found for astrometry:\n",
	prefs.nastrinstrustr, prefs.nastrinstrustr>1? "instruments":"instrument");
  for (i=0; i<prefs.nastrinstrustr; i++)
    {
    QPRINTF(OUTPUT, "\nInstrument A%-2d:\n", i+1);
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

