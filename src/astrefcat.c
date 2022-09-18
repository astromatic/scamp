/*
*				astrefcat.c
*
* Manage astrometric reference catalogs (query and load).
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
*	Last modified:		01/12/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "astrefcat.h"
#include "field.h"
#include "key.h"
#include "prefs.h"
#include "samples.h"
#include "url.h"

samplestruct	refsample;
keystruct       refkey[] = {
  {"X_WORLD", "Barycenter position along world x axis",
        &refsample.wcspos[0], H_FLOAT, T_DOUBLE, "%15e", "deg"},
  {"Y_WORLD", "Barycenter position along world y axis",
        &refsample.wcspos[1], H_FLOAT, T_DOUBLE, "%15e", "deg"},
  {"ERRA_WORLD", "World RMS position error along major axis",
        &refsample.wcsposerr[0], H_FLOAT, T_FLOAT, "%12e", "deg"},
  {"ERRB_WORLD", "World RMS position error along minor axis",
        &refsample.wcsposerr[1], H_FLOAT, T_FLOAT, "%12e", "deg"},
  {"PMALPHA_J2000", "RA component of proper motion vector",
        &refsample.wcsprop[0], H_FLOAT, T_FLOAT, "%12e", "mas/yr"},
  {"PMDELTA_J2000", "Declination component of proper motion vector",
        &refsample.wcsprop[1], H_FLOAT, T_FLOAT, "%12e", "mas/yr"},
  {"PMALPHAERR_J2000", "Proper motion uncertainty along RA",
        &refsample.wcsproperr[0], H_FLOAT, T_FLOAT, "%12e", "mas/yr"},
  {"PMDELTAERR_J2000", "Proper motion uncertainty along declination",
        &refsample.wcsproperr[1], H_FLOAT, T_FLOAT, "%12e", "mas/yr"},
  {"MAG", "Generic magnitude",
        &refsample.mag, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"MAGERR", "Generic magnitude RMS error",
        &refsample.magerr, H_FLOAT, T_FLOAT, "%8.4f", "mag"},
  {"OBSDATE", "Observation date",
        &refsample.epoch, H_FLOAT, T_DOUBLE, "%13.8f", "yr"},
  {""},
  };

astrefstruct	astrefcats[] = 
 {
  {"NONE", "", {""}, {""}, {""}, 0, 0},

  {"file", "", {""},
	{"1","2","3","4","5","6","7","8","9","10","11","12",""},
	{"1","2","3","4","5","6","7","8","9","10","11","12",""},
	12, 0},

  {"USNO-A2.0", "I/252", {"RAJ2000","DEJ2000","Epoch","Bmag","Rmag",""},
	{"Bmag", "Rmag",""}, {"Bj", "Rf",""},
	2, 0},

  {"USNO-B1.0", "I/284", {"Flags","RAJ2000","DEJ2000","e_RAJ2000","e_DEJ2000",
			"Epoch","pmRA","pmDE","e_pmRA","e_pmDE",
			"B1mag","B2mag","R1mag","R2mag","Imag",""},
	{"B1mag", "R1mag", "Imag",""}, {"Bj", "Rf", "In",""},
	3, 0},

  {"GSC-2.3", "I/305", {"Class","RAJ2000","DEJ2000","e_RAdeg","e_DEdeg",
			"Epoch","Umag","e_Umag","Bmag","e_Bmag","jmag","e_jmag",
			"Vmag","e_Vmag",""},
	{"Umag", "Bmag", "jmag", "Vmag", "Fmag", "Nmag", ""},
	{"U", "B", "Bj", "V", "Rf", "In", ""},
	6, 3},

  {"2MASS", "II/246", {"Cflg","RAJ2000","DEJ2000","errMaj","errMin","errPA",
		"JD","Jmag","e_Jmag","Hmag","e_Hmag","Kmag","e_Kmag",""},
	{"Jmag", "Hmag", "Kmag",""},  {"J", "H", "Ks",""}, 
	3, 0},

  {"DENIS-3", "B/denis", {"RAJ2000","DEJ2000", "ObsJD",
		"Imag","e_Imag","Jmag","e_Jmag","Kmag","e_Kmag",""},
	{"Imag", "Jmag", "Kmag", ""}, {"i", "J", "Ks", ""},
	3, 0},

  {"UCAC-4", "I/322A",  {"Na", "RAJ2000","DEJ2000","ePos","EpRA","EpDE",
		"pmRA","pmDE","e_pmRA","e_pmDE","f.mag","e_a.mag",""},
	{"f.mag", ""}, {"R", ""},
	1, 0},

  {"URAT-1", "I/329", {"Ns","RAJ2000","DEJ2000","sigm","Epoch",
		"pmRA","pmDE","e_pm","f.mag","e_f.mag",""},
	{"f.mag", ""}, {"R", ""},
	1, 0},

  {"SDSS-R9", "V/139", {"mode","Q","RA_ICRS","DE_ICRS","e_RA_ICRS","e_DE_ICRS",
		"ObsDate", "umag","e_umag","gmag","e_gmag","rmag","e_rmag",
		"imag","e_imag","zmag","e_zmag",""},
	{"umag", "gmag", "rmag", "imag", "zmag", ""},
	{"u", "g", "r", "i", "z", ""},
	5, 2},

  {"SDSS-R12", "V/147", {"mode","Q","RA_ICRS","DE_ICRS","e_RA_ICRS","e_DE_ICRS",
		"ObsDate", "umag","e_umag","gmag","e_gmag","rmag","e_rmag",
		"imag","e_imag","zmag","e_zmag",""},
	{"umag", "gmag", "rmag", "imag", "zmag", ""},
	{"u", "g", "r", "i", "z", ""},
	5, 2},

  {"NOMAD-1.0", "I/297", {"RAJ2000","DEJ2000","e_RAJ2000","e_DEJ2000",
		"pmRA","pmDE","e_pmRA","e_pmDE",
		"Bmag","Vmag","Rmag","Jmag","Hmag","Kmag",""},
	{"Bmag", "Vmag", "Rmag", "Jmag", "Hmag", "Kmag",""},
	{"B", "V", "R", "J", "H", "Ks", ""},
	6, 2},

  {"PPMX", "I/312", {"RAJ2000","DEJ2000","e_RAJ2000","e_DEJ2000",
		"pmRA","pmDE","e_pmRA","e_pmDE",
		"Bmag","e_Bmag","Vmag","e_Vmag","Rmag","e_Rmag",
		"Jmag","e_Jmag","Hmag","e_Hmag","Kmag","e_Kmag",""},
	{"Bmag", "Vmag", "Rmag", "Cmag", "Jmag", "Hmag", "Kmag", ""},
	{"B", "V", "R", "Rf", "J", "H", "Ks", ""},
	7, 3},

  {"CMC-15", "I/327", {"RAJ2000","DEJ2000","e_RAJ2000","e_DEJ2000","MJD-51263",
		"r'mag","e_r'mag","Jmag","Hmag","Ksmag",""},
	{"r'mag", "Jmag", "Hmag", "Ksmag",""}, {"r'", "J", "H", "Ks",""},
	4, 0,},

  {"TYCHO-2", "I/259/tyc2", {"pflag","RAmdeg","DEmdeg","e_RAmdeg","e_DEmdeg",
		"pmRA","pmDE","e_pmRA","e_pmDE",
		"BTmag","e_BTmag","VTmag","e_VTmag",""},
	{"BTmag", "VTmag",""}, {"BT", "VT",""},
	2, 1},

  {"IGSL", "I/324", {"RAJ2000","DEJ2000","e_RAJ2000","e_DEJ2000",
		"alphaEpoch","deltaEpoch","pmRA","pmDE","e_pmRA","e_pmDE",
		"magBJ","e_magBJ","magG","e_magG","magRF","e_magRF",
		"magGrvs","e_magGrvs",""},
	{"magBJ", "magG", "magRF", "magGrvs", ""},
	{"Bj", "G", "Rf", "Grvs", ""},
	7, 3},

  {"ALLWISE", "II/328", {"ccf","RAJ2000","DEJ2000","eeMaj","eeMin","eePA",
		"pmRA","pmDE","e_pmRA","e_pmDE",
		"Jmag","e_Jmag","Hmag","e_Hmag","Kmag","e_Kmag",
		"W1mag","e_W1mag","W2mag","e_W2mag","W3mag","e_W3mag",
		"W4mag","e_W4mag",""},
	{"Jmag", "Hmag", "Kmag", "W1mag", "W2mag", "W3mag", "W4mag", ""},
	{"J", "H", "Ks", "W1", "W2", "W3", "W4", ""}, 
	7, 3},

  {"GAIA-DR1", "I/337/gaia", {"Dup", "RA_ICRS","DE_ICRS","e_RA_ICRS","e_DE_ICRS",
		"Epoch","pmRA","pmDE","e_pmRA","e_pmDE",
		"<FG>","e_<FG>","<Gmag>",""},
	{"<Gmag>", ""},
	{"G", ""},
	1, 0},

  {"GAIA-DR2", "I/345/gaia2", {"Dup", "RA_ICRS","DE_ICRS","e_RA_ICRS","e_DE_ICRS",
		"Epoch","pmRA","pmDE","e_pmRA","e_pmDE",
		"Gmag","e_Gmag","BPmag","e_BPmag","RPmag","e_RPmag",""},
	{"Gmag", "BPmag","RPmag",""},
	{"G", "BP", "RP", ""},
	3, 0},

  {"GAIA-EDR3", "I/350/gaiaedr3", {"Dup", "RA_ICRS","DE_ICRS","e_RA_ICRS","e_DE_ICRS",
		"Epoch","pmRA","pmDE","e_pmRA","e_pmDE",
		"Gmag","e_Gmag","BPmag","e_BPmag","RPmag","e_RPmag",""},
	{"Gmag", "BPmag","RPmag",""},
	{"G", "BP", "RP", ""},
	3, 0},

  {"PANSTARRS-1", "II/349", {"Qual", "RAJ2000","DEJ2000","e_RAJ2000","e_DEJ2000",
		"Epoch", "gmag","e_gmag","rmag","e_rmag","imag","e_imag",
		"zmag","e_zmag","ymag","e_ymag",""},
	{"umag", "gmag", "rmag", "imag", "zmag", ""},
	{"g", "r", "i", "z", "y", ""},
	5, 2},

  {""}
 };

static char	*astref_strncpy(char *dest, char *src, int n);
static void	vizier_to_array(char *str, char (*cols)[COLUMN_SIZE]);

const char	astref_ctype[NAXIS][8]= {"RA---STG", "DEC--STG"};

const double	astref_crpix[NAXIS] = {8192.5, 8192.5},
		astref_cdelt[NAXIS] = {0.2*ARCSEC/DEG, 0.2*ARCSEC/DEG};

const int	astref_naxisn[NAXIS] = {16384, 16384};

/****** get_astreffield *********************************************************
PROTO   fieldstruct *get_astreffield(astrefenum refcat, double *wcspos,
				int lng, int lat, int naxis, double maxradius)
PURPOSE	Download reference catalog.
INPUT   Catalog name,
	Coordinate vector of the center,
	Longitude index,
	Latitude index,
	Number of axes (dimensions),
	Search radius (in degrees).
OUTPUT  Pointer to the reference field.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION	01/12/2020
*/
fieldstruct	*get_astreffield(astrefenum refcat, double *wcspos,
				int lng, int lat, int naxis, double maxradius)
  {
   astrefstruct	*astrefcat;
   fieldstruct	*field,*tfield;
   setstruct	*set;
   samplestruct	*sample;
   URL_FILE	*file;
   char		cols[MAX_COLUMN][COLUMN_SIZE],
		url[MAXCHAR], maglimcmd[128],
		str[MAXCHAR],
		*ctype[NAXIS],
		*sepoch[NAXIS], *sprop[NAXIS], *sproperr[NAXIS],
		*sflag, *smag, *smagerr,
		*colname, *col,
		*bandname, *vizierbandname, *catname,
		flag1,flag2, smode;
   double	poserr[NAXIS], prop[NAXIS], properr[NAXIS],
		mag[MAX_BAND], magerr[MAX_BAND], flux, fluxerr, epoch,
		alpha,delta, dist, poserra,poserrb,poserrtheta, cpt,spt;
   int		b,c,d,i,n,s,
		nsample,nsamplemax, nobs, mode, qual, band, nband, cindex;

/* One needs 2 angular coordinates here! */
  if (naxis<2)
    return NULL;

/* If these are not angular coordinates, file mode becomes mandatory */
  if (lng == lat && refcat!=ASTREFCAT_FILE)
    {
    warning("Cartesian coordinates found: ",
	"reference catalog switched to FILE mode");
    refcat = ASTREFCAT_FILE;
    }

  astrefcat = &astrefcats[(int)refcat];
  catname = astrefcat->name;
  nband = astrefcat->nband;
  if (!cistrcmp(prefs.astref_bandname, "DEFAULT", FIND_STRICT))
    band = astrefcat->defband;
  else if (!cistrcmp(prefs.astref_bandname, "BLUEST", FIND_STRICT))
    band = 0;
  else if (!cistrcmp(prefs.astref_bandname, "REDDEST", FIND_STRICT))
    band = nband-1;
  else if ((band=findkeys(prefs.astref_bandname,
		astrefcat->bandnames, FIND_STRICT))==RETURN_ERROR)
    {
    sprintf(str, "%s: no such band in astrometric reference catalog %s, use ",
	prefs.astref_bandname, catname);
    for (b=0; b<nband; b++)
      {
      if (b)
        strcat(str, b==nband-1? " or " : ", ");
      strcat(str, astrefcat->bandnames[b]);
      }
    error(EXIT_FAILURE, "*Error*: ", str);
    }
  astrefcat->band = band;
  bandname = astrefcat->bandname
	= astrefcat->bandnames[band];
  vizierbandname = astrefcat->vizierbandnames[band];

/* Call the right catalog */
  if (refcat==ASTREFCAT_FILE)
    {
    field = NULL;
    for (c=0; c<prefs.nastref_name; c++)
      {
      if ((tfield=load_astreffield(prefs.astref_name[c], wcspos, lng,lat,
		naxis, maxradius, band, prefs.astref_maglim)))
        {
        NFPRINTF(OUTPUT, "");
        QPRINTF(OUTPUT, " %d astrometric references loaded from %s\n",
		tfield->set[0]->nsample, tfield->rfilename);
        if (field)
          {
          union_samples(tfield->set[0]->sample, field->set[0],
		tfield->set[0]->nsample,
		ASTREF_ASSOCRADIUS, UNION_WCS);
          field->nsample = field->set[0]->nsample;
          end_field(tfield);
          }
        else field=tfield;
        }
      }
    if (!field)
        error(EXIT_FAILURE,"*Error*: No appropriate FITS-LDAC astrometric ",
			"reference catalog found");
    field->isrefcat = 1;
    return field;
    }

/* Prepare mag limit section of the command line */
  if (prefs.astref_maglim[0]>-99.0 || prefs.astref_maglim[1]<99.0)
    sprintf(maglimcmd, "%s=%f..%f&",
	vizierbandname, prefs.astref_maglim[0],prefs.astref_maglim[1]);
  else
    strcpy(maglimcmd, "");

/// Test all provided servers until one replies
  for (s=0; s<prefs.nref_server; s++) {
    colname = astrefcat->viziercolumns[0];
    sprintf(url,
	"http://%s/viz-bin/asu-tsv?&-mime=csv&-source=%s&-out.max=100000000"
	"&-out.meta=&%s-c=%.7f,%s%.7f&-c.rd=%.8g&-out=%s",
	prefs.ref_server[s],
	astrefcat->viziername,
	maglimcmd,
	wcspos[lng],
	wcspos[lat]>=0.0? "%2b" : "-",
	fabs(wcspos[lat]),
	maxradius,
	colname);
    while (*(colname += COLUMN_SIZE)) {
      strcat(url, ",");
      strcat(url, colname); 
    }

    sprintf(str,"Querying %s for %s astrometric reference stars...",
	prefs.ref_server[s],
	catname);
    NFPRINTF(OUTPUT, str);
    FPRINTF(OUTPUT, "Querying url %s \n", url);
    file = url_fopen(url, prefs.ref_timeout[s]);
    // Test and skip the first line
    *str = '\0';
    if (url_fgets(str, MAXCHAR, file) && *str == '#')
      break;
    else {
      url_fclose(file);
      if (s == prefs.nref_server-1) {
        if (*str)
          error(EXIT_FAILURE, "*Error*: No VizieR server at ",
		prefs.ref_server[s]);
        else
          error(EXIT_FAILURE, "*Error*: Connection timeout for ",
		prefs.ref_server[s]);
      } else {
        if (*str)
          warning("No VizieR server at ", prefs.ref_server[s]);
        else
          warning("Connection timeout for ", prefs.ref_server[s]);
      }
    }
  }

  set = init_set();
  prop[lng] = prop[lat] = properr[lng] = properr[lat] = 0.0;
  n = nsample = nsamplemax = 0;

/* Now examine each entry */
  while (url_fgets(str, MAXCHAR, file))
    if (*str != '#' && *str != '\n'  && *str != '-')
      {
      memset(cols, 0, MAX_COLUMN*COLUMN_SIZE);
      vizier_to_array(str, cols);
      cindex = 0;
      switch(refcat)
        {
        case ASTREFCAT_USNOA2:
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lat] = poserr[lng] = USNOA2_POSERR;
          sepoch[lng] = cols[cindex++];
          if (sepoch[lng][3] <= ' ')
            continue;
          epoch = atof(sepoch[lng]);
          mag[0] = atof(cols[cindex++]);
          mag[1] = atof(cols[cindex++]);
          magerr[0] = magerr[1] = USNOA2_BMAGERR;
          break;

        case ASTREFCAT_USNOB1:
/*-------- Avoid spikes */
          sflag = cols[cindex++];
          if (sflag[0]=='s' || sflag[1]=='s' || sflag[2]=='s')
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++]);
          if (poserr[lng] != 999.0)
            poserr[lng] *= MAS/DEG;
          else
            poserr[lng] = USNOB1_POSERR;
          poserr[lat] = atof(cols[cindex++]);
          if (poserr[lat] != 999.0)
            poserr[lat] *= MAS/DEG;
          else
            poserr[lat] = USNOB1_POSERR;
          epoch = atof(cols[cindex++]);
          prop[lng] = atof(cols[cindex++])*MAS/DEG;
          prop[lat] = atof(cols[cindex++])*MAS/DEG;
          properr[lng] = atof(cols[cindex++])*MAS/DEG;
          properr[lat] = atof(cols[cindex++])*MAS/DEG;
          if (properr[lng]!=0.0 || properr[lat]!=0.0)
            epoch = 2000.0;
          for (b=0; b<5; b++) {			/* 5, not nband!*/ 
            smag = cols[cindex++];
            if (smag[4] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = USNOB1_BMAGERR;
	    }
          }
/*-------- Merge B1 and B2, R1 and R2, shift In */
          if (mag[0] > 98.0)
            mag[0] = mag[1];
          if (mag[1] > 98.0)
            mag[1] = mag[0];
          mag[0] = (mag[0] + mag[1])*0.5;
          if (magerr[0] > 98.0)
            magerr[0] = magerr[1];
          if (magerr[1] > 98.0)
            magerr[1] = magerr[0];
          magerr[0] = (magerr[0] + magerr[1])*0.5;
          if (mag[2] > 98.0)
            mag[2] = mag[3];
          if (mag[3] > 98.0)
            mag[3] = mag[2];
          mag[1] = (mag[2] + mag[3])*0.5;
          if (magerr[2] > 98.0)
            magerr[2] = magerr[3];
          if (magerr[3] > 98.0)
            magerr[3] = magerr[2];
          magerr[1] = (magerr[2] + magerr[3])*0.5;
          mag[2] = mag[4];
          magerr[2] = magerr[4];
          break;

        case ASTREFCAT_GSC23:
          if (atoi(cols[cindex++])==5)
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*ARCSEC/DEG;
          poserr[lat] = atof(cols[cindex++])*ARCSEC/DEG;
          epoch = atof(cols[cindex++]);
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[4] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = atof(smagerr);
            }
          }
          break;

        case ASTREFCAT_2MASS:
/*-------- Avoid contaminated observations */
          sflag = cols[cindex++];
          if (sflag[0]!='0' ||sflag[1]!='0' || sflag[2]!='0')
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserra = atof(cols[cindex++]);
          poserrb = atof(cols[cindex++]);
          poserrtheta = atof(cols[cindex++]);
/*-------- Project uncertainties on alpha and delta axes */
          cpt = cos(poserrtheta*DEG);
          spt = sin(poserrtheta*DEG);
          poserr[lng] = sqrt(spt*spt*poserra*poserra+cpt*cpt*poserrb*poserrb)
				*ARCSEC/DEG;
          poserr[lat] = sqrt(cpt*cpt*poserra*poserra+spt*spt*poserrb*poserrb)
				*ARCSEC/DEG;
/*-------- Convert JDs to epoch */
          epoch = 2000.0 + (atof(cols[cindex++]) - JD2000)/365.25;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[4] <= ' ' || smagerr[4] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = atof(smagerr);
            }
          }
          break;

        case ASTREFCAT_DENIS3:
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lat] = poserr[lng] = DENIS3_POSERR;
/*-------- Convert JDs to epoch */
          epoch = 2000.0 + (atof(cols[cindex++]) - JD2000)/365.25;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[2] <= ' ' || smagerr[2] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = atof(smagerr);
            }
          }
          break;

        case ASTREFCAT_UCAC4:
/*-------- Avoid poor observations */
          if (atoi(cols[cindex++])<2)
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = poserr[lat] = atof(cols[cindex++]) * MAS/DEG;
          sepoch[lng] = cols[cindex++];
          sepoch[lat] = cols[cindex++];
          sprop[lng] = cols[cindex++];
          sprop[lat] = cols[cindex++];
          sproperr[lng] = cols[cindex++];
          sproperr[lat] = cols[cindex++];
          if (sprop[lng][6]<= ' ' || sprop[lat][6]<= ' ') {
            prop[lng] = prop[lat] = properr[lng] = properr[lat] = 0.0;
            epoch = 0.5*(atof(sepoch[lng]) + atof(sepoch[lat]));
          } else {
            prop[lng] = atof(sprop[lng])*MAS/DEG;
            prop[lat] = atof(sprop[lat])*MAS/DEG;
            properr[lng] = (sproperr[lng][2]==' '?
				1000.0 : atof(sproperr[lng])) * MAS/DEG;
            properr[lat] = (sproperr[lat][2]==' '?
				1000.0 : atof(sproperr[lat])) * MAS/DEG;
            epoch = 2000.0;
          }
          smag = cols[cindex++];
          smagerr = cols[cindex++];
          mag[0] = atof(smag);
          magerr[0] = (smagerr[2]==' ')? 0.9 : atof(smagerr);
          break;

        case ASTREFCAT_URAT1:
/*-------- Avoid poor observations */
          nobs = atoi(cols[cindex++]);
//          if (nobs<2)
//            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = poserr[lat] = atof(cols[cindex++]) * MAS/DEG;
          epoch = atof(cols[cindex++]);
          prop[lng] = atof(cols[cindex++]) * MAS/DEG;
          prop[lat] = atof(cols[cindex++]) * MAS/DEG;
          properr[lng] = properr[lat] = atof(cols[cindex++]) * MAS/DEG;
          mag[0] = atof(cols[cindex++]);
          magerr[0] = atof(cols[cindex++]);
          break;

        case ASTREFCAT_SDSSR9:
/*-------- Avoid missing or poor observations, and secondary detections */
          mode = atoi(cols[cindex++]);
          qual = atoi(cols[cindex++]);
          if (mode==2 || qual<2 || qual>3)
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*ARCSEC/DEG;
          poserr[lat] = atof(cols[cindex++])*ARCSEC/DEG;
          epoch = atof(cols[cindex++]);
          for (b=0; b<nband; b++) {
            mag[b] = atof(cols[cindex++]);
            magerr[b] = atof(cols[cindex++]);
          }
          break;

        case ASTREFCAT_SDSSR12:
/*-------- Avoid missing or poor observations, and secondary detections */
          mode = atoi(cols[cindex++]);
          qual = atoi(cols[cindex++]);
          if (mode==2 || qual<2 || qual>3)
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*ARCSEC/DEG;
          poserr[lat] = atof(cols[cindex++])*ARCSEC/DEG;
          epoch = atof(cols[cindex++]);
          for (b=0; b<nband; b++) {
            mag[b] = atof(cols[cindex++]);
            magerr[b] = atof(cols[cindex++]);
          }
          break;

        case ASTREFCAT_NOMAD1:
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*MAS/DEG;
          if (poserr[lng] < TINY)
            poserr[lng] = USNOB1_POSERR;
          poserr[lat] = atof(cols[cindex++])*MAS/DEG;
          if (poserr[lat] < TINY)
            poserr[lat] = USNOB1_POSERR;
          epoch = 2000.0;
          prop[lng] = atof(cols[cindex++]) * MAS/DEG;
          prop[lat] = atof(cols[cindex++]) * MAS/DEG;
          properr[lng] = atof(cols[cindex++])*MAS/DEG;
          properr[lat] = atof(cols[cindex++])*MAS/DEG;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[2] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = NOMAD1_MAGERR;
            }
          }
          break;

        case ASTREFCAT_PPMX:
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*MAS/DEG;
          poserr[lat] = atof(cols[cindex++])*MAS/DEG;
          epoch = 2000.0;
          prop[lng] = atof(cols[cindex++]) * MAS/DEG;
          prop[lat] = atof(cols[cindex++]) * MAS/DEG;
          properr[lng] = atof(cols[cindex++])*MAS/DEG;
          properr[lat] = atof(cols[cindex++])*MAS/DEG;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = (b!=2 && b!=3) ? cols[cindex++] : (char *)NULL;
            if (smag[2] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = smagerr? atof(smagerr)*0.001 : GSC_MAGERR;
            }
          }
          break;

        case ASTREFCAT_CMC15:
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*ARCSEC/DEG;
          if (poserr[lng] < TINY)
            poserr[lng] = CMC15_POSERR;
          poserr[lat] = atof(cols[cindex++])*ARCSEC/DEG;
          if (poserr[lat] < TINY)
            poserr[lat] = CMC15_POSERR;
          epoch = 2000.0 + (atof(cols[cindex++]) - (JD2000-2451263.5))/365.25;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = b==0 ? cols[cindex++] : (char *)NULL;
            if (smag[2] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = smagerr? atof(smagerr) : DEFAULT_MAGERR;
              if (magerr[b] == 0.0)
                magerr[b] = DEFAULT_MAGERR;
            }
          }
          break;

        case ASTREFCAT_TYCHO2:
/*-------- Reject (fainter) sources without mean positions and proper motions */
          sflag = cols[cindex++]; 
          if (sflag[0]=='X')
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*MAS/DEG;
          poserr[lat] = atof(cols[cindex++])*MAS/DEG;
          epoch = 2000.0;
          prop[lng] = atof(cols[cindex++]) * MAS/DEG;
          prop[lat] = atof(cols[cindex++]) * MAS/DEG;
          properr[lng] = atof(cols[cindex++])*MAS/DEG;
          properr[lat] = atof(cols[cindex++])*MAS/DEG;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[2] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = smagerr[2] <= ' '? atof(smagerr): DEFAULT_MAGERR;
            }
          }
          break;

        case ASTREFCAT_IGSL:
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*ARCSEC/DEG;
          poserr[lat] = atof(cols[cindex++])*ARCSEC/DEG;
          sepoch[lng] = cols[cindex++];
          sepoch[lat] = cols[cindex++];
          prop[lng] = atof(cols[cindex++]) * MAS/DEG;
          prop[lat] = atof(cols[cindex++]) * MAS/DEG;
          properr[lng] = atof(cols[cindex++])*MAS/DEG;
          properr[lat] = atof(cols[cindex++])*MAS/DEG;
          epoch = 2000.0;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[2] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = smagerr[2] <= ' '? atof(smagerr): DEFAULT_MAGERR;
            }
          }
          break;

        case ASTREFCAT_ALLWISE:
/*-------- Avoid contaminated observations */
          sflag = cols[cindex++];
          if (sflag[0]!='0' || sflag[1]!='0' || sflag[2]!='0' || sflag[3]!='0')
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserra = atof(cols[cindex++]);
          poserrb = atof(cols[cindex++]);
          poserrtheta = atof(cols[cindex++]);
/*-------- Project uncertainties on alpha and delta axes */
          cpt = cos(poserrtheta*DEG);
          spt = sin(poserrtheta*DEG);
          poserr[lng] = sqrt(spt*spt*poserra*poserra+cpt*cpt*poserrb*poserrb)
				*ARCSEC/DEG;
          poserr[lat] = sqrt(cpt*cpt*poserra*poserra+spt*spt*poserrb*poserrb)
				*ARCSEC/DEG;
          epoch = 2000.0;
          prop[lng] = atof(cols[cindex++]) * MAS/DEG;
          prop[lat] = atof(cols[cindex++]) * MAS/DEG;
          properr[lng] = atof(cols[cindex++]) * MAS/DEG;
          properr[lat] = atof(cols[cindex++]) * MAS/DEG;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[4] <= ' ' || smagerr[4] <= ' ') {
              mag[b] = magerr[b] = 99.0;
}
            else {
              mag[b] = atof(smag);
              magerr[b] = atof(smagerr);
            }
          }
          break;

        case ASTREFCAT_GAIADR1:
/*-------- Reject duplicated sources */
          sflag = cols[cindex++];
/*
          if (sflag[0]!='0')
            continue;
*/
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*MAS/DEG;
          poserr[lat] = atof(cols[cindex++])*MAS/DEG;
          epoch = atof(cols[cindex++]);
          sprop[lng] = cols[cindex++];
          sprop[lat] = cols[cindex++];
          sproperr[lng] = cols[cindex++];
          sproperr[lat] = cols[cindex++];
          if (sprop[lng][4]<= ' ' || sprop[lat][4]<= ' ') {
            prop[lng] = prop[lat] = properr[lng] = properr[lat] = 0.0;
          } else {
            prop[lng] = atof(sprop[lng])*MAS/DEG;
            prop[lat] = atof(sprop[lat])*MAS/DEG;
            properr[lng] = (sproperr[lng][1]==' '?
				1000.0 : atof(sproperr[lng])) * MAS/DEG;
            properr[lat] = (sproperr[lat][1]==' '?
				1000.0 : atof(sproperr[lat])) * MAS/DEG;
          }
          flux = atof(cols[cindex++]);
          fluxerr = atof(cols[cindex++]);
          mag[0] = atof(cols[cindex++]);
          if (flux <= 0.0)
            mag[0] = magerr[0] = 99.0;
          else
            magerr[0] = 1.0857 * fluxerr / flux;
          break;

        case ASTREFCAT_GAIADR2:
        case ASTREFCAT_GAIAEDR3:
/*-------- Reject duplicated sources */
          sflag = cols[cindex++];
/*
          if (sflag[0]!='0')
            continue;
*/
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*MAS/DEG;
          poserr[lat] = atof(cols[cindex++])*MAS/DEG;
          epoch = atof(cols[cindex++]);
          sprop[lng] = cols[cindex++];
          sprop[lat] = cols[cindex++];
          sproperr[lng] = cols[cindex++];
          sproperr[lat] = cols[cindex++];
          if (sprop[lng][4]<= ' ' || sprop[lat][4]<= ' ') {
            prop[lng] = prop[lat] = properr[lng] = properr[lat] = 0.0;
          } else {
            prop[lng] = atof(sprop[lng])*MAS/DEG;
            prop[lat] = atof(sprop[lat])*MAS/DEG;
            properr[lng] = (sproperr[lng][1]==' '?
				1000.0 : atof(sproperr[lng])) * MAS/DEG;
            properr[lat] = (sproperr[lat][1]==' '?
				1000.0 : atof(sproperr[lat])) * MAS/DEG;
          }
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[4] <= ' ' || smagerr[4] <= ' ') {
              mag[b] = magerr[b] = 99.0;
}
            else {
              mag[b] = atof(smag);
              magerr[b] = atof(smagerr);
            }
          }
          break;

        case ASTREFCAT_PANSTARRS1:
          if (atoi(cols[cindex++]) & 4 == 0)	// Test PS1 reliability
            continue;
          alpha = atof(cols[cindex++]);
          delta = atof(cols[cindex++]);
          poserr[lng] = atof(cols[cindex++])*ARCSEC/DEG;
          poserr[lat] = atof(cols[cindex++])*ARCSEC/DEG;
          epoch = 2000.0 + (atof(cols[cindex++]) - (JD2000 - 2400000.5))/365.25;
          for (b=0; b<nband; b++) {
            smag = cols[cindex++];
            smagerr = cols[cindex++];
            if (smag[2] <= ' ')
              mag[b] = magerr[b] = 99.0;
            else {
              mag[b] = atof(smag);
              magerr[b] = smagerr[2] <= ' '? atof(smagerr): DEFAULT_MAGERR;
            }
          }
          break;

        case ASTREFCAT_NONE:
        default:
          break;
        }

      if (!(n%10000))
        {
        sprintf(str,"%-.36s (%s band) : %d / %d references stored",
		catname,bandname, nsample,n);
        NFPRINTF(OUTPUT, str);
        }
/*------ Apply some selection over flags, fluxes,.. */
      if (1)
        {
/*------ Allocate memory for the first shipment */
        if (!set->sample)
          {
          nsample = 0;
          nsamplemax = LSAMPLE_DEFSIZE;
          malloc_samples(set, nsamplemax);
          }

/*------ Increase storage space to receive new candidates if needed */
        if (nsample>=nsamplemax)
          {
           int	nadd=(int)(1.62*nsamplemax);
          nsamplemax = nadd>nsamplemax?nadd:nsamplemax+1;
          realloc_samples(set, nsamplemax);
          }

        sample = set->sample + nsample;

        if (mag[band] < 98.0)
          {
          sample->mag = mag[band];
          sample->magerr = magerr[band];
          }
        else
          {
          for (b=0; b<nband; b++)
            if (mag[b] < 98.0)
              {
              sample->mag = mag[b];
              sample->magerr = magerr[b];
              break;
              }
/*-------- If no magnitude is suitable, drop this source */
          if (mag[b] > 98.0)
            continue;
          }
/*
        if (nband>1 && mag[0]<98.0 && mag[1]<98.0)
          sample->colour = mag[0] - mag[1];
        else
          sample->colour = 99.0;
*/
/*
printf("%10f %10f +/-%4.0f %4.0f  %4.0f %4.0f +/-%4.0f %4.0f  %7.2f  %5.2f +/- %4.2f\n",
alpha,delta,poserr[lng]*DEG/MAS,poserr[lat]*DEG/MAS,
prop[lng]*DEG/MAS,prop[lat]*DEG/MAS,properr[lng]*DEG/MAS,properr[lat]*DEG/MAS,
epoch, sample->mag, sample->magerr);
*/
        sample->flux = 0.0;
        sample->rawpos[lng] = sample->rawpos[lat] = 0.0;
        sample->rawposerr[lng] = sample->rawposerr[lat] = 0.0;
        sample->wcspos[lng] = alpha;
        sample->wcspos[lat] = delta;
        sample->wcsposerr[lng] = poserr[lng];
        sample->wcsposerr[lat] = poserr[lat];
        sample->wcsprop[lng] = prop[lng];
        sample->wcsprop[lat] = prop[lat];
        sample->wcsproperr[lng] = properr[lng];
        sample->wcsproperr[lat] = properr[lat];
        sample->epoch = epoch;
        sample->spread = sample->spreaderr = 0.0;
        sample->sexflags = 0;	/* SEx flags not relevant for ref. sources*/
        sample->scampflags = 0;
        sample->imaflags = 0;
        sample->set = set;
        nsample++;
        }
      n++;
      }

  url_fclose(file);

/* Escape here if no source returned */
  if (!nsample)
    return NULL;

  set->nsample = nsample;

/* Don't waste memory! */
  realloc_samples(set, nsample);

/* Now "create" the field */
  QCALLOC(field, fieldstruct, 1);
  QMALLOC(field->set, setstruct *, 1);
  field->set[0] = set;
  field->nset = 1;
/* This is a reference catalog */
  field->isrefcat = 1;
  field->astromlabel = field->photomlabel = -1;
  set->lng = field->lng = lng;
  set->lat = field->lat = lat;
  set->naxis = field->naxis = naxis;
  field->maxradius = set->radius = maxradius;
  set->field = field;
  for (d=0; d<naxis; d++)
    {
    set->wcspos[d] = field->meanwcspos[d] = wcspos[d];
    ctype[d] = (char *)astref_ctype[d];
    }

/* Create a dummy WCS structure */
  set->wcs = create_wcs(ctype, wcspos, (double *)astref_crpix,
			(double *)astref_cdelt,
			(int *)astref_naxisn, naxis);

  field->nsample = set->nsample;

  return field;
  }


/****** save_astreffield ******************************************************
PROTO   void save_astreffield(char *filename, fieldstruct *reffield)
PURPOSE	Save an astrometric reference catalog in (pseudo-)LDAC format.
INPUT   Catalog name,
	pointer to the reference field to save.
OUTPUT  -.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 22/01/2019
*/
void	save_astreffield(char *filename,  fieldstruct *reffield)
  {
   static char	imtabtemplate[][80] = {
"SIMPLE  =                    T / This is a FITS file",
"BITPIX  =                    8 / ",
"NAXIS   =                    2 / 2D data",
"NAXIS1  =                    1 / Number of rows",
"NAXIS2  =                    1 / Number of columns",
"EXTEND  =                    T / This file may contain FITS extensions",
"END                            "};
  catstruct	*cat;
  tabstruct	*asctab, *imtab, *objtab;
  keystruct	*key, *objkeys;
  setstruct	*set;
  samplestruct	*sample,
		objsample;
  char		*buf;
  long		dptr;
  int		i,k,n,s;

/* Create a new output catalog */
  cat = new_cat(1);
  init_cat(cat);
  strcpy(cat->filename, filename);
  if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
    error(EXIT_FAILURE, "*Error*: cannot open for writing ", filename);

/* Primary header */
  save_tab(cat, cat->tab);

  for (s=0; s<reffield->nset; s++)
    {
    set = reffield->set[s];
/*-- We create a dummy table (only used through its header) */
    QCALLOC(asctab, tabstruct, 1);
    asctab->headnblock = 1 + (sizeof(imtabtemplate)-1)/FBSIZE;
    QCALLOC(asctab->headbuf, char, asctab->headnblock*FBSIZE);
    memcpy(asctab->headbuf, imtabtemplate, sizeof(imtabtemplate));
    for (buf = asctab->headbuf, i=FBSIZE*asctab->headnblock; i--; buf++)
      if (!*buf)
        *buf = ' ';
    write_wcs(asctab, set->wcs);
/*-- (dummy) LDAC Image header */

    imtab = new_tab("LDAC_IMHEAD");
    key = new_key("Field Header Card");
    key->ptr = asctab->headbuf;
    asctab->headbuf = NULL;
    free_tab(asctab);
    key->naxis = 2;
    QMALLOC(key->naxisn, int, key->naxis);
    key->naxisn[0] = 80;
    key->naxisn[1] = 36;
    key->htype = H_STRING;
    key->ttype = T_STRING;
    key->nobj = 1;
    key->nbytes = 80*(fitsfind(key->ptr, "END     ")+1);
    add_key(key, imtab, 0);
    save_tab(cat, imtab);
    free_tab(imtab);

/*-- LDAC Object header */
    objtab = new_tab("LDAC_OBJECTS");
    objtab->cat = cat;
/*-- Set key pointers */
    QCALLOC(objkeys, keystruct, (sizeof(refkey) / sizeof(keystruct)));
    dptr = (long)((char *)&objsample - (char *)&refsample);
    for (k=0; refkey[k].name[0]; k++)
      {
      objkeys[k] = refkey[k];
      objkeys[k].ptr = (void *)((char *)objkeys[k].ptr + dptr); /* a trick */
      add_key(&objkeys[k],objtab, 0);
      }
    init_writeobj(cat, objtab, &buf);
    sample = set->sample;   
    for (n=set->nsample; n--;)
      {
      objsample = *(sample++);
// Save proper motions in mas/yr
      objsample.wcsprop[0] *= DEG / MAS;
      objsample.wcsprop[1] *= DEG / MAS;
      objsample.wcsproperr[0] *= DEG / MAS;
      objsample.wcsproperr[1] *= DEG / MAS;
      write_obj(objtab, buf);
      }
    end_writeobj(cat, objtab, buf);
    objtab->key = NULL;
    objtab->nkey = 0;
    free_tab(objtab);
    free(objkeys);
    }

  free_cat(&cat, 1);

  return;
  }

/****** load_astreffield ******************************************************
PROTO   fieldstruct *load_astreffield(char *filename, double *wcspos,
			int lng, int lat, int naxis, double maxradius, int band,
			double *maglim)
PURPOSE	Load a reference catalog in (pseudo-)LDAC format.
INPUT   Catalog name,
	pointer to the field center coordinates,
	longitude index,
	latitude index,
	number of axes (dimensions),
	search radius (in degrees),
	band index,
	magnitude range.
OUTPUT  Pointer to the reference catalog field structure.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 02/10/2009
*/
fieldstruct	*load_astreffield(char *filename, double *wcspos,
				int lng, int lat,
				int naxis, double maxradius, int band,
				double *maglim)
  {
   catstruct	*cat;
   tabstruct	*tab;
   fieldstruct	*field;
   setstruct	*set;
   char		str[MAXCHAR],
		*rfilename, *pstr, *pspath;
   int		i, n, nsample;


/* A short, "relative" version of the filename */
  if (!(rfilename = strrchr(filename, '/')))
    rfilename = filename;
  else
    rfilename++;

  sprintf(str,"Examining Catalog %s...", rfilename);
  NFPRINTF(OUTPUT, str);

/*-- Read input catalog */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: No such catalog: ", filename);
  QCALLOC(field, fieldstruct, 1);
  QMALLOC(field->set, setstruct *, 1);
  strcpy (field->filename, filename);
  field->rfilename = rfilename;

/* Create a file name with a "header" extension */
  strcpy(field->hfilename, filename);
  if (!(pstr = strrchr(field->hfilename, '.')))
    pstr = field->hfilename+strlen(field->hfilename);
  sprintf(pstr, "%s", prefs.ahead_suffix);

/* Extract the path from the filename */
#ifdef HAVE_GETENV
  pspath = getenv("PWD");
#else
  pspath = NULL;
#endif
  if (*field->filename == '/')
    strcpy(field->path, field->filename);
  else
    {
    strcpy(field->path, pspath? pspath: ".");
    if (*field->filename != '.' && (pstr = strchr(field->filename, '/')))
      {
      strcat(field->path, "/");
      strcat(field->path, pstr+1);
      }
    }
  if ((pstr = strrchr(field->path, '/')))
    *pstr = '\0';

/* Identify image headers in catalog  */
  tab = cat->tab;
 
  n = 0;

/* Find the object table */
  sprintf(str,"Loading Catalog %s...", rfilename);
  NFPRINTF(OUTPUT, str);
  tab = cat->tab;
  set = NULL;
  nsample = 0;
  n = 0;
  for (i=cat->ntab; i--; tab=tab->nexttab)
    if (!strcmp("LDAC_OBJECTS", tab->extname)
	|| !strcmp("OBJECTS", tab->extname))
    {
    if (field->nset>1)
      sprintf(str, "%s [%d]", rfilename, n+1);
    else
      strcpy(str, rfilename);
    set = read_astrefsamples(set,tab, str, wcspos,lng,lat,naxis,maxradius,band,
			maglim);
    nsample += set->nsample;
    n++;
    }

  field->nsample = nsample;
  free_cat(&cat, 1);

  if (!n)
    {
/*-- No source found: return a NULL pointer */
    end_field(field);
    return NULL;
    }

/* This is a reference catalog */
  field->astromlabel = field->photomlabel = -1;
  field->epoch = -1.0;
  field->set[0] = set;
  field->nset = 1;
  set->lng = field->lng = lng;
  set->lat = field->lat = lat;
  set->naxis = field->naxis = naxis;
  set->radius = maxradius;
  set->field = field;

  field->nsample = set->nsample;

  return field;
  }


/****** read_astrefsamples ****************************************************
PROTO	setstruct *read_astrefsamples(setstruct *set, tabstruct *tab,
				char *rfilename,
				double *wcspos, int lng, int lat, int naxis,
				double maxradius, int band, double *maglim)
PURPOSE	Read a set of astrometric reference samples.
INPUT	Set structure pointer,
	Pointer to the tab that contains the catalog,
        Reduced filename.
	Coordinate vector of the center,
	Longitude index,
	Latitude index,
	Number of axes (dimensions),
	Search radius (in degrees),
	band index,
	magnitude range.
OUTPUT  setstruct pointer (allocated if the input setstruct pointer is NULL).
NOTES   The filename is used for error messages only. Global preferences are
	used.
AUTHOR  E. Bertin (IAP)
VERSION 28/06/2020
*/
setstruct *read_astrefsamples(setstruct *set, tabstruct *tab, char *rfilename,
				double *wcspos, int lng, int lat, int naxis,
				double maxradius, int band, double *maglim)


  {
   tabstruct		*keytab;
   keystruct		*key;
   samplestruct		*sample;
   char			str[MAXCHAR];
   char			*buf;
   unsigned short	*flags,
			objflags;
   float		*xm,*ym, *xpm,*ypm, *xpmerr,*ypmerr,
			*mag, *magerr, *obsdate, *erra,*errb;
   double		*dxm, *dym, *dxpm,*dypm, *dxpmerr,*dypmerr,
			*dmag, *dmagerr, *dobsdate, *derra, *derrb,
			x,y, dx,dy,dfac, ea,eb, maxradius2, mmag;
   int			n, nsample,nsamplemax, nobj, maglimflag;

/* One needs 2 angular coordinates here! */
  dxm = dym = dxpm = dypm = dxpmerr = dypmerr = dmag = derra = derrb \
	= dmagerr = NULL;
  xm = ym = xpm = ypm = xpmerr = ypmerr = mag = erra = errb = magerr = NULL;
  maxradius2 = maxradius*maxradius;
  dfac = (lng!=lat)? cos(wcspos[lat]*DEG) : 1.0;

/* If a NULL pointer is provided, we allocate a new set */
  if (!set)
    {
    set = init_set();
    nsample = nsamplemax = 0;
    }
  else
    nsample = nsamplemax = set->nsample;

/* Init the single-row tab */
  keytab = init_readobj(tab, &buf);

  if (!(key = name_to_key(keytab, prefs.astrefcent_key[0])))
    {
    sprintf(str, "*Error*: %s parameter not found in catalog ",
		prefs.astrefcent_key[0]);
    error(EXIT_FAILURE, str, rfilename);
    }
  if (key->ttype == T_DOUBLE)
    dxm = (double *)key->ptr;
  else
    xm = (float *)key->ptr;

  nobj = key->nobj;

  if (!(key = name_to_key(keytab, prefs.astrefcent_key[1])))
    {
    sprintf(str, "*Error*: %s parameter not found in catalog ",
		prefs.astrefcent_key[1]);
    error(EXIT_FAILURE, str, rfilename);
    }
  if (key->ttype == T_DOUBLE)
    dym = (double *)key->ptr;
  else
    ym = (float *)key->ptr;

  if (!(key = name_to_key(keytab, prefs.astreferr_key[0])))
    {
    sprintf(str, "*Error*: %s parameter not found in catalog ",
		prefs.astreferr_key[0]);
    error(EXIT_FAILURE, str, rfilename);
    }
  if (key->ttype == T_DOUBLE)
    derra = (double *)key->ptr;
  else
    erra = (float *)key->ptr;

  if (!(key = name_to_key(keytab, prefs.astreferr_key[1])))
    {
    sprintf(str, "*Error*: %s parameter not found in catalog ",
		prefs.astreferr_key[1]);
    error(EXIT_FAILURE, str, rfilename);
    }
  if (key->ttype == T_DOUBLE)
    derrb = (double *)key->ptr;
  else
    errb = (float *)key->ptr;

  if ((key = name_to_key(keytab, prefs.astrefprop_key[0])))
    {
    if (key->ttype == T_DOUBLE)
      dxpm = (double *)key->ptr;
    else
      xpm = (float *)key->ptr;
    }

  if ((key = name_to_key(keytab, prefs.astrefprop_key[1])))
    {
    if (key->ttype == T_DOUBLE)
      dypm = (double *)key->ptr;
    else
      ypm = (float *)key->ptr;
    }

  if ((key = name_to_key(keytab, prefs.astrefproperr_key[0])))
    {
    if (key->ttype == T_DOUBLE)
      dxpmerr = (double *)key->ptr;
    else
      xpmerr = (float *)key->ptr;
    }

  if ((key = name_to_key(keytab, prefs.astrefproperr_key[1])))
    {
    if (key->ttype == T_DOUBLE)
      dypmerr = (double *)key->ptr;
    else
      ypmerr = (float *)key->ptr;
    }

  if (!(key = name_to_key(keytab, prefs.astrefmag_key)))
    {
    sprintf(str, "*Error*: %s parameter not found in catalog ",
		prefs.astrefmag_key);
    error(EXIT_FAILURE, str, rfilename);
    }
  if (key->ttype == T_DOUBLE)
    dmag = (double *)key->ptr;
  else
    mag = (float *)key->ptr;

  if (!(key = name_to_key(keytab, prefs.astrefmagerr_key)))
    {
    sprintf(str, "%s parameter not found in catalog ", prefs.astrefmagerr_key);
    warning(str, rfilename);
    }
  else
    {
    if (key->ttype == T_DOUBLE)
      dmagerr = (double *)key->ptr;
    else
      magerr = (float *)key->ptr;
    }

  if (!(key = name_to_key(keytab, prefs.astrefobsdate_key)))
    {
    sprintf(str, "%s parameter not found in catalog ", prefs.astrefobsdate_key);
    warning(str, rfilename);
    dobsdate = NULL;
    obsdate = NULL;
    }
  else
    {
    if (key->ttype == T_DOUBLE)
      dobsdate = (double *)key->ptr;
    else
      obsdate = (float *)key->ptr;
    }

/* Check that catalog contains enough bands if needed */
  if (band && (!key->naxis || band>=*key->naxisn))
    {
    sprintf(str, "*Error*: band #%d not found in catalog ", band+1);
    error(EXIT_FAILURE, str, rfilename);
    }

  if (!(key = name_to_key(keytab, "FLAGS")))
    warning("FLAGS parameter not found in catalog ", rfilename);
  flags = key? (unsigned short *)key->ptr : NULL;

/* Check magnitude limits only if needed */
  maglimflag = (maglim[0]>-99.0 || maglim[1]<99.0)? 1 : 0;

/* Now examine each vector of the shipment */
  for (n=0; nobj--; n++)
    {
    objflags = 0;
    read_obj(keytab,tab, buf);
    if (!(n%10000))
      {
      sprintf(str,"%-.36s: %d / %d references stored",
		rfilename,nsample,n);
      NFPRINTF(OUTPUT, str);
      }
/*---- Apply some selection over flags, fluxes... */
    mmag = mag? mag[band] : dmag[band];
    if (maglimflag && (mmag<maglim[0] || mmag>maglim[1]))
      continue;
    if (flags)
      {
      if (*flags & prefs.sexflags_mask)
        continue;
/*---- Mapping from SExtractor flags is straightforward */
      objflags = *flags;
      if (objflags & OBJ_SATUR)		/* A saturated object */
        set->nsaturated++;
      }
    x = fmod((xm? *xm : *dxm) +360.0,360.0);
    y = ym? *ym : *dym;
    dx = x-wcspos[lng];
    if (dx>180.0)
      dx -= 360.0;
    else if (dx<-180.0)
      dx += 360.0;
    dx *= dfac;
    dy = y-wcspos[lat];
    if (dx*dx+dy*dy > maxradius2)
      continue;
/*-- ... and check the integrity of the sample */
/*-- Allocate memory for the first shipment */
    if (!set->sample)
      {
      nsample = 0;
      nsamplemax = LSAMPLE_DEFSIZE;
      malloc_samples(set, nsamplemax);
      }

/*-- Increase storage space to receive new candidates if needed */
    if (nsample>=nsamplemax)
      {
       int	nadd=(int)(1.62*nsamplemax);
      nsamplemax = nadd>nsamplemax?nadd:nsamplemax+1;
      realloc_samples(set, nsamplemax);
      }

    sample = set->sample + nsample;
    sample->set = set;
    sample->sexflags = objflags;
    sample->mag = mag? mag[band] : dmag[band];
    sample->magerr = magerr? magerr[band] : (dmagerr? dmagerr[band] : 0.0);
    sample->epoch = dobsdate? *dobsdate : (obsdate? *obsdate : 0.0);
    sample->flux = 0.0;
    sample->wcspos[lng] = x;
    sample->wcspos[lat] = y;
    ea = erra? *erra : *derra;
    eb = errb? *errb : *derrb;
    sample->wcsposerr[lng] = sample->wcsposerr[lat] = sqrt(ea*ea+eb*eb);
    sample->wcsprop[lng] = (xpm ? *xpm : (dxpm ? *dxpm : 0.0)) * MAS / DEG;
    sample->wcsprop[lat] = (ypm ? *ypm : (dypm ? *dypm : 0.0)) * MAS / DEG;
    sample->wcsproperr[lng] = (xpmerr ? *xpmerr : (dxpmerr ? *dxpmerr : 0.0))
				* MAS / DEG;
    sample->wcsproperr[lat] = (ypmerr ? *ypmerr : (dypmerr ? *dypmerr : 0.0))
				* MAS / DEG;
/*-- In case of a contamination, position errors are easily doubled */
    if (flags && *flags>0)
      sample->wcsposerr[lng] = (sample->wcsposerr[lat] *= 2.0);
    nsample++;
    }

  end_readobj(keytab, tab, buf);

  set->nsample = nsample;

/* Don't waste memory! */
  if (nsample)
    realloc_samples(set, nsample);

  return set;
  }

/***i** astref_strncpy ****************************************************
PROTO	char *astref_strncpy(char *dest, char *src, int n)
PURPOSE	Copy a piece of string with max length n, and end it with '\0'.
INPUT	Destination string,
	input string,
	max number of characters.
OUTPUT  Pointer to dest.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 10/09/2012
*/
static char	*astref_strncpy(char *dest, char *src, int n)

  {
   char	*destt;
   int	i;

  destt = dest;
  for (i=n; i-- && *src != '\0';)
    *(destt++) = *(src++);
  *destt = '\0';

  return dest;
  }


/***i** vizier_to_array ****************************************************
PROTO	void vizier_to_array(char *str, char (*cols)[COLUMN_SIZE])
PURPOSE	Copy columns from the Vizier output to an array of strings
INPUT	Input Vizier string,
	column array.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION 19/02/2018
*/
static void	vizier_to_array(char *str, char (*cols)[COLUMN_SIZE])

  {
   char		*col;

  if (*str)
    do {
      col = *(cols++);
      while (*str && *str != ';')
        *(col++) = *(str++); 
      *col = '\0';
    } while (*(str++));

  return;
  }


