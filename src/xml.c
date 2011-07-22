/*
*				xml.c
*
* Handle XML metadata.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		22/07/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "globals.h"
#include "astrefcat.h"
#include "fits/fitscat.h"
#include "field.h"
#include "fgroup.h"
#include "key.h"
#include "prefs.h"
#include "samples.h"
#include "cplot.h"
#include "xml.h"

extern time_t		thetime,thetime2;	/* from makeit.c */
extern pkeystruct	key[];			/* from preflist.h */
extern char		keylist[][32];		/* from preflist.h */
 
fgroupstruct		**fgroups_xml;
fieldstruct		**fields_xml;
int			nfield_xml=0, ngroup_xml=0;
char			**astrinstrustr;
int			nastrinstrustr;


/****** init_xml ************************************************************
PROTO	void init_xml(fieldstruct **fields, int nfield,
			fgroupstruct **fgroups, int ngroup)
PURPOSE	Initialize a set of meta-data kept in memory before being written to the
	XML file
INPUT	Array of pointers to fields,
	number of fields,
	array of pointers to fgroups,
	number of groups.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/10/2006
 ***/
void	init_xml(fieldstruct **fields, int nfield,
			fgroupstruct **fgroups, int ngroup)
  {
   fields_xml = fields;
   nfield_xml = nfield;
   fgroups_xml = fgroups;
   ngroup_xml = ngroup;

  return;
  }


/****** end_xml ************************************************************
PROTO	void end_xml(void)
PURPOSE	Here only for consistency.
INPUT	-.
OUTPUT	.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	06/10/2006
 ***/
void	end_xml(void)
  {

  return;
  }


/****** write_xml ************************************************************
PROTO	int	write_xml(char *filename)
PURPOSE	Save meta-data to an XML file/stream.
INPUT	XML file name.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	01/06/2007
 ***/
int	write_xml(char *filename)
  {
   FILE		*file;

  if (!strcmp(filename, "STDOUT"))
    file = stdout;
  else if (!(file = fopen(prefs.xml_name, "w")))
    return RETURN_ERROR;

  write_xml_header(file);
/*
  write_vo_fields(file);

  fprintf(file, "   <DATA>\n");
  if (prefs.mergedcat_type == CAT_FITS_LDAC)
    fprintf(file,
	"   <FITS extnum=\"%2\"><STREAM href=\"%s%s\" /> </FITS>",
	prefs.mergedcat_name[0] == '/'? "file://" : "file:",
	prefs.mergedcat_name);
  fprintf(file, "   </DATA>\n");
  fprintf(file, "  </TABLE>\n");
*/

  write_xml_meta(file, (char *)NULL);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  fclose(file);

  return RETURN_OK;
  }


/****** write_xml_header ******************************************************
PROTO	int	write_xml_header(FILE *file)
PURPOSE	Save an XML-VOtable header to an XML file/stream
INPUT	file or stream pointer.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	Global preferences are used.
AUTHOR	E. Bertin (IAP)
VERSION	06/10/2006
 ***/
int	write_xml_header(FILE *file)
  {
   char		sysname[16];

  fprintf(file, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  fprintf(file, "<?xml-stylesheet type=\"text/xsl\" href=\"%s\"?>\n",
	prefs.xsl_name);
  fprintf(file, "<VOTABLE "
	"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
	"xsi:noNamespaceSchemaLocation="
	"\"http://www.ivoa.net/xml/VOTable/v1.1\">\n");
  fprintf(file, "<DESCRIPTION>produced by %s</DESCRIPTION>\n", BANNER);
  fprintf(file, "<!-- VOTable description at "
	"http://www.ivoa.net/Documents/latest/VOT.html -->\n");
  fprintf(file, "<RESOURCE ID=\"%s\" name=\"%s\">\n", BANNER,
		nfield_xml? fields_xml[0]->rfilename: BANNER);
  fprintf(file, " <DESCRIPTION>Data related to %s"
	"</DESCRIPTION>\n", BANNER);
  fprintf(file, " <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n");
  sprintf(sysname, "ICRS");

  fprintf(file, " <COOSYS ID=\"J2000\" equinox=\"J2000\""
	" epoch=\"2000.0\" system=\"ICRS\"/>\n");

  return RETURN_OK;
  }


/****** write_xml_meta ********************************************************
PROTO	int	write_xml_meta(FILE *file, char *msgerror)
PURPOSE	Save meta-data to an XML-VOTable file or stream
INPUT	Pointer to the output file (or stream),
	Pointer to an error msg (or NULL).
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP) C. Marmo (IAP)
VERSION	22/07/2011
 ***/
int	write_xml_meta(FILE *file, char *msgerror)
  {
   fgroupstruct		*fgroup;
   fieldstruct		*field;
   struct tm		*tm;
   double		deg2arcsec,deg2arcmin;
   char			*pspath,*psuser, *pshost, *str;
   int			*cp,
			d,f,f2,g,i,l,n,len, lng,lat, naxis;
#ifdef HAVE_PLPLOT
   char			plotfilename[MAXCHAR],
			*pstr;
   int			j,t, nplot, pnplot, pngindex, pngflag;
#endif

  QCALLOC(cp, int, prefs.ncplot_type);

/* Processing date and time if msg error present */
  if (msgerror)
    {
    thetime2 = time(NULL);
    tm = localtime(&thetime2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
        tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
        tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetime2, thetime);
    }

/* Username */
  psuser = pspath = pshost = NULL;
#ifdef HAVE_GETENV
  if (!(psuser=getenv("USERNAME")))	/* Cygwin,... */
    psuser = getenv("LOGNAME");		/* Linux,... */
  pspath = getenv("PWD");
  pshost = getenv("HOSTNAME");
#endif

  fprintf(file, " <RESOURCE ID=\"MetaData\" name=\"MetaData\">\n");
  fprintf(file, "  <DESCRIPTION>%s meta-data</DESCRIPTION>\n", BANNER);
  fprintf(file, "  <INFO name=\"QUERY_STATUS\" value=\"OK\" />\n");
  fprintf(file, "  <PARAM name=\"Software\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.title;meta.software\" value=\"%s\"/>\n",
	BANNER);
  fprintf(file, "  <PARAM name=\"Version\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.version;meta.software\" value=\"%s\"/>\n",
	MYVERSION);
  fprintf(file, "  <PARAM name=\"Soft_URL\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.ref.url;meta.software\" value=\"%s\"/>\n",
	WEBSITE);
  fprintf(file, "  <PARAM name=\"Soft_Auth\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.author;meta.software\" value=\"%s\"/>\n",
	"Emmanuel Bertin");
  fprintf(file, "  <PARAM name=\"Soft_Ref\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.bib.bibcode;meta.software\" value=\"%s\"/>\n",
	"2006ASPC..351..112B");
  fprintf(file, "  <PARAM name=\"NThreads\" datatype=\"int\""
	" ucd=\"meta.number;meta.software\" value=\"%d\"/>\n",
    	prefs.nthreads);
  fprintf(file, "  <PARAM name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.sdate_end);
  fprintf(file, "  <PARAM name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"time.event.end;meta.software\" value=\"%s\"/>\n",
	prefs.stime_end);
  fprintf(file, "  <PARAM name=\"Duration\" datatype=\"float\""
	" ucd=\"time.event;meta.software\" value=\"%.0f\" unit=\"s\"/>\n",
	prefs.time_diff);

  fprintf(file, "  <PARAM name=\"User\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	psuser);
  fprintf(file, "  <PARAM name=\"Host\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.curation\" value=\"%s\"/>\n",
	pshost);
  fprintf(file, "  <PARAM name=\"Path\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta.dataset\" value=\"%s\"/>\n",
	pspath);

  if (msgerror)
    {
    fprintf(file, "\n  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!! an Error occured"
	" !!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file,"  <PARAM name=\"Error_Msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\" value=\"%s\"/>\n", msgerror);
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n");
    fprintf(file, "  <!-- !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	"!!!!!!!!!!!!!!!!!!!! -->\n\n");
    }

  if (nfield_xml)
    {
    naxis = fields_xml[0]->naxis;
    lng = fields_xml[0]->lng;
    lat = fields_xml[0]->lat;
    }
  else
    naxis = lng = lat = 0;
  deg2arcsec = (lng!=lat) ? (DEG/ARCSEC) : 1.0;
  deg2arcmin = (lng!=lat) ? (DEG/ARCMIN) : 1.0;

/* Meta-data for each field (catalog) */
  fprintf(file, "  <TABLE ID=\"Fields\" name=\"Fields\">\n");
  fprintf(file, "   <DESCRIPTION>Data gathered by %s for every "
	" input catalog</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <!-- NFields may be 0"
	" if an error occurred early in the processing -->\n");
  fprintf(file, "   <PARAM name=\"NFields\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n", nfield_xml);

/* Test if PNG plots are being produced */
#ifdef HAVE_PLPLOT
  nplot = pnplot = pngflag = 0;
  for (j=0; j<prefs.ncplot_device; j++)
    if ((prefs.cplot_device[j] == CPLOT_PNG))
      {
      pngflag = 1;
      break;
      }

/* Check-plots */
  if (pngflag && (pngindex=cplot_check(CPLOT_ALLSKY)) != RETURN_ERROR)
    {
    strcpy(plotfilename, prefs.cplot_name[pngindex]);
    if (!(pstr = strrchr(plotfilename, '.')))
      pstr = plotfilename+strlen(plotfilename);
    sprintf(pstr, "_1.png");
    fprintf(file, "   <PARAM name=\"AllSkyPlot\" datatype=\"char\""
  	    " ucd=\"meta.id;meta.dataset\" value=\"%s\"/>\n",
	plotfilename);
    cp[nplot++] = pngindex;
    }
#endif
  fprintf(file, "   <FIELD name=\"Catalog_Number\" datatype=\"int\""
	" ucd=\"meta.record;meta.table;meta.file\"/>\n");
  fprintf(file, "   <FIELD name=\"Catalog_Name\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.table;meta.file\"/>\n");
  fprintf(file, "   <FIELD name=\"Image_Ident\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;obs.field\"/>\n");
  fprintf(file, "   <FIELD name=\"NExtensions\" datatype=\"int\""
        " ucd=\"meta.record\"/>\n");
  fprintf(file, "   <FIELD name=\"NAxis\" datatype=\"int\""
        " ucd=\"pos.wcs.naxis\"/>\n");
  fprintf(file, "   <FIELD name=\"Lng_Axis\" datatype=\"int\""
        " ucd=\"meta.id;pos.eq.ra\"/>\n");
  fprintf(file, "   <FIELD name=\"Lat_Axis\" datatype=\"int\""
        " ucd=\"meta.id;pos.eq.de\"/>\n");
  fprintf(file, "   <FIELD name=\"Ext_Header\" datatype=\"boolean\""
	" ucd=\"meta.code\"/>\n");
  fprintf(file, "   <FIELD name=\"NDetect\" datatype=\"int\""
        " ucd=\"meta.number;src\"/>\n");
  fprintf(file, "   <FIELD name=\"Group\" datatype=\"int\""
        " ucd=\"meta.id.parent;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"Astr_Instrum\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id.parent;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"Phot_Instrum\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id.parent;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"Photom_Flag\" datatype=\"boolean\""
	" ucd=\"meta.code;phot\"/>\n");
  fprintf(file, "   <FIELD name=\"Photom_Link\" datatype=\"boolean\""
	" ucd=\"meta.code;phot\"/>\n");
  fprintf(file, "   <FIELD name=\"Observation_Date\" datatype=\"double\""
	" ucd=\"time.epoch;obs.field\" unit=\"yr\"/>\n");
  fprintf(file, "   <FIELD name=\"Field_Coordinates\" datatype=\"double\""
	" arraysize=\"%d\" ucd=\"pos.eq;obs.image\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "deg":"pix");
  fprintf(file, "   <FIELD name=\"Pixel_Scale\" datatype=\"float\""
	" arraysize=\"%d\"  ucd=\"instr.pixel;obs.image;stat.mean\""
	" unit=\"%s\"/>\n", naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"Max_Radius\" datatype=\"float\""
	" ucd=\"phys.size.radius\" unit=\"%s\"/>\n", lng!=lat? "arcmin":"pix");
  fprintf(file, "   <FIELD name=\"ZeroPoint_Corr\" datatype=\"float\""
	" ucd=\"phot.mag;phot.calib;arith.zp\" unit=\"mag\"/>\n");
  if (prefs.match_flag)
    {
    fprintf(file, "   <!-- =========== MATCHing statistics =========== -->\n");
    fprintf(file, "   <FIELD name=\"DPixelScale\" datatype=\"float\""
	" ucd=\"instr.pixel;obs.image;arith.ratio\"/>\n");
    fprintf(file, "   <FIELD name=\"DPosAngle\" datatype=\"float\""
	" ucd=\"pos.posAng;obs.image;arith.diff\" unit=\"deg\"/>\n");
    fprintf(file, "   <FIELD name=\"AS_Contrast\" datatype=\"float\""
	" ucd=\"stat.correlation;arith.ratio\"/>\n");
    fprintf(file, "   <FIELD name=\"DX\" datatype=\"float\""
	" ucd=\"pos.eq;arith.diff\" unit=\"deg\"/>\n");
    fprintf(file, "   <FIELD name=\"DY\" datatype=\"float\""
	" ucd=\"pos.eq;arith.diff\" unit=\"deg\"/>\n");
    fprintf(file, "   <FIELD name=\"XY_Contrast\" datatype=\"float\""
	" ucd=\"stat.correlation;arith.ratio\"/>\n");
    fprintf(file, "   <FIELD name=\"Shear\" datatype=\"float\""
	" ucd=\"phys.size.axisRatio;obs.image\"/>\n");
    fprintf(file, "   <FIELD name=\"Shear_PosAngle\" datatype=\"float\""
	" ucd=\"pos.posAng;obs.image\" unit=\"deg\"/>\n");
    }

  fprintf(file,"   <!-- =========== Astrometric statistics =========== -->\n");
  fprintf(file, "   <FIELD name=\"Chi2_Internal\" datatype=\"float\""
	" ucd=\"stat.fit.chi2\"/>\n");
  fprintf(file, "   <FIELD name=\"NDeg_Internal\" datatype=\"int\""
	" ucd=\"stat.fit.dof\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Internal_HighSN\" datatype=\"float\""
	" ucd=\"stat.fit.chi2\"/>\n");
  fprintf(file, "   <FIELD name=\"NDeg_Internal_HighSN\" datatype=\"int\""
	" ucd=\"stat.fit.dof\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromOffset_Reference\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"arith.diff;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromSigma_Reference\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.stdev;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromCorr_Reference\" datatype=\"float\""
	" ucd=\"stat.correlation;pos.eq;obs.field\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Reference\" datatype=\"float\""
	" ucd=\"stat.fit.chi2\"/>\n");
  fprintf(file, "   <FIELD name=\"NDeg_Reference\" datatype=\"int\""
	" ucd=\"stat.fit.dof\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromOffset_Reference_HighSN\""
	" datatype=\"float\" arraysize=\"%d\""
	" ucd=\"arith.diff;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromSigma_Reference_HighSN\""
	" datatype=\"float\" arraysize=\"%d\""
	" ucd=\"stat.stdev;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromCorr_Reference_HighSN\""
	" datatype=\"float\" ucd=\"stat.correlation;pos.eq;obs.field\"/>\n");
  fprintf(file, "   <FIELD name=\"Chi2_Reference_HighSN\" datatype=\"float\""
	" ucd=\"stat.fit.chi2\"/>\n");
  fprintf(file, "   <FIELD name=\"NDeg_Reference_HighSN\" datatype=\"int\""
	" ucd=\"stat.fit.dof\"/>\n");
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (f=0; f<nfield_xml; f++)
    {
    field = fields_xml[f];
    fprintf(file, "    <TR>\n"
	"     <TD>%d</TD><TD>%s</TD><TD>%s</TD>\n"
	"     <TD>%d</TD><TD>%d</TD><TD>%d</TD><TD>%d</TD><TD>%c</TD>\n"
	"     <TD>%d</TD><TD>%d</TD><TD>A%d</TD><TD>P%d</TD><TD>%c</TD><TD>%c</TD>\n"
	"     <TD>%.9f</TD><TD>%.10g",
	field->fieldindex+1, field->rfilename, field->ident,
	field->nset,field->naxis,field->lng,field->lat,field->headflag?'T':'F',
	field->nsample, field->fgroup->no,
		field->astromlabel+1,field->photomlabel+1,
		field->photomflag==1? 'T':'F',field->photomflag==1? 'T':'F',
	field->epoch,
	field->meanwcspos[0]);
    for (d=1; d<field->naxis; d++)
      fprintf(file, " %.10g", field->meanwcspos[d]);
    fprintf(file, "</TD><TD>%.6g ", field->meanwcsscale[0]*deg2arcsec);
    for (d=1; d<field->naxis; d++)
      fprintf(file, " %.6g", field->meanwcsscale[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g</TD><TD>%.6g</TD>\n",
	field->maxradius*deg2arcmin, field->dmagzero);
    if (prefs.match_flag)
      fprintf(file, "     <TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n     "
	"<TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%.6g</TD>\n",
	field->match_dscale, field->match_dangle, field->match_asig,
	field->match_dlng, field->match_dlat,field->match_sig,
	field->match_shear, field->match_sangle);	
    fprintf(file, "     <TD>%.6g</TD><TD>%d</TD><TD>%.6g</TD><TD>%d</TD>\n",
	(double)field->chi2_int, (int)field->nchi2_int,
	(double)field->chi2_int_hsn, (int)field->nchi2_int_hsn);
    fprintf(file, "     <TD>%.6g", field->offset_ref[0]*deg2arcsec);
    for (d=1; d<field->naxis; d++)
      fprintf(file, " %.6g", field->offset_ref[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g", (double)field->sig_referr[0]*deg2arcsec);
    for (d=1; d<field->naxis; d++)
      fprintf(file, " %.6g", (double)field->sig_referr[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%d</TD><TD>%.6g",
	(double)field->sig_corr_ref, (double)field->chi2_ref,
	(int)field->nchi2_ref, field->offset_ref_hsn[0]*deg2arcsec);
    for (d=1; d<field->naxis; d++)
      fprintf(file, " %.6g", field->offset_ref_hsn[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g",(double)field->sig_referr_hsn[0]*deg2arcsec);
    for (d=1; d<field->naxis; d++)
      fprintf(file, " %.6g", (double)field->sig_referr_hsn[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%d</TD>\n    </TR>\n",
	(double)field->sig_corr_ref_hsn, (double)field->chi2_ref_hsn,
	(int)field->nchi2_ref_hsn);

    }
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");


  if (ngroup_xml)
    {
    naxis = fgroups_xml[0]->naxis;
    lng = fgroups_xml[0]->lng;
    lat = fgroups_xml[0]->lat;
    }
  else
    naxis = lng = lat = 0;
  deg2arcsec = (lng!=lat) ? (DEG/ARCSEC) : 1.0;
  deg2arcmin = (lng!=lat) ? (DEG/ARCMIN) : 1.0;

/* Meta-data for each group of fields */
  fprintf(file, "  <TABLE ID=\"FGroups\" name=\"FGroups\">\n");
  fprintf(file, "   <DESCRIPTION>Data gathered by %s for every "
	" group of fields found</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <PARAM name=\"NFGroups\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n", ngroup_xml);
  fprintf(file, "   <FIELD name=\"Name\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"Index\" datatype=\"int\""
	" ucd=\"meta.record;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NFields\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NAxis\" datatype=\"int\""
        " ucd=\"pos.wcs.naxis\"/>\n");
  fprintf(file, "   <FIELD name=\"Lng_Axis\" datatype=\"int\""
        " ucd=\"meta.id;pos.eq.ra\"/>\n");
  fprintf(file, "   <FIELD name=\"Lat_Axis\" datatype=\"int\""
        " ucd=\"meta.id;pos.eq.de\"/>\n");
  fprintf(file, "   <FIELD name=\"Field_Coordinates\" datatype=\"double\""
	" arraysize=\"%d\" ucd=\"pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "deg":"pix");
  fprintf(file, "   <FIELD name=\"Pixel_Scale\" datatype=\"float\""
	" arraysize=\"%d\"  ucd=\"instr.pixel;obs.field;stat.mean\""
	" unit=\"%s\"/>\n", naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"Max_Radius\" datatype=\"float\""
	" ucd=\"phys.size.radius;obs.field\" unit=\"%s\"/>\n",
	lng!=lat? "arcmin":"pix");
  fprintf(file,"   <!-- =========== Astrometric statistics =========== -->\n");
  fprintf(file, "   <FIELD name=\"AstRef_Catalog\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"AstRef_Band\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"instr.bandpass\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromSigma_Internal\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.stdev;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromCorr_Internal\" datatype=\"float\""
	" ucd=\"stat.correlation;pos.eq;obs.field\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromChi2_Internal\" datatype=\"float\""
	" ucd=\"stat.fit.chi2\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromNDets_Internal\" datatype=\"int\""
	" ucd=\"meta.number;src\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromSigma_Internal_HighSN\""
	" datatype=\"float\" arraysize=\"%d\""
	" ucd=\"stat.stdev;pos.eq;obs.field\""
	" unit=\"%s\"/>\n", naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromCorr_Internal_HighSN\""
	" datatype=\"float\" ucd=\"stat.correlation;pos.eq;obs.field\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromChi2_Internal_HighSN\""
	" datatype=\"float\" ucd=\"stat.fit.chi2\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromNDets_Internal_HighSN\""
	" datatype=\"int\" ucd=\"meta.number;src\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromOffset_Reference\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"arith.diff;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromSigma_Reference\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.stdev;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromCorr_Reference\" datatype=\"float\""
	" ucd=\"stat.correlation;pos.eq;obs.field\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromChi2_Reference\" datatype=\"float\""
	" ucd=\"stat.fit.chi2\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromNDets_Reference\" datatype=\"int\""
	" ucd=\"meta.number;src\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromOffset_Reference_HighSN\""
	" datatype=\"float\" arraysize=\"%d\""
	" ucd=\"arith.diff;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromSigma_Reference_HighSN\""
	" datatype=\"float\" arraysize=\"%d\""
	" ucd=\"stat.stDev;pos.eq;obs.field\" unit=\"%s\"/>\n",
	naxis, lng!=lat? "arcsec":"pix");
  fprintf(file, "   <FIELD name=\"AstromCorr_Reference_HighSN\""
	" datatype=\"float\" ucd=\"stat.correlation;pos.eq;obs.field\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromChi2_Reference_HighSN\""
	" datatype=\"float\" ucd=\"stat.fit.chi2\"/>\n");
  fprintf(file, "   <FIELD name=\"AstromNDets_Reference_HighSN\""
	" datatype=\"int\" ucd=\"meta.number;src\"/>\n");
  fprintf(file,"   <!-- =========== Photometric statistics =========== -->\n");
  fprintf(file, "   <PARAM name=\"NPhotInstru\" datatype=\"int\""
	" ucd=\"meta.number;meta.em\" value=\"%d\"/>\n", prefs.nphotinstrustr);
  fprintf(file, "   <FIELD name=\"PhotInstru_Name\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;instr.bandpass\"/>\n");
  fprintf(file, "   <FIELD name=\"PhotSigma_Internal\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.error;phot.mag\" unit=\"mag\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file, "   <FIELD name=\"PhotChi2_Internal\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.chi2;phot.mag\" unit=\"mag\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file, "   <FIELD name=\"PhotNDets_Internal\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.number;src\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file,"   <FIELD name=\"PhotSigma_Internal_HighSN\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.error;phot.mag\" unit=\"mag\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file,"   <FIELD name=\"PhotChi2_Internal_HighSN\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.chi2;phot.mag\" unit=\"mag\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file,"   <FIELD name=\"PhotNDets_Internal_HighSN\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.number;src\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file, "   <FIELD name=\"PhotSigma_Reference\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.error;phot.mag\" unit=\"mag\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file, "   <FIELD name=\"PhotChi2_Reference\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.chi2;phot.mag\" unit=\"mag\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file, "   <FIELD name=\"PhotNDets_Reference\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.number;src\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file,"   <FIELD name=\"PhotSigma_Reference_HighSN\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.error;phot.mag\" unit=\"mag\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file,"   <FIELD name=\"PhotChi2_Reference_HighSN\" datatype=\"float\""
	" arraysize=\"%d\" ucd=\"stat.chi2;phot.mag\" unit=\"mag\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file,"   <FIELD name=\"PhotNDets_Reference_HighSN\" datatype=\"int\""
	" arraysize=\"%d\" ucd=\"meta.number;src\"/>\n",
	prefs.nphotinstrustr);

/* Check-plots */
#ifdef HAVE_PLPLOT
  if (pngflag)
    {
    fprintf(file,
	"   <!-- =========== Check Plots for each group =========== -->\n");
    pnplot = nplot;
    if ((pngindex=cplot_check(CPLOT_FGROUPS)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"FgroupsPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_CHI2)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"Chi2Plot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_ADERROR1D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"IntErr1DimPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_ADERROR2D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"IntErr2DimPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_REFERROR1D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"RefErr1DimPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_REFERROR2D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"RefErr2DimPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_PHOTERROR)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"PhotErrPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_PHOTERRORVSMAG)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"PhotErrMagPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_PHOTZP)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"PhotZPPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_PHOTZP3D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"PhotZP3DPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_ASTRCOLSHIFT1D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"ColShiftPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_REFPROP)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"RefPropPlot\" datatype=\"char\""
  	    " arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    }
#endif

  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (g=0; g<ngroup_xml; g++)
    {
    fgroup = fgroups_xml[g];
    fprintf(file, "   <TR>\n     <TD>G%d</TD><TD>%d</TD><TD>%d</TD>\n"
	"     <TD>%d</TD><TD>%d</TD><TD>%d</TD>\n"
	"     <TD>%.10g",
	g+1, g+1, fgroup->nfield,
	fgroup->naxis, fgroup->lng, fgroup->lat,
	fgroup->meanwcspos[0]);
    for (d=1; d<fgroup->naxis; d++)
      fprintf(file, " %.10g", fgroup->meanwcspos[d]);
    fprintf(file, "</TD><TD>%.6g ", fgroup->meanwcsscale[0]*deg2arcsec);
    for (d=1; d<fgroup->naxis; d++)
      fprintf(file, " %.6g", fgroup->meanwcsscale[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g</TD>\n", fgroup->maxradius*deg2arcmin);
    fprintf(file, "     <TD>%s</TD><TD>%s</TD>\n",
	astrefcat[(int)prefs.astrefcat].name,
	astrefcat[(int)prefs.astrefcat].bandname ?
		astrefcat[(int)prefs.astrefcat].bandname : "");
    fprintf(file, "     <TD>%.6g", fgroup->sig_interr[0]*deg2arcsec);
    for (d=1; d<fgroup->naxis; d++)
      fprintf(file, " %.6g", fgroup->sig_interr[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%d</TD>\n     <TD>%.6g",
	fgroup->sig_corr_int, fgroup->chi2_int, fgroup->nintmatch,
	fgroup->sig_interr_hsn[0]*deg2arcsec);
    for (d=1; d<fgroup->naxis; d++)
      fprintf(file, " %.6g", fgroup->sig_interr_hsn[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%d</TD>\n     <TD>%.6g",
	fgroup->sig_corr_int_hsn, fgroup->chi2_int_hsn, fgroup->nintmatch_hsn,
	fgroup->offset_ref[0]*deg2arcsec);
    for (d=1; d<fgroup->naxis; d++)
      fprintf(file, " %.6g", fgroup->offset_ref[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g", fgroup->sig_referr[0]*deg2arcsec);
    for (d=1; d<fgroup->naxis; d++)
      fprintf(file, " %.6g", fgroup->sig_referr[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%d</TD>\n     <TD>%.6g",
	fgroup->sig_corr_ref, fgroup->chi2_ref, fgroup->nrefmatch,
	fgroup->offset_ref_hsn[0]*deg2arcsec);
    for (d=1; d<fgroup->naxis; d++)
      fprintf(file, " %.6g", fgroup->offset_ref_hsn[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g", fgroup->sig_referr_hsn[0]*deg2arcsec);
    for (d=1; d<fgroup->naxis; d++)
      fprintf(file, " %.6g", fgroup->sig_referr_hsn[d]*deg2arcsec);
    fprintf(file, "</TD><TD>%.6g</TD><TD>%.6g</TD><TD>%d</TD>\n",
	fgroup->sig_corr_ref_hsn, fgroup->chi2_ref_hsn, fgroup->nrefmatch_hsn);
    fprintf(file, "     <TD>P1");
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, ",P%d", i+1);
    fprintf(file, "</TD>\n     <TD>%.6g", fgroup->sig_intmagerr[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %.6g", fgroup->sig_intmagerr[i]);
    fprintf(file, "</TD><TD>%.6g", fgroup->chi2_intmag[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %.6g", fgroup->chi2_intmag[i]);
    fprintf(file, "</TD><TD>%d", fgroup->nintmagmatch[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %d", fgroup->nintmagmatch[i]);
    fprintf(file, "</TD>\n     <TD>%.6g", fgroup->sig_intmagerr_hsn[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %.6g", fgroup->sig_intmagerr_hsn[i]);
    fprintf(file, "</TD><TD>%.6g", fgroup->chi2_intmag_hsn[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %.6g", fgroup->chi2_intmag_hsn[i]);
    fprintf(file, "</TD><TD>%d", fgroup->nintmagmatch_hsn[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %d", fgroup->nintmagmatch_hsn[i]);
    fprintf(file, "</TD>\n     <TD>%.6g", fgroup->sig_refmagerr[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %.6g", fgroup->sig_refmagerr[i]);
    fprintf(file, "</TD><TD>%.6g", fgroup->chi2_refmag[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %.6g", fgroup->chi2_refmag[i]);
    fprintf(file, "</TD><TD>%d", fgroup->nrefmagmatch[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %d", fgroup->nrefmagmatch[i]);
    fprintf(file, "</TD>\n     <TD>%.6g", fgroup->sig_refmagerr_hsn[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %.6g", fgroup->sig_refmagerr_hsn[i]);
    fprintf(file, "</TD><TD>%.6g", fgroup->chi2_refmag_hsn[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %.6g", fgroup->chi2_refmag_hsn[i]);
    fprintf(file, "</TD><TD>%d", fgroup->nrefmagmatch_hsn[0]);
    for (i=1; i<prefs.nphotinstrustr; i++)
      fprintf(file, " %d", fgroup->nrefmagmatch_hsn[i]);

/*-- Checkplots */
#ifdef HAVE_PLPLOT
    if (pngflag)
      {
      for (t=pnplot; t<nplot; t++)
        {
        strcpy(plotfilename, prefs.cplot_name[cp[t]]);
        if (!(pstr = strrchr(plotfilename, '.')))
          pstr = plotfilename+strlen(plotfilename);
        sprintf(pstr, "_%d.png", g+1);
        fprintf(file, "</TD>\n     <TD>%s",plotfilename);
        }        
      }
#endif
    fprintf(file, "</TD>\n    </TR>\n");
    }
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Meta-data for each astrometric instrument */
  fprintf(file, "  <TABLE ID=\"Astrometric_Instruments\""
	" name=\"Astrometric_Instruments\">\n");
  fprintf(file, "   <DESCRIPTION>Data gathered by %s for every "
	" astrometric instrument identified</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <PARAM name=\"NAstromInstru\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	prefs.nastrinstrustr);
  fprintf(file, "   <FIELD name=\"Name\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"Index\" datatype=\"int\""
	" ucd=\"meta.record;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NFields\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NExtensions\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NKeys\" datatype=\"int\""
	" ucd=\"meta.number\"/>\n");
  fprintf(file, "   <FIELD name=\"Keys\" datatype=\"*\""
	" ucd=\"meta.note\"/>\n");

/* Check-plots */
#ifdef HAVE_PLPLOT
  if (pngflag)
    {
    pnplot = nplot;
    if ((pngindex=cplot_check(CPLOT_DISTORT)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"DistPlot\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_REFSYSMAP2D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"RefSysPlot\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_PIXERROR1D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"PixErr1DimPlot\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_SUBPIXERROR1D)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"SubPixErr1DimPlot\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    if ((pngindex=cplot_check(CPLOT_SHEAR)) != RETURN_ERROR)
      {
      fprintf(file, "   <FIELD name=\"ShearPlot\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
      cp[nplot++] = pngindex;
      }
    }
#endif

  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (i=0; i<prefs.nastrinstrustr; i++)
    {
    len = fitsfind(prefs.astrinstrustr[i], "END     ");
    f2 = 0;
    for (g=0;g<ngroup_xml; g++)
      for (f=0;f<fgroups_xml[g]->nfield;f++)
        if (fgroups_xml[g]->field[f]->astromlabel==i)
          f2++;
    fprintf(file, "    <TR>\n"
	"     <TD>A%d</TD><TD>%d</TD><TD>%d</TD><TD>%d</TD>\n"
	"     <TD>%d</TD><TD>%32.32s",
	i+1, i+1, f2, prefs.nastrinstruext[i],
	len, prefs.astrinstrustr[i]);
    for (l=1; l<len; l++)
      fprintf(file, ",%32.32s", prefs.astrinstrustr[i]+l*80);

/*-- Check-plots */
#ifdef HAVE_PLPLOT
    if (pngflag)
      {
      for (t=pnplot; t<nplot; t++)
        {
        strcpy(plotfilename, prefs.cplot_name[cp[t]]);
        if (!(pstr = strrchr(plotfilename, '.')))
            pstr = plotfilename+strlen(plotfilename);
        sprintf(pstr, "_%d.png", i+1);
        fprintf(file, "</TD>\n     <TD>%s",plotfilename);        
        }
      }
#endif
    fprintf(file, "</TD>\n    </TR>\n");
    }
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");


/* Meta-data for each photometric instrument */
  fprintf(file, "  <TABLE ID=\"Photometric_Instruments\""
	" name=\"Photometric_Instruments\">\n");
  fprintf(file, "   <DESCRIPTION>Data gathered by %s for every "
	" photometric instrument identified</DESCRIPTION>\n", BANNER);
  fprintf(file, "   <PARAM name=\"NPhotomInstru\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\" value=\"%d\"/>\n",
	prefs.nphotinstrustr);
  fprintf(file, "   <FIELD name=\"Name\" datatype=\"char\""
	" arraysize=\"*\" ucd=\"meta.id;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"Index\" datatype=\"int\""
	" ucd=\"meta.record;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"NFields\" datatype=\"int\""
	" ucd=\"meta.number;meta.dataset\"/>\n");
  fprintf(file, "   <FIELD name=\"MagZeroPoint_Output\" datatype=\"float\""
	" ucd=\"phot.mag;phot.calib;arith.zp\" unit=\"mag\"/>\n");
  fprintf(file, "   <FIELD name=\"NKeys\" datatype=\"int\""
	" ucd=\"meta.number\"/>\n");
  fprintf(file, "   <FIELD name=\"Keys\" datatype=\"*\""
	" ucd=\"meta.note\"/>\n");
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (i=0; i<prefs.nphotinstrustr; i++)
    {
    len = fitsfind(prefs.photinstrustr[i], "END     ");
    f2 = 0;
    for (g=0;g<ngroup_xml; g++)
      for (f=0;f<fgroups_xml[g]->nfield;f++)
        if (fgroups_xml[g]->field[f]->photomlabel==i)
          f2++;
    fprintf(file, "    <TR>\n"
	"     <TD>P%d</TD><TD>%d</TD><TD>%d</TD><TD>%.6g</TD>\n"
	"     <TD>%d</TD><TD>%32.32s",
	i+1, i+1, f2, prefs.magzero_out[i],
	len, prefs.photinstrustr[i]);
    for (l=1; l<len; l++)
      fprintf(file, ",%32.32s", prefs.photinstrustr[i]+l*80);
    fprintf(file, "</TD>\n    </TR>\n");
    }
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Warnings */
  fprintf(file, "  <TABLE ID=\"Warnings\" name=\"Warnings\">\n");
  fprintf(file,
	"   <DESCRIPTION>%s warnings (limited to the last %d)</DESCRIPTION>\n",
	BANNER, WARNING_NMAX);
  fprintf(file, "   <FIELD name=\"Date\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Time\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta;time.event.end\"/>\n");
  fprintf(file, "   <FIELD name=\"Msg\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"meta\"/>\n");
  fprintf(file, "   <DATA><TABLEDATA>\n");
  for (str = warning_history(); *str; str = warning_history())
    fprintf(file, "    <TR><TD>%10.10s</TD><TD>%8.8s</TD><TD>%s</TD></TR>\n",
	str, str+11, str+22);
  fprintf(file, "   </TABLEDATA></DATA>\n");
  fprintf(file, "  </TABLE>\n");

/* Configuration file */
  fprintf(file, "  <RESOURCE ID=\"Config\" name=\"Config\">\n");
  fprintf(file, "   <DESCRIPTION>%s configuration</DESCRIPTION>\n", BANNER);
  fprintf(file,
	"   <PARAM name=\"Command_Line\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param\" value=\"%s",
	prefs.command_line[0]);
  for (n=1; n<prefs.ncommand_line; n++)
    fprintf(file, " %s", prefs.command_line[n]);
  fprintf(file, "\"/>\n");
  fprintf(file,
	"   <PARAM name=\"Prefs_Name\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"obs.param;meta.file\" value=\"%s\"/>\n",
	prefs.prefs_name);

  free(cp);

  if (!msgerror)
    {
/*-- Field grouping */
    write_xmlconfigparam(file, "FGroup_Radius", "arcmin",
		"phot.calib;obs.param","%.6g");
/*-- Reference catalogs */
    write_xmlconfigparam(file, "Ref_Server", "", "meta.ref.url", "%s");
    write_xmlconfigparam(file, "Ref_Port", "", "meta.ref.url", "%d");
    write_xmlconfigparam(file, "CDSClient_Exec", "", "meta;meta.file", "%s");
    write_xmlconfigparam(file, "AstRef_Catalog", "",
		"meta.code;meta.file", "%s");
    write_xmlconfigparam(file, "AstRef_Band", "", "instr.bandpass", "%s");
    write_xmlconfigparam(file, "AstRefCat_Name", "",
		"meta.id;meta.file;meta.dataset", "%s");
    write_xmlconfigparam(file, "AstRefCent_Keys", "", "meta.id;pos.eq", "%s");
    write_xmlconfigparam(file,"AstRefErr_Keys", "",
		"meta.id;pos.errorEllipse","%s");
    write_xmlconfigparam(file,"AstRefMag_Key", "",
		"meta.id;phot.mag","%s");
    write_xmlconfigparam(file,"Save_RefCatalog", "",
		"meta.code","%c");
    write_xmlconfigparam(file, "RefOut_CatPath", "",
		"meta.id;meta.file","%s");
/*-- Merged output catalogs */
    write_xmlconfigparam(file, "MergedOutCat_Name", "",
		"meta.id;meta.file;meta.dataset", "%s");
    write_xmlconfigparam(file, "MergedOutCat_Type", "",
		"meta.code;meta.dataset", "%s");
/*-- Full output catalogs */
    write_xmlconfigparam(file, "FullOutCat_Name", "",
		"meta.id;meta.file;meta.dataset", "%s");
    write_xmlconfigparam(file, "FullOutCat_Type", "",
		"meta.code;meta.dataset", "%s");
/*-- Pattern matching */
    write_xmlconfigparam(file, "Match", "",
		"meta.code;obs.param","%c");
    write_xmlconfigparam(file, "Match_NMax", "",
		"meta.number;src;obs.param", "%d");
    write_xmlconfigparam(file, "PixScale_MaxErr", "",
		"instr.precision;instr.pixel;obs.param", "%.6g");
    write_xmlconfigparam(file, "PosAngle_MaxErr", "deg",
		"instr.precision;pos.posAng;obs.param", "%.6g");
    write_xmlconfigparam(file, "Position_MaxErr", "arcmin",
		"instr.precision;pos.eq;obs.param", "%.6g");
    write_xmlconfigparam(file, "Match_Resol", "arcsec",
		"stat.param;phys.angSize","%.6g");
    write_xmlconfigparam(file, "Match_Flipped", "",
		"instr.param;obs.param", "%c");
    write_xmlconfigparam(file, "Mosaic_Type", "",
		"instr.param;obs.param", "%s");
    write_xmlconfigparam(file, "FixFocalPlane_NMin", "",
		"meta.number;src;obs.param", "%d");
/*-- Cross-identification */
    write_xmlconfigparam(file,
		"CrossId_Radius", "arcsec",
		"pos.angDistance;stat.max", "%.6g");
/*-- Astrometric solution */
    write_xmlconfigparam(file, "Solve_Astrom", "", "meta.code;obs.param","%c");
    write_xmlconfigparam(file,"AstrInstru_Key", "",
		"meta.id;instr.setup;pos","%s");
    write_xmlconfigparam(file, "Stability_Type", "",
		"instr.param;obs.param", "%s");
    write_xmlconfigparam(file, "Centroid_Keys", "",
		"meta.id;pos.barycenter;src", "%s");
    write_xmlconfigparam(file, "CentroidErr_Keys", "",
		"meta.id;pos.errorEllipse;src", "%s");
    write_xmlconfigparam(file, "Distort_Keys", "",
		"meta.id;stat.fit.param;src", "%s");
    write_xmlconfigparam(file, "Distort_Groups", "",
		"meta.id;stat.fit.param", "%d");
    write_xmlconfigparam(file,"Distort_Degrees", "",
		"meta.id;stat.fit.param", "%d");
    write_xmlconfigparam(file,"Astref_Weight", "",
			 "meta.id;stat.fit.param", "%.6g");
    write_xmlconfigparam(file, "AstrClip_NSigma", "",
		"meta.id;stat.param;pos", "%.6g");
    write_xmlconfigparam(file, "Correct_ColourShifts", "",
		"meta.code;obs.param", "%c");
/*-- Photometric solution */
    write_xmlconfigparam(file, "Solve_Photom", "",
		"meta.code;phot.calib;obs.param", "%c");
    write_xmlconfigparam(file, "MagZero_Out", "mag",
		"arith.zp;phot.calib;obs.param", "%.6g");
    write_xmlconfigparam(file, "MagZero_IntErr", "mag",
		"stat.error;phot.mag;obs.param", "%.6g");
    write_xmlconfigparam(file, "MagZero_RefErr", "mag",
		"stat.error;phot.calib;obs.param","%.6g");
    write_xmlconfigparam(file, "PhotInstru_Key", "",
		"meta.id;instr.setup;phot.calib","%s");
    write_xmlconfigparam(file, "MagZero_Key", "",
		"meta.id;arith.zp;phot.calib","%s");
    write_xmlconfigparam(file, "ExpoTime_Key", "",
		"meta.id;time.duration.expo","%s");
    write_xmlconfigparam(file, "AirMass_Key", "",
		"meta.id;obs.airmass","%s");
    write_xmlconfigparam(file,"Extinct_Key", "",
		"meta.id;obs.atmos.extinction","%s");
    write_xmlconfigparam(file,"PhotomFlag_Key", "",
		"meta.id;obs.atmos", "%s");
    write_xmlconfigparam(file,"PhotFlux_Key", "",
		"meta.id;phot.flux;src", "%s");
    write_xmlconfigparam(file,"PhotFluxerr_Key", "",
		"meta.id;phot.flux;src", "%s");
    write_xmlconfigparam(file, "PhotClip_NSigma", "",
		"meta.id;stat.param;phot", "%d");
/*-- Check-plots */
    write_xmlconfigparam(file, "CheckPlot_CKey", "", "meta.id;meta.code", "%s");
    write_xmlconfigparam(file, "CheckPlot_Dev", "", "meta.code", "%s");
    write_xmlconfigparam(file, "CheckPlot_Res", "", "meta.number;meta", "%d");
    write_xmlconfigparam(file, "CheckPlot_AntiAlias", "", "meta.code", "%c");
    write_xmlconfigparam(file, "CheckPlot_Type", "", "meta.code", "%s");
    write_xmlconfigparam(file, "CheckPlot_Name", "", "meta.id;meta.file", "%s");
    write_xmlconfigparam(file, "CheckImage_Type", "", "meta.code", "%s");
    write_xmlconfigparam(file, "CheckImage_Name", "",
		"meta.id;meta.file;meta.fits", "%s");
/*-- Miscellaneous */
    write_xmlconfigparam(file, "SN_Thresholds", "",
		"stat.snr;phot.flux;obs.param", "%.6g");
    write_xmlconfigparam(file, "Flags_Mask", "", "meta.code.qual", "%d");
    write_xmlconfigparam(file, "WeightFlags_Mask", "", "meta.code.qual", "%d");
    write_xmlconfigparam(file, "ImaFlags_Mask", "", "meta.code.qual", "%d");
    write_xmlconfigparam(file, "AHeader_Global", "", "meta.id;meta.file","%s");
    write_xmlconfigparam(file, "AHeader_Suffix", "",
		"meta.id.part;meta.file","%s");
    write_xmlconfigparam(file, "Header_Suffix", "",
		"meta.id.part;meta.file","%s");
    write_xmlconfigparam(file, "Header_Type", "", "meta.code;meta.file","%s");
    write_xmlconfigparam(file, "Verbose_Type", "", "meta.code","%s");
    write_xmlconfigparam(file, "Write_XML", "", "meta.code","%s");
    write_xmlconfigparam(file, "NThreads", "",
		"meta.number;meta.software", "%d");
    }

  fprintf(file, "  </RESOURCE>\n");
  fprintf(file, " </RESOURCE>\n");

  return RETURN_OK;
  }


/****** write_xmlerror ******************************************************
PROTO	int	write_xmlerror(char *msgerror)
PURPOSE	Save meta-data to a simplified XML file in case of a catched error
INPUT	a character string.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	02/10/2007
 ***/
void	write_xmlerror(char *filename, char *msgerror)
  {
   FILE			*file;
   int			pipe_flag;

  pipe_flag = 0;
  if (!strcmp(filename, "STDOUT"))
    {
    file = stdout;
    pipe_flag = 1;
    }
  else if (!(file = fopen(filename, "w")))
    return;

  write_xml_header(file);
  write_xml_meta(file, msgerror);

  fprintf(file, "</RESOURCE>\n");
  fprintf(file, "</VOTABLE>\n");

  if (!pipe_flag)
    fclose(file);

  return;
  }


/****** write_xmlconfigparam **************************************************
PROTO	int write_xmlconfigparam(FILE *file, char *name, char *unit,
		char *ucd, char *format)
PURPOSE	Write to a VO-table the configuration parameters.
INPUT	Output stream (file) pointer,
	Name of the parameter keyword,
	unit,
	UCD string,
	printf() format to use in "value".
OUTPUT	RETURN_OK if the keyword exists, RETURN_ERROR otherwise.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	27/04/2010
 ***/
int	write_xmlconfigparam(FILE *file, char *name, char *unit,
		 char *ucd, char *format)
  {
   char		value[MAXCHAR], uunit[MAXCHAR];
   int		i,j,n;

  for (i=0; key[i].name[0] && cistrcmp(name, key[i].name, FIND_STRICT); i++);
  if (!key[i].name[0])
    return RETURN_ERROR;

  if (*unit)
    sprintf(uunit, " unit=\"%s\"", unit);
  else
    *uunit = '\0';
  switch(key[i].type)
    {
    case P_FLOAT:
      sprintf(value, format, *((double *)key[i].ptr));
      fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, uunit, ucd, value);
      break;
    case P_FLOATLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, format, ((double *)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, uunit, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, format, ((double *)key[i].ptr)[j]);
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, uunit, ucd);
      break;
    case P_INT:
      sprintf(value, format, *((int *)key[i].ptr));
      fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"int\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, uunit, ucd, value);
      break;
    case P_INTLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, format, ((int *)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"int\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, uunit, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, format, ((int *)key[i].ptr)[j]);
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\"%s datatype=\"double\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, uunit, ucd);
      break;
    case P_BOOL:
      sprintf(value, "%c", *((int *)key[i].ptr)? 'T':'F');
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_BOOLLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        sprintf(value, "%c", ((int *)key[i].ptr)[0]? 'T':'F');
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
		" arraysize=\"%d\" ucd=\"%s\" value=\"%s",
		name, n, ucd, value);
        for (j=1; j<n; j++)
          {
          sprintf(value, "%c", ((int *)key[i].ptr)[j]? 'T':'F');
          fprintf(file, " %s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"boolean\""
		" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    case P_STRING:
      strcpy(value, (char *)key[i].ptr);
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_STRINGLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        strcpy(value, ((char **)key[i].ptr)[0]);
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"%s",
		name, ucd, value);
        for (j=1; j<n; j++)
          {
          strcpy(value, ((char **)key[i].ptr)[j]);
          fprintf(file, ",%s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    case P_KEY:
      strcpy(value, key[i].keylist[*((int *)key[i].ptr)]);
      fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\" arraysize=\"*\""
	" ucd=\"%s\" value=\"%s\"/>\n",
	name, ucd, value);
      break;
    case P_KEYLIST:
      n = *(key[i].nlistptr);
      if (n)
        {
        strcpy(value, key[i].keylist[((int *)key[i].ptr)[0]]);
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"%s",
		name, ucd, value);
        for (j=1; j<n; j++)
          {
          strcpy(value, key[i].keylist[((int *)key[i].ptr)[j]]);
          fprintf(file, ",%s", value);
          }
        fprintf(file, "\"/>\n");
        }
      else
        fprintf(file, "   <PARAM name=\"%s\" datatype=\"char\""
		" arraysize=\"*\" ucd=\"%s\" value=\"\"/>\n",
		name, ucd);
      break;
    default:
        error(EXIT_FAILURE, "*Internal Error*: Type Unknown",
		" in write_xmlconfigparam()");
    }

  return RETURN_OK;
  }

