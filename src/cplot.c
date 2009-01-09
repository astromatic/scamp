/*
 				cplot.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	CPlot sublib
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:       Call a plotting library (PLPlot).
*
*	Last modify:	28/08/2008
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

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

struct {cplotdevenum device; char *devname; char *extension;}
		cplot_device[] = {{CPLOT_NULL, "null", ""},
		{CPLOT_XWIN, "xwin",":0"},
		{CPLOT_TK, "tk",""},
		{CPLOT_XTERM, "xterm",""},
		{CPLOT_PLMETA, "plmeta",".plm"},
		{CPLOT_PS, "ps", ".ps"},
		{CPLOT_PSC, "psc", ".ps"},
		{CPLOT_XFIG, "xfig", ".fig"},
		{CPLOT_LJIIP, "ljiip", ".lj"},
		{CPLOT_LJHPGL, "lj_hpgl",".hpg"},
		{CPLOT_IMP, "imp",".imp"},
		{CPLOT_PBM, "pbm",".pbm"},
		{CPLOT_PNG, "png",".png"},
		{CPLOT_JPEG, "jpeg",".jpg"},
		{CPLOT_PSTEX, "pstex", ".ps"},
		{CPLOT_AQT, "aqt", ""},
		{CPLOT_PDF, "pdf", ".pdf"},
		{CPLOT_SVG, "svg", ".svg"},
		{CPLOT_NULL,"",""}};

int	plotnum[CPLOT_NTYPES];
int	plotdev[CPLOT_NTYPES];
char	plotfilename[MAXCHAR];
int	plotaaflag;
 
/****** cplot_check ***********************************************************
PROTO	int cplot_check(cplotenum cplottype)
PURPOSE	Check that the specified check-plot type has been requested by user,
	by returning its index in CPLOT_TYPE keyword list, or RETURN_ERROR
	(-1) otherwise.
INPUT	Check-plot type.
OUTPUT	Index in CPLOT_TYPE keyword list, or RETURN_ERROR (-1) otherwise.
[5~[5~NOTES	Uses the global preferences.
5~AUTHOR	E. Bertin (IAP)
VERSION	18/10/2002
 ***/
int	 cplot_check(cplotenum cplottype)

  {
   int		i;

  for (i=0; i<prefs.ncplot_type; i++)
    if (cplottype == prefs.cplot_type[i])
      return i;

  return RETURN_ERROR;
  }


/****** cplot_init ***********************************************************
PROTO	int cplot_init(int nx, int ny, cplotenum cplottype)
PURPOSE	Initialize a check plot.
INPUT	Number of plots along the x axis,
	number of plots along the y axis,
	plot type,
	device number.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	19/01/2005
 ***/
int	cplot_init(int nx, int ny, cplotenum cplottype)
  {
   char		str[MAXCHAR],
		*pstr;
   int		j, num, cval, dev, argc;

/* Check that plot was requested */
  cval = cplot_check(cplottype);
/* return otherwise */
  if (cval == RETURN_ERROR)
    return RETURN_ERROR;
  dev = plotdev[cval]++;
/* Run convert to antialias the check-plot image (quick and dirty) */
  if (prefs.cplot_antialiasflag && dev>0
	&& (prefs.cplot_device[dev-1] == CPLOT_PNG
	|| prefs.cplot_device[dev-1] == CPLOT_JPEG))
    {
    sprintf(str, "convert -geometry \"%dx%d\" -antialias %s %s",
	prefs.cplot_res[0]? prefs.cplot_res[0] : CPLOT_DEFRESX,
	prefs.cplot_res[1]? prefs.cplot_res[1] : CPLOT_DEFRESY,
	plotfilename, plotfilename);
    system(str);
    }

  if (dev>=prefs.ncplot_device)
    return RETURN_ERROR;
  num = plotnum[cval]+1;

/* Provide the right extension to the output filename */
  strcpy(plotfilename, prefs.cplot_name[cval]);
  if (!(pstr = strrchr(plotfilename, '.')))
    pstr = plotfilename+strlen(plotfilename);

  for (j=0; *(cplot_device[j].devname); j++)
    if (prefs.cplot_device[dev]==cplot_device[j].device)
      break;

  if (cplot_device[j].device != CPLOT_XWIN)
    {
    sprintf(pstr, "_%d%s", num, cplot_device[j].extension);
    plsfnam(plotfilename);		/* file name */
    }
  plssub(nx, ny);
  plsdev(cplot_device[j].devname);

  plotaaflag = 0;
  if (cplot_device[j].device == CPLOT_PNG
	|| cplot_device[j].device == CPLOT_JPEG)
    {
/*-- Set custom resolutions */
    if (prefs.cplot_antialiasflag)
      {
/*---- Oversample for antialiasing */
      sprintf(str, "%dx%d",
	(prefs.cplot_res[0]? prefs.cplot_res[0] : CPLOT_DEFRESX)*CPLOT_AAFAC,
	(prefs.cplot_res[1]? prefs.cplot_res[1] : CPLOT_DEFRESY)*CPLOT_AAFAC);
      plsetopt("-geometry", str);
      plotaaflag = 1;
      }
    else if (prefs.cplot_res[0])
      {
/*---- No oversampling */
      sprintf(str, "%dx%d", prefs.cplot_res[0], prefs.cplot_res[1]);
      plsetopt("-geometry", str);
      }
    plsetopt("-drvopt","24bit");
    }
  else
    {
/*-- Small hack to reset driver options */
    argc = 0;
    plParseOpts(&argc, NULL, PL_PARSE_NOPROGRAM);
    }

  plfontld(1);
  plscolbg(255,255,255);	/* Force the background colour to white */
  plscol0(15, 0,0,0);		/* Force the foreground colour to black */
  plscol0(3, 0,200,0);		/* Force the green colour to darken */
  plscol0(7, 128,128,128);	/* Force the midground colour to grey */
  plscol0(8, 64,0,0);		/* Force the brown colour to darken */
  plinit();

  return RETURN_OK;
  }


/****** cplot_end ************************************************************
PROTO	int cplot_end(void)
PURPOSE	Terminate a CPlot (and save additional plots if required).
INPUT	-.
OUTPUT	RETURN_OK if everything went fine, RETURN_ERROR otherwise.
NOTES	.
AUTHOR	E. Bertin (IAP)
VERSION	03/01/2004
 ***/
int	cplot_end(cplotenum cplottype)
  {
   int		cval;

/* Check that plot was requested */
  cval = cplot_check(cplottype);
/* return otherwise */
  if (cval == RETURN_ERROR)
    return RETURN_ERROR;

  plotdev[cval] = 0;
  plotnum[cval]++;

/* Copy plot to other devices if required *
  pleop();
  plgstrm(&cur_strm);			* get current stream *
  for (i=1; i<prefs.ncplot_device; i++)
    {
    plmkstrm(&new_strm);		* create a new one * 
    strcpy(filename, prefs.cplot_name[cval]);
    if (!(pstr = strrchr(filename, '.')))
      pstr = filename+strlen(filename);
    for (j=0; *(cplot_device[j].devname); j++)
      if (prefs.cplot_device[i]==cplot_device[j].device)
        break;
    if (cplot_device[j].device != CPLOT_XWIN)
      {
      sprintf(pstr, "_%d%s", num, cplot_device[j].extension);
      plsfnam(filename);		* file name *
      }
    plsdev(cplot_device[j].devname);	* device *
    plcpstrm(cur_strm, 0);	* copy old stream parameters to new stream *
    plreplot();			* do the save by replaying the plot buffer *
    plend1();			* finish the device *
    }

  plsstrm(cur_strm);		* return to previous stream *
  plend();
*/

  return RETURN_OK;
  }


/****** cplot_degtosexal ******************************************************
PROTO	char	*cplot_degtosexal(char *str, double alpha, double step)
PURPOSE	Convert degrees to hh mm ss.xx coordinates in PLPLOT string format.
INPUT	Pointer to the output character string,
	alpha coordinate in degrees,
	step (precision) in degrees.
OUTPUT	Pointer to the output character string.
NOTES	At least 30 bytes must have been allocated for str.
AUTHOR	E. Bertin (IAP)
VERSION	20/07/2004
 ***/
char	*cplot_degtosexal(char *str, double alpha, double step)
  {
   int		hh, mm;
   double	dh, dm, ss;

  if (alpha<0.0 || alpha >360.0)
    alpha = fmod(alpha+360.0, 360.0);
  alpha += 1e-8;
  hh = (int)(dh = alpha/15.0);
  mm = (int)(dm = 60.0*(dh - hh));
  ss = 60.0*(dm - mm);
  if (step*DEG<0.999*15.0*ARCSEC)
    sprintf(str,"%02d#uh#d%02d#um#d%05.2f#us#d", hh, mm, ss);
  else if (step*DEG<0.999*15.0*ARCMIN)
    sprintf(str,"%02d#uh#d%02d#um#d%02.0f#us#d", hh, mm, ss);
  else if (step*DEG<0.999*15.0*DEG)
    sprintf(str,"%02d#uh#d%02d#um#d", hh, mm);
  else
    sprintf(str,"%02d#uh#d", hh);

  return str;
  }


/****** cplot_degtosexde *****************************************************
PROTO	char	*cplot_degtosexde(char *str, double delta, double step)
PURPOSE	Convert degrees to dd mm ss.xx coordinates in PLPLOT string format.
INPUT	Pointer to the output character string,
	delta coordinate in degrees,
	step (precision) in degrees.
OUTPUT	Pointer to the output character string.
NOTES	At least 30 bytes must have been allocated for str.
AUTHOR	E. Bertin (IAP)
VERSION	26/10/2005
 ***/
char	*cplot_degtosexde(char *str, double delta, double step)
  {
   char		sign;
   double	ddm, ds;
   int		dd, dm;

  sign = delta<1e-8?'-':'+';
  if (delta<1e-8)
    delta = -delta;
  delta += 1e-8;
  if (delta<-90.0)
    delta = -90.0;
  else if (delta>90.0)
    delta = 90.0;
  dd = (int)delta;
  dm = (int)(ddm = (60.0*(delta - dd)));
  ds = 60.0*(ddm - dm);
  if (step*DEG<0.999*ARCSEC)
    sprintf(str,"%c%02d#(2218)%02d#(2216)%05.2f#(2216)#(2216)", sign,dd,dm,ds);
  else if (step*DEG<0.999*ARCMIN)
    sprintf(str,"%c%02d#(2218)%02d#(2216)%2.0f#(2216)#(2216)", sign,dd,dm,ds);
  else if (step*DEG<0.999*DEG)
    sprintf(str,"%c%02d#(2218)%02d#(2216)", sign,dd,dm);
  else
    sprintf(str,"%c%02d#(2218)", sign,dd);

  return str;
  }

