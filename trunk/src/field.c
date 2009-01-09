 /*
				field.c

*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	Part of:	SCAMP
*
*	Author:		E.BERTIN (IAP)
*
*	Contents:	Handle catalogs.
*
*	Last modify:	12/12/2006
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <pthread.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "header.h"
#include "wcs/wcs.h"
#include "field.h"
#include "prefs.h"
#include "samples.h"
#ifdef USE_THREADS
#include "threads.h"
#endif

/*------------------- global variables for multithreading -------------------*/
#ifdef USE_THREADS
 pthread_t		*thread;
 pthread_mutex_t	instrumutex, readmutex, sortmutex;
 threads_gate_t		*pthread_startgate, *pthread_stopgate;
 fieldstruct		**pthread_fields;
 int			*pthread_fviewflag,
			pthread_endflag, pthread_nfield,
			pthread_findex, pthread_fviewindex;
#endif


/****** load_field ***********************************************************
PROTO   fieldstruct *load_field(char *filename)
PURPOSE Read catalog(s) and load field data.
INPUT   Character string that contains the file name.
OUTPUT  A pointer to the created field structure.
NOTES   Global preferences are used. The function is not reentrant because
	of static variables (prefs structure members are updated).
AUTHOR  E. Bertin (IAP)
VERSION 12/12/2006
*/
fieldstruct	*load_field(char *filename)
  {
   wcsstruct	*wcs;
   catstruct	*cat;
   tabstruct	*tab, *imatab;
   keystruct	*key;
   fieldstruct	*field;
   setstruct	**set;
   h_type	htype;
   t_type	ttype;
   char		str[MAXCHAR], label[72],
		*rfilename, *pstr, *astrombuf, *photombuf, *pspath;
   int		i, j, n, s, nsample, line;
   
/* A short, "relative" version of the filename */
  if (!(rfilename = strrchr(filename, '/')))
    rfilename = filename;
  else
    rfilename++;

  sprintf(str,"Examining Catalog %s", rfilename);
  NFPRINTF(OUTPUT, str);

/*-- Read input catalog */
  if (!(cat = read_cat(filename)))
    error(EXIT_FAILURE, "*Error*: No such catalog: ", filename);
  QCALLOC(field, fieldstruct, 1);
  QMALLOC(field->set, setstruct *, MAXSET);
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
  set = field->set;

/* Complete primary HDU first */
  field->headflag |= !read_aschead(prefs.ahead_global, 0, cat->tab);
  field->headflag |= !read_aschead(field->hfilename, 0, cat->tab);

/* Give a "colour" to the present field */
  field->cplot_colour = 15;
  fitsread(cat->tab->headbuf, prefs.cplot_colourkey, &field->cplot_colour,
		H_INT,T_LONG);
/* Put an astrometric label to the present field */
  field->astromlabel = 0;
/* Create a dummy FITS header to store all keyword values */
  QMALLOC(astrombuf, char, FBSIZE);
  memset(astrombuf, ' ', FBSIZE);
  strncpy(astrombuf, "END     ", 8);
  for (s=0; s<prefs.nastrinstru_key; s++)
    {
    fitsadd(astrombuf, prefs.astrinstru_key[s], "");
    if ((line=fitsfind(cat->tab->headbuf,prefs.astrinstru_key[s]))
	!= RETURN_ERROR)
      {
      fitspick(cat->tab->headbuf+line*80, str,(void *)label,
	       &htype,&ttype, str);
      fitswrite(astrombuf, prefs.astrinstru_key[s], label, htype, ttype);
      }
    }

/* Put a photometric label to the present field */
  field->photomlabel = 0;
/* Create dummy a FITS header to store all keyword values */
  QMALLOC(photombuf, char, FBSIZE);
  memset(photombuf, ' ', FBSIZE);
  strncpy(photombuf, "END     ", 8);
  for (s=0; s<prefs.nphotinstru_key; s++)
    {
    fitsadd(photombuf, prefs.photinstru_key[s], "");
    if ((line=fitsfind(cat->tab->headbuf,prefs.photinstru_key[s]))
	!= RETURN_ERROR)
      {
      fitspick(cat->tab->headbuf+line*80, str,(void *)label,
	       &htype,&ttype, str);
      fitswrite(photombuf, prefs.astrinstru_key[s], label, htype, ttype);
      }
    }
  n = 0;

/* Now scan other HDUs */
  for (i=cat->ntab; i--; tab=tab->nexttab)
    if ((!strcmp("LDAC_IMHEAD",tab->extname))
	&& (key=read_key(tab, "Field Header Card")))
      {
      set[n] = init_set();
/*---- Create a new table from scratch to hold the image header */
      imatab = new_tab("Image header");
      free(imatab->headbuf);
      imatab->headnblock = 1 + (key->nbytes-1)/FBSIZE;
      QCALLOC(imatab->headbuf, char, imatab->headnblock*FBSIZE);
      memcpy(imatab->headbuf, key->ptr, key->nbytes);
      imatab->cat = cat;
      readbasic_head(imatab);
      field->headflag |= !read_aschead(prefs.ahead_global, n, imatab);
      field->headflag |= !read_aschead(field->hfilename, n, imatab);
      if (!imatab->headbuf
	|| fitsread(imatab->headbuf, "OBJECT  ", field->ident,
		H_STRING,T_STRING)!= RETURN_OK)
        strcpy(field->ident, "no ident");
      set[n]->imatab = imatab;
      if (field->cplot_colour==15)
        fitsread(imatab->headbuf, prefs.cplot_colourkey, &field->cplot_colour,
		H_INT,T_LONG);
/*---- Try to read the astrometric label again */
      for (s=0; s<prefs.nastrinstru_key; s++)
        {
        fitsadd(astrombuf, prefs.astrinstru_key[s], "");
        if ((line=fitsfind(imatab->headbuf, prefs.astrinstru_key[s]))
	    != RETURN_ERROR)
          {
          fitspick(imatab->headbuf+line*80, str,(void *)label,
	       &htype,&ttype, str);
          fitswrite(astrombuf, prefs.astrinstru_key[s], label, htype, ttype);
          }
        }
/*---- Try to read the photometric label again */
      for (s=0; s<prefs.nphotinstru_key; s++)
        {
        fitsadd(photombuf, prefs.photinstru_key[s], "");
        if ((line=fitsfind(imatab->headbuf, prefs.photinstru_key[s]))
	    != RETURN_ERROR)
          {
          fitspick(imatab->headbuf+line*80, str,(void *)label,
	       &htype,&ttype, str);
          fitswrite(photombuf, prefs.photinstru_key[s], label, htype, ttype);
          }
        }
      n++;
      }

  field->nset = n;
  if (!n)
    error(EXIT_FAILURE,"*Error*: No SExtractor FITS-LDAC header found in ",
	rfilename);

  if (field->cplot_colour<0 || field->cplot_colour>15)
    warning("CHECKPLOT field colour out of range, defaulted to ", "15");

#ifdef USE_THREADS
  QPTHREAD_MUTEX_LOCK(&instrumutex);
#endif
/* Compare the dummy astrometric FITS header to the ones previously stored */
  for (j=0; j<prefs.nastrinstrustr; j++)
    if (!strncmp((const char *)prefs.astrinstrustr[j], astrombuf,
	80*prefs.nastrinstru_key))
      {
      field->astromlabel = j;
      break;
      }
  if (j>=prefs.nastrinstrustr)
    {
    QMEMCPY(astrombuf, prefs.astrinstrustr[prefs.nastrinstrustr], char,FBSIZE);
    field->astromlabel = prefs.nastrinstrustr++;
    }
  free(astrombuf);

/* Use the derived astrometric label index to associate the right */
/* mosaic and stability types to the present field */
  field->mosaic_type = prefs.mosaic_type[field->astromlabel]; 
  field->stability_type = prefs.stability_type[field->astromlabel]; 

/* Compare the dummy photometric FITS header to the ones previously stored */
  for (j=0; j<prefs.nphotinstrustr; j++)
    if (!strncmp((const char *)prefs.photinstrustr[j], photombuf,
	80*prefs.nphotinstru_key))
      {
      field->photomlabel = j;
      break;
      }
  if (j>=prefs.nphotinstrustr)
    {
    QMEMCPY(photombuf, prefs.photinstrustr[prefs.nphotinstrustr], char, FBSIZE);
    field->photomlabel = prefs.nphotinstrustr++;
    }
  free(photombuf);
#ifdef USE_THREADS
  QPTHREAD_MUTEX_UNLOCK(&instrumutex);
#endif

/* For every header the catalog contains */
  for (i=0; i<field->nset; i++)
    {
/*-- Manage WCS info */
    wcs = set[i]->wcs = read_wcs(set[i]->imatab);
    set[i]->lng = wcs->lng;
    set[i]->lat = wcs->lat;
/*-- Precess to 2000.0 if the equinox is different */
    if (fabs(wcs->equinox-2000.0)>0.001)
      {
      if (!i)
        {
        sprintf(str, "precessing EQUINOX %7.2f to %7.2f", wcs->equinox, 2000.0);
        NFPRINTF(OUTPUT, "");
        warning(str, "");
        }
      precess_wcs(wcs, wcs->equinox, 2000.0);
      }
/*-- Indicate what is the parent field */
    set[i]->field = field;
    }

/* Find the object table */
  sprintf(str,"Loading Catalog %s", rfilename);
  tab = cat->tab;
  nsample = 0;
  n = 0;
  for (i=cat->ntab; i--; tab=tab->nexttab)
    if (!strcmp("LDAC_OBJECTS", tab->extname)
	|| !strcmp("OBJECTS", tab->extname))
      {
      if (field->nset>1)
        sprintf(str, "%s [%d/%d]", rfilename, n+1, field->nset);
      else
        strcpy(str, rfilename);
      if (n>field->nset)
        {
        warning("Too many object catalogs in ", rfilename);
        break;
        }
      read_samples(set[n], tab, str);
      nsample += set[n]->nsample;
      free_tab(set[n]->imatab);
      set[n]->imatab = NULL;
      n++;
      }

  field->nsample = nsample;
  free_cat(&cat, 1);

  if (!n)
    {
    end_field(field);
    error(EXIT_FAILURE,"*Error*: No SExtractor FITS-LDAC catalog found in ",
	rfilename);
    }

  return field;
  }


/****** locate_field *********************************************************
PROTO   void locate_field(fieldstruct *field)
PURPOSE Compute basic field characteristics.
INPUT   Pointer to field structure.
OUTPUT  A pointer to the created field structure.
NOTES   Global preferences are used.
AUTHOR  E. Bertin (IAP)
VERSION 07/02/2005
*/
void	locate_field(fieldstruct *field)
  {
   setstruct		**pset, *set;
   wcsstruct		*wcs;
   double		*scale[NAXIS],*scalet[NAXIS],
			*wcsmean,
			cosalpha,sinalpha, sindelta, dist, maxradius, airmass;
   int			i, s, lat,lng, nset, naxis;

/* Some initializations */
  nset = field->nset;
  cosalpha = sinalpha = sindelta = 0.0;
  wcs = field->set[0]->wcs;
  naxis = field->naxis = wcs->naxis;
  wcsmean = field->meanwcspos;
  for (i=0; i<naxis; i++)
    {
    QMALLOC(scale[i], double, nset);
    scalet[i] = scale[i];
    wcsmean[i] = 0.0;
    }

/* Go through each set */
  pset = field->set;
  set = *(pset++);
  for (s=nset; s--; set=*(pset++))
    {
    wcs = set->wcs;
    lng = wcs->lng;
    lat = wcs->lat;
/*-- Locate set */
    locate_set(set);
    if (lat != lng)
      {
      cosalpha += cos(set->wcspos[lng]*DEG);
      sinalpha += sin(set->wcspos[lng]*DEG);
      sindelta += sin(set->wcspos[lat]*DEG);
      }
    for (i=0; i<naxis; i++)
      {
      if (lat==lng || (i!=lng && i!=lat))
        wcsmean[i] += set->wcspos[i];
      *(scalet[i]++) = set->wcsscale[i];
      }
    }


/* Now make the stats on each axis */
  lng = field->lng = field->set[0]->wcs->lng;
  lat = field->lat = field->set[0]->wcs->lat;
  field->epoch = field->set[0]->wcs->obsdate;
  for (i=0; i<naxis; i++)
    {
    if (lat!=lng && (i==lng))
      {
      wcsmean[i] = atan2(sinalpha/nset,cosalpha/nset)/DEG;
      wcsmean[i] = fmod(wcsmean[i]+360.0, 360.0);
      }
    else if (lat!=lng && (i==lat))
      wcsmean[i] = asin(sindelta/nset)/DEG;
    else
      wcsmean[i] /= nset;
    field->meanwcsscale[i] = dhmedian(scale[i], nset);
    }

/* Compute the radius of the field and mean airmass */
  airmass = maxradius = 0.0;
  pset = field->set;
  set=*(pset++);
  for (s=nset; s--; set=*(pset++))
    {
/*-- The distance is the distance to the center + the diagonal of the image */
    dist = wcs_dist(set->wcs, set->wcspos, field->meanwcspos) + set->radius;
    if (dist>maxradius)
      maxradius = dist;
    airmass += set->airmass;
    }

  field->maxradius = maxradius;
  field->airmass = airmass / nset;

/* Free memory */
  for (i=0; i<naxis; i++)
    free(scale[i]);

  return;
  }


/****** end_field ***********************************************************
PROTO   void end_field(fieldstruct *field)
PURPOSE Deallocate field data.
INPUT   Field pointer.
OUTPUT  -.
NOTES   -.
AUTHOR  E. Bertin (IAP)
VERSION	24/06/2004
*/
void	end_field(fieldstruct *field)
  {
   int	i;

  if (field->set)
    {
    for (i=0; i<field->nset; i++)
      if (field->set[i])
        end_set(field->set[i]);
    free(field->set);
    }
  free(field);

  return;
  }


/****** print_fieldinfo ******************************************************
PROTO	void print_fileinfo(fieldstruct *field)
PURPOSE	Print info about a field.
INPUT	Pointer to the field.
OUTPUT	-.
NOTES	-.
AUTHOR	E. Bertin (IAP)
VERSION	24/06/2004
 ***/
void	print_fieldinfo(fieldstruct *field)

  {
   tabstruct		*imatab;
   setstruct		*set;

  set = field->set[0];
  imatab = set->imatab;
  QPRINTF(OUTPUT, "%s:  \"%-19.19s\"  %s %3d set%s %7d detection%s\n",
        field->rfilename, field->ident,
        field->headflag? "EXTERN. HEADER" : "no ext. header",
	field->nset, field->nset>1 ? "s":"",
	field->nsample, field->nsample>1 ? "s":"");

  return;
  }


/****** dhmedian ******************************************************
PROTO	double   dhmedian(double *ra, int n)
PURPOSE	Compute the median of an array of doubles, using the Heapsort
	algorithm (based on Num.Rec algo.).
INPUT	Pointer to the array,
	Number of array elements.
OUTPUT	Median of the array.
NOTES	Warning: the order of input data is modified!.
AUTHOR	E. Bertin (IAP)
VERSION	22/07/2002
 ***/
double   dhmedian(double *ra, int n)

  {
   int		l, j, ir, i;
   double	rra;


  if (n<2)
    return *ra;
  ra--;
  for (l = ((ir=n)>>1)+1;;)
    {
    if (l>1)
      rra = ra[--l];
    else
      {
      rra = ra[ir];
      ra[ir] = ra[1];
      if (--ir == 1)
        {
        ra[1] = rra;
        return n&1? ra[n/2+1] : (ra[n/2]+ra[n/2+1])/2.0;
        }
      }
    for (j = (i=l)<<1; j <= ir;)
      {
      if (j < ir && ra[j] < ra[j+1])
        ++j;
      if (rra < ra[j])
        {
        ra[i] = ra[j];
        j += (i=j);
        }
      else
        j = ir + 1;
      }
    ra[i] = rra;
    }

/* (the 'return' is inside the loop!!) */
  }


#ifdef USE_THREADS

/****** pthread_load_field ***************************************************
PROTO   void *pthread_load_field(void *arg)
PURPOSE thread that takes care of reading catalogs.
INPUT   Pointer to the thread number.
OUTPUT  -.
NOTES   Relies on global variables.
AUTHOR  E. Bertin (IAP)
VERSION 27/09/2004
 ***/
void    *pthread_load_field(void *arg)
  {
   int findex, proc;

  findex = -1;
  proc = *((int *)arg);
  threads_gate_sync(pthread_startgate);
  while (!pthread_endflag)
    {
    QPTHREAD_MUTEX_LOCK(&readmutex);
    if (findex>-1)
/*---- Indicate that the field info is now suitable for viewing */
      pthread_fviewflag[findex] = 1;
    while (pthread_fviewindex<pthread_nfield
	&& pthread_fviewflag[pthread_fviewindex])
      print_fieldinfo(pthread_fields[pthread_fviewindex++]);
    if (pthread_findex<pthread_nfield)
      {
      findex = pthread_findex++;
      QPTHREAD_MUTEX_UNLOCK(&readmutex);
/*---- Load catalogs */
      pthread_fields[findex] = load_field(prefs.file_name[findex]);
/*---- Compute basic field astrometric features (center, field size,...) */
      locate_field(pthread_fields[findex]);
      }
    else
      {
      QPTHREAD_MUTEX_UNLOCK(&readmutex);
/*---- Wait for the input buffer to be updated */
      threads_gate_sync(pthread_stopgate);
/* ( Master thread process loads and saves new data here ) */
      threads_gate_sync(pthread_startgate);
      }
    }

  pthread_exit(NULL);

  return (void *)NULL;
  }

/****** pthread_load_fields ***************************************************
PROTO   void pthread_load_fields(fieldstruct **fields, int nfield)
PURPOSE Read catalogs in parallel using threads.
INPUT   Pointer to field structure pointers,
	number of fields.
OUTPUT  -.
NOTES   Relies on global variables.
AUTHOR  E. Bertin (IAP)
VERSION 27/10/2006
 ***/
void    pthread_load_fields(fieldstruct **fields, int nfield)
  {
   static pthread_attr_t	pthread_attr;
   int				*proc,
				p;

/* Number of active threads */
  nproc = prefs.nthreads;
  pthread_fields = fields;
  pthread_nfield = nfield;
  QCALLOC(pthread_fviewflag, int, nfield);
/* Set up multi-threading stuff */
  QMALLOC(proc, int, nproc);
  QMALLOC(thread, pthread_t, nproc);
  QPTHREAD_MUTEX_INIT(&readmutex, NULL);
  QPTHREAD_MUTEX_INIT(&instrumutex, NULL);
  QPTHREAD_MUTEX_INIT(&sortmutex, NULL);
  QPTHREAD_ATTR_INIT(&pthread_attr);
  QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
  pthread_startgate = threads_gate_init(nproc+1, NULL);
  pthread_stopgate = threads_gate_init(nproc+1, NULL);
/* Start the reading threads */
  for (p=0; p<nproc; p++)
    {
    proc[p] = p;
    QPTHREAD_CREATE(&thread[p], &pthread_attr, &pthread_load_field, &proc[p]);
    }
  QPTHREAD_MUTEX_LOCK(&readmutex);
  pthread_findex = pthread_fviewindex = 0;
  pthread_endflag = 0;
  QPTHREAD_MUTEX_UNLOCK(&readmutex);
/* Release threads!! */
  threads_gate_sync(pthread_startgate);
/* ( Slave threads process the current buffer data here ) */
  threads_gate_sync(pthread_stopgate);
  pthread_endflag = 1;
/* (Re-)activate existing threads... */
  threads_gate_sync(pthread_startgate);
/* ... and shutdown all threads */
  for (p=0; p<nproc; p++)
    QPTHREAD_JOIN(thread[p], NULL);
/* Clean up multi-threading stuff */
  threads_gate_end(pthread_startgate);
  threads_gate_end(pthread_stopgate);
  QPTHREAD_MUTEX_DESTROY(&readmutex);
  QPTHREAD_MUTEX_DESTROY(&instrumutex);
  QPTHREAD_ATTR_DESTROY(&pthread_attr);
  free(pthread_fviewflag);
  free(proc);
  free(thread);
  }


/****** pthread_end_fields ****************************************************
PROTO   void pthread_end_fields(fieldstruct **fields, int nfield)
PURPOSE Free structures and MUTEXes related to field parallel handling
INPUT   Pointer to field structure pointers,
	number of fields.
OUTPUT  -.
NOTES   Relies on global variables.
AUTHOR  E. Bertin (IAP)
VERSION 18/09/2006
 ***/
void    pthread_end_fields(fieldstruct **fields, int nfield)
  {
   int		f;

  QPTHREAD_MUTEX_DESTROY(&sortmutex);
  for (f=0; f<nfield; f++)
    end_field(fields[f]);
  free(fields);

  return;
  }

#endif

