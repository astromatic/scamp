/*
 *    field.c
 *
 * Manage catalogues (individual exposures).
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 * This file part of: SCAMP
 *
 * Copyright:  (C) 2002-2023 IAP/CNRS/SorbonneU/CEA/UParisSaclay
 *
 * License:  GNU General Public License
 *
 * SCAMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * SCAMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
 *
 * Last modified:  05/12/2023
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
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
pthread_t		*field_thread;
pthread_mutex_t		field_instrumutex, field_readmutex;
extern pthread_mutex_t	sample_sortmutex;
fieldstruct		**pthread_field_fields;
int			*pthread_field_fviewflag,
			pthread_field_endflag, pthread_field_nfield,
			pthread_field_findex, pthread_field_fviewindex;
#endif


/****** load_field ***********************************************************
PROTO   fieldstruct *load_field(char *filename, int fieldindex, char *hfilename)
PURPOSE Read catalog(s) and load field data.
INPUT   Character string that contains the file name.
OUTPUT  A pointer to the created field structure.
NOTES   Global preferences are used. The function is not reentrant because
        of static variables (prefs structure members are updated),
        FITS header filename (null=none).
AUTHOR  E. Bertin (CEA/AIM/UParisSaclay)
VERSION 09/04/2025
 */
fieldstruct *load_field(char *filename, int fieldindex, char *hfilename)
{
    wcsstruct *wcs;
    catstruct *cat;
    tabstruct *tab, *imatab;
    keystruct *key;
    fieldstruct *field;
    setstruct **set;
    h_type htype;
    t_type ttype;
    char str[MAXCHAR], label[72], keystr[16],
         *rfilename, *pstr, *pspath;
    int  d, i, j, n, s, nsample, line;

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
    field->fieldindex = fieldindex;
    field->isrefcat = 0;

    strcpy (field->filename, filename);
    field->rfilename = rfilename;

    /*-- Check if a "aheader" filename is provided */
    if (hfilename && *hfilename)
        strcpy(field->hfilename, hfilename);
    else { 
        /*-- Create a file name with a "aheader" extension */
        strcpy(field->hfilename, filename);
        if (!(pstr = strrchr(field->hfilename, '.')))
            pstr = field->hfilename+strlen(field->hfilename);
        sprintf(pstr, "%s", prefs.ahead_suffix);
    }

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
    /* Complete primary HDU first */
    field->headflag |= !read_aschead(prefs.ahead_global, 0, cat->tab);
    field->headflag |= !read_aschead(field->hfilename, 0, cat->tab);

    /* Give a "colour" to the present field */
    field->cplot_colour = 15;
    fitsread(cat->tab->headbuf, prefs.cplot_colourkey, &field->cplot_colour,
            H_INT, T_LONG, 0);

    /* Put an astrometric label to the present field */
    field->astromlabel = 0;
    /* Create a dummy FITS header to store all keyword values */
    QMALLOC(field->astrombuf, char, FBSIZE);
    memset(field->astrombuf, ' ', FBSIZE);
    strncpy(field->astrombuf, "END     ", 8);
    for (s=0; s<prefs.nastrinstru_key; s++) {
        fitsadd(field->astrombuf, prefs.astrinstru_key[s], "");
        if ((line=fitsfind(cat->tab->headbuf,prefs.astrinstru_key[s]))
                != RETURN_ERROR) {
            fitspick(cat->tab->headbuf+line*80, str,(void *)label,
                    &htype,&ttype, str);
            fitswrite(field->astrombuf, prefs.astrinstru_key[s], label,
                htype, ttype);
        }
    }

    /* Put a photometric label to the present field */
    field->photomlabel = 0;
    /* Create dummy a FITS header to store all keyword values */
    QMALLOC(field->photombuf, char, FBSIZE);
    memset(field->photombuf, ' ', FBSIZE);
    strncpy(field->photombuf, "END     ", 8);
    for (s=0; s<prefs.nphotinstru_key; s++) {
        fitsadd(field->photombuf, prefs.photinstru_key[s], "");
        if ((line=fitsfind(cat->tab->headbuf,prefs.photinstru_key[s]))
                != RETURN_ERROR) {
            fitspick(cat->tab->headbuf+line*80, str,(void *)label,
                    &htype,&ttype, str);
            fitswrite(field->photombuf, prefs.photinstru_key[s], label,
                htype, ttype);
        }
    }
    n = 0;

    /* Now scan other HDUs */
    tab = cat->tab;
    set = field->set;
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
                        H_STRING, T_STRING, MAXCHAR) != RETURN_OK)
                strcpy(field->ident, "no ident");
            set[n]->imatab = imatab;
            if (field->cplot_colour==15)
                fitsread(imatab->headbuf, prefs.cplot_colourkey,
                    &field->cplot_colour, H_INT, T_LONG, 0);
            /*---- Try to read the astrometric label again */
            for (s=0; s<prefs.nastrinstru_key; s++) {
                fitsadd(field->astrombuf, prefs.astrinstru_key[s], "");
                if ((line=fitsfind(imatab->headbuf, prefs.astrinstru_key[s]))
                        != RETURN_ERROR) {
                    fitspick(imatab->headbuf+line*80, str,(void *)label,
                        &htype,&ttype, str);
                    fitswrite(field->astrombuf, prefs.astrinstru_key[s], label,
                        htype, ttype);
                }
            }
            /*---- Try to read the photometric label again */
            for (s=0; s<prefs.nphotinstru_key; s++) {
                fitsadd(field->photombuf, prefs.photinstru_key[s], "");
                if ((line=fitsfind(imatab->headbuf, prefs.photinstru_key[s]))
                        != RETURN_ERROR)
                {
                    fitspick(imatab->headbuf+line*80, str,(void *)label,
                        &htype,&ttype, str);
                    fitswrite(field->photombuf, prefs.photinstru_key[s], label,
                        htype, ttype);
                }
            }
            n++;
        }

    field->nset = n;
    if (!n)
        error(EXIT_FAILURE,"*Error*: No SExtractor FITS-LDAC header found in ",
                rfilename);

    /* Save some memory */
    QREALLOC(field->set, setstruct *, field->nset);
    set = field->set;

    if (field->cplot_colour<0 || field->cplot_colour>15)
        warning("CHECKPLOT field colour out of range, defaulted to ", "15");

    field->projection_type = prefs.projection_type[field->fieldindex];

    /* For every header the catalog contains */
    for (i=0; i<field->nset; i++)
    {
        /*-- Set CTYPEs */
        if (field->projection_type != PROJECTION_SAME)
            for (d=0; d<NAXIS; d++)
            {
                sprintf(keystr, "CTYPE%1d  ", d+1);
                if (fitsread(set[i]->imatab->headbuf, keystr, str,
                	H_STRING,T_STRING, 16)
                        == RETURN_OK
                        && (pstr=strrchr(str, '-')))
                {
                    sprintf(pstr+1, field->projection_type==PROJECTION_TPV? "TPV":"TAN");
                    fitswrite(set[i]->imatab->headbuf, keystr, str, H_STRING,T_STRING);
                }
            }
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
        /*-- Force coordinate system to be ICRS */
        wcs->radecsys = RDSYS_ICRS;
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
            set[n]->setindex = n;
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


/****** label_field *********************************************************
  PROTO   void label_fields(fieldstruct *field)
  PURPOSE Set astrometric and photometric labels and related settings 
  INPUT   Pointer to field structure.
  OUTPUT  -.
  NOTES   Relies on global variables.
  AUTHOR  E. Bertin (CEA/AIM/UParisSaclay)
  VERSION 09/04/2025
 ***/
void    label_field(fieldstruct *field) {
   char str[MAXCHAR];
   int  j;

    // Compare the dummy astrometric FITS header to the ones previously stored
    for (j=0; j<prefs.nastrinstrustr; j++)
       if (!strncmp((const char *)prefs.astrinstrustr[j], field->astrombuf,
           80*prefs.nastrinstru_key) && field->nset == prefs.nastrinstruext[j]) {
            field->astromlabel = j;
            break;
        }
    if (j>=prefs.nastrinstrustr) {
        QMEMCPY(field->astrombuf, prefs.astrinstrustr[prefs.nastrinstrustr],
            char, FBSIZE);
        prefs.nastrinstruext[prefs.nastrinstrustr] = field->nset;
        field->astromlabel = prefs.nastrinstrustr++;
        if (prefs.nastrinstrustr > MAXASTRINSTRU)
        {
            sprintf(str, "%d", prefs.nastrinstrustr);
            error(EXIT_FAILURE,"*Error*: Too many astrometric instruments: ", str);
        }
    }
    // Compare the dummy photometric FITS header to the ones previously stored
    for (j=0; j<prefs.nphotinstrustr; j++)
        if (!strncmp((const char *)prefs.photinstrustr[j], field->photombuf,
                    80*prefs.nphotinstru_key)) {
            field->photomlabel = j;
            break;
        }
    if (j>=prefs.nphotinstrustr) {
        QMEMCPY(field->photombuf, prefs.photinstrustr[prefs.nphotinstrustr],
            char, FBSIZE);
        field->photomlabel = prefs.nphotinstrustr++;
        if (prefs.nphotinstrustr > MAXPHOTINSTRU) {
            sprintf(str, "%d", prefs.nphotinstrustr);
            error(EXIT_FAILURE,"*Error*: Too many photometric instruments: ",
                str);
        }
    }

    /* Use the derived astrometric label index to associate the right */
    /* mosaic and stability types to the present field */
    field->mosaic_type = prefs.mosaic_type[field->astromlabel]; 
    field->stability_type = prefs.stability_type[field->astromlabel]; 
}


/****** locate_field *********************************************************
  PROTO   void locate_field(fieldstruct *field)
  PURPOSE Compute basic field characteristics.
  INPUT   Pointer to field structure.
  OUTPUT  A pointer to the created field structure.
  NOTES   Global preferences are used.
  AUTHOR  E. Bertin (IAP)
  VERSION 11/11/2017
 */
void locate_field(fieldstruct *field)
{
    setstruct  **pset, *set;
    wcsstruct  *wcs;
    double  *scale[NAXIS],*scalet[NAXIS],
    *wcsmean,
    cosalpha,sinalpha, sindelta, dist, maxradius,
    airmass,airmassmin,airmassmax, epoch,epochmin,epochmax,
    expotime,expotimemin,expotimemax;
    int   i, s, lat,lng, nset, naxis, nairmass,nepoch,nexpotime;

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
    for (s=nset; s--;)
    {
        set = *(pset++);
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

    /* Compute the field radius, as well as airmass,epoch,exposure time stats */
    airmass = epoch = expotime = maxradius = 0.0;
    airmassmax = expotimemax = epochmax
        = -(airmassmin = expotimemin = epochmin = BIG);
    nairmass = nexpotime = nepoch = 0;
    pset = field->set;
    for (s=nset; s--;)
    {
        set=*(pset++);
        /*-- The distance is the distance to the center + the diagonal of the image */
        dist = wcs_dist(set->wcs, set->wcspos, field->meanwcspos) + set->radius;
        if (dist>maxradius)
            maxradius = dist;
        if (set->airmass != 0.0)
        {
            airmass += set->airmass;
            if (set->airmass < airmassmin)
                airmassmin = set->airmass;
            if (set->airmass > airmassmax)
                airmassmax = set->airmass;
            nairmass++;
        }
        if (set->epochmin!=0.0)
        {
            epoch += set->epoch;
            if (set->epochmin < epochmin)
                epochmin = set->epochmin;
            if (set->epochmax > epochmax)
                epochmax = set->epochmax;
            nepoch++;
        }
        if (set->expotime != 0.0)
        {
            expotime += set->expotime;
            if (set->expotime < expotimemin)
                expotimemin = set->expotime;
            if (set->expotime > expotimemax)
                expotimemax = set->expotime;
            nexpotime++;
        }
    }

    field->maxradius = maxradius;

    /* Update sample uncertainties */
    for (s=0; s<nset; s++)
        update_samples(field->set[s], maxradius);

    if ((nairmass))
    {
        field->airmass = airmass / nairmass;
        field->airmassmin = airmassmin;
        field->airmassmax = airmassmax;
    }
    else
        field->airmass = field->airmassmin = field->airmassmax = 0.0;
    if ((nepoch))
    {
        field->epoch = epoch / nepoch;
        field->epochmin = epochmin;
        field->epochmax = epochmax;
    }
    else
        field->epoch = field->epochmin = field->epochmax = 0.0;
    if ((nexpotime))
    {
        field->expotime = expotime / nexpotime;
        field->expotimemin = expotimemin;
        field->expotimemax = expotimemax;
    }
    else
        field->expotime = field->expotimemin = field->expotimemax = 0.0;

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
AUTHOR  E. Bertin (CEA/AIM/UParisSaclay)
VERSION 09/05/2025
*/
void end_field(fieldstruct *field) {
    int i;

    if (field->set) {
        for (i=0; i<field->nset; i++)
            if (field->set[i])
                end_set(field->set[i]);
        free(field->set);
    }
    free(field->astrombuf);
    free(field->photombuf);
    free(field);

    return;
}


/****** print_fieldinfo ******************************************************
  PROTO void print_fileinfo(fieldstruct *field)
  PURPOSE Print info about a field.
  INPUT Pointer to the field.
  OUTPUT -.
  NOTES -.
  AUTHOR E. Bertin (IAP)
  VERSION 24/06/2004
 ***/
void print_fieldinfo(fieldstruct *field)

{
    tabstruct  *imatab;
    setstruct  *set;

    set = field->set[0];
    imatab = set->imatab;
    QPRINTF(OUTPUT, "%s:  \"%-19.19s\"  %s %3d set%s %7d detection%s\n",
            field->rfilename, field->ident,
            field->headflag? "EXTERN. HEADER" : "no ext. header",
            field->nset, field->nset>1 ? "s":"",
            field->nsample, field->nsample>1 ? "s":"");

    return;
}


/****** get_field_meanepoch *************************************************
  PROTO	double get_fields_meanepoch(fieldstruct **fields, int nfield)
  PURPOSE	Return the average epoch of all input fields.
  INPUT Pointer to field structure pointers,
  	Number of fields.
  OUTPUT	Average epoch of input fields.
  NOTES -.
  AUTHOR E. Bertin (DAp/CEA/UParisSaclay)
  VERSION 05/12/2023
 ***/
double	get_fields_meanepoch(fieldstruct **fields, int nfield)

{
    double	epoch = 0.0;
    int	f;

    for (f=0; f<nfield; f++) {
        epoch += fields[f]->epoch;
    }
	
    return nfield ? epoch / (double)nfield : 0.0;
}


/****** dhmedian ************************************************************
  PROTO double   dhmedian(double *ra, int n)
  PURPOSE Compute the median of an array of doubles, using the Heapsort
  algorithm (based on Num.Rec algo.).
  INPUT Pointer to the array,
  Number of array elements.
  OUTPUT Median of the array.
  NOTES Warning: the order of input data is modified!.
  AUTHOR E. Bertin (IAP)
  VERSION 22/07/2002
 ***/
double   dhmedian(double *ra, int n)

{
    int  l, j, ir, i;
    double rra;


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
  VERSION 12/08/2020
 ***/
void    *pthread_load_field(void *arg)
{
   int findex, proc;

    findex = -1;
    proc = *((int *)arg);
    threads_gate_sync(pthread_startgate);
    while (!pthread_field_endflag) {
        QPTHREAD_MUTEX_LOCK(&field_readmutex);
        if (findex>-1)
            // Indicate that the field is now suitable for order-dependent
            // operations such as labeling or info display.
            pthread_field_fviewflag[findex] = 1;
        while (pthread_field_fviewindex<pthread_field_nfield
                && pthread_field_fviewflag[pthread_field_fviewindex]) {
            label_field(pthread_field_fields[pthread_field_fviewindex]);
            print_fieldinfo(pthread_field_fields[pthread_field_fviewindex++]);
        }
        if (pthread_field_findex<pthread_field_nfield) {
            findex = pthread_field_findex++;
            QPTHREAD_MUTEX_UNLOCK(&field_readmutex);
            // Load catalogs
            pthread_field_fields[findex] = load_field(prefs.file_name[findex], findex,
                    prefs.ahead_name[findex]);
            // Compute basic field astrometric features (center, field size,...)
            locate_field(pthread_field_fields[findex]);
        } else {
            QPTHREAD_MUTEX_UNLOCK(&field_readmutex);
            // Wait for the input buffer to be updated
            threads_gate_sync(pthread_stopgate);
            // ( Master thread process loads and saves new data here )
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
  VERSION 12/08/2020
 ***/
void    pthread_load_fields(fieldstruct **fields, int nfield)
{
    static pthread_attr_t pthread_attr;
    int    *proc,
           p;

    /* Number of active threads (must be limited on manycore systems) */
    nproc = prefs.nthreads;
    if (nproc>MAXNTHREADS_LOAD)
        nproc = MAXNTHREADS_LOAD;
    pthread_field_fields = fields;
    pthread_field_nfield = nfield;
    QCALLOC(pthread_field_fviewflag, int, nfield);
    /* Set up multi-threading stuff */
    QMALLOC(proc, int, nproc);
    QMALLOC(field_thread, pthread_t, nproc);
    QPTHREAD_MUTEX_INIT(&field_readmutex, NULL);
    QPTHREAD_MUTEX_INIT(&field_instrumutex, NULL);
    QPTHREAD_MUTEX_INIT(&sample_sortmutex, NULL);
    QPTHREAD_ATTR_INIT(&pthread_attr);
    QPTHREAD_ATTR_SETDETACHSTATE(&pthread_attr, PTHREAD_CREATE_JOINABLE);
    pthread_startgate = threads_gate_init(nproc+1, NULL);
    pthread_stopgate = threads_gate_init(nproc+1, NULL);
    /* Start the reading threads */
    for (p=0; p<nproc; p++)
    {
        proc[p] = p;
        QPTHREAD_CREATE(&field_thread[p], &pthread_attr, &pthread_load_field,
        		&proc[p]);
    }
    QPTHREAD_MUTEX_LOCK(&field_readmutex);
    pthread_field_findex = pthread_field_fviewindex = 0;
    pthread_field_endflag = 0;
    QPTHREAD_MUTEX_UNLOCK(&field_readmutex);
    /* Release threads!! */
    threads_gate_sync(pthread_startgate);
    /* ( Slave threads process the current buffer data here ) */
    threads_gate_sync(pthread_stopgate);
    pthread_field_endflag = 1;
    /* (Re-)activate existing threads... */
    threads_gate_sync(pthread_startgate);
    /* ... and shutdown all threads */
    for (p=0; p<nproc; p++)
        QPTHREAD_JOIN(field_thread[p], NULL);
    /* Clean up multi-threading stuff */
    threads_gate_end(pthread_startgate);
    threads_gate_end(pthread_stopgate);
    QPTHREAD_MUTEX_DESTROY(&field_readmutex);
    QPTHREAD_MUTEX_DESTROY(&field_instrumutex);
    QPTHREAD_ATTR_DESTROY(&pthread_attr);
    free(pthread_field_fviewflag);
    free(proc);
    free(field_thread);
}


/****** pthread_end_fields ****************************************************
  PROTO   void pthread_end_fields(fieldstruct **fields, int nfield)
  PURPOSE Free structures and MUTEXes related to field parallel handling
  INPUT   Pointer to field structure pointers,
  number of fields.
  OUTPUT  -.
  NOTES   Relies on global variables.
  AUTHOR  E. Bertin (IAP)
  VERSION 12/08/2020
 ***/
void    pthread_end_fields(fieldstruct **fields, int nfield)
{
    int  f;

    QPTHREAD_MUTEX_DESTROY(&sample_sortmutex);
    for (f=0; f<nfield; f++)
        end_field(fields[f]);
    free(fields);

    return;
}

#endif

