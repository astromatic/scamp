/*
 * makeit.c
 *
 * Main loop.
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 * This file part of: SCAMP
 *
 * Copyright: (C) 2002-2016 IAP/CNRS/UPMC
 *
 * License: GNU General Public License
 *
 * SCAMP is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * SCAMP is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
 *
 * Last modified: 20/10/2016
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_THREADS
#include <pthread.h>
#endif
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "define.h"
#include "globals.h"
#include "astrefcat.h"
#include "astrsolve.h"
#include "astrstats.h"
#include "catout.h"
#include "colour.h"
#include "cplot.h"
#include "dgeomap.h"
#include "fft.h"
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "header.h"
#include "match.h"
#include "merge.h"
#include "mosaic.h"
#include "photsolve.h"
#include "prefs.h"
#include "proper.h"
#ifdef USE_THREADS
#include "threads.h"
#endif
#include "wcs/wcs.h"
#include "xml.h"
#include "chealpixstore.h"
#include "crossid2.h"

time_t thetime, thetime2;
static PixelStore* new_pixstore(int, int, fieldstruct**, fieldstruct**);

/********************************** makeit ***********************************/
void makeit(void)
{
    static char filename[MAXCHAR], extension[MAXCHAR], str[MAXCHAR];
    fgroupstruct **fgroups;
    fieldstruct **fields, **reffields;
    setstruct *set;
    struct tm *tm;
    double alpha,delta;
    char *pstr,
         sign;
    int i,f,g, nfield, ngroup, nsample, nclip, hh,mm,dd,dm;

    /* Install error logging */
    error_installfunc(write_error);

    /* Processing start date and time */
    thetime = time(NULL);
    tm = localtime(&thetime);
    sprintf(prefs.sdate_start,"%04d-%02d-%02d",
            tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_start,"%02d:%02d:%02d",
            tm->tm_hour, tm->tm_min, tm->tm_sec);

    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT,
            "----- %s %s started on %s at %s with %d thread%s\n\n",
            BANNER,
            MYVERSION,
            prefs.sdate_start,
            prefs.stime_start,
            prefs.nthreads,
            prefs.nthreads>1? "s":"");

    nfield = prefs.nfile;

    /* End here if no filename has been provided */
    if (!nfield)
    {
        /*-- Processing end date and time */
        thetime2 = time(NULL);
        tm = localtime(&thetime2);
        sprintf(prefs.sdate_end,"%04d-%02d-%02d",
                tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
        sprintf(prefs.stime_end,"%02d:%02d:%02d",
                tm->tm_hour, tm->tm_min, tm->tm_sec);
        prefs.time_diff = difftime(thetime2, thetime);

        /*-- Write XML */
        if (prefs.xml_flag)
        {
            init_xml(NULL, 0, NULL, 0);
            write_xml(prefs.xml_name);
            end_xml();
        }
        return;
    }

    QMALLOC(fields, fieldstruct *, nfield);

    /*--------------------------- Read the Catalogs--------------------------- */
    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT, "----- %d inputs:\n", nfield);

#ifdef USE_THREADS
    pthread_load_fields(fields, nfield);
#else
    for (f=0; f<nfield; f++)
    {
        /*-- Load catalogs */
        fields[f] = load_field(prefs.file_name[f], f, prefs.ahead_name[f]);
        NFPRINTF(OUTPUT, "");

        /*-- Compute basic field astrometric features (center, field size,...) */
        locate_field(fields[f]);
        print_fieldinfo(fields[f]);
    }
#endif

    nsample = 0;
    for (f=0; f<nfield; f++)
        nsample += fields[f]->nsample;
    prefs.ndets = nsample;

    QPRINTF(OUTPUT, "\n----- %d detections loaded\n", nsample);

    /* Group fields on the sky */
    fgroups = group_fields(fields, nfield, &ngroup);
    NFPRINTF(OUTPUT, "");
    print_instruinfo();
    print_fgroupinfo(fgroups, ngroup);

    adjust_mosaic(fields, nfield);

#ifdef HAVE_PLPLOT
    /* Plot fields on the sky */
    cplot_allsky(fgroups, ngroup);
#endif

    /* One reference catalog per group */
    QCALLOC(reffields, fieldstruct *, ngroup);

    /* Read the reference catalogs */
    if (prefs.astrefcat != ASTREFCAT_NONE)
    {
        NFPRINTF(OUTPUT, "");
        QPRINTF(OUTPUT, "----- Reference catalogs:\n\n");
        for (g=0; g<ngroup; g++)
        {
            reffields[g] = get_astreffield(prefs.astrefcat,
                    fgroups[g]->meanwcspos, fgroups[g]->lng, fgroups[g]->lat,
                    fgroups[g]->naxis, fgroups[g]->maxradius+prefs.radius_maxerr);
            if (reffields[g])
            {
                NFPRINTF(OUTPUT, "");
                QPRINTF(OUTPUT, " Group %2d: %8d standard%s found in %s (%s band)\n",
                        g+1, reffields[g]->nsample, reffields[g]->nsample>1?"s":"",
                        astrefcats[prefs.astrefcat].name,
                        astrefcats[prefs.astrefcat].bandname);
                /*------ Save reference catalog (on request) */
                if (prefs.astrefcat != ASTREFCAT_FILE && prefs.outrefcat_flag)
                {
                    alpha = reffields[g]->meanwcspos[fgroups[g]->lng];
                    hh = (int)(alpha/15.0);
                    mm = (int)(60.0*(alpha/15.0 - hh));
                    delta = reffields[g]->meanwcspos[fgroups[g]->lat];
                    sign = delta<0.0?'-':'+';
                    delta = fabs(delta);
                    dd = (int)delta;
                    dm = (int)(60.0*(delta - dd));
                    sprintf(str, "%s/%s_%02d%02d%c%02d%02d_r%-.0f.cat",
                            prefs.outref_path,
                            astrefcats[prefs.astrefcat].name,
                            hh,mm,sign,dd,dm, reffields[g]->maxradius*DEG/ARCMIN);
                    save_astreffield(str, reffields[g]);
                }
            }
            else
            {
                sprintf(str, "No source found in reference catalog(s) for group %d; "
                        "wrong sky zone?", g+1);
                error(EXIT_FAILURE, "*Error*: ", str);
            }
        }
    }

    /* Find where fields are in the sky */
    if (prefs.match_flag && prefs.astrefcat != ASTREFCAT_NONE)
    {
        fft_init(prefs.nthreads/nfield);
        QPRINTF(OUTPUT, "\n----- Astrometric matching:\n\n");
#ifdef USE_THREADS
        pthread_match_fields(fgroups, reffields, ngroup);
#else
        for (g=0; g<ngroup; g++)
        {
            NFPRINTF(OUTPUT, "");
            QPRINTF(OUTPUT, " Group %2d: %8d standard%s in %s (band %s)\n",
                    g+1, reffields[g]->nsample, reffields[g]->nsample>1?"s":"",
                    astrefcats[prefs.astrefcat].name,
                    astrefcats[prefs.astrefcat].bandname);
            QIPRINTF(OUTPUT, " instruments pos.angle scale "
                    "cont. shift cont.");
            if (reffields[g])
                for (f=0; f<fgroups[g]->nfield; f++)
                {
                    match_field(fgroups[g]->field[f], reffields[g]);
                    print_matchinfo(fgroups[g]->field[f]);
                }
        }
#endif
        fft_end();
    }

    QPRINTF(OUTPUT, "\n");



    for (g=0; g<ngroup; g++)
    {
        /*-- Reproject all fields from a group to a common projection */
        reproj_fgroup(fgroups[g], reffields[g], 0);
        /*-- Perform cross-identifications across catalogs */
        sprintf(str, "Making preliminary cross-identifications in group %d", g+1);
        NFPRINTF(OUTPUT, str);
    }

    //debug_crossid(ngroup, fgroups);
    PixelStore *ps;
    ps = new_pixstore(nfield, ngroup, reffields, fields) ;
    CrossId_crossSamples(ps, prefs.crossid_radius);
    PixelStore_free(ps);

    fprintf(stderr, "end cross 1\n\n");
    if (prefs.solvastrom_flag)
    {
        /*-- Compute global astrometric solution: 1st iteration */
        astrsolve_fgroups(fgroups, ngroup);

        NFPRINTF(OUTPUT, "");
        QPRINTF(OUTPUT, " \n----- Astrometric clipping:\n\n");
        for (g=0; g<ngroup; g++)
        {
            /*---- Reproject all fields from a group to a common projection (update) */
            reproj_fgroup(fgroups[g], reffields[g], 0);
            /*---- Perform cross-identifications across catalogs */
            sprintf(str, "Making cross-identifications in group %d", g+1);
            NFPRINTF(OUTPUT, str);
        }

        ps = new_pixstore(nfield, ngroup, reffields, fields) ;
        CrossId_crossSamples(ps, prefs.crossid_radius);
        PixelStore_free(ps);

        for (g=0; g<ngroup; g++) {
            sprintf(str, "Computing astrometric stats for group %d", g+1);
            NFPRINTF(OUTPUT, str);
            astrstats_fgroup(fgroups[g], reffields[g], prefs.sn_thresh[1]);
            sprintf(str, "Astrometric clipping in group %d", g+1);
            NFPRINTF(OUTPUT, str);
            nclip = astrclip_fgroup(fgroups[g], reffields[g], prefs.astrclip_nsig);
            NFPRINTF(OUTPUT, "");
            QPRINTF(OUTPUT, " Group %2d: %d/%d detections removed\n",
                    g+1, nclip, fgroups[g]->nintmatch);
        }

        /*-- Compute global astrometric solution: 2nd iteration */
        astrsolve_fgroups(fgroups, ngroup);
    }

    /* Display internal astrometric stats */
    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT, " \n----- Astrometric stats (internal) :\n\n");
    QIPRINTF(OUTPUT,
            " All detections | "
            " High S/N ");
    QIPRINTF(OUTPUT,
            " dAXIS1 dAXIS2 chi2 ndets | "
            "dAXIS1 dAXIS2 chi2 ndets");
    for (g=0; g<ngroup; g++)
    {
        /*-- Reproject all fields from a group to a common projection (update) */
        reproj_fgroup(fgroups[g], reffields[g], 0);
        /*-- Perform cross-identifications across catalogs */
        sprintf(str, "Making cross-identifications in group %d", g+1);
        NFPRINTF(OUTPUT, str);
    }
    ps = new_pixstore(nfield, ngroup, reffields, fields) ;
    CrossId_crossSamples(ps, prefs.crossid_radius);
    PixelStore_free(ps);
    for (g=0; g<ngroup; g++)
    {
        astrstats_fgroup(fgroups[g], reffields[g], prefs.sn_thresh[1]);
        nclip = astrclip_fgroup(fgroups[g], reffields[g], prefs.astrclip_nsig);
        astrstats_fgroup(fgroups[g], reffields[g], prefs.sn_thresh[1]);
        if (fgroups[g]->nintmatch>0)
        {
            QPRINTF(OUTPUT, 
                    "Group %2d: %6.3g\" %6.3g\" %6.2g %7d %6.3g\" %6.3g\" %6.2g %7d\n",
                    g+1,
                    fgroups[g]->sig_interr[0]*DEG/ARCSEC,
                    fgroups[g]->sig_interr[1]*DEG/ARCSEC,
                    fgroups[g]->chi2_int, fgroups[g]->nintmatch,
                    fgroups[g]->sig_interr_hsn[0]*DEG/ARCSEC,
                    fgroups[g]->sig_interr_hsn[1]*DEG/ARCSEC,
                    fgroups[g]->chi2_int_hsn, fgroups[g]->nintmatch_hsn);
        }
    }

    /* Display external astrometric stats */
    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT, " \n----- Astrometric stats (external):\n\n");
    QIPRINTF(OUTPUT,
            " All detections | "
            " High S/N ");
    QIPRINTF(OUTPUT,
            " dAXIS1 dAXIS2 chi2 nstars | "
            "dAXIS1 dAXIS2 chi2 nstars");
    for (g=0; g<ngroup; g++)
    {
        QPRINTF(OUTPUT, 
                "Group %2d: %6.3g\" %6.3g\" %6.2g %7d %6.3g\" %6.3g\" %6.2g %7d\n",
                g+1,
                fgroups[g]->sig_referr[0]*DEG/ARCSEC,
                fgroups[g]->sig_referr[1]*DEG/ARCSEC,
                fgroups[g]->chi2_ref, fgroups[g]->nrefmatch,
                fgroups[g]->sig_referr_hsn[0]*DEG/ARCSEC,
                fgroups[g]->sig_referr_hsn[1]*DEG/ARCSEC,
                fgroups[g]->chi2_ref_hsn, fgroups[g]->nrefmatch_hsn);
    }

    if (prefs.solvphotom_flag)
    {
        /*-- Compute global photometric solution: 1st iteration */
        photsolve_fgroups(fgroups, ngroup);

        NFPRINTF(OUTPUT, "");
        QPRINTF(OUTPUT, " \n----- Photometric clipping:\n\n");
        for (g=0; g<ngroup; g++)
        {
            compmags_fgroup(fgroups[g]);
            for (i=0; i<prefs.nphotinstrustr; i++)
            {
                /*------ Compute photometric stats */
                sprintf(str, "Computing photometric stats for group %d / P%d",
                        g+1, i+1);
                NFPRINTF(OUTPUT, str);
                photstats_fgroup(fgroups[g], i, prefs.sn_thresh[1]);
                sprintf(str, "Photometric clipping in group %d / P%d", g+1, i+1);
                if (fgroups[g]->nintmagmatch[i]>0)
                {
                    nclip = photclip_fgroup(fgroups[g], i, prefs.photclip_nsig);
                    NFPRINTF(OUTPUT, "");
                    QPRINTF(OUTPUT, " Group %2d / P%-2d : %d/%d detections removed\n",
                            g+1, i+1, nclip, fgroups[g]->nintmagmatch[i]);
                }
            }
        }

        /*-- Compute global photometric solution: 2nd iteration */
        photsolve_fgroups(fgroups, ngroup);
    }

    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT, " \n----- Photometric stats (internal):\n\n");
    QIPRINTF(OUTPUT, " All detections | "
            " High S/N ");
    QIPRINTF(OUTPUT, " Instru mag RMS chi2 ndets | "
            "mag RMS chi2 ndets");
    for (g=0; g<ngroup; g++)
    {
        compmags_fgroup(fgroups[g]);
        for (i=0; i<prefs.nphotinstrustr; i++)
        {
            /*---- Compute photometric stats */
            photstats_fgroup(fgroups[g], i, prefs.sn_thresh[1]);
            nclip = photclip_fgroup(fgroups[g], i, prefs.photclip_nsig);
            photstats_fgroup(fgroups[g], i, prefs.sn_thresh[1]);
            if (fgroups[g]->nintmagmatch[i]>0)
            {
                QPRINTF(OUTPUT,
                        "Group %2d: P%-2d %7.3g %7.2g %7d %7.3g %7.2g %7d\n",
                        g+1, i+1,
                        fgroups[g]->sig_intmagerr[i], fgroups[g]->chi2_intmag[i],
                        fgroups[g]->nintmagmatch[i],
                        fgroups[g]->sig_intmagerr_hsn[i],
                        fgroups[g]->chi2_intmag_hsn[i],
                        fgroups[g]->nintmagmatch_hsn[i]);
            }
        }
    }

    NFPRINTF(OUTPUT, "");
    QPRINTF(OUTPUT, " \n----- Photometric stats (external):\n\n");
    QIPRINTF(OUTPUT, " All detections | "
            " High S/N ");
    QIPRINTF(OUTPUT, " Instru mag RMS chi2 nstars | "
            "mag RMS chi2 nstars");
    for (g=0; g<ngroup; g++)
        for (i=0; i<prefs.nphotinstrustr; i++)
        {
            if (fgroups[g]->nrefmagmatch[i]>0)
            {
                QPRINTF(OUTPUT,
                        "Group %2d: P%-2d %7.3g %7.2g %7d %7.3g %7.2g %7d\n",
                        g+1, i+1,
                        fgroups[g]->sig_refmagerr[i], fgroups[g]->chi2_refmag[i],
                        fgroups[g]->nrefmagmatch[i],
                        fgroups[g]->sig_refmagerr_hsn[i],
                        fgroups[g]->chi2_refmag_hsn[i],
                        fgroups[g]->nrefmagmatch_hsn[i]);
            }
        }

    QPRINTF(OUTPUT, "\n");

    NFPRINTF(OUTPUT, "Merging detections...");
    for (g=0; g<ngroup; g++)
        merge_fgroup(fgroups[g], reffields[g]);

    /* Compute colour indices */
    NFPRINTF(OUTPUT, "Computing global color indices");
    colour_fgroup(fgroups, ngroup);

    /* Correct colour shifts */
    if (prefs.colourshiftcorr_flag)
        /*-- Correct positions for colour-dependent effects */
        for (g=0; g<ngroup; g++)
        {
            sprintf(str, "Computing colour shifts in group %d", g+1);
            NFPRINTF(OUTPUT, str);
            astrcolshift_fgroup(fgroups[g], reffields[g]);
        }

    if (prefs.propmotioncorr_flag && prefs.solvastrom_flag)
    {
        /*-- Re-do Cross-ID to recover possibly fast moving objects */
        NFPRINTF(OUTPUT, "Pairing detections...");

        ps = new_pixstore(nfield, ngroup, reffields, fields) ;
        CrossId_crossSamples(ps, prefs.crossid_radius);
        PixelStore_free(ps);

        NFPRINTF(OUTPUT, "Merging detections...");
        for (g=0; g<ngroup; g++)
            merge_fgroup(fgroups[g], reffields[g]);
        for (g=0; g<ngroup; g++)
        {
            sprintf(str, "Computing proper motions in group %d", g+1);
            NFPRINTF(OUTPUT, str);
            astrprop_fgroup(fgroups[g]);
        }
        for (g=0; g<ngroup; g++)
            /*---- Reproject to a common projection while correcting for proper motions */
            reproj_fgroup(fgroups[g], reffields[g], 1);
        NFPRINTF(OUTPUT, "Pairing detections...");

        ps = new_pixstore(nfield, ngroup, reffields, fields) ;
        CrossId_crossSamples(ps, prefs.crossid_radius);
        PixelStore_free(ps);

        NFPRINTF(OUTPUT, "Merging detections...");
        for (g=0; g<ngroup; g++)
            merge_fgroup(fgroups[g], reffields[g]);
        /*-- Compute global astrometric solution: 3rd iteration */
        astrsolve_fgroups(fgroups, ngroup);
    }

    /* Compute final proper motions and parallaxes */
#ifdef HAVE_PLPLOT
    if (prefs.propmotion_flag || prefs.parallax_flag
            || cplot_check(CPLOT_REFPROP)!=RETURN_ERROR
            || cplot_check(CPLOT_ADPROP2D)!=RETURN_ERROR)
#else
        if (prefs.propmotion_flag || prefs.parallax_flag)
#endif
        {
            /*-- Re-do Cross-ID to recover possibly fast moving objects */
            NFPRINTF(OUTPUT, "Pairing detections...");

            ps = new_pixstore(nfield, ngroup, reffields, fields) ;
            CrossId_crossSamples(ps, prefs.crossid_radius);
            PixelStore_free(ps);
            NFPRINTF(OUTPUT, "Merging detections...");
            for (g=0; g<ngroup; g++)
                merge_fgroup(fgroups[g], reffields[g]);
            for (g=0; g<ngroup; g++)
            {
                sprintf(str, "Computing proper motions in group %d", g+1);
                NFPRINTF(OUTPUT, str);
                astrprop_fgroup(fgroups[g]);
            }
            for (g=0; g<ngroup; g++)
                /*---- Reproject all fields from a group to a common projection (update) */
                reproj_fgroup(fgroups[g], reffields[g], prefs.propmotioncorr_flag);
        }

    /* Save headers */
    NFPRINTF(OUTPUT, "Saving image headers...");
    for (f=0; f<nfield; f++) {
        /*---- Check if a header filename is provided */
        if (prefs.head_name[f] && *(prefs.head_name[f]))
            strcpy(filename, prefs.head_name[f]);
        else { 
            /*---- Create a file name with a "header" extension */
            strcpy(filename, fields[f]->filename);
            if (!(pstr = strrchr(filename, '.')))
                pstr = filename+strlen(filename);
            sprintf(pstr, "%s", prefs.head_suffix);
        }
        write_aschead(filename, fields[f]);
    }

#ifdef HAVE_PLPLOT

    /* Plot field and source positions */
    NFPRINTF(OUTPUT, "Generating group plots...");
    for (g=0; g<ngroup; g++)
        cplot_fgroup(fgroups[g], reffields[g]);
    for (g=0; g<ngroup; g++)
        cplot_astrepoch3d(fgroups[g]);
    /* Plot photometric relations */
    cplot_photom(fgroups, ngroup, reffields);
    for (i=0; i<prefs.nastrinstrustr; i++)
        cplot_shear(fgroups, ngroup, i);

    /* Plot astrometric errors in alpha and delta */
    NFPRINTF(OUTPUT, "Generating astrometric plots...");
    for (g=0; g<ngroup; g++)
        cplot_aderrhisto2d(fgroups[g], prefs.sn_thresh[1]);
    for (g=0; g<ngroup; g++)
        cplot_aderrhisto1d(fgroups[g], prefs.sn_thresh[1]);
    for (g=0; g<ngroup; g++)
        cplot_referrhisto2d(fgroups[g], reffields[g], prefs.sn_thresh[1]);
    for (g=0; g<ngroup; g++)
        cplot_referrhisto1d(fgroups[g], reffields[g], prefs.sn_thresh[1]);
    for (g=0; g<ngroup; g++)
        cplot_chi2(fgroups[g]);

    /* Plot sub-pixel astrometric error dependency */
    for (i=0; i<prefs.nastrinstrustr; i++)
        cplot_pixerrhisto1d(fgroups, ngroup, i, prefs.sn_thresh[1]);
    for (i=0; i<prefs.nastrinstrustr; i++)
        cplot_xpixerrhisto2d(fgroups, ngroup, i);
    for (i=0; i<prefs.nastrinstrustr; i++)
        cplot_ypixerrhisto2d(fgroups, ngroup, i);
    for (i=0; i<prefs.nastrinstrustr; i++)
        cplot_subpixerrhisto1d(fgroups, ngroup, i, prefs.sn_thresh[1]);

    /* Plot astrometric distortions */
    for (i=0; i<prefs.nastrinstrustr; i++)
        for (f=0; f<nfield; f++)
            if (fields[f]->astromlabel == i)
            {
                cplot_distort(fields[f]);
                break;
            }
    for (i=0; i<prefs.nastrinstrustr; i++)
        cplot_astintsysmap(fgroups, ngroup, i, prefs.sn_thresh[1]);
    for (i=0; i<prefs.nastrinstrustr; i++)
        cplot_astrefsysmap(fgroups, ngroup, i, prefs.sn_thresh[1]);

    NFPRINTF(OUTPUT, "Generating photometric plots...");
    for (g=0; g<ngroup; g++)
        cplot_photzp(fgroups[g]);
    for (g=0; g<ngroup; g++)
        cplot_photzp3d(fgroups[g]);
    for (g=0; g<ngroup; g++)
        cplot_photerrhisto(fgroups[g], reffields[g], prefs.sn_thresh[1]);
    for (g=0; g<ngroup; g++)
        cplot_photerrhistomag(fgroups[g], reffields[g], prefs.sn_thresh[1]);

    NFPRINTF(OUTPUT, "Generating color shift plots...");
    for (g=0; g<ngroup; g++)
        cplot_astrcolshift1d(fgroups[g], prefs.sn_thresh[1]);
    NFPRINTF(OUTPUT, "Generating proper-motion plots...");
    for (g=0; g<ngroup; g++)
        cplot_astrefprop(fgroups[g], reffields[g], prefs.sn_thresh[1]);
    for (g=0; g<ngroup; g++)
        cplot_adprophisto2d(fgroups[g], prefs.sn_thresh[1]);
#endif

    init_xml(fields, nfield, fgroups, ngroup);

    /* Processing end date and time */
    thetime2 = time(NULL);
    tm = localtime(&thetime2);
    sprintf(prefs.sdate_end,"%04d-%02d-%02d",
            tm->tm_year+1900, tm->tm_mon+1, tm->tm_mday);
    sprintf(prefs.stime_end,"%02d:%02d:%02d",
            tm->tm_hour, tm->tm_min, tm->tm_sec);
    prefs.time_diff = difftime(thetime2, thetime);

    /* Save merged catalogs */
    if (prefs.mergedcat_type != CAT_NONE)
    {
        for (g=0; g<ngroup; g++)
        {
            /*---- Write one catalog per field group */
            sprintf(str, "Saving merged catalog for group %d", g+1);
            NFPRINTF(OUTPUT, str);
            strcpy(filename, prefs.mergedcat_name);
            if (!(pstr = strrchr(filename, '.')))
            {
                pstr = filename+strlen(filename);
                extension[0] = (char)'\0';
            }
            else
                strcpy(extension, pstr);
            sprintf(pstr, "_%d%s", g+1, extension);
            writemergedcat_fgroup(filename, fgroups[g]);
        }
    }

    if (prefs.dgeomap_flag) {
        /*-- Compute and write differential geometry vector maps */
        NFPRINTF(OUTPUT, "Generating differential geometry vector maps ...");
        for (i=0; i<prefs.nastrinstrustr; i++) {
            strcpy(filename, prefs.dgeomap_name);
            if (!(pstr = strrchr(filename, '.'))) {
                pstr = filename+strlen(filename);
                extension[0] = (char)'\0';
            } else
                strcpy(extension, pstr);
            sprintf(pstr, "_%0d%s", i+1, extension);
            dgeomap_instru(fields, nfield, i, filename);
        }
    }

    /* Save full catalogs */
    if (prefs.fullcat_type != CAT_NONE)
    {
        for (g=0; g<ngroup; g++)
        {
            /*---- Write one catalog per field group */
            sprintf(str, "Saving full catalog for group %d", g+1);
            NFPRINTF(OUTPUT, str);
            strcpy(filename, prefs.fullcat_name);
            if (!(pstr = strrchr(filename, '.')))
            {
                pstr = filename+strlen(filename);
                extension[0] = (char)'\0';
            }
            else
                strcpy(extension, pstr);
            sprintf(pstr, "_%d%s", g+1, extension);
            writefullcat_fgroup(filename, fgroups[g]);
        }
    }

    /* Write XML */
    if (prefs.xml_flag)
        write_xml(prefs.xml_name);

    end_xml();

    /* Clean-up stuff */
    NFPRINTF(OUTPUT, "Cleaning up...");
    for (g=0; g<ngroup; g++)
    {
        end_fgroup(fgroups[g]);
        if (reffields[g])
            end_field(reffields[g]);
    }
    free(fgroups);
    free(reffields);

#ifdef USE_THREADS
    pthread_end_fields(fields, nfield);
#else
    for (f=0; f<nfield; f++)
        end_field(fields[f]);
    free(fields);
#endif

    return;
}


/****** write_error ********************************************************
  PROTO void write_error(char *msg1, char *msg2)
  PURPOSE Manage files in case of a catched error
  INPUT a character string,
  another character string
  OUTPUT RETURN_OK if everything went fine, RETURN_ERROR otherwise.
  NOTES -.
  AUTHOR E. Bertin (IAP)
  VERSION 02/10/2006
 ***/
void write_error(char *msg1, char *msg2)
{
    char error[MAXCHAR];

    sprintf(error, "%s%s", msg1,msg2);
    if (prefs.xml_flag)
        write_xmlerror(prefs.xml_name, error);
    end_xml();

    return;
}


static PixelStore*
new_pixstore(
        int nfield, 
        int ngroup, 
        fieldstruct **reffields, 
        fieldstruct **fields) 
{
    /* Initialize healpix values and stores */
    int64_t nsides = pow(2, prefs.healpix_resolution);
    PixelStore *ps = PixelStore_new(nsides);
    struct set *set;

    int i, f, g;
    for (i=0; i<nfield; i++) {
        for (f=0; f<fields[i]->nset; f++) {
            set = fields[i]->set[f];
            for (g=0; g < set->nsample;g++) {
                struct sample *s = &set->sample[g];
                s->id = g;
                PixelStore_add(ps, &set->sample[g]);
            }
        }
    }

    if (prefs.astrefcat != ASTREFCAT_NONE) {
        for (i=0; i<ngroup; i++) {
            for (f=0; f<reffields[i]->nset; f++) {
                set = reffields[i]->set[f];
                for (g=0; g < set->nsample; g++) {
                    struct sample *s = &set->sample[g];
                    s->id = g;
                    PixelStore_add(ps, &set->sample[g]);
                }
            }
        }
    }

    PixelStore_sort(ps);
    return ps;
}

