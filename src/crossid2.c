/*
 * Efficient catalogs cross matching functions.
 *
 * Copyright (C) 2017 University of Bordeaux. All right reserved.
 * Written by Emmanuel Bertin
 * Written by Sebastien Serre
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "crossid2.h"
#include "chealpix.h"
#include "chealpixstore.h"

#define NNEIGHBORS 8

static long cross_pixel(HealPixel*,PixelStore*,double);
static int crossmatch(struct sample*,struct sample*, double radius);
static double dist(double*,double*);


long
CrossId_crossSamples(
        PixelStore *pixstore,
        double  radius_arcsec)
{
    int i;
    PixelStore_updateSamplePos(pixstore);

    /* arcsec to radiant */
    double radius = radius_arcsec * TO_RAD;

    /* radiant to vector */
    double va[3], vb[3];
    ang2vec(0,0,va);
    ang2vec(HALFPI - radius,0,vb);

    /* vector to euclidean distance */
    radius = dist(va, vb);

    long nmatches = 0;
    for (i=0; i<pixstore->npixels; i++) {
        HealPixel *pix = PixelStore_get(pixstore, pixstore->pixelids[i]);
        nmatches += cross_pixel(pix, pixstore, radius);
    }

    fprintf(stderr,
            "Crossmatch end: %li matches for all pixels!\n", nmatches);

    return nmatches;

}


/* cross all samples from one pixel to himself, and all neighbors pixel 
   samples. */
static long
cross_pixel(
        HealPixel *pix, 
        PixelStore *store, 
        double radius)
{

    long nbmatches = 0;

    /*
     * Iterate over HealPixel structure which old sample structures
     * belonging to him.
     */
    int j, k, l;
    int field_index_start, field_index_stop;

    struct sample *current_spl;
    double distance, old_distance;
    bool rematch_test_sample;

    for (j=0; j<pix->nsamples; j++) {
        current_spl = pix->samples[j];

        /*
         * Eliminate previous fields, we only match with upper fields
         */
        PixelStore_getHigherFields(pix, current_spl, 
                &field_index_start, &field_index_stop);

        /*
         * First cross match with samples of the pixel between them
         */

        /* 
         * TODO break after the first field that contains a match, but test
         * all samples for this field to get the best. Then stop. 
         */
        for(k=field_index_start; k<field_index_stop; k++)
            if (crossmatch(current_spl, pix->samples[k], radius))
                break;

        /*
         * Then iterate against neighbors pixels
         */
        HealPixel *test_pixel;
        for (k=0; k<NNEIGHBORS; k++) {

            /* Allready crossed by an neighbor ? */
            if (pix->tneighbors[k] == true)
                continue;

            /*
             * Does the pixel exists? It may be a neighbor of current pixel,
             * but not be initialized because it does not contains
             * any samples.
             */
            test_pixel = pix->pneighbors[k];
            if (test_pixel == NULL)
                continue;

            /*
             * Eliminate previous to current fields, we only match with upper 
             * fields.
             */
            PixelStore_getHigherFields(test_pixel, current_spl, 
                    &field_index_start, &field_index_stop);

            for (l=field_index_start; l<field_index_stop; l++)
                if(crossmatch(current_spl, test_pixel->samples[l], radius))
                    break;

        }
    }

    return nbmatches;

}


static int
crossmatch(
        struct sample *a, 
        struct sample *b, 
        double radius)
{

    double distance = dist(a->vector, b->vector);

    if (distance > radius) 
        return 0;

    if (a->nextsamp)
        if (dist(a->vector, a->nextsamp->vector) < distance) 
            return 0;

    if (b->prevsamp)
        if (dist(b->vector, b->prevsamp->vector) < distance) 
            return 0;

    if (a->nextsamp)
        a->nextsamp->prevsamp = NULL;

    if (b->prevsamp)
        b->prevsamp->nextsamp = NULL;

    a->nextsamp = b;
    b->prevsamp = a;

    return 1;
}


static double
dist(double *va, double *vb)
{
    double x = va[0] - vb[0];
    double y = va[1] - vb[1];
    double z = va[2] - vb[2];

    return x*x + y*y + z*z;
}
