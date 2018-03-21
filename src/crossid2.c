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
#include "field.h"

#define NNEIGHBORS 8

static long cross_pixel(HealPixel*,PixelStore*,double);
static int crossmatch(struct sample*,struct sample*, double, PixelStore*);
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

struct field_id {
    int epoch;
    int fieldindex;
    int havematch;
};
/* cross all samples from one pixel to himself, and all neighbors pixel 
   samples. */
static long
cross_pixel(
        HealPixel *pix, 
        PixelStore *store, 
        double radius)
{
    long nbmatches = 0;
    int j, k, l;
    int field_index_start, field_index_stop;

    struct sample *current_spl, *test_spl;
    double distance, old_distance;
    bool rematch_test_sample;
    struct field_id current_field;

    /*
     * Mark other pixels as done so they will not cross identify with me later.
     */
    HealPixel *test_pix;
    for (j=0; j<NNEIGHBORS; j++) {
        test_pix = pix->pneighbors[j];
        if (test_pix) {
            for (k=0; k<NNEIGHBORS; k++) {
                if (test_pix->neighbors[l] == pix->id) {
                    test_pix->tneighbors[l] = true;
                    break;
                }
            }
        }
    }

    /*
     * Then iterate over our pixel samples
     */
    for (j=0; j<pix->nsamples; j++) {
        current_spl = pix->samples[j];

        /*
         * First cross match with samples of the pixel between them.
         */

        /*
         * Eliminate previous fields, we only match with upper fields
         * TODO this can be optimized, by using this function once for a set
         * of samples of the same field.
         */
        PixelStore_getHigherFields(pix, current_spl, 
                &field_index_start, &field_index_stop);

        /*
         * All this to iterate over all samples of a field, and stop at
         * the next field if there is a match. A match in a field must not
         * stop iteration for the same field.
         */


        /* get the next field id */
        current_field.epoch      = 0;
        current_field.fieldindex = 0;
        current_field.havematch  = 0;

        for(k=field_index_start; k<field_index_stop; k++) {
            test_spl = pix->samples[k];

            if (    current_field.epoch      != test_spl->epoch || 
                    current_field.fieldindex != test_spl->set->field->fieldindex) 
            {
                /* we are changing of field */
                if (current_field.havematch) { 
                    /* if previous field have a match break */
                    break;
                } else {
                    /* else change field and continue */
                    current_field.epoch = test_spl->epoch;
                    current_field.fieldindex = test_spl->set->field->fieldindex;
                }
            } 

            current_field.havematch = crossmatch(
                    current_spl, test_spl, radius, store);
        }


        /*
         * Then iterate against neighbors pixels
         * TODO this can be optimized, by using this loop once for a set
         * of samples.
         */
        for (k=0; k<NNEIGHBORS; k++) {

            /* Allready crossed by an neighbor ? */
            if (pix->tneighbors[k] == true)
                continue;

            /*
             * Does the pixel exists? It may be a neighbor of current pixel,
             * but not be initialized because it does not contains
             * any samples.
             */
            test_pix = pix->pneighbors[k];
            if (test_pix == NULL)
                continue;

            /*
             * Eliminate previous to current fields, we only match with upper 
             * fields.
             */
            PixelStore_getHigherFields(test_pix, current_spl, 
                    &field_index_start, &field_index_stop);

            /*
             * All this to iterate over all samples of a field, and stop at
             * the next field if there is a match. A match in a field must not
             * stop iteration for the same field 
             */
            
            /* initialize current field to a false value */
            current_field.epoch = 0;
            current_field.fieldindex = 0;
            current_field.havematch = 0;

            for (l=field_index_start; l<field_index_stop; l++) {
                test_spl = test_pix->samples[l];

                if (    current_field.epoch      != test_spl->epoch || 
                        current_field.fieldindex != test_spl->set->field->fieldindex) 
                {
                    /* we are changing of field */
                    if (current_field.havematch) { 
                        /* if previous field have a match break */
                        break;
                    } else {
                        /* else change field and continue */
                        current_field.epoch = test_spl->epoch;
                        current_field.fieldindex = test_spl->set->field->fieldindex;
                    }
                } 

                current_field.havematch = crossmatch(
                        current_spl, test_spl, radius, store);
            }
        }
    }

    return nbmatches;

}


static void
relink_sample(struct sample *sample, PixelStore *store) 
{
    fprintf(stderr, "\n\nRELINK SAMPLE\n\n");
    sample->nextsamp = NULL;
    HealPixel *pix = PixelStore_getPixelFromSample(store, sample);
}


static int
crossmatch(
        struct sample *a, 
        struct sample *b, 
        double radius,
        PixelStore *store)
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
        relink_sample(b->prevsamp, store);

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
