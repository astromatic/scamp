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

#include "chealpix.h"
#include "chealpixstore.h"
#include "crossid2.h"
#include "field.h"

#define NNEIGHBORS 8

typedef enum {CROSS_UP, CROSS_DOWN} cross_direction;
struct field_id {
    int fieldindex;
    int havematch;
};

static void cross_sample(samplestruct*, HealPixel*, PixelStore*, double, cross_direction);
static int cross_time_closest_sample(samplestruct*,samplestruct*,
        double,PixelStore*, struct field_id*);
static int crossmatch(samplestruct*,samplestruct*, double, PixelStore*);

static inline double
dist(double *va, double *vb)
{
    double x = va[0] - vb[0];
    double y = va[1] - vb[1];
    double z = va[2] - vb[2];

    return x*x + y*y + z*z;
}

void
CrossId_crossSamples(
        PixelStore *pixstore,
        double  radius_arcsec)
{

    /* sort samples */
    PixelStore_sort(pixstore);

    /*
     * Convert radius in arcsec then in 3d euclidean distance
     * arcsec -> rad -> vector -> distance
     */
    double radius_radian = radius_arcsec / 3600 * TO_RAD;
    double va[3], vb[3];
    ang2vec(0, 0, va);
    ang2vec(radius_radian, 0, vb);
    double maxdist = dist(va, vb);

    /*
     * Iterate over all created Healpix pixels. They contains a least one
     * sample
     */

    int i;
    for (i=0; i<pixstore->npixels; i++) {
        HealPixel *pix = PixelStore_get(pixstore, pixstore->pixelids[i]);
        int j;
        for (j=0; j<pix->nsamples; j++) {
            cross_sample(pix->samples[j], pix, pixstore, maxdist, CROSS_UP);
        }
    }
}

/* TODO factorize */
static void
cross_sample(
        samplestruct *current_spl,
        HealPixel *pix,
        PixelStore *store,
        double maxdist,
        cross_direction direction)
{

    int field_index_pivot;
    struct field_id fi;

    /*
     * Eliminate previous fields, we only match with upper fields
     * TODO this can be optimized, by using this function once for a set
     * of samples of the same field.
     */
    int i;
    fi.fieldindex = -1;
    fi.havematch  = 0;
    switch (direction) {
        case CROSS_UP:
            field_index_pivot = PixelStore_getHigherFields(pix, current_spl);
            for (i=field_index_pivot; i<pix->nsamples; i++) 
            {
                if (cross_time_closest_sample(
                            current_spl, pix->samples[i],
                            maxdist, store, &fi) == 1)
                    break;
            }
            break;
        case CROSS_DOWN:
            field_index_pivot = PixelStore_getLowerFields(pix, current_spl);
            for (i=field_index_pivot; i>-1; i--) 
            {
                if (cross_time_closest_sample(
                            pix->samples[i], current_spl,
                            maxdist, store, &fi) == 1)
                    break;
            }

            break;
    }

    /*
     * Then iterate against neighbors pixels
     * TODO this can be optimized, by using this loop once for a set
     * of samples.
     */
    for (i=0; i<NNEIGHBORS; i++) {

        /*
         * Does the pixel exists? It may be a neighbor of current pixel,
         * but not be initialized because it does not contains
         * any samples. TODO opti: get pixel once.
         */
        HealPixel *test_pix = PixelStore_get(store, pix->neighbors[i]);
        if (test_pix == NULL)
            continue;

        /*
         * Eliminate previous to current fields, we only match with upper
         * fields (or the oposite).
         */
        int j;
        fi.fieldindex = -1;
        fi.havematch = 0;
        switch (direction) {
            case CROSS_UP:
                field_index_pivot = PixelStore_getHigherFields(test_pix, current_spl);
                for (j=field_index_pivot; j<test_pix->nsamples; j++) {
                    if (cross_time_closest_sample(current_spl, test_pix->samples[j],
                                maxdist, store, &fi) == 1)
                        break;
                }
                break;
            case CROSS_DOWN:
                field_index_pivot = PixelStore_getLowerFields(test_pix, current_spl);
                for (j=field_index_pivot; j>-1; j--) {
                    if (cross_time_closest_sample(test_pix->samples[j], current_spl,
                                maxdist, store, &fi) == 1)
                        break;
                }
                break;
        }
    }
}



static int
cross_time_closest_sample(
        samplestruct *current_spl,
        samplestruct *test_spl,
        double maxdist,
        PixelStore *store,
        struct field_id *current_field)
{

    if (current_field->fieldindex != test_spl->set->field->fieldindex) {
        /* We are changing of field */

        if (current_field->havematch) {
            /* if previous field have a match break */
            return 1;
        } 

        /* else change field and continue */
        current_field->fieldindex = test_spl->set->field->fieldindex;
    }

    current_field->havematch = crossmatch(current_spl, test_spl, maxdist, store);

    return 0;
}


static int
crossmatch(
        samplestruct *a,
        samplestruct *b,
        double maxdist,
        PixelStore *store)
{

    /*
     * "a" and "b" comme allready sorted by epoch, with "a" being older than
     * b. Calling this function in the opposite order is a bug.
     */
    /*
    if (PixelStore_compare(a, b) > 0)
        exit(1);
    */


    /*
     * First if distance exceed the maximum distance, stop here.
     */
    double distance = dist(a->vector, b->vector);
    if (distance > maxdist)
        return 0;

    /*
     * "a" allready have a nextsamp. So we have to determinate what action to
     * take, even before comparing distances...
     */
    if (a->nextsamp) {
        /*
         * Compare the new candiate a->nextsamp and a->nextsamp. 
         */
        int a_cmp = PixelStore_compare(b, a->nextsamp);

        /* 
         * If b is farther in time than a->nextsamp, stop here.
         *
         * If "a" have a valid ->nextsamp closer in time than "b" 
         * we will not compare "a"<->"b" distance and keep the 
         * actual nextsamp for "a".
         */
        if (a_cmp > 0)
            return 0;

        /*
         * If "a->nextsamp" belong to the same field layer than "b". So we
         * compare distances, and go further if "b" is a better match.
         */
        if (a_cmp == 0) {
            double a_next_dist = dist(a->vector, a->nextsamp->vector);
            if (a_next_dist <= distance)
                return 0;
        } 
        /* Else whatever the distance is, "a" will want to link to "b" */
    }

    /*
     * The same as above, but the opposite (search in descending order)
     */
    if (b->prevsamp) {

        int b_cmp = PixelStore_compare(a, b->prevsamp);

        if (b_cmp < 0)
            return 0;

        if (b_cmp == 0) {
            double b_prev_dist = dist(b->vector, b->prevsamp->vector);
            if (b_prev_dist <= distance)
                return 0;
        }
        /* Else whatever the distance is, "b" will want to link to "a" */
    }


    /* 
     * Here, we want to link a->b. But we also want previous linked samples
     * to search another link.
     *
     * The order of these next liens IS important.
     */

    /* If an old link exists, take it and reset their link information */
    samplestruct *relink_up    = b->prevsamp;
    samplestruct *relink_down  = a->nextsamp;

    if (relink_down)
        relink_down->prevsamp = NULL;
    if (relink_up)
        relink_up->nextsamp = NULL;

    a->nextsamp = b;
    b->prevsamp = a;




    /* 
     * We have now reached a stable state, with a and b linked, and possibly
     * old links pointing to NULL samples.
     */




    /*We want to recurse for each of them
     * and try to find a new link. There are many cases that need this 
     * recursion.
     *
     * One of the trickiest:
     * - let's take 4 samples "w" "x" "y" "z".
     * - "w" match "x", succeed
     * - "y" match "x", fail because "w" is closer to "x" than "y"
     * - "z" match "w", succeed
     * - "x" would match "y", but the comparison is allready done, and "y" 
     * and "x" will stay unlinked.
     *
     * So recursion on unlink is required.
     */
    if (relink_up) {
        HealPixel *pix = PixelStore_getPixelFromSample(store, relink_up); 
        cross_sample(relink_up, pix, store, maxdist, CROSS_UP);
    }

    if (relink_down) {
        HealPixel *pix = PixelStore_getPixelFromSample(store, relink_down);
        cross_sample(relink_down, pix, store, maxdist, CROSS_DOWN);
    }

    /* alleluia */
    return 1;
}


