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

static int ds = 0;

typedef enum {TIME_CLOSEST, RAW_CLOSEST} closest_algo;
typedef enum {CROSS_UP, CROSS_DOWN} cross_direction;
struct field_id {
    int fieldindex;
    int havematch;
};

static closest_algo closest = TIME_CLOSEST;
static void cross_sample(struct sample*, HealPixel*, PixelStore*, double, cross_direction);
static int cross_time_closest_sample(struct sample*,struct sample*,
        double,PixelStore*, struct field_id*);
static int crossmatch(struct sample*,struct sample*, double, PixelStore*);
static void relink_sample_up(struct sample*, PixelStore*, double);
static void relink_sample_down(struct sample*, PixelStore*, double);
static double dist(double*,double*);


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

static void
cross_sample(
        struct sample *current_spl,
        HealPixel *pix,
        PixelStore *store,
        double maxdist,
        cross_direction direction)
{

    if (direction == CROSS_DOWN) {
        fprintf(stderr, "not implemented\n");
        exit(1);
    }

    int field_index_pivot;
    struct field_id fi;

    /*
     * Eliminate previous fields, we only match with upper fields
     * TODO this can be optimized, by using this function once for a set
     * of samples of the same field.
     */
    switch (direction) {
        case CROSS_UP:
            field_index_pivot = PixelStore_getHigherFields(pix, current_spl);
            break;
        case CROSS_DOWN:
            field_index_pivot = PixelStore_getLowerFields(pix, current_spl);
            break;
    }

    int i;
    switch (closest) {
        case TIME_CLOSEST:
            fi.fieldindex = -1;
            fi.havematch = 0;
            if (direction == CROSS_UP)
                for (i=field_index_pivot; i<pix->nsamples; i++) 
                {
                    if (cross_time_closest_sample(current_spl, pix->samples[i],
                                maxdist, store, &fi) == 1)
                    {
                        break;
                    }
                }
            else
                for (i=0; i>field_index_pivot; i++) 
                {
                    if (cross_time_closest_sample(pix->samples[i], current_spl,
                                maxdist, store, &fi) == 1)
                    {
                        break;
                    }
                }

            break;
        case RAW_CLOSEST:
            if (direction == CROSS_UP)
                for (i=field_index_pivot; i<pix->nsamples; i++)
                    crossmatch(current_spl, pix->samples[i], maxdist, store);
            else
                for (i=0; i<field_index_pivot; i++)
                    crossmatch(pix->samples[i], current_spl, maxdist, store);
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
        switch (direction) {
            case CROSS_UP:
                field_index_pivot = PixelStore_getHigherFields(test_pix, current_spl);
                break;
            case CROSS_DOWN:
                field_index_pivot = PixelStore_getLowerFields(test_pix, current_spl);
                break;
        }

        int j;
        switch (closest) {
            case TIME_CLOSEST:
                fi.fieldindex = -1;
                fi.havematch = 0;
                if (direction == CROSS_UP)
                    for (j=field_index_pivot; j<test_pix->nsamples; j++) 
                    {
                        if (cross_time_closest_sample(current_spl, test_pix->samples[j],
                                    maxdist, store, &fi) == 1)
                        {
                            break;
                        }
                    }
                else
                    for (i=0; i>field_index_pivot; i++) 
                    {
                        if (cross_time_closest_sample(test_pix->samples[j], current_spl,
                                    maxdist, store, &fi) == 1)
                        {
                            break;
                        }
                    }

                break;
            case RAW_CLOSEST:
                if (direction == CROSS_UP)
                    for (j=field_index_pivot; j<test_pix->nsamples; j++)
                        crossmatch(current_spl, test_pix->samples[j], maxdist, store);
                else
                    for (j=0; j<field_index_pivot; j++)
                        crossmatch(test_pix->samples[j], current_spl, maxdist, store);
                break;
        }
    }
}



static int
cross_time_closest_sample(
        struct sample *current_spl,
        struct sample *test_spl,
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
        struct sample *a,
        struct sample *b,
        double maxdist,
        PixelStore *store)
{

    double distance = dist(a->vector, b->vector);
    struct sample *relink_up    = NULL;
    struct sample *relink_down  = NULL;

    if (distance > maxdist)
        return 0;

    if (a->nextsamp) {
        double a_next_dist = dist(a->vector, a->nextsamp->vector);
        if (a_next_dist < distance)
            return 0;
        else if (a_next_dist == distance)
            if (PixelStore_compare(b, a->nextsamp) >= 0)
                return 0;
    }

    if (b->prevsamp) {
        double b_prev_dist = dist(b->vector, b->prevsamp->vector);
        if (b_prev_dist < distance) {
            return 0;
        } else if (b_prev_dist == distance)
            if (PixelStore_compare(a, b->prevsamp) <= 0)
                return 0;
    }

    relink_up   = b->prevsamp;
    relink_down = a->nextsamp;

    a->nextsamp = b;
    b->prevsamp = a;

    if (relink_up)
        relink_sample_up(relink_up, store, maxdist);
    if (relink_down)
        relink_sample_down(relink_down, store, maxdist);

    return 1;
}

static void
relink_sample_down(struct sample *sample, PixelStore *store, double maxdist){
    sample->prevsamp = NULL;
    HealPixel *pix = PixelStore_getPixelFromSample(store, sample);
    cross_sample(sample, pix, store, maxdist, CROSS_DOWN);
}

static void
relink_sample_up(struct sample *sample, PixelStore *store, double maxdist)
{
    sample->nextsamp = NULL;
    HealPixel *pix = PixelStore_getPixelFromSample(store, sample);
    cross_sample(sample, pix, store, maxdist, CROSS_UP);
}

static double
dist(double *va, double *vb)
{
    double x = va[0] - vb[0];
    double y = va[1] - vb[1];
    double z = va[2] - vb[2];

    return x*x + y*y + z*z;
}
