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
#ifdef USE_THREADS
#include <pthread.h>
#endif

#include "crossid2.h"
#include "chealpix.h"
#include "chealpixstore.h"
#include "field.h"
#define NNEIGHBORS 8

typedef enum {TIME_CLOSEST, RAW_CLOSEST} closest_algo;
struct field_id {
    int fieldindex;
    int havematch;
};
#ifdef USE_THREADS
struct thread_args {
    PixelStore  *store;
    int64_t     *pixelindex;
    int         npixs;
    double      radius;
};
static pthread_mutex_t GLOBAL_MUTEX = PTHREAD_MUTEX_INITIALIZER;
#endif /* USE_THREADS */

static closest_algo closest = TIME_CLOSEST;
static void* pthread_cross_pixel(void*);
static void cross_pixel(HealPixel*,PixelStore*,double);
static void cross_sample(struct sample*, HealPixel*, PixelStore*, bool, double);
static int cross_time_closest_sample(
        struct sample*,struct sample*,double,PixelStore*, struct field_id*);
static int crossmatch(struct sample*,struct sample*, double, PixelStore*);
static void relink_sample(struct sample*, PixelStore*, double);
static double dist(double*,double*);


void
CrossId_crossSamples(
        PixelStore *pixstore,
        double  radius_arcsec)
{
    /*
     * First build samples vectors. These are built from the lon and lat
     * sample values that are modified by astrometric resolution.
     */
    PixelStore_updateSamplePos(pixstore);

    /*
     * Convert radius in arcsec then in 3d euclidean distance
     * arcsec -> rad -> vector -> distance
     */
    double radius = radius_arcsec * TO_RAD;
    double va[3], vb[3];
    ang2vec(0, 0, va);
    ang2vec(HALFPI - radius, 0, vb);
    double radius_dist = dist(va, vb);

    /*
     * Iterate over all created Healpix pixels. They contains a least one
     * sample
     */

    int i;
#ifdef USE_THREADS
#define NTHREADS 28
    pthread_t *threads;
    struct thread_args *args;
    int *npixs;
    QMALLOC(threads, pthread_t, NTHREADS);
    QMALLOC(args, struct thread_args, NTHREADS);
    QMALLOC(npixs, int, NTHREADS);
    int np = pixstore->npixels / NTHREADS;
    for (i=0; i<NTHREADS; i++)
        npixs[i] = np;
    npixs[0] += pixstore->npixels % NTHREADS;
    int64_t *pixelindex = pixstore->pixelids;
    for (i=0; i<NTHREADS; i++) {
        struct thread_args *arg = &args[i];
        arg->store = pixstore;
        arg->radius = radius;
        arg->pixelindex = pixelindex;
        arg->npixs  = npixs[i];
        pthread_create(&threads[i], NULL, pthread_cross_pixel, arg);
        pixelindex += npixs[i];
    }

    for (i=0; i<NTHREADS; i++)
        pthread_join(threads[i], NULL);
    free(threads);
    free(args);
    free(npixs);
#else
    for (i=0; i<pixstore->npixels; i++)
        cross_pixel(
            PixelStore_get(pixstore, pixstore->pixelids[i]),
            pixstore,
            radius_dist);

#endif /* USE_THREADS */
}

static void*
pthread_cross_pixel(void *args)
{
    struct thread_args *ta = (struct thread_args*) args;
    int i;
    for (i=0; i<ta->npixs; i++) {
        HealPixel *pix = PixelStore_get(ta->store, ta->pixelindex[i]);
        cross_pixel(pix, ta->store, ta->radius);
    }
}


static void
cross_pixel(
        HealPixel *pix,
        PixelStore *store,
        double radius)
{

#ifdef USE_THREADS
    pthread_mutex_lock(&GLOBAL_MUTEX);
    pthread_mutex_lock(&pix->mutex);
#endif /* USE_THREADS */
    /*
     * Mark other pixels as done. They will then not cross identify with
     * this pixel later.
     */
    HealPixel *test_pix;
    int i, j;
    for (i=0; i<NNEIGHBORS; i++) {
        test_pix = pix->pneighbors[i];
        if (test_pix) {
            for (j=0; j<NNEIGHBORS; j++) {
                if (test_pix->neighbors[j] == pix->id) {
                    test_pix->tneighbors[j] = true;
                    break;
                }
            }
        }
    }
#ifdef USE_THREADS
    pthread_mutex_unlock(&GLOBAL_MUTEX);
#endif /* USE_THREADS */
    /*
     * For each sample of this pixel, find the best match. There are two
     * method for this:
     * - RAW_CLOSEST: find the closest sample in all upper fields.
     * - TIME_CLOSEST: for each upper fields starting from the lower, find a
     *   closest match . If it is found, stop iterating upper fields.
     */
    for (i=0; i<pix->nsamples; i++)
        cross_sample(pix->samples[i], pix, store, false, radius);
#ifdef USE_THREADS
    pthread_mutex_unlock(&pix->mutex);
#endif /* USE_THREADS */

}


static void
cross_sample(
        struct sample *current_spl,
        HealPixel *pix,
        PixelStore *store,
        bool force_neighbor_cross,
        double radius)
{
    int field_index_start;
    struct field_id fi;

    /*
     * Eliminate previous fields, we only match with upper fields
     * TODO this can be optimized, by using this function once for a set
     * of samples of the same field.
     */
    field_index_start = PixelStore_getHigherFields(pix, current_spl);

    int i;
    switch (closest) {
        case TIME_CLOSEST:
            fi.fieldindex = fi.havematch = 0;
            for (i=field_index_start; i<pix->nsamples; i++) {
                if (cross_time_closest_sample(current_spl, pix->samples[i],
                            radius, store, &fi) == 1)
                    break;
            }
            break;
        case RAW_CLOSEST:
            for (i=field_index_start; i<pix->nsamples; i++)
                crossmatch(current_spl, pix->samples[i], radius, store);
            break;
    }

    /*
     * Then iterate against neighbors pixels
     * TODO this can be optimized, by using this loop once for a set
     * of samples.
     */
    for (i=0; i<NNEIGHBORS; i++) {
        /* Allready crossed by an neighbor ? */
        if (force_neighbor_cross == false)
            if (pix->tneighbors[i] == true)
                continue;

        /*
         * Does the pixel exists? It may be a neighbor of current pixel,
         * but not be initialized because it does not contains
         * any samples.
         */
        HealPixel *test_pix = pix->pneighbors[i];
        if (test_pix == NULL)
            continue;

        /*
         * Eliminate previous to current fields, we only match with upper
         * fields.
         */
        field_index_start = PixelStore_getHigherFields(test_pix, current_spl);

        pthread_mutex_lock(&test_pix->mutex);
        int j;
        switch (closest) {
            case TIME_CLOSEST:
                fi.fieldindex = fi.havematch = 0;
                for (j=field_index_start; j<test_pix->nsamples; j++) {
                    if (cross_time_closest_sample(
                                current_spl, test_pix->samples[j],
                                radius, store, &fi) == 1)
                        break;
                }
                break;
            case RAW_CLOSEST:
                for (j=field_index_start; j<test_pix->nsamples; j++)
                    crossmatch(current_spl, test_pix->samples[j], radius, store);
                break;
        }
        pthread_mutex_unlock(&test_pix->mutex);

    }
}


static int
cross_time_closest_sample(
        struct sample *current_spl,
        struct sample *test_spl,
        double radius,
        PixelStore *store,
        struct field_id *current_field) {

        if (current_field->fieldindex != test_spl->set->field->fieldindex)
        {
            /* we are changing of field */
            if (current_field->havematch) {
                /* if previous field have a match break */
                return 0;
            } else {
                /* else change field and continue */
                current_field->fieldindex = test_spl->set->field->fieldindex;
            }
        }

        current_field->havematch = crossmatch(
                current_spl, test_spl, radius, store);

    return 1;
}


static int
crossmatch(
        struct sample *a,
        struct sample *b,
        double radius,
        PixelStore *store)
{

    double distance = dist(a->vector, b->vector);
    struct sample *relink = NULL;

    if (distance > radius)
        return 0;

    if (a->nextsamp && dist(a->vector, a->nextsamp->vector) <= distance)
        return 0;

    if (b->prevsamp && dist(b->vector, b->prevsamp->vector) <= distance)
        return 0;

    if (a->nextsamp)
        a->nextsamp->prevsamp = NULL;

    relink = b->prevsamp;

    a->nextsamp = b;
    b->prevsamp = a;

    if (relink)
        relink_sample(relink, store, radius);

    return 1;
}


static void
relink_sample(struct sample *sample, PixelStore *store, double radius)
{
    sample->nextsamp = NULL;
    HealPixel *pix = PixelStore_getPixelFromSample(store, sample);
    //cross_sample(sample, pix, store, true, radius);
}


/* up to one month of time distance would give a distance of 0 */
static double TIME_DISTANCE_LIMIT = 60 * 60 * 24 * 7 * 4;
static double TIME_DISTANCE_WEIGHT = 0.5;
static double
dist2(struct sample *a, struct sample *b)
{
    double x = a->vector[0] - b->vector[0];
    double y = a->vector[1] - b->vector[1];
    double z = a->vector[2] - b->vector[2];
    double t = abs(a->set->field->epoch - b->set->field->epoch) / TIME_DISTANCE_LIMIT * TIME_DISTANCE_WEIGHT;

    return x*x + y*y + z*z + t;

    
}

static double
dist(double *va, double *vb)
{
    double x = va[0] - vb[0];
    double y = va[1] - vb[1];
    double z = va[2] - vb[2];

    return x*x + y*y + z*z;
}
