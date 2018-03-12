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
#include <omp.h>
#include <time.h>
#include <pthread.h>

#include "crossid2.h"
#include "chealpixstore.h"

static void crossmatch(struct sample*,struct sample*);
static long cross_pixel(HealPixel*,PixelStore*,double);

static long ntestmatches;

static pthread_mutex_t CMUTEX = PTHREAD_MUTEX_INITIALIZER;

#define NNEIGHBORS 8

struct thread_args {
    PixelStore  *store;
    int64_t  *pixelindex;
    int   npixs;
    double   radius;
    int   *result;
};


static void*
pthread_cross_pixel(void *args) 
{
    struct thread_args *ta = (struct thread_args*) args;

    int i;
    int nmatches = 0;
    for (i=0; i<ta->npixs; i++) {
        HealPixel *pix = PixelStore_get(ta->store, ta->pixelindex[i]);
        nmatches += cross_pixel(pix, ta->store, ta->radius);
    }

    *(ta->result) = nmatches;
    return NULL;
}


long
Crossmatch_crossSamples(
        PixelStore *pixstore,
        double  radius_arcsec,
        int   nthreads)
{
    int i;

    /* arcsec to radiant */
    double radius = radius_arcsec / 3600 * TO_RAD;
    PixelStore_setMaxRadius(pixstore, radius);


    /* allocate mem */
    pthread_t *threads; 
    struct thread_args *args; 
    int *results;    
    int *npixs;
    QMALLOC(threads, pthread_t, nthreads);
    QMALLOC(args, struct thread_args, nthreads);
    QMALLOC(results, int, nthreads);
    QMALLOC(npixs, int, nthreads);


    /* distribute work between threads */
    int np = pixstore->npixels / nthreads;
    for (i=0; i<nthreads; i++)
        npixs[i] = np;
    npixs[0] += pixstore->npixels % nthreads;


    /* start threads */
    int64_t *pixelindex = pixstore->pixelids;
    for (i=0; i<nthreads; i++) {

        /* construct thread argument structure */
        struct thread_args *arg = &args[i];
        arg->store   = pixstore;
        arg->radius  = radius;
        arg->pixelindex = pixelindex;
        arg->npixs   = npixs[i];
        arg->result  = &results[i];

        /* launch! */
        pthread_create(&threads[i], NULL, pthread_cross_pixel, arg);

        /* increment pixelindex for next thread */
        pixelindex += npixs[i];

    }


    /* wait for joins */
    for (i=0; i<nthreads; i++)
        pthread_join(threads[i], NULL);


    /* reduce */
    long nmatches = 0;
    for (i=0; i<nthreads; i++)
        nmatches += results[i];


    /* cleanup */
    free(threads);
    free(args);
    free(npixs);
    free(results);


    fprintf(stderr,
            "Crossmatch end: %li matches for all pixels!\n", nmatches);

    return nmatches;

}

/**
 * This function prevent dead locks.
 *
 * THIS is the trickiest part of this file. pneighbors is are pointers,
 * tneighbors are boolean values indicating that the cross identification
 * between two pixels is done.
 *
 * When a pixel want to cross his neigbhors, he first execute this, setting
 * all neighbors "tneigbhors" value for himself to true. Neigbhors will then
 * state that they do not need to cross. This is done only if my 
 * "tneighbor" value for the foreign pixel is not true. If my "tneighbor" is
 * true, this would indicate that the neighbor has allready reserved the 
 * cross identification. I must not touch his "tnieghbor" then, or he will 
 * not cross it two, and the pixels would never been cross identified.
 *
 * Called by cross_pixel, to notify neighbors that I handle myself for the
 * rest of the run. So do not cross with me.
 *
 */
static void
set_reserve_cross(HealPixel *a) 
{
    int i, j;
    HealPixel *b;

    pthread_mutex_lock(&CMUTEX);

    for (i=0; i<NNEIGHBORS; i++) {
        b = a->pneighbors[i];
        /* 
         * test if b is not null and if my self test for this neighbor
         * is positive wich signifies that this "b" has reserved the
         * crossmatch from him to me ("a") 
         */
        if (b && a->tneighbors[i] != true) {
            for (j=0; j<NNEIGHBORS; j++) {
                if (a == b->pneighbors[j]) {
                    b->tneighbors[j] = true;
                    break;
                }
            }
        }
    }

    pthread_mutex_unlock(&CMUTEX);

}


static long
cross_pixel(HealPixel *pix, PixelStore *store, double radius) 
{

    set_reserve_cross(pix);
    pthread_mutex_lock(&pix->mutex);

    long nbmatches = 0;

    /*
     * Iterate over HealPixel structure which old sample structures
     * belonging to him.
     */
    long j, k, l;

    struct sample *current_spl;
    struct sample *test_spl;


    for (j=0; j<pix->nsamples; j++) {
        current_spl = pix->samples[j];

        /*
         * First cross match with samples of the pixel between them
         */
        for(k=0; k<j; k++) {
            test_spl = pix->samples[k];

            if (current_spl->set->field == test_spl->set->field)
                continue;

            /*
            if (abs(current_spl->col - test_spl->col) > radius)
                continue;
                */

            crossmatch(current_spl, test_spl);

        }

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
             * Ok, then lock test pixel and iterate over samples.
             */
            pthread_mutex_lock(&test_pixel->mutex);

            for (l=0; l<test_pixel->nsamples; l++) {
                test_spl = test_pixel->samples[l];

                if (current_spl->set->field == test_spl->set->field)
                    continue;

                /*
                if (abs(current_spl->col - test_spl->col) > radius)
                    continue;
                    */

                crossmatch(current_spl, test_spl);

            }

            pthread_mutex_unlock(&test_pixel->mutex);
        }

        /*
           if (current_spl->bestMatch != NULL) {
           nbmatches++;
           }
         */
    }

    pthread_mutex_unlock(&pix->mutex);

    return nbmatches;

}


int
get_iterate_count() 
{
    return ntestmatches;
}


static inline double
dist(double *va, double *vb) 
{
    double x = va[0] - vb[0];
    double y = va[1] - vb[1];
    double z = va[2] - vb[2];

    /* TODO return x*x + y*y + z*z; */
    return sqrt(x*x + y*y + z*z);
}

static void
crossmatch(struct sample *current_spl, struct sample *test_spl) 
{
    ntestmatches++;
    /*
     * Get distance between samples
     */
    double distance = dist(current_spl->vector, test_spl->vector);

    /*
     * If distance is less than previous (or initial) update
     * best_distance and set test_spl to current_spl.bestMatch
     * TODO XXX another lock here?
     if (distance < current_spl->bestMatchDistance) {
     current_spl->bestMatch = test_spl;   / XXX false shared ! /
     current_spl->bestMatchDistance = distance; / XXX false shared ! /
     }

     if (distance < test_spl->bestMatchDistance) {
     test_spl->bestMatch = current_spl;
     test_spl->bestMatchDistance = distance;
     }
     */

}
