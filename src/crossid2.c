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
#include "chealpix.h"
#include "chealpixstore.h"

#define NNEIGHBORS 8

static void crossmatch(struct sample*,struct sample*, double radius);
static long cross_pixel(HealPixel*,PixelStore*,double);


static pthread_mutex_t CMUTEX = PTHREAD_MUTEX_INITIALIZER;

enum join_direction {RIGHT_JOIN, LEFT_JOIN};
struct thread_args {
    PixelStore  *store;
    int64_t  *pixelindex;
    int       npixs;
    double    radius;
    int      *result;
};


static inline double
dist(double *va, double *vb)
{
    double x = va[0] - vb[0];
    double y = va[1] - vb[1];
    double z = va[2] - vb[2];

    return x*x + y*y + z*z;
}




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
CrossId_crossSamples(
        PixelStore *pixstore,
        double  radius_arcsec,
        int   nthreads)
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
 * THIS is the trickiest part of this file (crossmatch is nice two). 
 * pneighbors are pointers,
 * tneighbors are boolean values indicating that the cross identification
 * between two pixels must be considered as done.
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


/* compare two samples by epoch first, the most important, then by field 
   pointer. */
static int
cmp_samples(
        struct sample *a, 
        struct sample *b) 
{
    if (!a) return 1;
    if (!b) return -1;
    if (a->epoch == b->epoch)
        return a->set->field > b->set->field ? 1 : a->set->field == b->set->field ? 0 : - 1;
    return a->epoch > b->epoch ? 1 : -1;
}


/* from a linked sample "s" get a sample that belong to the field "f" */
static inline struct sample*
get_field_sample(
        struct sample *s, 
        struct field *f)
{
    /* TODO optimize with a kind of bsearch(3) */
    while (s->prevsamp)
        s = s->prevsamp;

    do {
        if (s->set->field == f)
            return s;
    } while (s = s->nextsamp);

    return NULL;
}


/* 
 Here we join left and right:

For left as L and right as R, numbers in variables represents fields sorted by
epoch, uppercast var name, represent our sample:
L = {l1, l2, l3,  L4,  l5, l6, l7, l8, l13} 
R = {r5, r6, r9,  R10,  r11, r12}

With LEFT_JOIN it will look like:
{l1, l2, l3,  L4,  l5, l6, l7, l8,  R10,  r11, r12}

With RIGHT_JOIN it will look like:
{l1, l2, l3,  L4,  r5, r6, r9,  R10,  r11, r12}

*/
static void
join_samples(
        struct sample *left, 
        struct sample *right, 
        enum join_direction join) 
{

    /* get head of both linked samples */
    struct sample *left_head = left;
    while (left_head->prevsamp)
        left_head = left_head->prevsamp;

    struct sample *right_head = right;
    while (right_head->prevsamp)
        right_head = right_head->prevsamp;

    /* get the total number of samples */
    struct sample *c;
    int nsamples = 0;

    c = left_head;
    do { nsamples++; } while (c = c->nextsamp);
    c = right_head;
    do { nsamples++; } while (c = c->nextsamp);

    /* allocate a temporary array for nsamples, + 1 to allways have a last
       NULL sample for final linking */
    struct sample **out = calloc(nsamples + 1, sizeof(struct sample*));

    /* first take everythink from left_head to left */
    int i = 0;
    while (left_head != left) {
        out[i++] = left_head;
        left_head = left_head->nextsamp;
    }
    out[i++] = left_head;
    left_head = left_head->nextsamp;


    struct sample *right_remaining;
    switch(join) {
        case LEFT_JOIN:
            /* from there, take any left_head until right is reached */
            while (cmp_samples(left_head, right) < 0) {
                out[i++] = left_head;
                left_head = left_head->nextsamp;
            }
            right_remaining = right;
            break;

        case RIGHT_JOIN:
            /* from there, take any left_head until right_head is reached */
            while (cmp_samples(left_head, right_head) < 0) {
                out[i++] = left_head;
                left_head = left_head->nextsamp;
            }
            right_remaining = right_head;
            while (right_head) {
                out[i++] = right_head;
                right_head = right_head->nextsamp;
            }
            break;
    }

    /* consume remaining right samples */
    while (right_remaining) {
        out[i++] = right_remaining;
        right_remaining = right_remaining->nextsamp;
    }

    /* Relink all array of samples, they are correctly ordered */
    struct sample *prev = NULL;
    struct sample *current;
    int j;
    for (j=0;j<i; j++) {
        struct sample *current = out[j];
        current->prevsamp = prev;
        current->nextsamp = out[j+1];
        prev = current;
    }

    /* end WTF! */

    /* TODO maybe, avoid using malloc, and have to go to links head and count 
       every link of samples. XXX This will make the code unreadable. but 
       certainly faster. */

    free(out);
}


/* compare two samples and maybe join them */
static inline void
crossmatch(
        struct sample *a, 
        struct sample *b, 
        double radius)
{

    /* Sample of the same field, no need to go further */
    if (a->set->field == b->set->field)
        return;

    /* Get distance between two samples */
    double a_b_distance = dist(a->vector, b->vector);

    /* If distance exceeds the max radius, end here */
    if (a_b_distance > radius)
        return;

    /* get best match from field b for a and terminate if the current match
     is better ... */
    struct sample *a_best_bfield = get_field_sample(a, b->set->field);
    if (a_best_bfield)
        if (a_b_distance > dist(a->vector, a_best_bfield->vector))
            return;

    /* ... the opposite */
    struct sample *b_best_afield = get_field_sample(b, a->set->field);
    if (b_best_afield)
        if (a_b_distance > dist(b->vector, b_best_afield->vector))
            return;

    /* join samples */
    if (cmp_samples(a,b) > 0) {
        join_samples(b, a, RIGHT_JOIN);
    } else {
        join_samples(a, b, LEFT_JOIN);
    }

    return;
}


/* cross all samples from one pixel to himself, and all neighbors pixel 
   samples. */
static long
cross_pixel(
        HealPixel *pix, 
        PixelStore *store, 
        double radius)
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

            /*
            if (abs(current_spl->col - test_spl->col) > radius)
                continue;
             */

            crossmatch(current_spl, test_spl, radius);

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

                /*
                if (abs(current_spl->col - test_spl->col) > radius)
                    continue;
                    */

                crossmatch(current_spl, test_spl, radius);

            }

            pthread_mutex_unlock(&test_pixel->mutex);
        }

    }

    pthread_mutex_unlock(&pix->mutex);

    return nbmatches;

}
