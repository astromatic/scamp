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

static void crossmatch(struct sample*,struct sample*, double radius);
static long cross_pixel(HealPixel*,PixelStore*,double);

static long ntestmatches;

static pthread_mutex_t CMUTEX = PTHREAD_MUTEX_INITIALIZER;

static inline double
dist(double *va, double *vb)
{
    double x = va[0] - vb[0];
    double y = va[1] - vb[1];
    double z = va[2] - vb[2];

    return x*x + y*y + z*z;
}


#define NNEIGHBORS 8

struct thread_args {
    PixelStore  *store;
    int64_t  *pixelindex;
    int       npixs;
    double    radius;
    int      *result;
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
CrossId_crossSamples(
        PixelStore *pixstore,
        double  radius_arcsec,
        int   nthreads)
{
    int i;
    PixelStore_sort(pixstore);

    /* arcsec to radiant */
    double radius = radius_arcsec / 3600 * TO_RAD;

    /* radiant to vector */
    double va[3], vb[3];
    ang2vec(0,0,va);
    ang2vec(radius,0,vb);

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

static struct sample* get_sample_link_head(struct sample *s) {
    if (s->prevsamp)
        return get_sample_link_head(s->prevsamp);
    return s;
}
static struct sample* get_field_member(struct sample *s, struct field *f) {
    if (s == NULL)
        return NULL;
    if (s->set->field == f)
        return s;
    return get_field_member(s->nextsamp, f);
}
static void relink_sample_b(struct sample *ahead, struct sample *oldb, struct sample *newb) {
    if (ahead == oldb) {
        ahead->prevsamp->nextsamp = newb;
        newb->prevsamp = ahead->prevsamp;
        ahead->prevsamp = NULL;
        return;
    }
    relink_sample_b(ahead->nextsamp, oldb, newb);
}
static void relink_sample_a(struct sample *bhead, struct sample *olda, struct sample *newa) {
    if (bhead == olda) {
        olda->nextsamp->prevsamp = newa;
        newa->nextsamp = bhead->nextsamp;
        bhead->nextsamp = NULL;
    }
    relink_sample_a(bhead->nextsamp, olda, newa);
}
static void insert_sample_b(struct sample *head, struct sample *s) {
    if (head->epoch > s->epoch) {
        head->prevsamp->nextsamp = s;
        s->prevsamp = head->prevsamp;
        head->prevsamp = NULL;
        return;
    }
    if (head->nextsamp)
        return insert_sample_b(head->nextsamp, s);

    head->nextsamp = s;
    s->prevsamp = head;
}

static void insert_sample_a(struct sample *head, struct sample *s) {
    if (s->epoch > head->epoch) {
       s->nextsamp = head;
       if (head->prevsamp) {
            head->prevsamp->nextsamp = NULL;
       }
       head->prevsamp = s;
    }
    if (head->nextsamp)
        return insert_sample_a(head->nextsamp, s);

    head->prevsamp = s;
    s->nextsamp = head;

}
static inline void
crossmatch(struct sample *a, struct sample *b, double radius)
{

    /* Sample of the same field */
    if (a->set->field == b->set->field)
        return;

    /* Get distance between two samples */
    ntestmatches++;
    double a_b_distance = dist(a->vector, b->vector);

    /* If distance exceeds the limit, end here */
    if (a_b_distance > radius) {
        fprintf(stderr, "over radius\n");
        return;
    }

    /* Maybe switch a and b to assert that link will be from 
     * a->next to b->prev */
    if (a->epoch > b->epoch) {
        struct sample *a_tmp = a;
        a = b;
        b = a_tmp;
    }

    fprintf(stderr, "link done 1!\n");
    /* 
     * Get the current best distance from "a" to a sample from the field "b"
     * (if it exists). If current bestdist_a_to_fb exist and is lower than 
     * a_b_distance, return. 
     */
    struct sample *a_head  = get_sample_link_head(a);
    struct sample *b_from_a_old = get_field_member(a, b->set->field);
    if (b_from_a_old)
        if (a_b_distance > dist(b_from_a_old->vector, a->vector))
            return;

    fprintf(stderr, "link done 2!\n");
    /* 
     * A would want to link to "b", but maybe "b" have a better match to 
     * field "a" (in an other, allready crossed pixel). If current 
     * bestdist_b_to_fa exist and is lower than  a_b_distance return. 
     */
    struct sample *b_head  = get_sample_link_head(b);
    struct sample *a_from_b_old = get_field_member(b, a->set->field);
    if (a_from_b_old)
        if (a_b_distance > dist(a_from_b_old->vector, b->vector))
            return;

    fprintf(stderr, "link done 3!\n");
    if (b_from_a_old) {
        fprintf(stderr, "link done 3.1!\n");
        //relink_sample_b(a_head, b_from_a_old, b);
    } else {
        fprintf(stderr, "link done 3.2!\n");
        //insert_sample_b(a_head, b);
    }

    fprintf(stderr, "link done 4!\n");
    if (a_from_b_old) {
        fprintf(stderr, "link done 4.1!\n");
        //relink_sample_a(b_head, a_from_b_old, a);
    } else {
        fprintf(stderr, "link done 4.2!\n");
        //insert_sample_a(b_head, a);
    }

    fprintf(stderr, "link done!\n");
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


int
get_iterate_count()
{
    return ntestmatches;
}
