/**
 *
 * \file        crossid.c
 * \brief       Efficient catalogs cross matching functions.
 * \author      Sébastien Serre
 * \date        07/05/2018
 *
 * \copyright   Copyright (C) 2017 University of Bordeaux. All right reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * \details
 *
 * This file is in charge of linking samples to their best match, traversing
 * all fields in epoch order. It uses PixelStore as sample container, and
 * use it to iterate over all pixels. The most important functions in this
 * file are the static functions "cross_sample" and "crossmatch".
 *
 * I've found that two fields and four samples are enought to represent all
 * cases I have identified and have to deal with. The problem is that the
 * order in wich matches are done, can leave unliked samples. These example
 * are high level, ignoring pixels, just talking about distance between samples
 * and ordering in time.
 *
 * Let's take this as a first work case with field1 an field2 sorted by epoch
 * (with field1 older) and four samples "a", "b",  "y" and "z":
 *
 * \code{.unparsed}
 * field2 ---y----x--------------
 * field1 ------a--b-------------
 * \endcode
 *
 * \par Case 1 (the fine case)
 * - "b" link "x" as his closest match,
 * - "a" link "y" as his closest match, (would have taken "x", but "x" have a
 * better  match with "b")
 * - we have iterated all field1 samples, this end here for this field.
 *
 * Everything is fine.
 *
 * \par Case 2 (a problem)
 * - "a" link "x" as his closest match, (better than "y"), then
 * - "b" link "x" as his closest match, (unlinking a->next bestmach)
 * - we have iterated all field1 samples, this end here for this field.
 *
 * In this case, "a" will never match "y" where he should.
 *
 * So when "b" is linking himself to "x", it must tell "a" (the previous "x"
 * match) to search a new match in upper fields. This is done in the new
 * algorythm as "CROSS_UP" option. In this case, "CROSS_UP" would just
 * add a new step after unlinking any sample->next match. Example:
 *
 * - "a" link "x" as his closest match, (better than "y"), then
 * - "b" link "x" as his closest match, (unlinking a->next bestmach)
 * - CROSS_UP "a": "a" link "y" as his closest match (and fail for "x"
 * because "b" is a better match for "x")
 * - we have iterated all field1 samples, this end here for this field.
 *
 * Fine.
 *
 * \par Case 3 (CROSS_UP is not enought)
 *
 * Another test case:
 * \code{.unparsed}
 * field2 -----y--x--------------
 * field1 ------a----b-----------
 * \endcode
 *
 * We use relink_up:
 * - "b" link "x" as his closest match, then
 * - "a" link "x" as his closest match, (unlinking b->next),
 * - relink_up(b) found no valid match (or any other "bad" match).
 * - "a" link "y" as his closest match, (it is better than "x") and unlink
 * him from x->prev.
 * - we have iterated all field1 samples, this end here for this field.
 *
 * In this case, "b" and "x" are not matching, when they should. Plus the
 * order in wich the matches have been done have deleted the relation between
 * "x" and "b", "relink_up" does not sove this.
 *
 * So, when "a" link to "y" and unlink from "x", it left "x" prevsamples
 * empty. We need to tell "x" to find a new match back in time. This is done
 * with the "CROSS_DOWN" logic in the new code.
 *
 * \par Case 4 (fine again)
 *
 * Here whith the use of "relink_up" and "relink_down", we have the expected
 * result in all cases:
 *
 * - "b" link "x" as his closest match, then
 * - "a" link "x" as his closest match, (unlink b->next)
 * - relink_up "b": found no valid match (or any other "bad" match),
 * - "a" link "y" as his closest match, (unlink x->prev)
 * - relink_down "x": found "b" as his closest match
 * - we have iterated all field1 samples, this end here for this field.
 *
 * Everything is fine, "a" matches "y" and "b" matches "x".
 *
 * Note that these steps can recurse, because they might unlink
 * other samples. So a simple crossmatch can recurse many times
 * change and relink many samples until it reaches a stable state.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "chealpix.h"
#include "chealpixstore.h"
#include "crossid.h"
#include "field.h"
#include "prefs.h"


/* local types */
#define NNEIGHBORS 8
typedef enum {CROSS_UP, CROSS_DOWN} cross_direction;
struct field_id {
    int fieldindex;
    int havematch;
};


/* forward declarations */
static void cross_sample(samplestruct*, HealPixel*, PixelStore*, double, cross_direction);
static void try_match(cross_direction, samplestruct*, HealPixel*, PixelStore*, double);
static int  crossmatch(samplestruct*,samplestruct*, double, PixelStore*);
static void move_refsamples(PixelStore *pixstore);




/**
 * \details
 * This is the function used to compare samples distance. The normal euclidean
 * ditance would be sqrt(x*x + y*y + z*z), (x*x+y*y+z*z) but this one is
 * enought to compare and is faster.
 */
static inline double
dist(double *va, double *vb)
{
    double x = va[0] - vb[0];
    double y = va[1] - vb[1];
    double z = va[2] - vb[2];

    return x*x + y*y + z*z;
}

/**
 * \brief cross identify and link samples
 * \param field of groups
 * \param reffield reference field
 * \param cross-matching radius in degrees
 */
void
CrossId_run(
        fgroupstruct *fgroup,
        fieldstruct *reffield,
	double radius)
{

    fieldstruct **fields;
    int		nfield;

    nfield = fgroup->nfield;
    fields = fgroup->field;

    PixelStore ps;

    /* Define healpix resolution from the match radius:
     * 180.0 is the base ring size in degrees, nsides must be a
     * power of two that will build a pixel definition that fullfill the need
     * of rings.
     *
     * This get the nearest next power of two from the ideal total rings
     * required. [https://en.wikipedia.org/wiki/Logarithm#Change_of_base]
     * We set a 2.1e6 limit for the number of rings, as it corresponds to
     * nsides = 2^29
     */

    double nrings = 180.0 / radius;
    if (nrings > 2.1e6)
      nrings = 2.1e6;

    int64_t nsides_pow = (int64_t)ceil(log(nrings / 4.0 + 1.0) / log(2));

    /* minus 1 nsides power, to be sure to not loss any matches */
    if (nsides_pow>0)
      nsides_pow--;

    int64_t nsides = 2L << nsides_pow;

    PixelStore_new(nsides, &ps);

    struct set *set;
    int i, j, k;

    for (i=0; i<nfield; i++) {
        for (j=0; j<fields[i]->nset; j++) {
            set = fields[i]->set[j];
            for (k=0; k < set->nsample;k++) {
                struct sample *s = &set->sample[k];
                PixelStore_add(&ps, &set->sample[k]);
            }
        }
    }

    if (prefs.astrefcat != ASTREFCAT_NONE) {
        for (j=0; j<reffield->nset; j++) {
            set = reffield->set[j];
                for (k=0; k < set->nsample; k++) {
                    struct sample *s = &set->sample[k];
                    PixelStore_add(&ps, &set->sample[k]);
                }
        }
    }

    /* sort samples */
    PixelStore_sort(&ps);

    /*
     * Convert radius in arcsec then in 3d euclidean distance
     * arcsec -> rad -> vector -> distance
     */
    double radius_radian = radius * TO_RAD;
    double va[3], vb[3];
    ang2vec(0, 0, va);
    ang2vec(radius_radian, 0, vb);
    double maxdist = dist(va, vb);

    /*
     * Iterate over all created Healpix pixels. They contains a least one
     * sample.
     */
    for (i=0; i<ps.npixels; i++) {
        HealPixel *pix = PixelStore_get(&ps, ps.pixelids[i]);
        for (j=0; j<pix->nsamples; j++) {
            cross_sample(pix->samples[j], pix, &ps, maxdist, CROSS_DOWN);
        }
    }

    move_refsamples(&ps);
    PixelStore_free(&ps);
}



/**
 * \brief Take a sample and link it with his best match in "cross_direction"
 * \param current_spl The sample we want to link
 * \param pix The sample pixel
 * \param store The pixel store (to get neighbors)
 * \param maxdist The maximum radius distance
 * \param direction Match back or forward in time
 */
static void
cross_sample(
        samplestruct *current_spl,
        HealPixel *pix,
        PixelStore *store,
        double maxdist,
        cross_direction direction)
{

    try_match(direction, current_spl, pix, store, maxdist);

    int i;
    for (i=0; i<NNEIGHBORS; i++) {

        HealPixel *test_pix = PixelStore_get(store, pix->neighbors[i]);

        if (test_pix != NULL)
            try_match(direction, current_spl, test_pix, store, maxdist);

    }
}


#define A_HAVE_CLOSEST 1
#define B_HAVE_CLOSEST 2
static void
try_match(
        cross_direction direction,
        samplestruct *current_spl,
        HealPixel    *pix,
        PixelStore   *store,
        double       maxdist)
{
    int i, pivot;

    /* TODO optimize, when a match is done on a field for current_pix, do
       not test farther fields */
    if (direction == CROSS_UP) {
        pivot = PixelStore_getHigherFields(pix, current_spl);
        for (i=pivot; i<pix->nsamples; i++) {
            samplestruct *test_spl = pix->samples[i];
            if (crossmatch(current_spl, test_spl, maxdist, store) == A_HAVE_CLOSEST)
                break;
        }
    } else {
        pivot = PixelStore_getLowerFields(pix, current_spl);
        for (i=pivot; i>-1; i--) {
            samplestruct *test_spl = pix->samples[i];
            if (crossmatch(test_spl, current_spl, maxdist, store) == B_HAVE_CLOSEST)
                break;
        }
    }
}

/**
 * \brief Crossmatch two samples
 * \param a the first sample (older than b)
 * \param b the seconds sample (newer than a)
 * \param maxidist maximum radius
 * \param store our pixel store
 * \details This function can call "cross_sample" and then recurse. See the
 * documentation in this header file, and the code comments to understand
 * the how and why of his behaviour.
 *
 */
static int
crossmatch(
        samplestruct *a,
        samplestruct *b,
        double maxdist,
        PixelStore *store)
{

    double distance = dist(a->vector, b->vector);
    if (distance > maxdist)
        return 0;

    if (a->nextsamp) {
        switch (PixelStore_compare(a->nextsamp, b)) {
            case -1:
                return A_HAVE_CLOSEST;
            case 0:
                if (dist(a->nextsamp->vector, a->vector) <= distance)
                    return 0;
        }
    }

    if (b->prevsamp) {
        switch (PixelStore_compare(a, b->prevsamp)) {
            case -1:
                return B_HAVE_CLOSEST;
            case 0:
                if (dist(b->vector, b->prevsamp->vector) <= distance)
                    return 0;
        }
    }

    /*
     * Here, we want to link a->b. But we also want previous linked samples
     * to search another link.
     */

    /* If an old link exists, take it and reset their link information */
    samplestruct *relink_up    = b->prevsamp;
    samplestruct *relink_down  = a->nextsamp;

    /* unlink old matches */
    if (relink_down)
        relink_down->prevsamp = NULL;
    if (relink_up)
        relink_up->nextsamp = NULL;

    /* do the new match linkage */
    a->nextsamp = b;
    b->prevsamp = a;

    /*
     * We have now reached a stable state, with a and b linked, and possibly
     * old matches linkss pointing to NULL samples. This is where the old crossid
     * algorithm was ending.
     *
     * We now want to recurse for each unlinked samples so they can find a new
     * best match. There are many cases that needed this recursion (see doc in
     * this header file).
     */
    if (relink_up) {
        HealPixel *pix = PixelStore_getPixelFromSample(store, relink_up);
        cross_sample(relink_up, pix, store, maxdist, CROSS_UP);
    }

    if (relink_down) {
        HealPixel *pix = PixelStore_getPixelFromSample(store, relink_down);
        cross_sample(relink_down, pix, store, maxdist, CROSS_DOWN);
    }

    return 0;
}

/**
 * \brief Move every refsample at the begining of the linked samples.
 * \details TODO XXX: This should take special action when two refsamples are
 * in the same list. Here the first is set, and the others are removed.
 */
static void
move_refsamples(PixelStore *pixstore)
{
    int i;
    for (i=0; i<pixstore->npixels; i++) {
        HealPixel *pix = PixelStore_get(pixstore, pixstore->pixelids[i]);
        int j;
        for (j=0; j<pix->nsamples; j++) {

            samplestruct *spl = pix->samples[j];

            if (spl->set->field->isrefcat) {
                if (spl->prevsamp == NULL) {
                    /* no prevsamp, I am allready in first position, or have
                       no linked samples */
                    continue;
                }


                samplestruct *prev  = spl->prevsamp;
                if (spl->nextsamp) {
                    /* There is both a spl->prev and a spl->next, link them */
                    samplestruct *next = spl->nextsamp;
                    prev->nextsamp = next;
                    next->prevsamp = prev;
                } else {
                    prev->nextsamp = NULL;
                }

                spl->prevsamp = spl->nextsamp = NULL;

                /* rewind to the first sample */
                samplestruct *first = prev;
                while (first->prevsamp)
                    first = first->prevsamp;


                /* if there is allready a refsample first, discard
                   TODO XXX: take appropriate action for this, in case two ref
                   are in the match. This can append with two ref samples close
                   from each other and having différent epoc value. */
                if (first->set->field->isrefcat)
                    continue;

                /* then link our spl on top of the samples link */
                spl->nextsamp = first;
                first->prevsamp = spl;
            }

        }
    }
}

