/**
 *
 * \file        crossid2.c
 * \brief       Efficient catalogs cross matching functions.
 * \author      Sébastien Serre
 * \date        11/04/2018
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
 * Here whith the use of "relink_up" and "CROSS_DOWN", we have the expected 
 * result in all cases:
 *
 * - "b" link "x" as his closest match, then
 * - "a" link "x" as his closest match, (unlink b->next)
 * - relink_up "b": found no valid match (or any other "bad" match),
 * - "a" link "y" as his closest match, (unlink x->prev)
 * - CROSS_DOWN "x": found "b" as his closest match
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
#include "crossid2.h"
#include "field.h"


/* local types */
#define NNEIGHBORS 8
typedef enum {CROSS_UP, CROSS_DOWN} cross_direction;
struct field_id {
    int fieldindex;
    int havematch;
};

/* forward declarations */
static void cross_sample(samplestruct*, HealPixel*, PixelStore*, double, cross_direction);
static int cross_time_closest_sample(samplestruct*,samplestruct*,
        double,PixelStore*, struct field_id*);
static int crossmatch(samplestruct*,samplestruct*, double, PixelStore*);


/**
 * \details
 * This is the function used to compare samples distance. The normal euclidean
 * ditance would be sqrt(x*x + y*y + z*z), but this one is enought to compare
 * and is faster.
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
 * \brief Cross match and link all samples
 * \param pixstore The pixel store (PixelStore)
 * \param radius_arcsec The max radius
 * \details
 * This is where it would start to implement multi threading. 
 */
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
     * sample.
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

    int field_index_pivot;
    struct field_id fi;

    /*
     * Iterate over current_spl pixel samples.
     *
     * If i want to CROSS_UP, eliminate all samples older or bellonging to the
     * field of current_spl, then iterate from the closer samples in time,
     * and stop if a match is found in a field. "cross_time_closest" take
     * care of taking the best sample for a field.
     */
    int i;
    fi.fieldindex = -1;
    fi.havematch  = 0;
    switch (direction) {
        case CROSS_UP:
            field_index_pivot = PixelStore_getHigherFields(pix, current_spl);
            for (i=field_index_pivot; i<pix->nsamples; i++) {
                if (cross_time_closest_sample(
                            current_spl, pix->samples[i],
                            maxdist, store, &fi) == 1)
                    break;
            }
            break;
        case CROSS_DOWN:
            field_index_pivot = PixelStore_getLowerFields(pix, current_spl);
            for (i=field_index_pivot; i>-1; i--) {
                if (cross_time_closest_sample(
                            pix->samples[i], current_spl,
                            maxdist, store, &fi) == 1)
                    break;
            }

            break;
    }

    /*
     * Then iterate against neighbors pixels
     */
    for (i=0; i<NNEIGHBORS; i++) {

        /*
         * Does the pixel exists? It may be a neighbor of current pixel,
         * but not be initialized because it does not contains samples.
         */
        HealPixel *test_pix = PixelStore_get(store, pix->neighbors[i]);
        if (test_pix == NULL)
            continue;

        /*
         * Again only match a portion of the pix->sample, and stop if a field
         * contains a match.
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


/**
 * \brief Called from a list iterator over samples, return 1 when a match is found.
 * \param current_spl The sample we want to link
 * \param test_spl The sample we test
 * \param maxdist maximum radius
 * \param store the pixel store
 * \param current_field the structure used to keep track of the last match/field.
 * \return 0 when no matches are found, 1 where there is.
 * \details This function means is to iterate samples from different fields and
 * to stop when a match is found AND we are changing field, wich means that we
 * have the previous field best sample match. This function must be called
 * on a sample array sorted in the desired time order.
 */
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

    /*
     * Might set a->nextsamp equal to b, and b->prevsamp equal to a.
     *
     * So, "a" and "b" must comme allready sorted by epoch, 
     * with "a" being older tha "b". Calling this function in the opposite 
     * order is a bug.
     */
    //assert(PixelStore_compare(a, b) < 0);


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
         * Compare the new candiate and the old
         */
        switch (PixelStore_compare(b, a->nextsamp)) {
            case 1:
                /* 
                 * b is from a field farther in time from a than a->nextsamp. 
                 * Stop here.
                 */
                return 0;
            case 0:
                /*
                 * "a->nextsamp" belong to the same field layer than "b". So we
                 * compare distances, and go further if "b" is a better match.
                 */
                if (dist(a->vector, a->nextsamp->vector) <= distance)
                    return 0;
        }

        /* Else whatever the distance is, "a" will want to link to "b" */

    }

    /*
     * The same as above, but the opposite (search in descending order)
     */
    if (b->prevsamp) {

        switch (PixelStore_compare(a, b->prevsamp)) {
            case -1:
                return 0;
            case 0:
                if (dist(b->vector, b->prevsamp->vector) <= distance)
                    return 0;
        }
        /* Else whatever the distance is, "b" will want to link to "a" */
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

    /* alleluia */
    return 1;
}


