/**
 *
 * \file        chealpixstore.c
 * \brief       Healpix pixels storage mechanism.
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
 * This file is divided in four parts:
 * - public interface
 * - static functions
 * - AVL tree implementations,
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "chealpix.h"
#include "chealpixstore.h"
#include "field.h"

#define PIXELIDS_BASE_SIZE 1000 /* initial malloc or PixelStore->pixelids */
#define SPL_BASE_SIZE 10        /* initial malloc of HealPixel->samples */

/* forward AVL tree function declarations */
typedef struct pixel_avl pixel_avl;
struct pixel_avl {
    HealPixel pixel;        /* The samples being stored. Key is at pixel.id */
    pixel_avl *pBefore;     /* Other elements less than pixel.id */
    pixel_avl *pAfter;      /* Other elements greater than pixel.id */
    pixel_avl *pUp;         /* Parent element */
    short int height;       /* Height of this node.  Leaf==1 */
    short int imbalance;    /* Height difference between pBefore and pAfter */
};
static pixel_avl *pixelAvlSearch(pixel_avl*, const long);
static pixel_avl *pixelAvlInsert(pixel_avl **, pixel_avl *);
static void pixelAvlFree(pixel_avl*);

/* forward utility declaration */
static int cmp_samples(const void* a, const void *b);
static void pixelAvlSort(pixel_avl*);
static void pixelAvlPrint(pixel_avl*);


/**
 * \brief create a new pixel store 
 * \param nsides the healpix nsides resolution parameter
 * \param store
 * \details
 * The store must be released with PixelStore_free when no more used.
 */
void
PixelStore_new(int64_t nsides, PixelStore *store)
{
    store->pixels = NULL;
    store->nsides = nsides;
    store->npixels = 0;
    QMALLOC(store->pixelids, int64_t, PIXELIDS_BASE_SIZE);
    store->pixelids_size = PIXELIDS_BASE_SIZE;
}


/**
 * \brief free everything 
 * \param store
 */
void
PixelStore_free(PixelStore* store)
{
    pixelAvlFree((pixel_avl*) store->pixels);
    free(store->pixelids);
}


/**
 * \brief add a sample to a pixel store
 * \param store
 * \param spl the sample
 */
void
PixelStore_add(
        PixelStore      *store,
        samplestruct   *spl)
{

    spl->nextsamp = spl->prevsamp = NULL;

    double lon = spl->wcspos[spl->set->field->lng] * TO_RAD;
    double col = HALFPI - (spl->wcspos[spl->set->field->lat] * TO_RAD);
    ang2vec(col, lon, spl->vector);

    int64_t pixnum;
    ang2pix_nest64(store->nsides, col, lon, &pixnum);

    /* search for the pixel */
    pixel_avl *avlpix =
        pixelAvlSearch((pixel_avl*) store->pixels, pixnum);

    if (pixnum == 0)
        fprintf(stderr, "have pixnum 0\n");

    if (!avlpix) { // no such pixel

        int i;
        QCALLOC(avlpix, pixel_avl, 1);
        avlpix->pixel.id = pixnum;
        QCALLOC(avlpix->pixel.samples, samplestruct*, SPL_BASE_SIZE);
        avlpix->pixel.size = SPL_BASE_SIZE;
        avlpix->pixel.nsamples = 0;

        neighbours_nest64(store->nsides, pixnum, avlpix->pixel.neighbors);

        pixelAvlInsert((pixel_avl**) &store->pixels, avlpix);

        if (store->pixelids_size == store->npixels) {
            QREALLOC(store->pixelids,
                 long, store->pixelids_size * 2);
            store->pixelids_size *= 2;
        }
        store->pixelids[store->npixels] = pixnum;
        store->npixels++;

    }

    HealPixel *pix = &avlpix->pixel;
    if (pix->nsamples == pix->size) {
        QREALLOC(pix->samples, samplestruct*, pix->size * 2);
        pix->size *= 2;
    }

    pix->samples[pix->nsamples] = spl;
    pix->nsamples++;
}


/**
 * \brief get a piwel by index 
 * \param store
 * \param key
 */
HealPixel*
PixelStore_get(
        PixelStore *store,
        int64_t  key)
{
    pixel_avl *match_avl;

    match_avl = pixelAvlSearch((pixel_avl*) store->pixels, key);
    if (!match_avl)
        return (HealPixel*) NULL;
    return &match_avl->pixel;

}


/**
 * \brief return the pixel to wich the sample belong to
 * \param store
 * \param sample
 */
HealPixel*
PixelStore_getPixelFromSample(
        PixelStore      *store,
        samplestruct    *sample)
{
    double lon, col;
    vec2ang(sample->vector, &col, &lon);

    int64_t pixnum;
    ang2pix_nest64(store->nsides, col, lon, &pixnum);

    return PixelStore_get(store, pixnum);
}


/**
 * \brief sort samples by epoch and fields if epoch are equals
 * \param store
 */
void
PixelStore_sort(PixelStore* store)
{
    pixelAvlSort((pixel_avl*) store->pixels);
}


/**
 * \brief return an index value to the first next test sample
 * \param pix
 * \param pivot
 * \details
 * return the index at wich wen can start crossmatching for the sample "pivot"
 * avoiding unmatchables samples that have a timestamp lesser than the one
 * of the sample.
 */
int
PixelStore_getHigherFields(
        HealPixel       *pix,
        samplestruct    *pivot)
{
    int i;
    for (i=0; i<pix->nsamples; i++) {
        if (cmp_samples(&pivot, &pix->samples[i]) < 0)
            break;
    }
    return i;
    /*
       TODO test and use this optimized version.
    int max = pix->nsamples;
    int min = 0;
    int i;

    while (min < max)
    {
        i = (min + max) / 2;
        if (cmp_samples(&pivot, &pix->samples[i]) < 0)
            max = i;
        else
            min = i + 1;
    }

    return max;
    */
}


/**
 * \brief return an index value to the previous next test sample
 * \param pix
 * \param pivot
 * \details
 * return the index at wich wen can start crossmatching for the sample "pivot"
 * avoiding unmatchables samples that have a timestamp gretter than the one
 * of the sample.
 * Iterating with this value, must be done in reverse order:
 * for (i=ret; i>-1; i--)
 */
int
PixelStore_getLowerFields(
        HealPixel       *pix,
        samplestruct    *pivot)
{
    int i;
    for (i=0; i<pix->nsamples; i++) {
        if (cmp_samples(&pivot, &pix->samples[i]) <= 0)
            break;
    }
    return --i;
    /*
       TODO test and use this optimized version.
    int max = pix->nsamples;
    int min = 0;
    int i;

    while (min < max)
    {
        i = (min + max) / 2;
        if (cmp_samples(&pivot, &pix->samples[i]) < 0)
            max = i;
        else
            min = i + 1;
    }

    return max;
    */
}


/**
 * \brief print store values
 * \param store
 */
void
PixelStore_print(PixelStore* store)
{
    pixelAvlPrint((pixel_avl*) store->pixels);
}


/**
 * \brief compare two samples
 * \param a
 * \param b
 * \return 1 if a is newer than b, -1 if he is older, 0 if they comme from the same field.
 * \details
 * This compare uses the same function than the one used to sort samples.
 *
 * If epoch are equals, samples will be compared on field indexes.
 * If a and b come from the same field, they will the be considered equal.
 */
int
PixelStore_compare(samplestruct *a, samplestruct *b)
{
    return cmp_samples(&a, &b);
}



/* print */
static void pixelAvlPrint(pixel_avl *pix) {
    if (pix == NULL)
        return;
    pixelAvlPrint(pix->pAfter);
    pixelAvlPrint(pix->pBefore);

    fprintf(stderr, "pixel num: %li\n", pix->pixel.id);
}

/* free */
static void pixelAvlFree(pixel_avl *pix) {
    if (pix == NULL)
        return;

    pixelAvlFree(pix->pAfter);
    pixelAvlFree(pix->pBefore);

    free(pix->pixel.samples);
    free(pix);
}

/* compare */
static int cmp_samples(const void* a, const void *b) {
    samplestruct *sa = * (samplestruct**) a;
    samplestruct *sb = * (samplestruct**) b;
    struct field *fa = sa->set->field;
    struct field *fb = sb->set->field;

    if (fa->epoch < fb->epoch)
        return -1;
    else if (fa->epoch > fb->epoch)
        return 1;
    else
        return fa->fieldindex > fb->fieldindex ? 1 :
                    fa->fieldindex < fb->fieldindex ? -1 : 0;
}

/* sort */
static void pixelAvlSort(pixel_avl *pix) {
    if (pix == NULL)
        return;
    pixelAvlSort(pix->pAfter);
    pixelAvlSort(pix->pBefore);
    qsort(pix->pixel.samples,
            pix->pixel.nsamples,
            sizeof(samplestruct*),
            cmp_samples);

}



/*****************************************************************************
 * 1 AVL Tree implementation
 *
 * From the SQLite source code ext/misc/amatch.c and ext/misc/closure.c with
 * the following notice:
 *
 * The author disclaims copyright to this source code.  In place of
 * a legal notice, here is a blessing:
 *
 * May you do good and not evil.
 * May you find forgiveness for yourself and forgive others.
 * May you share freely, never taking more than you give.
 */

static int pixel_cmp(long a, long b) {return a < b ? -1 : (a > b ? 1 : 0);}

/* Recompute the amatch_avl.height and amatch_avl.imbalance fields for p.
 ** Assume that the children of p have correct heights.
 */
static void amatchAvlRecomputeHeight(pixel_avl *p) {
    short int hBefore = p->pBefore ? p->pBefore->height : 0;
    short int hAfter = p->pAfter ? p->pAfter->height : 0;
    p->imbalance = hBefore - hAfter; /* -: pAfter higher.  +: pBefore higher */
    p->height = (hBefore > hAfter ? hBefore : hAfter) + 1;
}

/*
 **  P    B
 ** / \     / \
 **   B   Z ==>  X   P
 **  / \      / \
 ** X   Y    Y   Z
 **
 */
static pixel_avl *amatchAvlRotateBefore(pixel_avl *pP) {
    pixel_avl *pB = pP->pBefore;
    pixel_avl *pY = pB->pAfter;
    pB->pUp = pP->pUp;
    pB->pAfter = pP;
    pP->pUp = pB;
    pP->pBefore = pY;
    if (pY)
        pY->pUp = pP;
    amatchAvlRecomputeHeight(pP);
    amatchAvlRecomputeHeight(pB);
    return pB;
}

/*
 **  P    A
 ** / \     / \
 **   X   A ==>  P   Z
 **   / \    / \
 **  Y   Z  X   Y
 **
 */
static pixel_avl *amatchAvlRotateAfter(pixel_avl *pP) {
    pixel_avl *pA = pP->pAfter;
    pixel_avl *pY = pA->pBefore;
    pA->pUp = pP->pUp;
    pA->pBefore = pP;
    pP->pUp = pA;
    pP->pAfter = pY;
    if (pY)
        pY->pUp = pP;
    amatchAvlRecomputeHeight(pP);
    amatchAvlRecomputeHeight(pA);
    return pA;
}

/*
 ** Return a pointer to the pBefore or pAfter pointer in the parent
 ** of p that points to p.  Or if p is the root node, return pp.
 */
static pixel_avl **pixelAvrFromPtr(pixel_avl *p, pixel_avl **pp) {
    pixel_avl *pUp = p->pUp;
    if (pUp == 0)
        return pp;
    if (pUp->pAfter == p)
        return &pUp->pAfter;
    return &pUp->pBefore;
}

/*
 ** Rebalance all nodes starting with p and working up to the root.
 ** Return the new root.
 */
static pixel_avl *pixelAvlBalance(pixel_avl *p) {
    pixel_avl *pTop = p;
    pixel_avl **pp;
    while (p) {
        amatchAvlRecomputeHeight(p);
        if (p->imbalance >= 2) {
            pixel_avl *pB = p->pBefore;
            if (pB->imbalance < 0)
                p->pBefore = amatchAvlRotateAfter(pB);
            pp = pixelAvrFromPtr(p, &p);
            p = *pp = amatchAvlRotateBefore(p);
        } else if (p->imbalance <= (-2)) {
            pixel_avl *pA = p->pAfter;
            if (pA->imbalance > 0)
                p->pAfter = amatchAvlRotateBefore(pA);
            pp = pixelAvrFromPtr(p, &p);
            p = *pp = amatchAvlRotateAfter(p);
        }
        pTop = p;
        p = p->pUp;
    }
    return pTop;
}

/* Search the tree rooted at p for an entry with zId.  Return a pointer
 ** to the entry or return NULL.
 */
static pixel_avl *pixelAvlSearch(pixel_avl *p, const long zId) {
    int c;
    while (p && (c = pixel_cmp(zId, p->pixel.id)) != 0) {
        p = (c < 0) ? p->pBefore : p->pAfter;
    }
    return p;
}

/* Insert a new node pNew.  Return NULL on success.  If the key is not
 ** unique, then do not perform the insert but instead leave pNew unchanged
 ** and return a pointer to an existing node with the same key.
 */
static pixel_avl *pixelAvlInsert(pixel_avl **ppHead, pixel_avl *pNew) {
    int c;
    pixel_avl *p = *ppHead;
    if (p == 0) {
        p = pNew;
        pNew->pUp = 0;
    } else {
        while (p) {
            c = pixel_cmp(pNew->pixel.id, p->pixel.id);
            if (c < 0) {
                if (p->pBefore) {
                    p = p->pBefore;
                } else {
                    p->pBefore = pNew;
                    pNew->pUp = p;
                    break;
                }
            } else if (c > 0) {
                if (p->pAfter) {
                    p = p->pAfter;
                } else {
                    p->pAfter = pNew;
                    pNew->pUp = p;
                    break;
                }
            } else {
                return p;
            }
        }
    }
    pNew->pBefore = 0;
    pNew->pAfter = 0;
    pNew->height = 1;
    pNew->imbalance = 0;
    *ppHead = pixelAvlBalance(p);
    return 0;
}


#if 0 /* NOT USED */
/* Remove node pOld from the tree.  pOld must be an element of the tree or
 ** the AVL tree will become corrupt.
 */
static void amatchAvlRemove(amatch_avl **ppHead, amatch_avl *pOld){
    amatch_avl **ppParent;
    amatch_avl *pBalance = 0;
    /* assert( amatchAvlSearch(*ppHead, pOld->zKey)==pOld ); */
    ppParent = amatchAvlFromPtr(pOld, ppHead);
    if( pOld->pBefore==0 && pOld->pAfter==0 ){
        *ppParent = 0;
        pBalance = pOld->pUp;
    }else if( pOld->pBefore && pOld->pAfter ){
        amatch_avl *pX, *pY;
        pX = amatchAvlFirst(pOld->pAfter);
        *amatchAvlFromPtr(pX, 0) = pX->pAfter;
        if( pX->pAfter ) pX->pAfter->pUp = pX->pUp;
        pBalance = pX->pUp;
        pX->pAfter = pOld->pAfter;
        if( pX->pAfter ){
            pX->pAfter->pUp = pX;
        }else{
            assert( pBalance==pOld );
            pBalance = pX;
        }
        pX->pBefore = pY = pOld->pBefore;
        if( pY ) pY->pUp = pX;
        pX->pUp = pOld->pUp;
        *ppParent = pX;
    }else if( pOld->pBefore==0 ){
        *ppParent = pBalance = pOld->pAfter;
        pBalance->pUp = pOld->pUp;
    }else if( pOld->pAfter==0 ){
        *ppParent = pBalance = pOld->pBefore;
        pBalance->pUp = pOld->pUp;
    }
    *ppHead = amatchAvlBalance(pBalance);
    pOld->pUp = 0;
    pOld->pBefore = 0;
    pOld->pAfter = 0;
    /* assert( amatchAvlIntegrity(*ppHead) ); */
    /* assert( amatchAvlIntegrity2(*ppHead) ); */
}
#endif
