/*
 * Healpix pixels storage mechanism.
 *
 * This file is divided in four parts:
 * - 1 AVL tree implementations,
 * - 2 static functions
 * - 3 public interface
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

#include <stdlib.h>
#include <stdio.h>

#include "chealpixstore.h"
#include "chealpix.h"

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
typedef struct pixel_avl pixel_avl;
struct pixel_avl {
    HealPixel pixel; /* The samples being stored. Key is at pixel.id */
    pixel_avl *pBefore; /* Other elements less than pixel.id */
    pixel_avl *pAfter; /* Other elements greater than pixel.id */
    pixel_avl *pUp; /* Parent element */
    short int height; /* Height of this node.  Leaf==1 */
    short int imbalance; /* Height difference between pBefore and pAfter */
};

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
pixel_avl *pixelAvlSearch(pixel_avl *p, const long zId) {
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
pixel_avl *pixelAvlInsert(pixel_avl **ppHead, pixel_avl *pNew) {
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

void pixelAvlFree(pixel_avl *pix) {
    if (pix == NULL)
        return;

    pixelAvlFree(pix->pAfter);
    pixelAvlFree(pix->pBefore);

    free(pix->pixel.samples);
    free(pix);
}

void pixelAvlLink(pixel_avl *root, pixel_avl *leaf) {
    if (leaf->pAfter)
        pixelAvlLink(root, leaf->pAfter);
    if (leaf->pBefore)
        pixelAvlLink(root, leaf->pBefore);

    int i;
    for (i=0; i<8; i++) {
        pixel_avl *f = pixelAvlSearch(root, leaf->pixel.neighbors[i]);
        leaf->pixel.pneighbors[i] = &f->pixel;
    }
}

void fix_pixel_neighbors(pixel_avl *root) {
    pixelAvlLink(root, root);
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
/**
 * AVL related functions end
 ******************************************************************************/



/******************************************************************************
 * 2 PRIVATE FUNCTIONS
 */
#define SPL_BASE_SIZE 50
    static void
insert_sample_into_avltree_store(
        PixelStore *store, 
        struct sample spl, 
        struct sample **ext) 
{

    /* search for the pixel */
    /*
    pixel_avl *avlpix =
        pixelAvlSearch((pixel_avl*) store->pixels, spl.pix_nest);

    if (!avlpix) { // no such pixel

        int i;

        avlpix = CALLOC(1, sizeof(pixel_avl));
        avlpix->pixel.id = spl.pix_nest;
        avlpix->pixel.samples = CALLOC(SPL_BASE_SIZE, sizeof(struct sample));
        avlpix->pixel.ext = CALLOC(SPL_BASE_SIZE, sizeof(struct sample***));
        avlpix->pixel.nsamples = 0;
        avlpix->pixel.size = SPL_BASE_SIZE;
        pthread_mutex_init(&avlpix->pixel.mutex, NULL);

        for (i=0;i<8;i++)
            avlpix->pixel.tneighbors[i] = false;
        neighbours_nest64(store->nsides, spl.pix_nest, avlpix->pixel.neighbors);

        pixelAvlInsert((pixel_avl**) &store->pixels, avlpix);

        if (store->pixelids_size == store->npixels) {
            store->pixelids = REALLOC(store->pixelids, 
                    sizeof(long) * store->pixelids_size * 2);
            store->pixelids_size *= 2;
        }
        store->pixelids[store->npixels] = spl.pix_nest;
        store->npixels++;

    }

    HealPixel *pix = &avlpix->pixel;
    if (pix->nsamples == pix->size) {
        pix->samples = REALLOC(pix->samples, sizeof(struct sample) * pix->size * 2);
        pix->ext  = REALLOC(pix->ext, sizeof(struct sample***)  * pix->size * 2);
        int i;
        for (i=0; i<pix->size; i++) {
            struct sample **ext2  =  pix->ext[i];
            *ext2    = &pix->samples[i];
        }
        pix->size *= 2;
    }

    pix->samples[pix->nsamples] = spl;
    *ext = &pix->samples[pix->nsamples];
    pix->ext[pix->nsamples] = ext;
    //printf("%li %li\n", pix->samples[pix->nsamples].id, (*pix->ext[pix->nsamples])->id);
    pix->nsamples++;
        */

}

#define PIXELIDS_BASE_SIZE 1000
static PixelStore*
new_store(int64_t nsides) {

    PixelStore *store;
    QMALLOC(store, PixelStore, 1);

    store->pixels = NULL;
    store->nsides = nsides;
    store->npixels = 0;
    QMALLOC(store->pixelids, int64_t, PIXELIDS_BASE_SIZE);
    store->pixelids_size = PIXELIDS_BASE_SIZE;

    return store;
}


/**
 * PRIVATE FUNCTIONS END
 ******************************************************************************/



/******************************************************************************
 * PUBLIC FUNCTIONS
 */


/*
   PixelStore*
   PixelStore_new(Field *fields, int nfields, int64_t nsides) {

   PixelStore *store = new_store();

   Field field;
   Set set;
   Sample *spl;

   int i, j, k;
   for (i = 0; i < nfields; i++) {
   field = fields[i];

   for (j = 0; j < field.nsets; j++) {
   set = field.sets[j];

   for (k = 0; k < set.nsamples; k++) {

   spl = &set.samples[k];
   spl->bestMatch = NULL;

   ang2pix_nest64(nsides,spl->col, spl->lon, &spl->pix_nest);
   ang2vec(spl->col, spl->lon, spl->vector);
   insert_sample_into_avltree_store(store, spl, nsides);

   }
   }
   }

   fix_pixel_neighbors(store->pixels);
   return store;
   }
 */

    PixelStore*
PixelStore_new(int64_t nsides) 
{
    return new_store(nsides);
}

    void
PixelStore_add(
        PixelStore  *store, 
        struct sample spl, 
        struct sample **ext)
{
    spl.prevsamp = spl.nextsamp = NULL;
    /* XXX TODO get colatitude and longitude for healpix */
    /*
       ang2pix_nest64(store->nsides, spl.col, spl.lon, &spl.pix_nest);
       ang2vec(spl.col, spl.lon, spl.vector);
     */
    insert_sample_into_avltree_store(store, spl, ext);
}


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

    void
pixelAvlSetMaxRadius(
        pixel_avl *leaf, 
        double   radius) 
{
    if (leaf->pAfter != NULL)
        pixelAvlSetMaxRadius(leaf->pAfter, radius);
    if (leaf->pBefore != NULL)
        pixelAvlSetMaxRadius(leaf->pBefore, radius);

    HealPixel *p = &leaf->pixel;
    int i;
    /*
    for (i=0; i<p->nsamples; i++)
        (&p->samples[i])->bestMatchDistance = radius;
        */

}


    void
PixelStore_setMaxRadius(
        PixelStore *store, 
        double   radius) 
{
    pixel_avl *root = store->pixels;

    /*
     * get the euclidean distance for this radius
     */
    double va[3], vb[3], euclidean_dist;
    ang2vec(0,0,va);
    ang2vec(radius,0, vb);
    euclidean_dist = euclidean_distance(va,vb);

    pixelAvlSetMaxRadius(root, euclidean_dist);

}


    void
PixelStore_free(PixelStore* store) 
{
    pixelAvlFree((pixel_avl*) store->pixels);
    free(store->pixelids);
    free(store);

}
/**
 * PUBLIC FUNCTIONS END
 ******************************************************************************/

