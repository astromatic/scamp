/**
 *
 * \file        chealpixstore.h
 * \brief       Helpix pixels storarge
 * \author      Sébastien Serre
 * \date        11/04/2018
 *
 * \copyright   Copyright (C) 2017 University of Bordeaux. All right reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#ifndef _CHEALPIXSTORE_H_
#define _CHEALPIXSTORE_H_

#include <stdbool.h>

#include "define.h"
#include "fgroup.h"
#include "globals.h"
#include "samples.h"

#define SC_PI 3.141592653589793238462643383279502884197                            
#define SC_TWOPI 6.283185307179586476925286766559005768394                         
#define HALFPI 1.570796326794896619231321691639751442099                        
#define SC_INV_HALFPI 0.6366197723675813430755350534900574                         
#define TO_RAD 0.0174532925199432957692369076848861271344

typedef struct HealPixel HealPixel;
struct HealPixel {
    int64_t id;                 /* healpix id */
    samplestruct **samples;    /* our samples */
    int nsamples;               /* number of samples belonging to this pixel */
    int size;                   /* for reallocation if required */
    int64_t neighbors[8];       /* Neighbors indexes */
};

typedef struct PixelStore {
    int64_t     nsides;
    void        *pixels;

    /* These are used to iterate over pixels */
    long        npixels;
    int64_t     *pixelids;
    int         pixelids_size; /* PRIVATE, for re allocation if required */

} PixelStore;

extern void
PixelStore_new(int64_t nsides, PixelStore *ps);

extern void
PixelStore_free(PixelStore *store);

extern void
PixelStore_add(PixelStore *store, samplestruct *spl);

extern HealPixel*
PixelStore_get(PixelStore *store, int64_t key);

extern HealPixel*
PixelStore_getPixelFromSample(PixelStore *store, samplestruct *sample);

extern int
PixelStore_getHigherFields(HealPixel *pix, samplestruct *pivot);

extern int
PixelStore_getLowerFields(HealPixel *pix, samplestruct *pivot);

extern void
PixelStore_sort(PixelStore *store);

extern int
PixelStore_compare(samplestruct *a, samplestruct *b);

#endif /* _CHEALPIXSTORE_H_ */
