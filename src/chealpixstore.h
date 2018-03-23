/*
 * Healpix pixels storage mechanism.
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


#ifndef _CHEALPIXSTORE_H_
#define _CHEALPIXSTORE_H_

#include <stdbool.h>
#include <pthread.h>

#include "define.h"
#include "globals.h"
#include "samples.h"

#define SC_PI 3.141592653589793238462643383279502884197                            
#define SC_TWOPI 6.283185307179586476925286766559005768394                         
#define HALFPI 1.570796326794896619231321691639751442099                        
#define SC_INV_HALFPI 0.6366197723675813430755350534900574                         
#define TO_RAD 0.0174532925199432957692369076848861271344

typedef struct HealPixel HealPixel;
/**
 * HealPixel store pointers to every samples of all fields, belonging to a
 * common healpix pixel.
 */
struct HealPixel {

    int64_t id;            /* healpix id */
    struct sample **samples;    /* our samples */
    int nsamples;       /* number of samples belonging to this pixel */
    int size;           /* for reallocation if required */
    int64_t neighbors[8];  /* Neighbors indexes */
    HealPixel *pneighbors[8];
    bool tneighbors[8]; /* check if neighbors have allready been matched */
	pthread_mutex_t mutex;

};

typedef struct PixelStore {
    int64_t     nsides;
    void        *pixels; /* our opaque data */

    /* These are used to iterate over pixels */
    long        npixels;
    int64_t     *pixelids;
    int         pixelids_size; /* PRIVATE, for re allocation if required */

} PixelStore;


extern PixelStore*
PixelStore_new(int64_t nsides);

/* 
 * Store "spl" in "store" in pixel id "key". Set "ext to contains a pointer to
 * our newly created sample. In case of internal reallocation, ext is 
 * updated automatically.
 */
extern void
PixelStore_add(PixelStore *store, struct sample *spl);

extern HealPixel*
PixelStore_get(PixelStore *store, int64_t key);

extern void
PixelStore_free(PixelStore *store);

extern void
PixelStore_setMaxRadius(PixelStore *store, double radius);

extern void
PixelStore_updateSamplePos(PixelStore *store);

int
PixelStore_getHigherFields(HealPixel *pix, struct sample *pivot);

HealPixel*
PixelStore_getPixelFromSample(PixelStore *store, struct sample *sample);
extern void
PixelStore_sort(PixelStore *store);

extern int
PixelStore_compare(struct sample *a, struct sample *b);
#endif /* _CHEALPIXSTORE_H_ */
