/**
 * test chealpix store
 *
 * Avl tree implementation is safe (come from SQLite) , and bugs on it would 
 * crash all other tests.
 *
 * Here we test the two functions "getHigherFields" and "getLowerFields" with
 * the following samples:
 *
 * s0->field1->epoch 0
 * s1->field1->epoch 0
 * s2->field2->epoch 0
 * s3->field2->epoch 0
 * s4->field3->epoch 10
 * s5->field3->epoch 10
 *
 */
 

#include <assert.h>
#include <stdio.h>

#include "../src/samples.h"
#include "../src/chealpixstore.h"

char		gstr[MAXCHAR];

int main(int argc, char **argv)
{
    PixelStore ps;
    PixelStore_new(2, &ps);

    /* prepare fake field structures */
    fieldstruct f1;
    f1.epoch = 0;
    f1.fieldindex = 0;
    f1.isrefcat = 0;
    f1.lat = 0;
    f1.lng = 1;
    setstruct s1;
    s1.field = &f1;

    fieldstruct f2;
    f2.epoch = 0;
    f2.fieldindex = 1;
    f2.isrefcat = 0;
    f2.lat = 0;
    f2.lng = 1;
    setstruct s2;
    s2.field = &f2;

    fieldstruct f3;
    f3.epoch = 1;
    f3.fieldindex = 1;
    f3.isrefcat = 0;
    f3.lat = 0;
    f3.lng = 1;
    setstruct s3;
    s3.field = &f3;


    /* prepare our samples */
    samplestruct spls[6];
    samplestruct *sp;

    sp = &spls[0];
    sp->wcspos[0] = sp->wcspos[1] = 0.0;
    sp->set = &s1;
    sp->epoch = 0;
    PixelStore_add(&ps, sp);

    sp = &spls[1];
    sp->wcspos[0] = sp->wcspos[1] = 0.0;
    sp->set = &s1;
    sp->epoch = 0;
    PixelStore_add(&ps, sp);

    sp = &spls[2];
    sp->wcspos[0] = sp->wcspos[1] = 0.0;
    sp->set = &s2;
    sp->epoch = 0;
    PixelStore_add(&ps, sp);

    sp = &spls[3];
    sp->wcspos[0] = sp->wcspos[1] = 0.0;
    sp->set = &s2;
    sp->epoch = 0;
    PixelStore_add(&ps, sp);

    sp = &spls[4];
    sp->wcspos[0] = sp->wcspos[1] = 0.0;
    sp->set = &s3;
    sp->epoch = 10;
    PixelStore_add(&ps, sp);

    sp = &spls[5];
    sp->wcspos[0] = sp->wcspos[1] = 0.0;
    sp->set = &s3;
    sp->epoch = 10;
    PixelStore_add(&ps, sp);

    /* get the only pixel used here */
    HealPixel *pix = PixelStore_getPixelFromSample(&ps, &spls[0]);


    int index;

    /*
     * testing getHigerFields
     */

    /* getHigherFields for s0 must return index 2 */
    index = PixelStore_getHigherFields(pix, &spls[0]);
    assert(index == 2);

    /* getHigherFields for s1 must return index 2 */
    index = PixelStore_getHigherFields(pix, &spls[1]);
    assert(index == 2);

    /* getHigherFields for s2 must return index 4 */
    index = PixelStore_getHigherFields(pix, &spls[2]);
    assert(index == 4);

    /* getHigherFields for s3 must return index 4 */
    index = PixelStore_getHigherFields(pix, &spls[3]);
    assert(index == 4);

    /* getHigherFields for s4 must return index 6 */
    index = PixelStore_getHigherFields(pix, &spls[4]);
    assert(index == 6);

    /* getHigherFields for s5 must return index 6 */
    index = PixelStore_getHigherFields(pix, &spls[5]);
    assert(index == 6);




    /*
     * testing getLowerFields
     */

    /* getLowerFields for s0 must return -1 */
    index = PixelStore_getLowerFields(pix, &spls[0]);
    assert(index == -1);

    /* getLowerFields for s1 must return -1 */
    index = PixelStore_getLowerFields(pix, &spls[1]);
    assert(index == -1);

    /* getLowerFields for s2 must return 1 */
    index = PixelStore_getLowerFields(pix, &spls[2]);
    assert(index == 1);

    /* getLowerFields for s3 must return 1 */
    index = PixelStore_getLowerFields(pix, &spls[3]);
    assert(index == 1);

    /* getLowerFields for s4 must return 3 */
    index = PixelStore_getLowerFields(pix, &spls[4]);
    assert(index == 3);

    /* getLowerFields for s5 must return 3 */
    index = PixelStore_getLowerFields(pix, &spls[5]);
    assert(index == 3);

    /* end */
    PixelStore_free(&ps);

    return 0;
}
