/**
 * test crossid
 *
 */
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../src/crossid.h"
#include "../src/field.h"
#include "../src/prefs.h"

prefstruct	prefs;
char		gstr[MAXCHAR];

static char testfile[] = "extra/744331p.ascii";

static int
realloc_samples_x(setstruct *set, int splsize)
{
    set->sample = realloc(set->sample, sizeof(samplestruct) * splsize * 2);
    return splsize * 2;
}

static fieldstruct*
load_ascii_field(int id)
{
    fieldstruct *field = malloc(sizeof(fieldstruct));
    field->lng = 0;
    field->lat = 1;
    field->fieldindex = id;
    setstruct *set = malloc(sizeof(setstruct));
    int splsize = 1000;
    set->setindex = 0;
    set->nsample = 0;
    set->sample = malloc(sizeof(samplestruct) * splsize);
    field->nset = 1;
    field->set = malloc(sizeof(void*));
    field->set[0] = set;
    set->field = field;

    FILE *cat = fopen(testfile, "r");
    int index;
    double x, y;
    while ((fscanf(cat, "%i %lf %lf\n", &index, &x, &y) == 3)) {
        if (splsize == set->nsample)
            splsize = realloc_samples_x(set, splsize);
        samplestruct *spl = &set->sample[set->nsample++];
        spl->wcspos[0] = x;
        spl->wcspos[1] = y;
        spl->set = set;
        spl->epoch = 1.0;
    }
    return field;

}

void
free_ascii_field(fieldstruct *field)
{
    free(field->set[0]->sample);
    free(field->set[0]);
    free(field->set);
    free(field);
}

void
compare_fields(int n, fieldstruct **fields)
{
    fieldstruct *a = fields[0];
    fieldstruct *b = fields[1];

    int i;
    samplestruct **linkedspl = malloc(sizeof(void*) * n);
    for (i=0; i<fields[0]->set[0]->nsample; i++) {
        int j;
        for (j=0; j<n; j++) {
            linkedspl[j] = &fields[j]->set[0]->sample[i];
        }
        for (j=1; j<n; j++) {
            assert(abs(linkedspl[0]->wcspos[0] - linkedspl[j]->wcspos[0]) < FLT_EPSILON);
            assert(abs(linkedspl[0]->wcspos[1] - linkedspl[j]->wcspos[1]) < FLT_EPSILON);
        }

        assert(linkedspl[0]->prevsamp == NULL);
        assert(linkedspl[n-1]->nextsamp == NULL);

        samplestruct *prev = linkedspl[0];
        for (j=1; j<n; j++) {
            samplestruct *next = linkedspl[j];
            assert(prev->nextsamp == next);
            assert(next->prevsamp == prev);
            prev = next;
        }
    }
}

void
move_samples(char axe, double arcsec, fieldstruct *field)
{
    int i;
    double deg = arcsec / 3600;
    for (i=0; i<field->set[0]->nsample; i++) {
        samplestruct *sp = &field->set[0]->sample[i];
        if (axe == 'x') {
            sp->wcspos[0] += deg;
        } else {
            sp->wcspos[1] += deg;
        }
    }
}

int
main(int argc, char **argv)
{
    prefs.astrefcat = ASTREFCAT_NONE;
    fieldstruct *fields[10];
    fieldstruct **rfields;

    fgroupstruct fg;
    fg.field = fields;

    int i;
    for (i=0; i<10; i++)
        fields[i] = load_ascii_field(i);

    for (i=2; i<=10; i++) {
        fg.nfield = i;
        CrossId_run(&fg, NULL, 2.0 * ARCSEC/DEG);
        compare_fields(i, fields);
    }

    for (i=0; i<10; i++) {
        free_ascii_field(fields[i]);
    }


    exit(0);
}
