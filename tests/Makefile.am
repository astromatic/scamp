# Test Makefile for SCAMP
AUTOMAKE_OPTIONS = parallel-tests
check_PROGRAMS = \
                 test_chealpix \
                 test_chealpixstore \
                 test_crossid_single_catalog \
                 test_crossid_single_catalog_moving
TESTS		= test_chealpix test_chealpixstore test_crossid_single_catalog test_crossid_single_catalog_moving test_integration01.py
EXTRA_DIST	= test_integration01.conf extra/megacam.ahead extra/GAIA-DR1_1334+3754_r46.cat extra/7*p.cat

ALL_SCAMP_SOURCES = $(top_builddir)/src/*c $(top_builddir)/src/*h

test_chealpix_SOURCES = test_chealpix.c $(ALL_SCAMP_SOURCES)
test_chealpix_LDADD = $(top_builddir)/src/chealpix.o 

test_chealpixstore_SOURCES = test_chealpixstore.c $(ALL_SCAMP_SOURCES)
test_chealpixstore_LDADD = $(top_builddir)/src/chealpix.o $(top_builddir)/src/chealpixstore.o

test_crossid_single_catalog_SOURCES = test_crossid_single_catalog.c $(ALL_SCAMP_SOURCES)
test_crossid_single_catalog_LDADD = $(top_builddir)/src/chealpix.o \
                     $(top_builddir)/src/chealpixstore.o \
                     $(top_builddir)/src/crossid.o \
                     $(top_builddir)/src/field.o \
                     $(top_builddir)/src/threads.o \
                     $(top_builddir)/src/header.o \
                     $(top_builddir)/src/fitswcs.o \
                     $(top_builddir)/src/samples.o \
                     $(top_builddir)/src/wcs/libwcs_c.a \
                     $(top_builddir)/src/fits/libfits.a 

test_crossid_single_catalog_moving_SOURCES = test_crossid_single_catalog_moving.c $(ALL_SCAMP_SOURCES)
test_crossid_single_catalog_moving_LDADD = $(top_builddir)/src/chealpix.o \
                     $(top_builddir)/src/chealpixstore.o \
                     $(top_builddir)/src/crossid.o \
                     $(top_builddir)/src/field.o \
                     $(top_builddir)/src/threads.o \
                     $(top_builddir)/src/header.o \
                     $(top_builddir)/src/fitswcs.o \
                     $(top_builddir)/src/samples.o \
                     $(top_builddir)/src/wcs/libwcs_c.a \
                     $(top_builddir)/src/fits/libfits.a 


distclean-local:
		-rm *.head *.png *.xml
