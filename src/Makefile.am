# Program Makefile for SCAMP
# Copyright (C) 2002-2009 Emmanuel Bertin.
SUBDIRS		= fits wcs

if USE_PLPLOT
CPLOTSOURCE	= astrcplot.c cplot.c photcplot.c
endif

bin_PROGRAMS	= scamp
scamp_SOURCES	= astrefcat.c astrsolve.c catout.c check.c \
		  $(CPLOTSOURCE) crossid.c fft.c fgroup.c field.c fitswcs.c \
		  header.c main.c makeit.c match.c misc.c mosaic.c \
		  photsolve.c prefs.c proper.c samples.c threads.c xml.c \
		  astrefcat.h astrsolve.h cathead.h catout.h check.h cplot.h \
		  crossid.h define.h fft.h fgroup.h field.h fitswcs.h \
		  globals.h header.h key.h match.h misc.h mosaic.h \
		  photsolve.h preflist.h prefs.h proper.h samples.h threads.h \
		  types.h wcscelsys.h xml.h
scamp_LDADD	= $(top_builddir)/src/fits/libfits.a \
		  $(top_builddir)/src/wcs/libwcs_c.a
DATE=`date +"%Y-%m-%d"`

