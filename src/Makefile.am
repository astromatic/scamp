#
#				Makefile.am
#
# src Makefile.am. Process this file with automake to generate a Makefile
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SCAMP
#
#	Copyright:		(C) 2002-2015 IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	SCAMP is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SCAMP is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		02/06/2015
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBDIRS		= fits wcs

if USE_PLPLOT
CPLOTSOURCE	= astrcplot.c cplot.c photcplot.c
endif

bin_PROGRAMS	= scamp
scamp_SOURCES	= astrefcat.c astrsolve.c astrstats.c catout.c check.c \
		  colour.c $(CPLOTSOURCE) crossid.c dgeomap.c fft.c fgroup.c \
		  field.c fitswcs.c header.c main.c makeit.c match.c merge.c \
		  misc.c mosaic.c photsolve.c pqueue.c prefs.c proper.c \
		  quadtree.c samples.c threads.c url.c xml.c \
		  astrefcat.h astrsolve.h astrstats.h cathead.h catout.h \
		  check.h colour.h cplot.h crossid.h define.h dgeomap.h fft.h \
		  fgroup.h field.h fitswcs.h globals.h header.h key.h match.h \
		  merge.h misc.h mosaic.h photsolve.h preflist.h prefs.h \
		  proper.h pqueue.h quadtree.h samples.h threads.h types.h \
		  url.h wcscelsys.h xml.h chealpix.h chealpix.c chealpixstore.h \
          chealpixstore.c crossid2.h crossid2.c
scamp_LDADD	= $(srcdir)/fits/libfits.a \
		  $(srcdir)/wcs/libwcs_c.a
DATE=`date +"%Y-%m-%d"`

