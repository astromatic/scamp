/*
*				cathead.h
*
* Merged catalogue headers.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SCAMP is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
* 	(at your option) any later version.
*	SCAMP is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		05/02/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _CATOUT_H_
#include "catout.h"
#endif

/* Merged output catalog fields */

mergedsamplestruct	refmergedsample;
keystruct		refmergedkey[] = {
  {"NPOS", "Number of overlapping positions",
	&refmergedsample.npos, H_INT, T_LONG,
	"%4d", "", "meta.number", ""},
  {"X_WORLD", "Barycenter position along world x axis",
	&refmergedsample.wcspos[0], H_FLOAT, T_DOUBLE,
	"%13.9f", "deg", "pos.eq.ra;stat.mean", "deg"},
  {"Y_WORLD", "Barycenter position along world y axis",
	&refmergedsample.wcspos[1], H_FLOAT, T_DOUBLE,
	"%13.9f", "deg", "pos.eq.de;stat.mean", "deg"},
  {"ERRA_WORLD", "RMS position error along major world axis",
	&refmergedsample.wcsposerr[0], H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.error;stat.max;pos.errorEllipse;meta.main", "deg"},
  {"ERRB_WORLD", "RMS position error along minor world axis",
	&refmergedsample.wcsposerr[1], H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.error;stat.min;pos.errorEllipse;meta.main", "deg"},
  {"ERRTHETA_WORLD", "Error ellipse pos. angle (CCW/world-x)",
	&refmergedsample.wcspostheta, H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.error;stat.min;pos.errorEllipse;meta.main", "deg"},
  {"DISPX_WORLD", "RMS dispersion of pos along x world axis",
	&refmergedsample.wcsposdisp[0], H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.stdev;stat.max;pos.errorEllipse;meta.main", "deg"},
  {"DISPY_WORLD", "RMS dispersion of pos along y world axis",
	&refmergedsample.wcsposdisp[1], H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.stdev;stat.min;pos.errorEllipse;meta.main", "deg"},
  {"PMX_WORLD", "Proper motion along world x axis",
	&refmergedsample.wcsprop[0], H_FLOAT, T_FLOAT,
	"%10e", "mas/yr", "pos.pm;pos.eq.ra;stat.fit", "mas/yr"},
  {"PMY_WORLD", "Proper motion along world y axis",
	&refmergedsample.wcsprop[1], H_FLOAT, T_FLOAT,
	"%10e", "mas/yr", "pos.pm;pos.eq.de;stat.fit", "mas/yr"},
  {"PMXERR_WORLD", "P.motion uncertainty along world x axis",
	&refmergedsample.wcsproperr[0], H_FLOAT, T_FLOAT,
	"%10e", "mas/yr", "stat.error;pos.pm;pos.eq.ra;stat.fit", "mas/yr"},
  {"PMYERR_WORLD", "P.motion uncertainty along world y axis",
	&refmergedsample.wcsproperr[1], H_FLOAT, T_FLOAT,
	"%10e", "mas/yr", "stat.error;pos.pm;pos.eq.de;stat.fit", "mas/yr"},
  {"PARALLAX_WORLD", "Trigonometric parallax",
	&refmergedsample.wcsparal, H_FLOAT, T_FLOAT,
	"%10e", "mas", "pos.parallax.trig;stat.fit", "mas"},
  {"PARALLAXERR_WORLD", "Trignonometric parallax uncertainty",
	&refmergedsample.wcsproperr, H_FLOAT, T_FLOAT,
	"%10e", "mas", "stat.error;pos.parallax.trig;stat.fit", "mas"},
  {"EPOCH", "Mean epoch",
	&refmergedsample.epoch, H_FLOAT, T_FLOAT,
	"%15.10f", "yr", "time.epoch;stat.mean", "yr"},
  {"EPOCH_MIN", "Minimum epoch",
	&refmergedsample.epochmin, H_FLOAT, T_FLOAT,
	"%15.10f", "yr", "time.epoch;stat.min", "yr"},
  {"EPOCH_MAX", "Maximum epoch",
	&refmergedsample.epochmax, H_FLOAT, T_FLOAT,
	"%15.10f", "yr", "time.epoch;stat.max", "yr"},
  {"NMAG", "Number of overlaps for each band",
	&refmergedsample.nmag, H_INT, T_LONG,
	"%4d", "", "meta.number", "",
	1, &refmergedsample.nband},
  {"MAG", "Magnitude for each band",
	&refmergedsample.mag, H_FLOAT, T_FLOAT,
	"%9.5f", "mag", "phot.mag", "mag",
	1, &refmergedsample.nband},
  {"MAGERR", "RMS mag error estimate for each band",
	&refmergedsample.magerr, H_FLOAT, T_FLOAT,
	"%9.5f", "mag", "stat.error;phot.mag", "mag",
	1, &refmergedsample.nband},
  {"MAG_DISP", "RMS mag dispersion for each band",
	&refmergedsample.magdisp, H_FLOAT, T_FLOAT,
	"%8.4f", "mag", "stat.stdev;phot.mag", "mag",
	1, &refmergedsample.nband},
  {"FLAGS", "Extraction flags",
	&refmergedsample.sexflags, H_INT, T_SHORT,
	"%3d", "", "meta.code.qual", ""},
  {""},
  };


/* Skycat header */

const char	skycathead[] = "QueryResult\n\n"
	"# Config entry for original catalog server:\n"
	"serv_type: catalog\n"
	"long_name: SExtractor catalog\n"
	"short_name: SExCat\n"
	"symbol: id diamond %4.1f\n"
	"# End config entry\n\n"
	"id\tra\tdec\tmag";

const char	skycattail[] = "[EOD]";

