/*
*				cathead.h
*
* Merged and full catalogue headers.
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
*	Last modified:		22/07/2011
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _CATOUT_H_
#include "catout.h"
#endif

/* Merged output catalog fields */

mergedsamplestruct	refmergedsample;
keystruct		refmergedkey[] = {
  {"SOURCE_NUMBER", "Source index",
	&refmergedsample.sourceindex, H_INT, T_LONG,
	"%10d", "", "meta.number", ""},
  {"NPOS_TOTAL", "Total number of overlapping detections",
	&refmergedsample.npos_tot, H_INT, T_LONG,
	"%4d", "", "meta.number", ""},
  {"NPOS_OK", "Number of unsaturated, uncropped overlapping detections",
	&refmergedsample.npos_ok, H_INT, T_LONG,
	"%4d", "", "meta.number", ""},
  {"ALPHA_J2000", "Position along right ascension",
	&refmergedsample.wcspos[0], H_FLOAT, T_DOUBLE,
	"%13.9f", "deg", "pos.eq.ra;stat.mean", "deg"},
  {"DELTA_J2000", "Position along declination",
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
  {"DISPALPHA_J2000", "RMS dispersion of pos along right ascension",
	&refmergedsample.wcsposdisp[0], H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.stdev;stat.max;pos.errorEllipse;meta.main", "deg"},
  {"DISPDELTA_J2000", "RMS dispersion of pos along declination",
	&refmergedsample.wcsposdisp[1], H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.stdev;stat.min;pos.errorEllipse;meta.main", "deg"},
  {"PMALPHA_J2000", "Proper motion along right ascension",
	&refmergedsample.wcsprop[0], H_FLOAT, T_FLOAT,
	"%10e", "mas/yr", "pos.pm;pos.eq.ra;stat.fit", "mas/yr"},
  {"PMDELTA_J2000", "Proper motion along declination",
	&refmergedsample.wcsprop[1], H_FLOAT, T_FLOAT,
	"%10e", "mas/yr", "pos.pm;pos.eq.de;stat.fit", "mas/yr"},
  {"PMALPHAERR_J2000", "P.motion uncertainty along right ascension",
	&refmergedsample.wcsproperr[0], H_FLOAT, T_FLOAT,
	"%10e", "mas/yr", "stat.error;pos.pm;pos.eq.ra;stat.fit", "mas/yr"},
  {"PMDELTAERR_J2000", "P.motion uncertainty along declination",
	&refmergedsample.wcsproperr[1], H_FLOAT, T_FLOAT,
	"%10e", "mas/yr", "stat.error;pos.pm;pos.eq.de;stat.fit", "mas/yr"},
  {"PARALLAX_WORLD", "Trigonometric parallax",
	&refmergedsample.wcsparal, H_FLOAT, T_FLOAT,
	"%10e", "mas", "pos.parallax.trig;stat.fit", "mas"},
  {"PARALLAXERR_WORLD", "Trignonometric parallax uncertainty",
	&refmergedsample.wcsparalerr, H_FLOAT, T_FLOAT,
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
  {"COLOR", "Color index",
	&refmergedsample.colour, H_FLOAT, T_FLOAT,
	"%8.4f", "", "phot.color", "mag"},
  {"FLAGS_EXTRACTION", "Extraction flags",
	&refmergedsample.sexflags, H_INT, T_SHORT,
	"%3d", "", "meta.code.qual", ""},
  {"FLAGS_SCAMP", "Calibration flags",
	&refmergedsample.scampflags, H_INT, T_SHORT,
	"%3d", "", "meta.code.qual", ""},
  {""},
  };

/* Full output catalog fields */

fullsamplestruct	reffullsample;
keystruct		reffullkey[] = {
  {"SOURCE_NUMBER", "Source index",
	&reffullsample.sourceindex, H_INT, T_LONG,
	"%10d", "", "meta.number", ""},
   {"CATALOG_NUMBER", "File index",
	&reffullsample.fieldindex, H_INT, T_LONG,
	"%7d", "", "meta.number", ""},
   {"EXTENSION", "Extension index",
	&reffullsample.setindex, H_INT, T_SHORT,
	"%5d", "", "meta.number", ""},
   {"ASTR_INSTRUM", "Astrometric instrument index",
	&reffullsample.astrinstruindex, H_INT, T_SHORT,
	"%5d", "", "meta.number", ""},
   {"PHOT_INSTRUM", "Photometric instrument index",
	&reffullsample.photinstruindex, H_INT, T_SHORT,
	"%5d", "", "meta.number", ""},
   {"X_IMAGE", "Position along x image axis",
	&reffullsample.rawpos[0], H_FLOAT, T_DOUBLE,
	"%11.4f", "pixel", "pos.cartesian.x", "pix"},
   {"Y_IMAGE", "Position along y image axis",
	&reffullsample.rawpos[1], H_FLOAT, T_DOUBLE,
	"%11.4f", "pixel", "pos.cartesian.y", "pix"},
   {"ERRA_IMAGE", "RMS position error along major axis",
	&reffullsample.rawposerr[0], H_FLOAT, T_FLOAT,
	"%9.5f", "pixel", "stat.error;stat.max;pos.errorEllipse;meta.main",
	"pix"},
   {"ERRB_IMAGE", "RMS position error along minor axis",
	&reffullsample.rawposerr[1], H_FLOAT, T_FLOAT,
	"%9.5f", "pixel", "stat.error;stat.min;pos.errorEllipse;meta.main",
	"pix"},
   {"ERRTHETA_IMAGE", "Error ellipse pos. angle (CCW/world-x)",
	&reffullsample.wcspostheta, H_FLOAT, T_FLOAT,
	"%6.2f", "deg", "stat.error;stat.min;pos.errorEllipse;meta.main", "deg"},
   {"ALPHA_J2000", "Position along right ascension",
	&reffullsample.wcspos[0], H_FLOAT, T_DOUBLE,
	"%13.9f", "deg", "pos.eq.ra;stat.mean", "deg"},
   {"DELTA_J2000", "Position along declination",
	&reffullsample.wcspos[1], H_FLOAT, T_DOUBLE,
	"%13.9f", "deg", "pos.eq.de;stat.mean", "deg"},
   {"ERRA_WORLD", "RMS position error along major world axis",
	&reffullsample.wcsposerr[0], H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.error;stat.max;pos.errorEllipse;meta.main", "deg"},
   {"ERRB_WORLD", "RMS position error along minor world axis",
	&reffullsample.wcsposerr[1], H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.error;stat.min;pos.errorEllipse;meta.main", "deg"},
   {"ERRTHETA_WORLD", "Error ellipse pos. angle (CCW/world-x)",
	&reffullsample.wcspostheta, H_FLOAT, T_FLOAT,
	"%10e", "deg", "stat.error;stat.min;pos.errorEllipse;meta.main", "deg"},
  {"EPOCH", "Epoch",
	&reffullsample.epoch, H_FLOAT, T_FLOAT,
	"%15.10f", "yr", "time.epoch;stat.mean", "yr"},
  {"MAG", "Magnitude in the current band",
	&reffullsample.mag, H_FLOAT, T_FLOAT,
	"%9.5f", "mag", "phot.mag", "mag"},
  {"MAGERR", "RMS mag error estimate in the current band",
	&reffullsample.magerr, H_FLOAT, T_FLOAT,
	"%9.5f", "mag", "stat.error;phot.mag", "mag"},
  {"FLAGS_EXTRACTION", "Extraction flags",
	&reffullsample.sexflags, H_INT, T_SHORT,
	"%3d", "", "meta.code.qual", ""},
  {"FLAGS_SCAMP", "Calibration flags",
	&reffullsample.scampflags, H_INT, T_SHORT,
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

