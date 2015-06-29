/**
* @file		dgeomap.c
* @brief	Manage differential geometry maps
* @date		20/01/2015
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 1993-2015 IAP/CNRS/UPMC
*
*	License:		GNU General Public License
*
*	SCAMP is free software: you can redistribute it and/or modify
*	it under the terms of the GNU General Public License as published by
*	the Free Software Foundation, either version 3 of the License, or
*	(at your option) any later version.
*	SCAMP is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "define.h"
#include "globals.h"
#include "dgeomap.h"
#include "fgroup.h"
#include "field.h"
#include "fits/fitscat.h"
#include "fitswcs.h"
#include "misc.h"
#include "prefs.h"
#include "quadtree.h"
#include "samples.h"


/****** dgeomap_instru ****************************************************//**
Write differential geometry maps for a given astrometric instrument
@param[in] fields	Pointer to an array of field pointers
@param[in] nfield	Number of fields
@param[out] 		RETURN_OK if everything went fine,
			RETURN_ERROR otherwise.
@author 		E. Bertin (IAP)
@date			20/01/2015
 ***/
int	dgeomap_instru(fieldstruct **fields, int nfield, int instru,
			char *filename) {

   fieldstruct		*field, *field0;
   setstruct		*set, *set0;
   samplestruct		*samp, *samp2, *sampb, *sampn;
   wcsstruct		*wcs0;
   catstruct		*cat;
   tabstruct		*tab;
   dgeopointstruct	*dgeoneighbours[DGEOMAP_NNEIGHBOURMAX],
			*dgeopoint;
   treestruct		*quadtree;
   double		rawpos2[NAXIS],
			xmin[NAXIS], xmax[NAXIS],
			var, w, dgeox, dgeoy, dgeow, sdgeow, hsdgeow, dstepin;
   float		*line[2],
			*dgeoin, *dgeoout, *line0, *line1, *pix, *pix0, *pix1,
			*pixin, *pixout,
			val0, val1, stepout, coeff, coeff0;
   short		sexflagmask;
   unsigned int		imaflagmask;
   char			str[64];
   int			d,f,g,i,j,k,n,p,s, naxis, npixinx, npixiny, npixin,
			npixoutx, npixouty, npixout, nset,
			ngeopoint, ngeopointmax, nneighbour, stepin, nlineout;

// SExtractor and image flags
  sexflagmask = (short)prefs.astr_sexflagsmask;
  imaflagmask = prefs.astr_imaflagsmask;

  nset = 0;
  for (f=0; f<nfield; f++) {
    field = fields[f];
//-- Pick up a suitable field with the right astrometric instrument
    if (field->astromlabel == instru) {
      field0 = field;
      nset = field0->nset;
      set0 = field0->set[0];
      if (!set0->wcs || set0->wcs->naxis<2)
        return RETURN_ERROR;
      break;
    }
  }

  if (!nset)
    return RETURN_ERROR;

// Initialize FITS output
  cat = new_cat(1);
  init_cat(cat);
  sprintf(cat->filename, "%s", filename);
  if (open_cat(cat, WRITE_ONLY) != RETURN_OK)
    error(EXIT_FAILURE,"*Error*: cannot open for writing ", filename);
  if (nset>1) {
    addkeywordto_head(cat->tab, "NEXTEND ", "Number of extensions");
    fitswrite(cat->tab->headbuf, "NEXTEND ", &nset, H_INT, T_LONG);
    save_head(cat, cat->tab);
  }

// Compute the 2D astrometric residual histogram for every set
  for (s=0; s<nset; s++) {
    sprintf(str,
	"Generating differential geometry map for instrument A%d, set %d/%d",
	instru+1, s+1, nset);
    NFPRINTF(OUTPUT, str);
    set0 = field0->set[s];
    wcs0 = set0->wcs;
    if (!wcs0 || wcs0->naxis<2)
      return RETURN_ERROR;

    xmin[0] = xmin[1] = -16.0;
    xmax[0] = wcs0->naxisn[0] + 16.0;
    xmax[1] = wcs0->naxisn[1] + 16.0;

    dgeopoint = NULL;
    ngeopointmax = ngeopoint = 0;
    for (f=0; f<nfield; f++) {
      field = fields[f];
      if (field->astromlabel == instru) {	// Only for astrometric instru
        set = field->set[s];
        samp = set->sample;
        for (n=set->nsample; n--; samp++) {
          if ((samp->sexflags & sexflagmask)	// Filter out flagged detections
		|| (samp->imaflags & imaflagmask))
            continue;
          if (samp->rawpos[0] < xmin[0] || samp->rawpos[0] >= xmax[0]
		|| samp->rawpos[1] < xmin[1] || samp->rawpos[1] >= xmax[1])
            continue;
          dgeow = dgeox = dgeoy = 0.0; 
//-------- Explore forward and backward directions
          sampn = sampb = samp;
          while ((sampn && (samp2 = sampn = sampn->nextsamp))
		|| ((samp2 = sampb = sampb->prevsamp)
			&& sampb->set->field->astromlabel>=0)) {
            if ((samp2->sexflags & sexflagmask)
		|| (samp2->imaflags & imaflagmask))
              continue;
            var = samp->wcsposerr[0]*samp->wcsposerr[0]
		+ samp2->wcsposerr[0]*samp2->wcsposerr[0];
            wcs_to_raw(set->wcs, samp2->wcspos, rawpos2);
//---------- Weighted sum of contributions from detections
            if (var > TINY) {
              w = 1.0 / var;
              dgeow += w;
              dgeox += (samp->rawpos[0] - rawpos2[0]) * w;
              dgeoy += (samp->rawpos[1] - rawpos2[1]) * w;
            }
          }
          if (dgeow > TINY) {
            if (ngeopoint >= ngeopointmax) {
//------------ (Re-)allocate memory for storing data points if necessary
              ngeopointmax += DGEOMAP_MEMINC;
              if (dgeopoint)
                QREALLOC(dgeopoint, dgeopointstruct, ngeopointmax)
              else
                QMALLOC(dgeopoint, dgeopointstruct, ngeopointmax)
            }
//---------- Store position, average difference vector and weight
            dgeopoint[ngeopoint].pos[0] = samp->rawpos[0];
            dgeopoint[ngeopoint].pos[1] = samp->rawpos[1];
            dgeopoint[ngeopoint].dpos[0] = dgeox / dgeow;
            dgeopoint[ngeopoint].dpos[1] = dgeoy / dgeow;
            dgeopoint[ngeopoint].weight = dgeow;
            ngeopoint++;
          }
        }
      }
    }

//-- Setup both the low (input) and high (output) maps
    dstepin = (double)(stepin = prefs.dgeomap_step);
    npixinx = ((npixoutx = wcs0->naxisn[0]) - 1) / stepin + 1;
    npixiny = ((npixouty = wcs0->naxisn[1]) - 1) / stepin + 1;
    npixin = npixinx * npixiny;
    npixout = npixoutx * npixouty;

    QCALLOC(dgeoout, float, 2*npixout);
    if (prefs.dgeomap_step > 1) {
      QCALLOC(dgeoin, float, 2*npixin);
    } else
      dgeoin = dgeoout;

//-- Create quad-tree for accelerating the computation of nearest neighbours
    quadtree = tree_newtree(xmin, xmax);
//-- Fill-in quad-tree
    for (i=0; i<ngeopoint; i++)
      tree_insertleaf(quadtree, dgeopoint[i].pos, dgeopoint + i);
    rawpos2[1] = 1.0;
    for (j=0; j<npixiny; j++, rawpos2[1] += dstepin) {
      rawpos2[0] = 1.0;
      p = j * npixinx;
      for (i=npixinx; i--; p++, rawpos2[0] += dstepin) {
//------ Find the nearest neighbours at current position
        nneighbour = tree_knn(quadtree, rawpos2, prefs.dgeomap_nnearest,
			(void **)dgeoneighbours);
        if (!nneighbour)
          continue;
        qsort(dgeoneighbours, nneighbour, sizeof(dgeopointstruct *),
		dgeomap_compdx);
        sdgeow = dgeow = 0.0;
        for (n=0; n<nneighbour; n++)
          sdgeow += dgeoneighbours[n]->weight;
        hsdgeow = 0.5 * sdgeow;
        for (n=0; n<nneighbour; n++)
          if ((dgeow += dgeoneighbours[n]->weight) > hsdgeow) {
            break;
          }
        if (n == 0)
          dgeoin[p] = (float)dgeoneighbours[0]->dpos[0];
        else if (n == nneighbour)
          dgeoin[p] = (float)dgeoneighbours[n-1]->dpos[0];
        else if ((w = dgeoneighbours[n-1]->weight + dgeoneighbours[n]->weight)
		> TINY)
          dgeoin[p] = (float)dgeoneighbours[n-1]->dpos[0]
		+ (sdgeow - 2.0*(dgeow - dgeoneighbours[n]->weight)
		+ dgeoneighbours[n-1]->weight) / w
		* (dgeoneighbours[n]->dpos[0] - dgeoneighbours[n-1]->dpos[0]);

        qsort(dgeoneighbours, nneighbour, sizeof(dgeopointstruct *),
		dgeomap_compdy);
        dgeow = 0.0;
        for (n=0; n<nneighbour; n++)
          if ((dgeow += dgeoneighbours[n]->weight) > hsdgeow) {
            break;
          }
        if (n == 0)
          dgeoin[p + npixin] = (float)dgeoneighbours[0]->dpos[1];
        else if (n == nneighbour)
          dgeoin[p + npixin] = (float)dgeoneighbours[n-1]->dpos[1];
        else if ((w = dgeoneighbours[n-1]->weight + dgeoneighbours[n]->weight)
		> TINY)
          dgeoin[p + npixin] = (float)dgeoneighbours[n]->dpos[1]
		+ (sdgeow - 2.0*(dgeow - dgeoneighbours[n]->weight)
		+ dgeoneighbours[n-1]->weight) / w
		* (dgeoneighbours[n]->dpos[1] - dgeoneighbours[n-1]->dpos[1]);
      }
    }

    if (stepin > 1) {
//---- Fill output map by interpolating input map
      stepout = 1.0 / stepin;
      QMALLOC(line[0], float, stepin * npixinx);
      QMALLOC(line[1], float, stepin * npixinx);
      pixin = dgeoin;
      pixout = dgeoout;
      for (d=0; d<2; d++) {		// 2 dimensions
        nlineout = 0;
        for (j=0; j<npixiny;) {
//-------- First interpolate current line and the next
          if (!j) {
            pix = line0 = line[j%2];
            for (i=npixinx; i--;) {
              if (i) {
                val0 = *pixin;
                val1 = *(++pixin);
                coeff = 0.0;
              } else {
                coeff = 1.0;
                pixin++;
              }
              for (k = stepin; k--; coeff += stepout)
                *(pix++) = val0 * (1.0 - coeff) + val1 * coeff;
            }
          }
          if (++j < npixiny) {
            pix = line1 = line[j%2];
            for (i=npixinx; i--;) {
              if (i) {
                val0 = *pixin;
                val1 = *(++pixin);
                coeff = 0.0;
              } else {
                coeff = 1.0;
                pixin++;
              }
              for (k = stepin; k--; coeff += stepout)
                *(pix++) = val0 * (1.0 - coeff) + val1 * coeff;
            }
            line0 = line[(j+1)%2];
            coeff = 0.0;
          } else {
            coeff = 1.0;	// Interpolate using the penultimate line
          }
          for (k = stepin; k--; coeff += stepout) {
            if (++nlineout > npixouty)
              break;
            coeff0 = 1.0 - coeff;
            pix0 = line0;
            pix1 = line1;
            for (i=npixoutx; i--;)
              *(pixout++) = *(pix0++) * coeff0 + *(pix1++) * coeff;
          }
        }
      }
      free(line[0]);
      free(line[1]);
      free(dgeoin);
    }

//-- Save the current set as a FITS extension
    cat->tab->bodybuf = NULL;	// Prevent deallocating old body with init_cat()
    init_cat(cat);		// Cleanout previous tab
    tab = cat->tab;
    tab->bitpix =  BP_FLOAT;
    tab->bytepix = t_size[T_FLOAT];
    write_wcs(tab, wcs0);
    tab->naxis = 3;       /* This is an image */
    QREALLOC(tab->naxisn, int, tab->naxis);
    tab->naxisn[0] = npixoutx;
    tab->naxisn[1] = npixouty;
    tab->naxisn[2] = 2;
    tab->tabsize = tab->bytepix*tab->naxisn[0]*tab->naxisn[1]*tab->naxisn[2];
    tab->bodybuf = (void *)dgeoout;
    update_head(tab);
    if (nset>1)
      ext_head(tab);
    addkeywordto_head(cat->tab, "CTYPE3  ", "Coordinate type");
    fitswrite(tab->headbuf, "CTYPE3  ", "COMPLEX ", H_STRING, T_STRING);
    save_tab(cat, tab);

    tree_endtree(quadtree);
    free(dgeopoint);
    free(dgeoout);
  }

  tab->bodybuf = NULL;// prevent free_cheat() from freeing the dgeomap ptr
  free_cat(&cat,1);

  return RETURN_OK;
}


/****** dgeomap_compdx ****************************************************//**
qsort() comparison function for shifts in x
@param[in] dgeopoint1	Pointer to the first dgeopoint
@param[in] dgeopoint2	Pointer to the second dgeopoint
@param[out] 		<0 if dx1>dx2, >0 if dx1<dx2, 0 otherwise
@author 		E. Bertin (IAP)
@date			13/01/2015
 ***/
int	dgeomap_compdx(const void *dgeopoint1, const void *dgeopoint2) {

   double	dx;

  dx = (*(dgeopointstruct **)dgeopoint1)->dpos[0]
	- (*(dgeopointstruct **)dgeopoint2)->dpos[0];

  return dx>0.0 ? 1: (dx<0.0? -1 : 0);
}


/****** dgeomap_compdy ****************************************************//**
qsort() comparison function for shifts in y
@param[in] dgeopoint1	Pointer to the first dgeopoint
@param[in] dgeopoint2	Pointer to the second dgeopoint
@param[out] 		<0 if dx1>dx2, >0 if dx1<dx2, 0 otherwise
@author 		E. Bertin (IAP)
@date			13/01/2015
 ***/
int	dgeomap_compdy(const void *dgeopoint1, const void *dgeopoint2) {

   double	dy;

  dy = (*(dgeopointstruct **)dgeopoint1)->dpos[1]
	- (*(dgeopointstruct **)dgeopoint2)->dpos[1];

  return dy>0.0 ? 1: (dy<0.0? -1 : 0);
}


