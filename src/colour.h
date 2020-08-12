/*
*				colour.h
*
* Include file for colour.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2008-2020 IAP/CNRS/SorbonneU
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
*	Last modified:		12/08/2020
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/


#ifndef _COLOUR_H_
#define _COLOUR_H_

/*--------------------------------- constants -------------------------------*/
#define		PCA_NITER	200	/* Max nb of iter. in colour_findpc() */
#define		PCA_CONVEPS	1e-6	/* colour_findpc() converg. criterion */

/*--------------------------------- typedefs --------------------------------*/


/*--------------------------- structure definitions -------------------------*/


/*-------------------------------- protos -----------------------------------*/

extern double	colour_findpc(double *covmat, double *vec, int nmat);

extern void	colour_fgroup(fgroupstruct **fgroups, int ngroup);
#endif // _COLOUR_H_
