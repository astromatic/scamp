/*
*				header.h
*
* Include file for header.c.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
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
*	Last modified:		10/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifndef _HEADER_H_
#define _HEADER_H_

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include "fits/fitscat.h"
#include "field.h"
#include "fitswcs.h"

/*----------------------------- Internal constants --------------------------*/
/*--------------------------------- typedefs --------------------------------*/
/*------------------------------- functions ---------------------------------*/

extern int		read_aschead(char *filename, int frameno,
				tabstruct *tab),
			write_aschead(char *filename, fieldstruct *field);
#endif // _HEADER_H_
