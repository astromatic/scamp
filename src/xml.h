/*
*				xml.h
*
* Include file for xml.c
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

#ifndef _XML_H_
#define _XML_H_

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include "field.h"
#include "fgroup.h"

/*----------------------------- Internal constants --------------------------*/
#ifndef XSL_URL
#define	XSL_URL	"."
#endif

/*--------------------------------- typedefs --------------------------------*/
/*------------------------------- functions ---------------------------------*/

extern int	write_xml(char *filename),
		write_xml_header(FILE *file),
		write_xml_meta(FILE *file, char *msgerror),
		write_xmlconfigparam(FILE *file, char *name, char *unit,
				char *ucd, char *format);
extern void	end_xml(void),
		init_xml(fieldstruct **fields, int nfield,
				fgroupstruct **fgroups, int ngroup),
		write_xmlerror(char *filename, char *msgerror);

#endif // _XML_H_
