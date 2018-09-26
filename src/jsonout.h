/*
 *               jsonout.h
 *
 * Generate scamp metadata
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 *
 *   This file part of:  SCAMP
 *
 *   Copyright:      (C) 2002-2018 Emmanuel Bertin -- IAP/CNRS/UPMC
 *
 *   License:        GNU General Public License
 *
 *   SCAMP is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *   SCAMP is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *   You should have received a copy of the GNU General Public License
 *   along with SCAMP. If not, see <http://www.gnu.org/licenses/>.
 *
 *   Last modified:      13/03/2018
 *
 *%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
#ifndef _JSONOUT_H_
#define _JSONOUT_H_

#include "fgroup.h"
#include "field.h"
#include <json.h>

/****** JsonOut_set_data ******************************************************
  PROTO void JsonOut_set_data(fieldstruct**,int,fgroupstruct**,int)
  PURPOSE Initialize fields and groups for json output.
  INPUT an array of fieldstruct*
  INPUT number of fields
  INPUT an array of fgroupstruct*
  INPUT number of fgroups
  OUTPUT
  NOTES -.
  AUTHOR    E. Bertin (IAP)
  VERSION   13/03/2018
 ***/
void JsonOut_set_data(fieldstruct**,int,fgroupstruct**,int);

/****** JsonOut_write *********************************************************
  PROTO void JsonOut_write()
  PURPOSE Generate json metadata file
  OUTPUT
  NOTES Requires a call to JsonOut_set_data
  AUTHOR    E. Bertin (IAP)
  VERSION   13/03/2018
 ***/
void JsonOut_write();

#endif /* _JSONOUT_H_ */
