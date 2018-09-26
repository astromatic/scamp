/*
 *               htmlout.h
 *
 * Generate scamp html output
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

#ifndef _HTMLOUT_H_
#define _HTMLOUT_H_

/****** HtmlOut_write *********************************************************
  PROTO void HtmlOut_write()
  PURPOSE Generate html output file
  OUTPUT -
  NOTES Requires that JsonOut_write as bean called, and the html template 
  readable.
  AUTHOR    E. Bertin (IAP)
  VERSION   13/03/2018
 ***/
void HtmlOut_write();

#endif /* _HTMLOUT_H_ */
