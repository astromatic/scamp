/**
 *
 * \file        crossid.h
 * \author      Emmanuel Bertin
 * \author      Sébastien Serre
 * \date        7/05/2018
 *
 * \copyright   Copyright (C) 2017 University of Bordeaux. All right reserved.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 */

#ifndef __CROSSID2_H__
#define __CROSSID2_H__

#include "fgroup.h"
#include "field.h"

extern void
CrossId_run(fgroupstruct *fgroup,  fieldstruct *reffield, double radius);


#endif /* __CROSSID2_H__ */
