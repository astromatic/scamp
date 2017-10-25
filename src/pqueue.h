/**
* @file		pqueue.h
* @brief	Include file for pqueue.c.
* @date		27/11/2014
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2014 IAP/CNRS/UPMC
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

#ifndef _PQUEUE_H_
#define _PQUEUE_H_

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif


//-------------------------------- constants ----------------------------------

#define	PQUEUE_MEMINC	64		/// Memory allocation increment

//--------------------------------- typedefs ----------------------------------

typedef struct pnode {
   void		*element;		/// Pointer to element;
   int		p;			/// Priority index (p-- when priority++)
}	pnodestruct;

typedef struct pqueue {
   struct pnode		*node;		/// Pointer to array of nodes
   int			nnode;		/// Number of active nodes
   int			nnodemax;	/// Number of allocated nodes
}	pqueuestruct;

//------------------------------ Prototypes -----------------------------------

pqueuestruct	*pqueue_new(void);

void		*pqueue_pullelement(pqueuestruct *queue);

int		pqueue_insertelement(pqueuestruct *queue,
				void *element, int priority);

void		pqueue_end(pqueuestruct *queue);

#endif // _PQUEUE_H_
