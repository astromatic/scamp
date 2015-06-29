/**
* @file		pqueue.c
* @brief	Manage priority queues
* @date		02/12/2014
* @copyright
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 1993-2014 IAP/CNRS/UPMC
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

#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "globals.h"
#include "pqueue.h"


/****** pqueue_insertelement **********************************************//**
Insert an element in a priority queue
@param[in] queue	Pointer to the priority queue
@param[in] element	Pointer to the element to be inserted
@param[in] priority	Priority index (the lower, the higher the priority)
@param[out] 		RETURN_OK if everything went fine,
			RETURN_ERROR otherwise.
@author 		E. Bertin (IAP)
@date			02/12/2014
 ***/
int	pqueue_insertelement(pqueuestruct *queue, void *element, int priority) {
   pnodestruct	*node, *cnode, *fnode,
		tnode;
   int		pos;

  if (queue->nnode >= queue->nnodemax) {	// (Re-)allocate mem if needed
    queue->nnodemax += PQUEUE_MEMINC;
    if (queue->node) {
      if (!(queue->node = (pnodestruct *)realloc((void *)queue->node,
				queue->nnodemax * sizeof(pnodestruct))))
        return RETURN_ERROR;
    }
    else if (!(queue->node = (pnodestruct *)malloc(queue->nnodemax
				* sizeof(pnodestruct))))
      return RETURN_ERROR;

  }

  pos = queue->nnode++;
  node = queue->node;
  cnode = node + pos;
  cnode->element = element;		// Add new element at the last position
  cnode->p = priority;
// Navigate up through the tree ...
// ... swapping node content until priority index exceeds that of parent node
  while (pos				// test and update level
	&& cnode->p < (fnode = node + (pos = (pos-1)/2))->p) {
    tnode = *fnode;			// Swap node content
    *fnode = *cnode;
    *cnode = tnode;
    cnode = fnode;			// Up one level
  }

  return RETURN_OK;
}


/****** pqueue_pullelement ************************************************//**
Retrieve the element with highest priority (lowest value)
@param[in] queue	Pointer to the priority queue
@param[out] 		Pointer to element with highest priority.
@author 		E. Bertin (IAP)
@date			27/11/2014
 ***/
void	*pqueue_pullelement(pqueuestruct *queue) {
   pnodestruct	*node, *cnode, *fnode,
		tnode;
   void		*element;
   int		pos, nnode;

  node = queue->node;
  if (queue->nnode)
    element = node[0].element;				// Store element
  else
    return NULL;					// Exit if no element

  nnode = --queue->nnode;				// New number of nodes
  node[pos = 0] = node[nnode];				// Move last elem to top
  cnode = node;
  while ((pos = pos * 2 + 1) < nnode) {			// Test and update level
    fnode = node + pos;
    if ((pos + 1) < nnode && (fnode + 1)->p < fnode->p) {
      fnode++;			// Select child with the lowest priority index
      pos++;
    }
    if (cnode->p < fnode->p)		// Exit if priority index low enough
      break;
    tnode = *fnode;					// Swap node content
    *fnode = *cnode;
    *cnode = tnode;
    cnode = fnode;					// Down one level
  }

  return element;
}


/****** pqueue_new *******************************************************//**
Initialize a new priority queue
@param[out] 		Pointer to the new priority queue, or NULL if memory
			allocation failed.
@author 		E. Bertin (IAP)
@date			26/11/2014
 ***/
pqueuestruct	*pqueue_new(void) {


  return (pqueuestruct *)calloc((size_t)1, sizeof(pqueuestruct));
}


/****** pqueue_end *******************************************************//**
Remove a node in a (quad-)tree and return the parent node
@param[in] 		Pointer to the queue
@author 		E. Bertin (IAP)
@date			26/11/2014
 ***/
void	pqueue_end(pqueuestruct *queue) {

  free(queue->node);
  free(queue);

  return;
}


