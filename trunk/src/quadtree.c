/**
* @file		tree.c
* @brief	Manage (quad-)trees
* @date		03/12/2014
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "define.h"
#include "globals.h"
#include "pqueue.h"
#include "quadtree.h"

/****** tree_insertleaf ***************************************************//**
Insert a leaf in a (quad-)tree
@param[in] tree		Pointer to the (quad-)tree
@param[in] coords	Pointer to x,y coordinates
@param[in] leaf		(void) pointer to the leaf to be inserted
@param[out] 		RETURN_OK if everything went fine,
			RETURN_ERROR otherwise.
@author 		E. Bertin (IAP)
@date			03/12/2014
 ***/
int	tree_insertleaf(treestruct *tree, double *coords, void *leaf) {

   treenodestruct	**pnode, **pnodeo,
			*node, *nodeo;
   void			*leafo;
   double		dval;
   unsigned int		index[TREE_NDIM],
			indexo[TREE_NDIM];
   int			d, depth, eflag;

  for (d=0; d<TREE_NDIM; d++) {
    dval = coords[d] * tree->indexscale[d] + tree->indexoffset[d];
    if (dval < 0.0 || dval > (double)TREE_MAXBINS)
      return RETURN_ERROR;
    index[d] = (int)dval;
  }

  depth = TREE_MAXDEPTH;
  leafo = NULL;
  pnode = &tree->node;
  node = NULL;
  while (depth) {
    if (!*pnode) {					// Node is empty
      if (!(*pnode = tree_newnode(node)))		// Create new node
        return RETURN_ERROR;
      node = *pnode;
      if (leafo)
        for (d=0; d<TREE_NDIM; d++)	// non-leaf new node index refers to ...
          node->index[d] = (index[d] >> depth) << depth;	// 1st corner
      else {
        node->leaf = leaf;
        for (d=0; d<TREE_NDIM; d++)
          node->index[d] = index[d];
        break;
      }
      
    } else
      node = *pnode;
    if (node->leaf) {			// Former leaf
      eflag = 0;

      for (d=0; d<TREE_NDIM; d++)
        eflag += (int)(index[d] == (indexo[d] = node->index[d]));
      if (eflag == TREE_NDIM)		// Exit if leaf with same coords exists
        return RETURN_ERROR;
      for (d=0; d<TREE_NDIM; d++)
        node->index[d] = (indexo[d] >> depth) << depth;// Truncate leaf index
      leafo = node->leaf;
      node->leaf = NULL;		// No longer a leaf node
    }
    depth--;
    pnode = node->child;
    for (d=0; d<TREE_NDIM; d++)		// Compute next node index
      pnode += ((index[d] >> depth) & 1) << d;
    if (leafo) {
      pnodeo = node->child;
      for (d=0; d<TREE_NDIM; d++)	// Compute next node index for old leaf
        pnodeo += ((indexo[d] >> depth) & 1) << d;
      if (pnode != pnodeo) {		// Child nodes can be used as leaves
        if (!(nodeo = *pnodeo = tree_newnode(node)))
          return RETURN_ERROR;
        nodeo->leaf = leafo;
        leafo = NULL;
        for (d=0; d<TREE_NDIM; d++)	// Copy old full leaf index
          nodeo->index[d] = indexo[d];
      }
    }
  }

  return RETURN_OK;
}


/****** tree_newnode ******************************************************//**
Initialize a new node in a (quad-)tree
@param[out] 		Pointer to the new node, or NULL if memory allocation
			failed.
@author 		E. Bertin (IAP)
@date			03/12/2014
 ***/
treenodestruct	*tree_newnode(treenodestruct *parentnode) {

   treenodestruct	*newnode;

  newnode = (treenodestruct *)calloc((size_t)1, sizeof(treenodestruct));
  newnode->parent = parentnode;
  newnode->depth = parentnode? parentnode->depth - 1 : TREE_MAXDEPTH;

  return newnode;
}


/****** tree_endnode ******************************************************//**
Remove a node in a (quad-)tree and return the parent node
@param[in] 		Pointer to the node
@param[out] 		Pointer to the parent node
@author 		E. Bertin (IAP)
@date			20/11/2014
 ***/
treenodestruct	*tree_endnode(treenodestruct *node) {

   treenodestruct *parent;

  parent = node->parent;
  free(node);

  return parent;
}


/****** tree_newtree ******************************************************//**
Initialize a new (quad-)tree
@param[in] 		Pointer to an array of minimum coordinates
@param[in] 		Pointer to an array of maximum coordinates
@param[out] 		Pointer to the new tree, or NULL if memory allocation
			failed.
@author 		E. Bertin (IAP)
@date			24/11/2014
 ***/
treestruct	*tree_newtree(double *min, double *max) {

   treestruct	*tree;
   double	cdiff, cdiffmax;
   int		d;

  if (!(tree = (treestruct *)calloc((size_t)1, sizeof(treestruct))))
    return NULL;

  cdiffmax = 0.0;
  for (d=0; d<TREE_NDIM; d++) {
    if ((cdiff = max[d] - min[d]) <= TINY)	// Improper coordinate range
      return NULL;
    if (cdiff > cdiffmax)
      cdiffmax = cdiff;
  }

  for (d=0; d<TREE_NDIM; d++)
    tree->indexoffset[d] = - min[d]
			* (tree->indexscale[d] = TREE_MAXBINS / cdiffmax);

   return tree;
}


/****** tree_endtree ******************************************************//**
End a (quad-)tree
@param[in] 		Pointer to the tree
@author 		E. Bertin (IAP)
@date			02/12/2014
 ***/
void	tree_endtree(treestruct *tree) {

   treenodestruct	*node, *nodeo;
   int			nchild[TREE_MAXDEPTH] = {TREE_NCHILD},
			cdepth = 0;

  nodeo = node = tree->node;
  while (cdepth >= 0) {				// Navigate the tree
    if (!nchild[cdepth]--
	|| !node || node->leaf) {		// Last child or leaf node
      node = node? tree_endnode(node) : nodeo;	// Free that node
      cdepth--;					// Up one level
    } else {
      nodeo = node;				// Backup current node pointer
      node = node->child[nchild[cdepth]];	// Goto next child
      nchild[++cdepth] = TREE_NCHILD;		// Reinit. child counter below
    }
  }

  free(tree);					// The tree structure itself

  return;
}


/****** tree_knn **********************************************************//**
Return the k nearest neighbours
@param[in] tree		Pointer to the (quad-)tree
@param[in] coords	Pointer to x,y coordinates
@param[in] k		Number of neighbours
@param[in] knodes	Pointer to an allocated array of at least k pointers
			that will point to the nearest leaf elements.
@param[out]		Actual number of retrieved leaf elements
@author 		E. Bertin (IAP)
@date			16/12/2014
 ***/
int	tree_knn(treestruct *tree, double *coords, int k, void **kleaf) {

   pqueuestruct		*queue;
   treenodestruct	*child, *node;
   double		dval;
   unsigned int		index[TREE_NDIM];
   int			c, d, kk, depth;

  if (!tree->node)				// Return 0 for an empty tree
    return 0;

  for (d=0; d<TREE_NDIM; d++) {
    dval = coords[d] * tree->indexscale[d] + tree->indexoffset[d];
    if (dval < 0.0 || dval > (double)TREE_MAXBINS)
      return RETURN_ERROR;
    index[d] = (int)dval;
  }

  queue = pqueue_new();				// Initialize priority queue
  pqueue_insertelement(queue, tree->node, tree_distance(tree->node, index));
  kk = 0;
  while (node = pqueue_pullelement(queue))  {	// Get closest node
    if (node->leaf) {				// If the node is a leaf ...
      kleaf[kk++] = node->leaf;			// ... it is a closest neighbour
      if (kk == k)				// exit if kth closest neighbour
        break;
    } else
      for (c=0; c<TREE_NCHILD; c++)
        if ((child = node->child[c]))
          pqueue_insertelement(queue, child, tree_distance(child, index));
  }

  pqueue_end(queue);				// Terminate priority queue

  return kk;
}


/****** tree_distance *****************************************************//**
Return Euclidean distance between a node at depth depth and a coordinate index
@param[in] node		Pointer to the (quad-)tree node
@param[in] index	Pointer to the coordinate index
@@param[out]		Euclidean distance (truncated to integer)
@author 		E. Bertin (IAP)
@date			14/01/2015
 ***/
int	tree_distance(treenodestruct *node, unsigned int *index) {

   double		dist, dinc, dval;
   int			d;

  dist = 0.0;
  if (node->leaf)			// Leaf node: we use the full leaf index
    for (d=0; d<TREE_NDIM; d++) {
      dval = (double)index[d] - (double)node->index[d];
      dist += dval * dval;
    }
  else {				// Non-leaf node: first corner index
    dinc = (double)(1 << node->depth);
    for (d=0; d<TREE_NDIM; d++) {
      dval = (double)index[d] - (double)node->index[d];
      if (dval >= 0.0) {
        if (dval < dinc)
          dval = 0.0;
        else if (dval >= dinc)		// Find the closest border
          dval -= dinc;
      }
      dist += dval*dval;
    }
  }

  return (int)sqrt(dist);
}
