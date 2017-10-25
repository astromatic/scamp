/**
* @file		quadtree.h
* @brief	Include file for quadtree.c.
* @date		10/12/2014
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

#ifndef _QUADTREE_H_
#define _QUADTREE_H_

#ifdef HAVE_CONFIG_H
#include        "config.h"
#endif

//-------------------------------- constants ----------------------------------

#define	TREE_NDIM	2		/// Dimension of the tree (e.g. 2D, 3D)
#define	TREE_NCHILD	(1<<TREE_NDIM)	/// Number of childs per node
#define	TREE_MAXDEPTH	30		/// Maximum tree depth
#define	TREE_MAXBINS	(1<<TREE_MAXDEPTH)	/// Max number of tree bins
#define	TREE_MEMINC	1024		/// Memory allocation increment

//--------------------------------- typedefs ----------------------------------

typedef struct treenode {
   struct treenode	*parent;		/// Pointer to parent node
   struct treenode	*child[TREE_NCHILD];	/// Pointers to child nodes
   void			*leaf;			/// Pointer to leaf
   unsigned int		index[TREE_NDIM];	/// Current leaf index
   int			depth;			/// Depth in the tree
}	treenodestruct;

typedef struct tree {
   double	indexscale[TREE_NDIM];		/// Coordinate scaling factor
   double	indexoffset[TREE_NDIM];		/// Coordinate offset
   struct treenode	*node;			/// Pointer to top node
}	treestruct;

//------------------------------ Prototypes -----------------------------------

treestruct	*tree_newtree(double *min, double *max);

treenodestruct	*tree_endnode(treenodestruct *node),
		*tree_newnode(treenodestruct *parentnode);

int		tree_distance(treenodestruct *node, unsigned int *index),
		tree_insertleaf(treestruct *tree, double *coords, void *leaf),
		tree_knn(treestruct *tree, double *coords, int k, void **kleaf);

void		tree_endtree(treestruct *tree);

#endif // _QUADTREE_H_
