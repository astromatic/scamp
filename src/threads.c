/*
*				threads.c
*
* Manage POSIX threads.
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	AstrOmatic software
*
*	Copyright:		(C) 2002-2010 Emmanuel Bertin -- IAP/CNRS/UPMC
*				(C) Mark Hays (original tutorial)
*
*	License:		GNU General Public License
*
*	AstrOmatic software is free software: you can redistribute it and/or
*	modify it under the terms of the GNU General Public License as
*	published by the Free Software Foundation, either version 3 of the
*	License, or (at your option) any later version.
*	AstrOmatic software is distributed in the hope that it will be useful,
*	but WITHOUT ANY WARRANTY; without even the implied warranty of
*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*	GNU General Public License for more details.
*	You should have received a copy of the GNU General Public License
*	along with AstrOmatic software.
*	If not, see <http://www.gnu.org/licenses/>.
*
*	Last modified:		13/10/2010
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#ifdef HAVE_CONFIG_H
#include	"config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_THREADS
#include "threads.h"
#endif

#include "define.h"
#include "types.h"
#include "globals.h"
#include "fits/fitscat.h"

#ifdef USE_THREADS

/******* threads_gate_init ***************************************************
PROTO	threads_gate_t *threads_gate_init(int nthreads, void (*func)(void))
PURPOSE	Create a new gate.
INPUT	Thread gate structure pointer,
	Void pointer to a function to exec prior to gate opening (can be NULL).
OUTPUT	-.
NOTES   From Mark Hays' POSIX tutorial.
	See e.g. http://www.cs.ualberta.ca/~paullu/C681/mark.hays.threads.html.
AUTHOR  E. Bertin (IAP)
VERSION 03/07/2003
 ***/
threads_gate_t *threads_gate_init(int nthreads, void (*func)(void))
  {
   threads_gate_t	*gate;

  QMALLOC(gate, threads_gate_t, 1);
  gate->ngate = 0;
  gate->nthreads = nthreads;
  gate->func = func;
  QPTHREAD_MUTEX_INIT(&gate->mutex, NULL);
  QPTHREAD_MUTEX_INIT(&gate->block, NULL);
  QPTHREAD_COND_INIT(&gate->condvar, NULL);
  QPTHREAD_COND_INIT(&gate->last, NULL);

  return gate;
  }


/******* threads_gate_end ****************************************************
PROTO	void threads_gate_end(threads_gate_t *gate)
PURPOSE	Destroy an existing gate.
INPUT	Thread gate structure pointer.
OUTPUT	-.
NOTES   From Mark Hays' POSIX tutorial.
	See e.g. http://www.cs.ualberta.ca/~paullu/C681/mark.hays.threads.html.
AUTHOR  E. Bertin (IAP)
VERSION 04/01/2010
 ***/
void threads_gate_end(threads_gate_t *gate)

  {
  gate->ngate=gate->nthreads=0;
  QPTHREAD_MUTEX_DESTROY(&gate->mutex);
  QPTHREAD_MUTEX_DESTROY(&gate->block);
  QPTHREAD_COND_DESTROY(&gate->condvar);
  QPTHREAD_COND_DESTROY(&gate->last);
  free(gate);

  return;
  }


/******* threads_gate_sync ****************************************************
PROTO	void threads_gate_sync(threads_gate_t *gate, void *func(void))
PURPOSE	Synchronize a number of POSIX threads.
INPUT	Thread gate structure pointer.
OUTPUT	-.
NOTES   From Mark Hays' POSIX tutorial.
	See e.g. http://www.cs.ualberta.ca/~paullu/C681/mark.hays.threads.html.
AUTHOR  E. Bertin (IAP)
VERSION 03/07/2003
 ***/
void threads_gate_sync(threads_gate_t *gate)

  {
  if (gate->nthreads<2)
    return;		/* trivial case */
  QPTHREAD_MUTEX_LOCK(&gate->block);		/* lock the block -- new */
						/* threads sleep here */
  QPTHREAD_MUTEX_LOCK(&gate->mutex);		/* lock the mutex */
  if (++(gate->ngate) < gate->nthreads)		/* are we the last one in? */
    {
    QPTHREAD_MUTEX_UNLOCK(&gate->block);	/* no, unlock block and */
    QPTHREAD_COND_WAIT(&gate->condvar, &gate->mutex);	/* go to sleep */
    }
  else						/* yes, we're last */
    {
    if (gate->func)
      (*gate->func)();				/* execute function */
    QPTHREAD_COND_BROADCAST(&gate->condvar);	/* wake everyone up and */
    QPTHREAD_COND_WAIT(&gate->last, &gate->mutex);/* sleep til all awake */
    QPTHREAD_MUTEX_UNLOCK(&gate->block);	/* release the block */
    }
  if (--(gate->ngate)==1)			/* next to last one out? */
    QPTHREAD_COND_BROADCAST(&gate->last);	/* yes, wake up last one */
  QPTHREAD_MUTEX_UNLOCK(&gate->mutex);		/* release the mutex */

  return;
  }

#endif
