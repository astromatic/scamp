/*
*				url.c
*
* Manage HTTP communications using libcurl (adapted from fopen.c by V.Sanders).
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
*
*	This file part of:	SCAMP
*
*	Copyright:		(C) 2002-2018 IAP/CNRS/SorbonneU
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
*	This file incorporates work covered by the copyright and permission
*	notice below.
*
*	Last modified:		05/09/2018
*
*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*****************************************************************************
 *
 * This example source code introduces a c library buffered I/O interface to
 * URL reads it supports fopen(), fread(), fgets(), feof(), fclose(),
 * rewind(). Supported functions have identical prototypes to their normal c
 * lib namesakes and are preceaded by url_ .
 *
 * Using this code you can replace your program's fopen() with url_fopen()
 * and fread() with url_fread() and it become possible to read remote streams
 * instead of (only) local files. Local files (ie those that can be directly
 * fopened) will drop back to using the underlying clib implementations
 *
 * Copyright (c) 2003 Simtec Electronics
 *
 * Re-implemented by Vincent Sanders <vince@kyllikki.org> with extensive
 * reference to original curl example code
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * This example requires libcurl 7.9.7 or later.
 */


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <sys/time.h>

#include CURL_H

#include "url.h"

/* we use a global one for convenience */
CURLM *multi_handle;

/****** write_callback ****************************************************//**
Curl calls this routine to get more data
 ***/
static size_t write_callback(char *buffer,
                             size_t size,
                             size_t nitems,
                             void *userp)
{
  char *newbuff;
  size_t rembuff;

  URL_FILE *url = (URL_FILE *)userp;
  size *= nitems;

  rembuff=url->buffer_len - url->buffer_pos; /* remaining space in buffer */

  if(size > rembuff) {
    /* not enough space in buffer */
    newbuff=realloc(url->buffer,url->buffer_len + (size - rembuff));
    if(newbuff==NULL) {
      fprintf(stderr,"callback buffer grow failed\n");
      size=rembuff;
    }
    else {
      /* realloc suceeded increase buffer size*/
      url->buffer_len+=size - rembuff;
      url->buffer=newbuff;
    }
  }

  memcpy(&url->buffer[url->buffer_pos], buffer, size);
  url->buffer_pos += size;

  return size;
}

/****** fill_buffer ****************************************************//**
Use to attempt to fill the read buffer up to requested number of bytes
 ***/
static int fill_buffer(URL_FILE *file, size_t want)
{
  fd_set fdread;
  fd_set fdwrite;
  fd_set fdexcep;
  struct timeval timeout;
  int rc;
  CURLMcode mc; /* curl_multi_fdset() return code */

  /* only attempt to fill buffer if transactions still running and buffer
   * doesnt exceed required size already
   */
  if((!file->still_running) || (file->buffer_pos > want))
    return 0;

  /* attempt to fill buffer */
  do {
    int maxfd = -1;
    long curl_timeo = -1;

    FD_ZERO(&fdread);
    FD_ZERO(&fdwrite);
    FD_ZERO(&fdexcep);

    /* set a suitable timeout to fail on */
    timeout.tv_sec = 60; /* 1 minute */
    timeout.tv_usec = 0;

    curl_multi_timeout(multi_handle, &curl_timeo);
    if(curl_timeo >= 0) {
      timeout.tv_sec = curl_timeo / 1000;
      if(timeout.tv_sec > 1)
        timeout.tv_sec = 1;
      else
        timeout.tv_usec = (curl_timeo % 1000) * 1000;
    }

    /* get file descriptors from the transfers */
    mc = curl_multi_fdset(multi_handle, &fdread, &fdwrite, &fdexcep, &maxfd);

    if(mc != CURLM_OK)
    {
      fprintf(stderr, "curl_multi_fdset() failed, code %d.\n", mc);
      break;
    }

    /* On success the value of maxfd is guaranteed to be >= -1. We call
       select(maxfd + 1, ...); specially in case of (maxfd == -1) there are
       no fds ready yet so we call select(0, ...) --or Sleep() on Windows--
       to sleep 100ms, which is the minimum suggested value in the
       curl_multi_fdset() doc. */

    if(maxfd == -1) {
#ifdef _WIN32
      Sleep(100);
      rc = 0;
#else
      /* Portable sleep for platforms other than Windows. */
      struct timeval wait = { 0, 100 * 1000 }; /* 100ms */
      rc = select(0, NULL, NULL, NULL, &wait);
#endif
    }
    else {
      /* Note that on some platforms 'timeout' may be modified by select().
         If you need access to the original value save a copy beforehand. */
      rc = select(maxfd+1, &fdread, &fdwrite, &fdexcep, &timeout);
    }

    switch(rc) {
    case -1:
      /* select error */
      break;

    case 0:
    default:
      /* timeout or readable/writable sockets */
      curl_multi_perform(multi_handle, &file->still_running);
      break;
    }
  } while(file->still_running && (file->buffer_pos < want));
  return 1;
}

/****** use_buffer ****************************************************//**
Use to remove want bytes from the front of a files buffer
 ***/
static int use_buffer(URL_FILE *file,int want)
{
  /* sort out buffer */
  if((file->buffer_pos - want) <=0) {
    /* ditch buffer - write will recreate */
    free(file->buffer);
    file->buffer=NULL;
    file->buffer_pos=0;
    file->buffer_len=0;
  }
  else {
    /* move rest down make it available for later */
    memmove(file->buffer,
            &file->buffer[want],
            (file->buffer_pos - want));

    file->buffer_pos -= want;
  }
  return 0;
}

/****** url_fopen ****************************************************//**
Check for URL and set up handle
 ***/
URL_FILE *url_fopen(const char *url, double timeout)
{
  URL_FILE *file;

  file = malloc(sizeof(URL_FILE));
  if(!file)
    return NULL;

  memset(file, 0, sizeof(URL_FILE));
  file->curl = curl_easy_init();

  curl_easy_setopt(file->curl, CURLOPT_URL, url);
  curl_easy_setopt(file->curl, CURLOPT_WRITEDATA, file);
  curl_easy_setopt(file->curl, CURLOPT_VERBOSE, 0L);
  curl_easy_setopt(file->curl, CURLOPT_FOLLOWLOCATION, 1L);
  curl_easy_setopt(file->curl, CURLOPT_WRITEFUNCTION, write_callback);
  if (timeout > 0)
    curl_easy_setopt(file->curl, CURLOPT_CONNECTTIMEOUT_MS,
	(long)(timeout * 1000.0));

  if(!multi_handle)
    multi_handle = curl_multi_init();

  curl_multi_add_handle(multi_handle, file->curl);

  /* lets start the fetch */
  curl_multi_perform(multi_handle, &file->still_running);

  if((file->buffer_pos == 0) && (!file->still_running)) {
    /* if still_running is 0 now, we should return NULL */

    /* make sure the easy handle is not in the multi handle anymore */
    curl_multi_remove_handle(multi_handle, file->curl);

    /* cleanup */
    curl_easy_cleanup(file->curl);

    free(file);

    file = NULL;
  }

  return file;
}

/****** url_fclose ****************************************************//**
Removes handle and cleanup
 ***/
int url_fclose(URL_FILE *file)
{
  int ret=0;/* default is good return */

  /* make sure the easy handle is not in the multi handle anymore */
  curl_multi_remove_handle(multi_handle, file->curl);

  /* cleanup */
  curl_easy_cleanup(file->curl);

  free(file->buffer);/* free any allocated buffer space */
  free(file);

  return ret;
}

/****** url_feof ****************************************************//**
Check if the end of URL data has been reached
 ***/
int url_feof(URL_FILE *file)
{
  int ret=0;

  if((file->buffer_pos == 0) && (!file->still_running))
    ret = 1;

  return ret;
}

/****** url_fread ****************************************************//**
Read URL data
 ***/
size_t url_fread(void *ptr, size_t size, size_t nmemb, URL_FILE *file)
{
  size_t want;

  want = nmemb * size;

  fill_buffer(file,want);

  /* check if theres data in the buffer - if not fill_buffer()
   * either errored or EOF */
  if(!file->buffer_pos)
    return 0;

  /* ensure only available data is considered */
  if(file->buffer_pos < want)
    want = file->buffer_pos;

  /* xfer data to caller */
  memcpy(ptr, file->buffer, want);

  use_buffer(file,want);

  return want / size;     /* number of items */

}

/****** url_fgets ****************************************************//**
Read URL ASCII lines
 ***/
char *url_fgets(char *ptr, size_t size, URL_FILE *file)
{
  size_t want = size - 1;/* always need to leave room for zero termination */
  size_t loop;

  fill_buffer(file,want);

  /* check if theres data in the buffer - if not fill either errored or EOF */
  if(!file->buffer_pos)
    return NULL;

  /* ensure only available data is considered */
  if(file->buffer_pos < want)
    want = file->buffer_pos;

  /*buffer contains data */
  /* look for newline or eof */
  for(loop=0;loop < want;loop++) {
    if(file->buffer[loop] == '\n') {
      want=loop+1;/* include newline */
      break;
    }
  }

  /* xfer data to caller */
  memcpy(ptr, file->buffer, want);
  ptr[want]=0;/* allways null terminate */

  use_buffer(file,want);

  return ptr;/*success */
}

/****** url_rewind ****************************************************//**
Set the stream position indicator back to the beginning of URL data
 ***/
void url_rewind(URL_FILE *file)
{
  /* halt transaction */
  curl_multi_remove_handle(multi_handle, file->curl);

  /* restart */
  curl_multi_add_handle(multi_handle, file->curl);

  /* ditch buffer - write will recreate - resets stream pos*/
  free(file->buffer);
  file->buffer=NULL;
  file->buffer_pos=0;
  file->buffer_len=0;

  return;
}

