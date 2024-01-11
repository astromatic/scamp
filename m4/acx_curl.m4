dnl
dnl				acx_curl.m4
dnl
dnl Figure out if the cURL library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2003-2023 Emmanuel Bertin -- IAP/CFHT/CNRS/UPMC
dnl
dnl	License:		GNU General Public License
dnl
dnl	AstrOmatic software is free software: you can redistribute it and/or
dnl	modify it under the terms of the GNU General Public License as
dnl	published by the Free Software Foundation, either version 3 of the
dnl	License, or (at your option) any later version.
dnl	AstrOmatic software is distributed in the hope that it will be useful,
dnl	but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl	GNU General Public License for more details.
dnl	You should have received a copy of the GNU General Public License
dnl	along with AstrOmatic software.
dnl	If not, see <http://www.gnu.org/licenses/>.
dnl
dnl	Last modified:		15/02/2023
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_CURL([CURL_LIBDIR, CURL_INCDIR,
dnl                      [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the cURL library and header files
dnl are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$CURL_LIBS $LIBS"
dnl        CFLAGS="$CFLAGS $CURL_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if libcurl
dnl is found (HAVE_CURL is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_CURL], [
AC_REQUIRE([AC_CANONICAL_HOST])

CURL_LIBS=""
OLIBS="$LIBS"
LIBS=""
acx_curl_ok=yes
if test x$2 = x && test x$1 = x; then
  AC_CHECK_PROG(acx_curlconfig_ok, [curl-config], [yes], [no])
  if test x$acx_curlconfig_ok = xyes; then
    [CURL_CFLAGS=` --cflags`]
    [CURL_LIBS=`curl-config --libs`]
    AC_DEFINE(CURL_H, "curl/curl.h", [cURL header filename.])
  fi
else
  acx_curlconfig_ok=no
fi
if test x$acx_curlconfig_ok = xno; then
  if test x$2 = x; then
    AC_CHECK_HEADER(curl/curl.h, [acx_curlhead_ok=yes], [acx_curlhead_ok=no])
    if test x$acx_curlhead_ok = xyes; then
      AC_DEFINE(CURL_H, "curl/curl.h", [cURL header filename.])
    else
      AC_CHECK_HEADER(curl.h, [acx_curlhead_ok=yes], [acx_curlhead_ok=no])
      if test x$acx_curlhead_ok = xyes; then
        AC_DEFINE(CURL_H, "curl.h", [cURL header filename.])
      else
        acx_curl_ok=no
      fi
    fi
  else
    AC_CHECK_HEADER($2/curl/curl.h,
		[acx_curlhead_ok=yes], [acx_curlhead_ok=no])
    if test x$acx_curlhead_ok = xyes; then
      AC_DEFINE(CURL_H, "curl/curl.h", [cURL header filename.])
      [CURL_CFLAGS="-I$2"]
    else
      AC_CHECK_HEADER($2/curl.h,
		[acx_curlhead_ok=yes], [acx_curlhead_ok=no])
      if test x$acx_curlhead_ok = xyes; then
        AC_DEFINE(CURL_H, "curl.h", [cURL header filename.])
        [CURL_CFLAGS="-I$2"]
      else
        acx_curl_ok=no
      fi
    fi
  fi
  if test x$1 = x; then
    AC_SEARCH_LIBS(curl_easy_init, [curl],
		[CURL_LIBS="$ac_cv_search_curl_easy_init"],
		[acx_curl_ok=no])
  else
    AC_SEARCH_LIBS(curl_easy_init, [curl],
		[CURL_LIBS="-L$1 $ac_cv_search_curl_easy_init"],
		[acx_curl_ok=no], [-L$1])
  fi
fi

LIBS="$OLIBS"
if test x$acx_curl_ok = xyes; then
  AC_SUBST(CURL_CFLAGS)
  AC_SUBST(CURL_LIBS)
  AC_DEFINE(HAVE_CURL,1,
	[Define if you have the cURL libraries and header files.])
  $3
else
  $4
fi

])dnl ACX_CURL

