dnl
dnl				acx_lbfgs.m4
dnl
dnl Figure out if the L-BFGS library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2003-2021 IAP/CNRS/SorbonneU
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
dnl	Last modified:		14/04/2021
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_LBFGS([LBFGS_LIBDIR, LBFGS_INCDIR,
dnl                      [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the L-BFGS library and header files
dnl are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$LBFGS_LIBS $LIBS"
dnl        CFLAGS="$CFLAGS $LBFGS_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if liblbfgs
dnl is found (HAVE_LBFGS is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_LBFGS], [
AC_REQUIRE([AC_CANONICAL_HOST])

LBFGS_LIBS=""
OLIBS="$LIBS"
LIBS=""
acx_lbfgs_ok=yes
acx_lbfgsconfig_ok=no
if test x$2 = x && test x$1 = x; then
  AC_CHECK_PROG(acx_lbfgsconfig_ok, [lbfgs-config], [yes], [no])
  if test x"$acx_lbfgsconfig_ok" = xyes; then
    [LBFGS_CFLAGS=`lbfgs-config --cflags`]
    [LBFGS_LIBS=`lbfgs-config --libs`]
    AC_DEFINE(LBFGS_H, "lbfgs.h", [L-BFGS header filename.])
  fi
fi
if test x$acx_lbfgsconfig_ok = xno; then
  if test x$2 = x; then
    AC_CHECK_HEADER(lbfgs.h, [acx_lbfgshead_ok=yes], [acx_lbfgshead_ok=no])
    if test x$acx_lbfgshead_ok = xyes; then
      AC_DEFINE(LBFGS_H, "lbfgs.h", [L-BFGS header filename.])
    else
      AC_CHECK_HEADER(lbfgs/lbfgs.h, [acx_lbfgshead_ok=yes], [acx_lbfgshead_ok=no])
      if test x$acx_lbfgshead_ok = xyes; then
        AC_DEFINE(LBFGS_H, "lbfgs/lbfgs.h", [L-BFGS header filename.])
      else
        acx_lbfgs_ok=no
      fi
    fi
  else
    AC_CHECK_HEADER($2/lbfgs.h,
		[acx_lbfgshead_ok=yes], [acx_lbfgshead_ok=no])
    if test x$acx_lbfgshead_ok = xyes; then
      AC_DEFINE(LBFGS_H, "lbfgs.h", [L-BFGS header filename.])
      [LBFGS_CFLAGS="-I$2"]
    else
      AC_CHECK_HEADER($2/lbfgs/lbfgs.h,
		[acx_lbfgshead_ok=yes], [acx_lbfgshead_ok=no])
      if test x$acx_lbfgshead_ok = xyes; then
        AC_DEFINE(LBFGS_H, "lbfgs/lbfgs.h", [L-BFGS header filename.])
        [LBFGS_CFLAGS="-I$2"]
      else
        acx_lbfgs_ok=no
      fi
    fi
  fi
  if test x$1 = x; then
    AC_SEARCH_LIBS(lbfgs, [lbfgs],
		[LBFGS_LIBS="$ac_cv_search_lbfgs"],
		[acx_lbfgs_ok=no])
  else
    AC_SEARCH_LIBS(lbfgs, [lbfgs],
		[LBFGS_LIBS="-L$1 $ac_cv_search_lbfgs"],
		[acx_lbfgs_ok=no], [-L$1])
  fi
fi

LIBS="$OLIBS"
if test x$acx_lbfgs_ok = xyes; then
  AC_SUBST(LBFGS_CFLAGS)
  AC_SUBST(LBFGS_LIBS)
  AC_DEFINE(HAVE_LBFGS,1,
	[Define if you have the L-BFGS libraries and header files.])
  $3
else
  $4
fi

])dnl ACX_LBFGS

