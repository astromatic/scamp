dnl
dnl				acx_atlas.m4
dnl
dnl Figure out if the ATLAS library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2016 IAP/CNRS/UPMC
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
dnl	Last modified:		03/10/2016
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_LAPACKE([LAPACKE_LIBSDIR, LAPACKE_INCDIR,
dnl                   BLAS_LIBS, LAPACKE64_FLAG, 
dnl                   [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$LAPACKE_LIBS $LIBS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if LAPACKe
dnl is found (HAVE_LAPACKE is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_LAPACKE], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

acx_lapacke_ok=no
if test x$2 = x; then
  AC_CHECK_HEADER(
    [lapacke/lapacke.h],
    [
      [acx_lapacke_ok=yes]
      AC_DEFINE_UNQUOTED(LAPACKE_H, ["lapacke/lapacke.h"],
	[LAPACKe header filename.])
    ],
    [AC_CHECK_HEADER(
      [lapacke.h],
      [
        [acx_lapacke_ok=yes]
        AC_DEFINE_UNQUOTED(LAPACKE_H, ["lapacke.h"],
		[LAPACKe header filename.])
      ],
      [LAPACKE_ERROR="LAPACKe include files not found!"]
    )]
  )
else
  AC_CHECK_HEADER(
    [$2/lapacke.h],
    [
      [acx_lapacke_ok=yes]
      AC_DEFINE_UNQUOTED(LAPACKE_H, ["$2/lapacke.h"], [LAPACKe header filename.])
    ],
    [AC_CHECK_HEADER(
      [$2/include/lapacke.h],
      [
        [acx_lapacke_ok=yes]
        AC_DEFINE_UNQUOTED(LAPACKE_H, ["$2/include/lapacke.h"],
		[LAPACKe header filename.])
      ],
      [LAPACKE_ERROR="LAPACKe include files not found in "$2"!"]
    )]
  )
fi

dnl -------------------------
dnl Search library files
dnl -------------------------

if test x$acx_lapacke_ok = xyes; then
  acx_lapacke_ok=no
  OLIBS="$LIBS"
  LIBS=""
  if test x$4 = xyes; then
    lapacke_suffix="64"
  else
    lapacke_suffix=""
  fi
  if test x$1 = x; then
    lapacke_libopt="-llapack$lapacke_suffix $3"
  else
    lapacke_libopt="-L$1 -llapack$lapacke_suffix $3"
  fi

  AC_SEARCH_LIBS(
        LAPACKE_dpotrf, [lapacke],
        [acx_lapacke_ok=yes],
        [LAPACKE_ERROR="LAPACK"$lapacke_suffix"/LAPACKe library files not found!"],
        $lapacke_libopt
  )
  LIBS="$OLIBS"
fi

dnl -------------------------------------------------------------------------
dnl Finally execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND
dnl -------------------------------------------------------------------------

if test x"$acx_lapacke_ok" = xyes; then
  AC_DEFINE(HAVE_LAPACKE,1,
	[Define if you have the LAPACKe libraries and header files.])
  LAPACKE_LIBS="-llapacke $lapacke_libopt $ac_cv_search_clapack_dpotrf"
  AC_SUBST(LAPACKE_CFLAGS)
  AC_SUBST(LAPACKE_LIBS)
  AC_SUBST(LAPACKE_WARN)
  $5
else
  AC_SUBST(LAPACKE_ERROR)
  $6
fi

])dnl ACX_LAPACKE

