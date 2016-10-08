dnl
dnl				acx_lapacke.m4
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
dnl	Last modified:		06/10/2016
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

dnl -----------------------
dnl Search main header file
dnl -----------------------

LAPACKE_ERROR=""
if test x$2 = x; then
  [acx_lapacke_incdir="lapacke/"]
  AC_CHECK_HEADER(
    ${acx_lapacke_incdir}lapacke.h,,
    [
      [acx_lapacke_incdir=""]
      AC_CHECK_HEADER(
        [lapacke.h],,
        [LAPACKE_ERROR="LAPACKe header file not found!"]
      )
    ]
  )
else
  acx_lapacke_incdir="$2/"
  AC_CHECK_HEADER(
    [${acx_lapacke_incdir}lapacke.h],,
    [
      [acx_lapacke_incdir="$2/include/"]
      AC_CHECK_HEADER(
        [${acx_lapacke_incdir}lapacke.h],,
        [LAPACKE_ERROR="LAPACKe header file not found in "$2"!"]
    )]
  )
fi

if test x$LAPACKE_ERROR = x; then
  AC_DEFINE_UNQUOTED(LAPACKE_H, "${acx_lapacke_incdir}lapacke.h",
	[LAPACKe header filename.])

dnl ---------------------------
dnl Search "config" header file
dnl ---------------------------
  AC_CHECK_HEADER(
    ${acx_lapacke_incdir}lapacke_config.h,
    AC_DEFINE(HAVE_LAPACK_CONFIG_H, 1,
      [Define if you have the LAPACKe config header file.]
    )
  )

dnl -------------------------
dnl Search library files
dnl -------------------------

  OLIBS="$LIBS"
  LIBS=""
  if test x$4 = xyes; then
    acx_lapacke_suffix="64"
    LAPACKE_CFLAGS="-DLAPACK_ILP64"
  else
    acx_lapacke_suffix=""
    LAPACKE_CFLAGS=""
  fi
  if test x$1 = x; then
    acx_lapacke_libopt="-llapack$acx_lapacke_suffix $3"
  else
    acx_lapacke_libopt="-L$1 -llapack$acx_lapacke_suffix $3"
  fi
  AC_SEARCH_LIBS(
        LAPACKE_dpotrf, [lapacke],,
        [LAPACKE_ERROR="LAPACK"$acx_lapacke_suffix"/LAPACKe library files not found!"],
        [$acx_lapacke_libopt]
  )
  LIBS="$OLIBS"
fi

dnl -------------------------------------------------------------------------
dnl Finally execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND
dnl -------------------------------------------------------------------------

if test x$LAPACKE_ERROR = x; then
  AC_DEFINE(HAVE_LAPACKE,1,
	[Define if you have the LAPACKe libraries and header files.])
  LAPACKE_LIBS="-llapacke $acx_lapacke_libopt $ac_cv_search_clapack_dpotrf"
  AC_SUBST(LAPACKE_CFLAGS)
  AC_SUBST(LAPACKE_LDFLAGS, "")
  AC_SUBST(LAPACKE_LIBS)
  AC_SUBST(LAPACKE_WARN)
  $5
else
  AC_SUBST(LAPACKE_ERROR)
  $6
fi

])dnl ACX_LAPACKE

