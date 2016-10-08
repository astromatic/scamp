dnl
dnl				acx_blas.m4
dnl
dnl Figure out if the (C)BLAS library and header files are installed.
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
dnl	Last modified:		08/10/2016
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_BLAS([BLAS_LIBSDIR, BLAS_INCDIR, BLASP_FLAG,
dnl                  BLAS64_FLAG, [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$BLAS_LIBS $LIBS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if BLAS
dnl is found (HAVE_BLAS is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_BLAS], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

BLAS_ERROR=""
if test x$2 = x; then
dnl We give our preference to OpenBLAS
  [acx_blas_incdir="openblas/"]
  AC_CHECK_HEADER(
    ${acx_blas_incdir}cblas.h,,
    [
      [acx_blas_incdir=""]
      AC_CHECK_HEADER(
        [cblas.h],,
        [BLAS_ERROR="BLAS header file not found!"]
      )
    ]
  )
else
  acx_blas_incdir="$2/"
  AC_CHECK_HEADER(
    [${acx_blas_incdir}cblas.h],,
    [
      [acx_blas_incdir="$2/include/"]
      AC_CHECK_HEADER(
        [${acx_blas_incdir}cblas.h],,
        [BLAS_ERROR="BLAS header file not found in "$2"!"]
    )]
  )
fi

if test x$BLAS_ERROR = x; then
  AC_DEFINE_UNQUOTED(BLAS_H, "${acx_blas_incdir}cblas.h", [BLAS header filename.])

dnl -------------------------
dnl Search BLAS library files
dnl -------------------------

  OLIBS="$LIBS"
  LIBS=""
  if test x$4 = xyes; then
    acx_blas_suffix="64"
    BLAS_CFLAGS="-DOPENBLAS_USE64BITINT"
  else
    acx_blas_suffix=""
    BLAS_CFLAGS=""
  fi
  if test x$1 = x; then
    acx_blas_libopt=""
  else
    acx_blas_libopt="-L$1"
  fi
  if test x$3 == xyes; then
    AC_SEARCH_LIBS(
      cblas_dgemm, ["openblasp"$blas_suffix],
      AC_DEFINE(HAVE_OPENBLASP,1,
		[Define if you have the OpenBLAS parallel libraries and header files.]),
      [AC_SEARCH_LIBS(
        cblas_dgemm, ["openblas"$blas_suffix],
	[BLAS_WARN="parallel OpenBLAS"$blas_suffix" not found, reverting to scalar OpenBLAS"$blas_suffix"!"],
        [AC_SEARCH_LIBS(
          cblas_dgemm, ["blas"$blas_suffix],
	  [BLAS_WARN="parallel OpenBLAS"$blas_suffix" not found, reverting to scalar BLAS"$blas_suffix"!"],
          [BLAS_ERROR="CBLAS"$blas_suffix" library files not found!"],
          $blas_libopt
        )],
        $blas_libopt
      )],
      $blas_libopt
    )
  else
    AC_SEARCH_LIBS(
      cblas_dgemm, ["openblas"$blas_suffix],,
      [AC_SEARCH_LIBS(
        cblas_dgemm, ["blas"$blas_suffix],
	[BLAS_WARN="OpenBLAS"$blas_suffix" not found, reverting to regular BLAS"$blas_suffix"!"],
        [BLAS_ERROR="CBLAS"$blas_suffix" library files not found!"],
        $blas_libopt
      )],
      $blas_libopt
    )
  fi
  LIBS="$OLIBS"
fi

dnl -------------------------------------------------------------------------
dnl Finally execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND
dnl -------------------------------------------------------------------------

if test x$BLAS_ERROR = x; then
  AC_DEFINE(HAVE_BLAS,1, [Define if you have the BLAS libraries and header files.])
  BLAS_LIBS="$blas_libopt $ac_cv_search_cblas_dgemm"
  AC_SUBST(BLAS_CFLAGS)
  AC_SUBST(BLAS_LDFLAGS, "")
  AC_SUBST(BLAS_LIBS)
  AC_SUBST(BLAS_WARN)
  $5
else
  AC_SUBST(BLAS_ERROR)
  $6
fi

])dnl ACX_BLAS

