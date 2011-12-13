dnl
dnl				acx_mkl_lapacke.m4
dnl
dnl Figure out if the INTEL MKL LAPACKe library and header files are installed.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2003-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
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
dnl	Last modified:		13/12/2011
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_MKL_LAPACKE([LAPACKE_LIBDIR, LAPACKE_INCDIR,
dnl			PARALLEL_FLAG,
dnl                     [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$MKL_LAPACKE_LIBS $LIBS"
dnl
dnl You may wish to use these variables in your default CFLAGS:
dnl
dnl        CFLAGS="$CFLAGS $MKL_LAPACKE_CFLAGS"
dnl
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if LAPACKe is found
dnl (HAVE_LAPACKE and possibly HAVE_LAPACKE_MP are defined first), and
dnl ACTION-IF-NOT-FOUND is a list of commands to run it if it is not found.

AC_DEFUN([ACX_MKL_LAPACKE], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

acx_lapacke_ok=no
if test x$2 = x; then
  if test x$1 = x; then
    AC_CHECK_HEADERS([mkl_lapacke.h],[acx_lapacke_ok=yes])
    if test x$acx_lapacke_ok = xyes; then
      AC_DEFINE(LAPACKE_H, "mkl_lapacke.h", [LAPACKe header filename.])
    else
      AC_CHECK_HEADERS(
		[$MKLROOT/include/mkl_lapacke.h], [acx_lapacke_ok=yes])
      if test x$acx_lapacke_ok = xyes; then
        AC_DEFINE_UNQUOTED(LAPACKE_H, "$MKLROOT/include/mkl_lapacke.h",
		[LAPACKe header filename.])
      else
        MKL_LAPACKE_ERROR="INTEL MKL LAPacke include files not found!"
      fi
    fi
  else
    AC_CHECK_HEADERS([$1/include/mkl_lapacke.h], [acx_lapacke_ok=yes])
    if test x$acx_lapacke_ok = xyes; then
      AC_DEFINE_UNQUOTED(LAPACKE_H, "$1/include/mkl_lapacke.h",
		[LAPACKe header filename.])
    else
      AC_CHECK_HEADERS([$1/../include/mkl_lapacke.h], [acx_lapacke_ok=yes])
      if test x$acx_lapacke_ok = xyes; then
        AC_DEFINE_UNQUOTED(LAPACKE_H, "$1/../include/mkl_lapacke.h",
		[LAPACKe header filename.])
      else
        AC_CHECK_HEADERS([$1/../../include/mkl_lapacke.h], [acx_lapacke_ok=yes])
        if test x$acx_lapacke_ok = xyes; then
          AC_DEFINE_UNQUOTED(LAPACKE_H, "$1/../../include/mkl_lapacke.h",
		[LAPACKe header filename.])
        else
          AC_CHECK_HEADERS([mkl_lapacke.h],[acx_lapacke_ok=yes])
          if test x$acx_lapacke_ok = xyes; then
            AC_DEFINE_UNQUOTED(LAPACKE_H, "mkl_lapacke.h",
		[LAPACKe header filename.])
          else
            MKL_LAPACKE_ERROR="INTEL MKL LAPacke include files not found!"
          fi
        fi
      fi
    fi
  fi
else
  AC_CHECK_HEADERS([$2/include/mkl_lapacke.h], [acx_lapacke_ok=yes])
  if test x$acx_lapacke_ok = xyes; then
    AC_DEFINE_UNQUOTED(LAPACKE_H, "$2/include/mkl_lapacke.h",
	[LAPACKe header filename.])
  else
    AC_CHECK_HEADERS([$2/mkl_lapacke.h], [acx_lapacke_ok=yes])
    if test x$acx_lapacke_ok = xyes; then
      AC_DEFINE_UNQUOTED(LAPACKE_H, "$2/mkl_lapacke.h",
	[LAPACKe header filename.])
    else
      MKL_LAPACKE_ERROR="INTEL MKL LAPacke include files not found in $2!"
    fi
  fi
fi

dnl --------------------
dnl Search library files
dnl --------------------

OLDFLAGS="$LDFLAGS"
OLIBS="$LIBS"
LIBS=""
if test x$acx_lapacke_ok = xyes; then
dnl --------------------
dnl Try to find INTEL architecture (Intel 64 or ia32)
dnl --------------------
  if icc -V 2>&1 | grep -i "Intel(R) 64" > /dev/null 2>&1; then
    mkl_root="$MKLROOT/intel64"
    mkl_cflags="-DMKL_ILP64"
    mkl_lapacke_lib="mkl_intel_ilp64"
  elif icc -V 2>&1 | grep -i "Intel(R)" > /dev/null 2>&1; then
    mkl_root="$MKLROOT/ia32"
    mkl_cflags=""
    mkl_lapacke_lib="mkl_intel"
  else
    mkl_root=""
    MKL_LAPACKE_WARN="INTEL compiler not detected"
    AC_SUBST(MKL_LAPACKE_WARN)
  fi
  if test x$3 = xyes; then
dnl Check the parallel version of the MKL:
    if test x$1 = x; then
      unset ac_cv_lib_"$mkl_lapacke_lib"_LAPACKE_dpotrf
      AC_CHECK_LIB($mkl_lapacke_lib, [LAPACKE_dpotrf],, [acx_lapacke_ok=no],
		[-lmkl_intel_thread -lmkl_core -openmp -lpthread])
      if test x$acx_lapacke_ok = xyes; then
        MKL_LAPACKE_LIBS="-l$mkl_lapacke_lib -lmkl_intel_thread -lmkl_core -openmp"
        AC_DEFINE(HAVE_LAPACKE_MP,1, [Define if you have the parallel LAPACKe libraries.])
      else
        unset ac_cv_lib_"$mkl_lapacke_lib"_LAPACKE_dpotrf
        acx_lapacke_ok=yes
        ACX_SEARCH_LIBDIR(MKL_LAPACKE_LIBS, $mkl_root/lib, $mkl_lapacke_lib,
			[LAPACKE_dpotrf],,[acx_lapacke_ok=no],
			[-lmkl_intel_thread -lmkl_core -openmp -lpthread])
        if test x$acx_lapacke_ok = xyes; then
          AC_DEFINE(HAVE_LAPACKE_MP,1, [Define if you have the parallel LAPACKe libraries.])
        else
          MKL_LAPACKE_ERROR="INTEL MKL parallel LAPACKe library files not found at usual locations!"
        fi
      fi
    else
      unset ac_cv_lib_"$mkl_lapacke_lib"_LAPACKE_dpotrf
      ACX_SEARCH_LIBDIR(MKL_LAPACKE_LIBS, $1/lib $1, $mkl_lapacke_lib,
			[LAPACKE_dpotrf],,[acx_lapacke_ok=no],
			[-lmkl_intel_thread -lmkl_core -openmp -lpthread])
      if test x$acx_lapacke_ok = xyes; then
        AC_DEFINE(HAVE_LAPACKE_MP,1, [Define if you have the parallel LAPACKe libraries.])
      else
        MKL_LAPACKE_ERROR="INTEL MKL parallel LAPACKe library files not found at the provided location!"
      fi
    fi
  else
dnl Check the serial version of the MKL:
    if test x$1 = x; then
      unset ac_cv_lib_"$mkl_lapacke_lib"_LAPACKE_dpotrf
      AC_CHECK_LIB($mkl_lapacke_lib, [LAPACKE_dpotrf],, [acx_lapacke_ok=no],
		[-lmkl_sequential -lmkl_core])
      if test x$acx_lapacke_ok = xyes; then
        MKL_LAPACKE_LIBS="-l$mkl_lapacke_lib -lmkl_sequential -lmkl_core"
      else
        unset ac_cv_lib_"$mkl_lapacke_lib"_LAPACKE_dpotrf
        acx_lapacke_ok=yes
        ACX_SEARCH_LIBDIR(MKL_LAPACKE_LIBS, $mkl_root/lib, $mkl_lapacke_lib,
			[LAPACKE_dpotrf],,[acx_lapacke_ok=no],
			[-lmkl_sequential -lmkl_core])
        if test x$acx_lapacke_ok = xno; then
          MKL_LAPACKE_ERROR="INTEL MKL serial LAPACKe library files not found at usual locations!"
        fi
      fi
      unset ac_cv_lib_"$mkl_lapacke_lib"_LAPACKE_dpotrf
      ACX_SEARCH_LIBDIR(MKL_LAPACKE_LIBS, $1/lib $1, $mkl_lapacke_lib,
			[LAPACKE_dpotrf],,[acx_lapacke_ok=no],
			[-lmkl_sequential -lmkl_core])
      if test x$acx_lapacke_ok = xno; then
          MKL_LAPACKE_ERROR="INTEL MKL parallel LAPACKe library files not found at usual locations!"
      fi
    fi
  fi
fi

LDFLAGS="$OLDFLAGS"
LIBS="$OLIBS"
if test x$acx_lapacke_ok = xyes; then
  AC_DEFINE(HAVE_LAPACKE,1, [Define if you have the LAPACKe libraries.])
  AC_SUBST(MKL_LAPACKE_CFLAGS, $mkl_cflags)
  AC_SUBST(MKL_LAPACKE_LIBS)
  $4
else
  AC_SUBST(MKL_LAPACKE_ERROR)
  $5
fi
])dnl ACX_MKL_LAPACKE
