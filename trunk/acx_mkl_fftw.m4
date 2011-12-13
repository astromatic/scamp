dnl
dnl				acx_mkl_fftw.m4
dnl
dnl Figure out if the INTEL MKL FFTW library and header files are installed.
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
dnl @synopsis ACX_MKL_FFTW([FFTW_LIBDIR, FFTW_INCDIR,
dnl			PARALLEL_FLAG, FLOAT_FLAG
dnl                     [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$MKL_FFTW_LIBS $LIBS"
dnl
dnl You may wish to use these variables in your default CFLAGS:
dnl
dnl        CFLAGS="$CFLAGS $MKL_FFTW_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if FFTW is found
dnl (HAVE_FFTW and possibly HAVE_FFTW_MP are defined first), and
dnl ACTION-IF-NOT-FOUND is a list of commands to run it if it is not found.

AC_DEFUN([ACX_MKL_FFTW], [
AC_REQUIRE([AC_CANONICAL_HOST])

dnl --------------------
dnl Search include files
dnl --------------------

acx_fftw_ok=no
if test x$2 = x; then
  if test x$1 = x; then
    AC_CHECK_HEADERS([fftw3_mkl.h],[acx_fftw_ok=yes])
    if test x$acx_fftw_ok = xyes; then
      AC_DEFINE(FFTW_H, "fftw3_mkl.h", [FFTW header filename.])
    else
      AC_CHECK_HEADERS([fftw/fftw3_mkl.h],[acx_fftw_ok=yes])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(FFTW_H, "fftw/fftw3_mkl.h", [FFTW header filename.])
      else
        AC_CHECK_HEADERS(
		[$MKLROOT/include/fftw/fftw3_mkl.h], [acx_fftw_ok=yes])
        if test x$acx_fftw_ok = xyes; then
          AC_DEFINE_UNQUOTED(FFTW_H, "$MKLROOT/include/fftw/fftw3_mkl.h",
		[FFTW header filename.])
        else
          MKL_FFTW_ERROR="INTEL MKL FFTW include files not found!"
        fi
      fi
    fi
  else
    AC_CHECK_HEADERS([$1/include/fftw/fftw3_mkl.h], [acx_fftw_ok=yes])
    if test x$acx_fftw_ok = xyes; then
      AC_DEFINE_UNQUOTED(FFTW_H, "$1/include/fftw/fftw3_mkl.h",
		[FFTW header filename.])
    else
      AC_CHECK_HEADERS([$1/../include/fftw/fftw3_mkl.h], [acx_fftw_ok=yes])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE_UNQUOTED(FFTW_H, "$1/../include/fftw/fftw3_mkl.h",
		[FFTW header filename.])
      else
        AC_CHECK_HEADERS([$1/../../include/fftw/fftw3_mkl.h],
		[acx_fftw_ok=yes])
        if test x$acx_ffw_ok = xyes; then
          AC_DEFINE_UNQUOTED(FFTW_H, "$1/../../include/fftw/fftw3_mkl.h",
		[FFTW header filename.])
        else
          AC_CHECK_HEADERS([fftw3_mkl.h],[acx_fftw_ok=yes])
          if test x$acx_fftw_ok = xyes; then
            AC_DEFINE_UNQUOTED(FFTW_H, "fftw3_mkl.h",
		[FFTW header filename.])
          else
            MKL_FFTW_ERROR="INTEL MKL FFTW include files not found!"
          fi
        fi
      fi
    fi
  fi
else
  AC_CHECK_HEADERS([$2/include/fftw/fftw3_mkl.h], [acx_fftw_ok=yes])
  if test x$acx_fftw_ok = xyes; then
    AC_DEFINE_UNQUOTED(FFTW_H, "$2/include/fftw/fftw3_mkl.h",
	[FFTW header filename.])
  else
    AC_CHECK_HEADERS([$2/fftw/fftw3_mkl.h], [acx_fftw_ok=yes])
    if test x$acx_fftw_ok = xyes; then
      AC_DEFINE_UNQUOTED(FFTW_H, "$2/fftw/fftw3_mkl.h",
	[FFTW header filename.])
    else
      AC_CHECK_HEADERS([$2/fftw3_mkl.h], [acx_fftw_ok=yes])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE_UNQUOTED(FFTW_H, "$2/fftw3_mkl.h", [FFTW header filename.])
      else
        MKL_FFTW_ERROR="INTEL MKL FFTW include files not found in $2!"
      fi
    fi
  fi
fi

dnl --------------------
dnl Search library files
dnl --------------------

OLIBS="$LIBS"
LIBS=""
if test x$acx_fftw_ok = xyes; then
dnl Try to find INTEL architecture (Intel 64 or ia32)
  if icc -V 2>&1 | grep -i "Intel(R) 64" > /dev/null 2>&1; then
    mkl_root="$MKLROOT/intel64"
    mkl_cflags="-DMKL_ILP64"
    mkl_fftw_lib="mkl_intel_ilp64"
  elif icc -V 2>&1 | grep -i "Intel(R)" > /dev/null 2>&1; then
    mkl_root="$MKLROOT/ia32"
    mkl_cflags=""
    mkl_fftw_lib="mkl_intel"
  else
    mkl_root=""
    MKL_FFTW_WARN="INTEL compiler not detected"
    AC_SUBST(MKL_FFTW_WARN)
  fi
  if test x$3 = xyes; then
dnl Check the parallel version of the MKL:
dnl Set precision
    if test x$4 = xyes; then
      mkl_func=fftwf_init_threads
    else
      mkl_func=fftw_init_threads
    fi
dnl Check if the function is in the library
    if test x$1 = x; then
      unset ac_cv_lib_"$mkl_fftw_lib"_"$mkl_func"
      AC_CHECK_LIB($mkl_fftw_lib, $mkl_func,, [acx_fftw_ok=no],
		[-lmkl_intel_thread -lmkl_core -openmp -lpthread])
      if test x$acx_fftw_ok = xyes; then
        MKL_FFTW_LIBS="-l$mkl_fftw_lib -lmkl_intel_thread -lmkl_core -openmp"
        AC_DEFINE(HAVE_FFTW_MP,1, [Define if you have the parallel FFTW libraries.])
      else
        unset ac_cv_lib_"$mkl_fftw_lib"_"$mkl_func"
        acx_fftw_ok=yes
        ACX_SEARCH_LIBDIR(MKL_FFTW_LIBS, $mkl_root/lib, $mkl_fftw_lib,
		$mkl_func,, [acx_fftw_ok=no],
		[-lmkl_intel_thread -lmkl_core -openmp -lpthread])
        if test x$acx_fftw_ok = xyes; then
          AC_DEFINE(HAVE_FFTW_MP,1, [Define if you have the parallel FFTW libraries.])
        else
          MKL_FFTW_ERROR="INTEL MKL parallel FFTW library files not found at usual locations!"
        fi
      fi
    else
      unset ac_cv_lib_"$mkl_fftw_lib"_"$mkl_func"
      ACX_SEARCH_LIBDIR(MKL_FFTW_LIBS, $1/lib $1, $mkl_fftw_lib,
		$mkl_func,, [acx_fftw_ok=no],
		[-lmkl_intel_thread -lmkl_core -openmp -lpthread])
      if test x$acx_fftw_ok = xyes; then
        AC_DEFINE(HAVE_FFTW_MP,1, [Define if you have the parallel FFTW libraries.])
      else
          MKL_FFTW_ERROR="INTEL MKL parallel FFTW library files not found at the provided location!"
      fi
    fi
  else
dnl Check the serial version of the MKL:
dnl Set precision
    if test x$4 = xyes; then
      mkl_func=fftwf_execute
    else
      mkl_func=fftw_execute
    fi
dnl Check if the function is in the library
    if test x$1 = x; then
      unset ac_cv_lib_"$mkl_fftw_lib"_"$mkl_func"
      AC_CHECK_LIB($mkl_fftw_lib, $mkl_func,, [acx_fftw_ok=no],
		[-lmkl_sequential -lmkl_core])
      if test x$acx_fftw_ok = xyes; then
        MKL_FFTW_LIBS="-l$mkl_fftw_lib -lmkl_sequential -lmkl_core"
      else
        unset ac_cv_lib_"$mkl_fftw_lib"_"$mkl_func"
        acx_fftw_ok=yes
        ACX_SEARCH_LIBDIR(MKL_FFTW_LIBS, $mkl_root/lib, $mkl_fftw_lib,
		$mkl_func,, [acx_fftw_ok=no],
		[-lmkl_sequential -lmkl_core])
        if test x$acx_fftw_ok = xno; then
          MKL_FFTW_ERROR="INTEL MKL serial FFTW library files not found at usual locations!"
        fi
      fi
    else
      unset ac_cv_lib_"$mkl_fftw_lib"_"$mkl_func"
      ACX_SEARCH_LIBDIR(MKL_FFTW_LIBS, $1/lib $1, $mkl_fftw_lib,
		$mkl_func,, [acx_fftw_ok=no],
		[-lmkl_sequential -lmkl_core])
      if test x$acx_fftw_ok = xno; then
        MKL_FFTW_ERROR="INTEL MKL parallel FFTW library files not found at the provided location!"
      fi
    fi
  fi
fi

LIBS="$OLIBS"
if test x$acx_fftw_ok = xyes; then
  AC_DEFINE(HAVE_FFTW,1, [Define if you have the FFTW libraries.])
  AC_SUBST(MKL_LAPACKE_CFLAGS, $mkl_cflags)
  AC_SUBST(MKL_FFTW_LIBS)
  $5
else
  AC_SUBST(MKL_FFTW_ERROR)
  $6
fi
])dnl ACX_MKL_FFTW
