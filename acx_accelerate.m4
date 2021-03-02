dnl
dnl				acx_accelerate.m4
dnl
dnl Set up options for using the Mac OSX Accelerate Framework.
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl	This file part of:	AstrOmatic software
dnl
dnl	Copyright:		(C) 2002-2020 IAP/CNRS/SorbonneU
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
dnl	Last modified:		12/08/2020
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_ACCEL([MKL_DIR, ILP64_FLAG, STATIC_FLAG, CONV_LIBS])
dnl
dnl This macro sets the MKL_CFLAGS, MKL_LDFLAGS and MKL_LIBS variables to
dnl for compiling and linking with INTEL's MKL. A coma-separated list of
dnl convenience libraries may be included in the linked group for static linking.
dnl You may wish to use these variables in your default CFLAGS:
dnl
dnl        CFLAGS="$CFLAGS $MKL_CFLAGS"
dnl
dnl You may wish to use these variables in your default LDFLAGS:
dnl
dnl        LDFLAGS="$LDFLAGS $MKL_LDLAGS"
dnl
dnl You may wish to use these variables in your default LIBS:
dnl
dnl        LIBS="$LIBS $MKL_LIBS"
dnl

AC_DEFUN([ACX_ACCEL], [
AC_REQUIRE([AC_CANONICAL_HOST])


dnl --------------------------
dnl Exit if host is not MacOSX
dnl --------------------------
case $host_os in
  darwin* ) ;;
  *)
    ACCEL_WARN="Accelerate only available on Mac OSX"
    AC_SUBST(ACCEL_WARN)
    exit
esac

acx_accelerate_ok=no
AC_CHECK_HEADERS([Accelerate/Accelerate.h], [acx_accelerate_ok=yes])
AC_SUBST(ACCEL_LIBS, "-framework Accelerate")
AC_SUBST(MKL_LDFLAGS, "")

dnl --------------------
dnl Set internal flags
dnl --------------------

AC_DEFINE(HAVE_ACCELERATE,1, [Define if you have the Accelerate libraries.])
AC_DEFINE(HAVE_LAPACK,1, [Define if you have the LAPACK libraries.])

dnl --------------------
dnl Set include files
dnl --------------------

AC_DEFINE(ACCELERATE_H, "Accelerate/Accelerate.h", [Accelerate header filename.])
AC_DEFINE(LAPACK_H, "Accelerate/Accelerate.h", [LAPACK header filename.])

])dnl ACX_ACCEL

