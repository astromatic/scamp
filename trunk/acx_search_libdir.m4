dnl                             acx_search_libdir.m4
dnl
dnl Search a library in a list of directories
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl     This file part of:      AstrOmatic software
dnl
dnl     Copyright:              (C) 2011 Emmanuel Bertin -- IAP/CNRS/UPMC
dnl                             (C) 2008 Duncan Simpson (original version)
dnl
dnl     Licenses:               GPL (this version)
dnl                             "Copying and distribution of this file, with or
dnl                             without modification, are permitted in any
dnl                             medium without royalty provided the copyright
dnl                             notice and this notice are preserved. This file
dnl                             is offered as-is, without any warranty."
dnl                             (original script)
dnl
dnl     AstrOmatic software is free software: you can redistribute it and/or
dnl     modify it under the terms of the GNU General Public License as
dnl     published by the Free Software Foundation, either version 3 of the
dnl     License, or (at your option) any later version.
dnl     AstrOmatic software is distributed in the hope that it will be useful,
dnl     but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl     GNU General Public License for more details.
dnl     You should have received a copy of the GNU General Public License
dnl     along with AstrOmatic software.
dnl     If not, see <http://www.gnu.org/licenses/>.
dnl
dnl     Last modified:          29/12/2011
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dnl
dnl @synopsis ACX_SEARCH_LIBDIR(VAR, EMPTYFLAG, LIBDIRS, LIBRARY, FUNCTION,
dnl			ACTION-IF-FOUND, ACTION-IF-NOT-FOUND,
dnl			EXTRA_LIBRARIES)
dnl
dnl This macro is a modification of AX_EXT_HAVE_LIB by Duncan Simpson
dnl <dps@simpson.demon.co.uk>

AC_DEFUN([ACX_SEARCH_LIBDIR], [
AC_REQUIRE([AC_CANONICAL_HOST])

libdir_ldflags=${LDFLAGS}
old_libs=$LIBS
libdir_found="no"
$1=""
if test x$2 = xyes; then
  AC_CHECK_LIB($4,$5,
	[libdir_found="yes"
	$1="-l$4 $8"],
	[],$8)
fi
for dir in $3; do
  if test $libdir_found = no && test -d ${dir}; then
    libdir_haslib_cvdir=`echo $dir | $as_tr_sh`
    AC_CACHE_CHECK([for $4 library with -L$dir],
	[libdir_cv${libdir_haslib_cvdir}_haslib_$4],
	[libdir_save_LIBS=$LIBS
	libdir_save_ldflags=${LDFLAGS}
	LIBS="-l$4 $8"
	LDFLAGS="-L$dir ${libdir_save_ldflags}"
	AC_TRY_LINK_FUNC([$5],
		[eval "libdir_cv${libdir_haslib_cvdir}_haslib_$4"="yes"],
		[eval "libdir_cv${libdir_haslib_cvdir}_haslib_$4"="no"])
	LIBS=$libdir_save_LIBS
	LDFLAGS=$libdir_save_ldflags])
    if eval `echo 'test x${'libdir_cv${libdir_haslib_cvdir}_haslib_$4'}' = "xyes"`; then
      $1="-L${dir} -l$4 $8"
      libdir_found="yes"
    fi
  fi
done
LDFLAGS="$libdir_ldflags"
LIBS="$old_libs"
if test x$libdir_found = xyes; then :
  $6
else :
  $7
fi
])
