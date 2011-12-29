dnl                             acx_search_incdir.m4
dnl
dnl Search an include file in a list of directories
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
dnl @synopsis ACX_SEARCH_INCDIR(VAR, EMPTYFLAG, INCDIRS, INCLUDE,
dnl			ACTION-IF-FOUND, ACTION-IF-NOT-FOUND)
dnl
dnl This macro is a modification of AX_EXT_HAVE_LIB by Duncan Simpson
dnl <dps@simpson.demon.co.uk>

AC_DEFUN([ACX_SEARCH_INCDIR], [
AC_REQUIRE([AC_CANONICAL_HOST])
AC_LANG_PUSH(C)
  hdr=`echo $4 | $as_tr_sh`
  incdir_found="no"
  $1=""
  if test x$2=xyes; then
    AC_CHECK_HEADER($4,[incdir_found="yes"])
  fi
  for incdir in $3; do
    if test "x${incdir_found}" = "xno"; then
      incdir_hashdr_cvdir=`echo $incdir | $as_tr_sh`
      AC_CACHE_CHECK([for $4 with -I$incdir],
	[incdir_cv${incdir_hashdr_cvdir}_hashdr_${hdr}],
	[incdir_have_hdr_save_cflags=${CFLAGS}
	CFLAGS="${CFLAGS} -I${incdir}"
	AC_COMPILE_IFELSE(
	[AC_LANG_PROGRAM([#include <$4>])],
	[incdir_found="yes"; eval "incdir_cv${incdir_hashdr_cvdir}_hashdr_${hdr}"="yes"],
	[incdir_found="no"; eval "incdir_cv${incdir_hashdr_cvdir}_hashdr_${hdr}"="no"])
	CFLAGS=$incdir_have_hdr_save_cflags])
      if eval `echo 'test x${'incdir_cv${incdir_hashdr_cvdir}_hashdr_${hdr}'}' = "xyes"`; then
        incdir_found="yes";
        $1="${incdir}"
        hdr=`echo $4 | $as_tr_cpp`
        AC_DEFINE_UNQUOTED(HAVE_${hdr}, 1,
		[Define this if you have the $1 header])
      fi
    fi
  done
AC_LANG_POP
if test x$incdir_found = xyes; then :
  $5
else :
  $6
fi
])

