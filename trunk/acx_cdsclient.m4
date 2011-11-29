dnl
dnl				acx_cdsclient.m4
dnl
dnl Figure out if the CDSclient package is installed.
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
dnl	Last modified:		29/11/2011
dnl
dnl %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dnl @synopsis ACX_CDSCLIENT([CDSCLIENT,[ACTION-IF-FOUND
dnl                         [, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the CDSClient programs are installed.
dnl ACTION-IF-FOUND is a list of shell commands to run if CDSClient
dnl is found (HAVE_CDSCLIENT is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.

AC_DEFUN([ACX_CDSCLIENT], [
AC_REQUIRE([AC_CANONICAL_HOST])

# We only need to check the availability of the "aclient_cgi" executable
if test x$1 = x; then
  AC_CHECK_PROG(acx_cdsclient_ok, [aclient_cgi], [yes], [no])
else
  AC_CHECK_PROG(acx_cdsclient_ok, [aclient_cgi], [yes], [no], $1)
fi

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_cdsclient_ok" = xyes; then
        AC_DEFINE(HAVE_CDSCLIENT,1,
        [Define if you have the CDSclient tools installed.])
        $2
else
        $3
fi

AC_SUBST(CDSCLIENT)

])dnl ACX_CDSCLIENT
