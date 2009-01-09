dnl @synopsis ACX_CDSCLIENT([CDSCLIENT,[ACTION-IF-FOUND
dnl                         [, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the CDSClient programs are installed.
dnl ACTION-IF-FOUND is a list of shell commands to run if CDSClient
dnl is found (HAVE_CDSCLIENT is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.
dnl
dnl @version $Id: acx_cdsclient.m4,v 1.0 2004/06/02 21:30:17 bertin Exp $
dnl @author Emmanuel Bertin <bertin@iap.fr>

AC_DEFUN([ACX_CDSCLIENT], [
AC_REQUIRE([AC_CANONICAL_HOST])

# We only need to check the availability of the "aclient" executable
AC_CHECK_PROG(acx_cdsclient_ok, $1, [yes], [no])

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
