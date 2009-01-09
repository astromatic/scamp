dnl @synopsis ACX_PLPLOT([PLPLOT_LIBDIR, PLPLOT_INCDIR,
dnl                      [ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]]])
dnl This macro figures out if the PlPlot library and header files
dnl are installed.
dnl You may wish to use these variables in your default LIBS and CFLAGS:
dnl
dnl        LIBS="$PLPLOT_LIBS $LIBS"
dnl        CFLAGS="$CFLAGS $PLPLOT_CFLAGS"
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if PlPlot
dnl is found (HAVE_PLPLOT is defined first), and ACTION-IF-NOT-FOUND
dnl is a list of commands to run it if it is not found.
dnl
dnl @version $Id: acx_plplot.m4,v 1.0 2008/08/28 21:30:17 bertin Exp $
dnl @author Emmanuel Bertin <bertin@iap.fr>

AC_DEFUN([ACX_PLPLOT], [
AC_REQUIRE([AC_CANONICAL_HOST])

PLPLOT_LIBS=""
OLIBS="$LIBS"
LIBS=""

acx_plplot_ok=yes
acx_plplotpkg_ok=no
if test x$2 = x && test x$1 = x; then
  AC_MSG_CHECKING([for PLPlot pkg-config info])
  if pkg-config --exists plplotd; then
    AC_MSG_RESULT([yes])
    [PLPLOT_CFLAGS=`pkg-config --cflags plplotd`]
    [PLPLOT_LIBS=`pkg-config --libs plplotd`]
    AC_DEFINE(PLPLOT_H, "plplot.h", [PLPlot header filename.])
    AC_DEFINE(PLPLOTP_H, "plplotP.h", [PLPlot private header filename.])
    acx_plplotpkg_ok=yes
  else
    AC_MSG_RESULT([no])
  fi
fi
if test x$acx_plplotpkg_ok = xno; then
  if test x$2 = x; then
    AC_CHECK_HEADER(plplot.h, [acx_plplothead_ok=yes], [acx_plplothead_ok=no])
    if test x$acx_plplothead_ok = xyes; then
      AC_DEFINE(PLPLOT_H, "plplot.h", [PLPlot header filename.])
      AC_DEFINE(PLPLOTP_H, "plplotP.h", [PLPlot private header filename.])
    else
      AC_CHECK_HEADER(plplot/plplot.h,
		[acx_plplothead_ok=yes], [acx_plplothead_ok=no])
      if test x$acx_plplothead_ok = xyes; then
        AC_DEFINE(PLPLOT_H, "plplot/plplot.h", [PLPlot header filename.])
        AC_DEFINE(PLPLOTP_H, "plplot/plplotP.h",
		[PLPlot private header filename.])
      else
        acx_plplot_ok=no
      fi
    fi
  else
    AC_CHECK_HEADER($2/plplot.h,
		[acx_plplothead_ok=yes], [acx_plplothead_ok=no])
    if test x$acx_plplothead_ok = xyes; then
      AC_DEFINE(PLPLOT_H, "plplot.h", [PLPlot header filename.])
      AC_DEFINE(PLPLOTP_H, "plplotP.h", [PLPlot private header filename.])
     [PLPLOT_CFLAGS="-I$2"]
    else
      acx_plplot_ok=no
    fi
  fi
  if test x$1 = x; then
    AC_CHECK_LIB(plplotd, c_plinit,, [acx_plplot_ok=no])
    [PLPLOT_LIBS="-lplplotd"]
  else
    AC_CHECK_LIB(plplotd, c_plinit,, [acx_plplot_ok=no], [-L$1])
    [PLPLOT_LIBS="-L$1 -lplplotd"]
  fi
fi

LIBS="$OLIBS"
if test x$acx_plplot_ok = xyes; then
  AC_SUBST(PLPLOT_CFLAGS)
  AC_SUBST(PLPLOT_LIBS)
  AC_DEFINE(HAVE_PLPLOT,1,
	[Define if you have the PLPlot libraries and header files.])
  $3
else
  $4
fi

])dnl ACX_PLPLOT

