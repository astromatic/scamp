#
#				configure.ac
#
# Process this file with autoconf to generate a configure script.
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#	This file part of:	SCAMP
#
#	Copyright:		(C) 2002-2011 Emmanuel Bertin -- IAP/CNRS/UPMC
#
#	License:		GNU General Public License
#
#	SCAMP is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
# 	(at your option) any later version.
#	SCAMP is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	You should have received a copy of the GNU General Public License
#	along with SCAMP.  If not, see <http://www.gnu.org/licenses/>.
#
#	Last modified:		17/12/2011
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# First, disable the annoying config.cache
define([AC_CACHE_LOAD],)
define([AC_CACHE_SAVE],)

# This is your standard Bertin source code...
AC_INIT([scamp],[1.11.1],[bertin@iap.fr])
AC_CONFIG_SRCDIR(src/makeit.c)
AC_CONFIG_AUX_DIR(autoconf)
AC_CONFIG_HEADERS(config.h)
AM_INIT_AUTOMAKE
date=`date +%Y-%m-%d`
date2=`date +"%a %b %d %Y"`
date3=`date +"%B %Y"`
AC_DEFINE_UNQUOTED(DATE, "$date", [Archive creation date])
AC_SUBST(PACKAGER, "Emmanuel Bertin")
AC_SUBST(DATE2, "$date2")
AC_SUBST(DATE3, "$date3")

# Include macros
sinclude(acx_search_libdir.m4)
sinclude(acx_atlas.m4)
sinclude(acx_cdsclient.m4)
sinclude(acx_fftw.m4)
sinclude(acx_mkl_fftw.m4)
sinclude(acx_mkl_lapacke.m4)
sinclude(acx_plplot.m4)
sinclude(acx_prog_cc_optim.m4)
sinclude(acx_pthread.m4)
sinclude(acx_urbi_resolve_dir.m4)

# Display pakage and version number
AC_MSG_RESULT([******** Configuring:  $PACKAGE_NAME $VERSION ($date) ********])

AC_MSG_CHECKING([if compilation flags should be set automatically])
AC_ARG_ENABLE(auto-flags,
	[AS_HELP_STRING([--enable-auto-flags],
	[Let the configure script choose the compilation flags \
(off by default)])],
	use_autoflag="yes"
	AC_MSG_RESULT([yes]),
	use_autoflag="no"
	AC_MSG_RESULT([no]))

# Provide special option for the Linux Intel C compiler
AC_MSG_CHECKING([whether we should use the INTEL compiler])
AC_ARG_ENABLE(icc,
	[AS_HELP_STRING([--enable-icc],
	[Use the Intel compiler (off by default)])],
        use_icc="yes"
        CC="icc"
	AC_MSG_RESULT([yes]),
        use_icc="no"
	AC_MSG_RESULT([no]))

# Provide special options for INTEL MKL
# We force the use of icc
AC_MSG_CHECKING([whether we should use INTEL's MKL])
AC_ARG_ENABLE(mkl,
	[AS_HELP_STRING([--enable-mkl],
	[Use INTEL's MKL for solvers and FFTs (off by default)])],
        use_mkl="yes"
        use_icc="yes"
        CC="icc"
	AC_MSG_RESULT([yes]),
        use_mkl="no"
	AC_MSG_RESULT([no]))

# Checks for programs.
AC_LANG(C)

AC_SEARCH_LIBS([strerror],[cposix])
if test "$use_autoflag" = "yes"; then
  CFLAGS=""
  LDFLAGS=""
  ACX_PROG_CC_OPTIM
fi
LT_INIT
AC_PROG_INSTALL

# Checks for libraries.
AC_CHECK_LIB(m, sin)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdlib.h string.h unistd.h sys/mman.h dlfcn.h])
if test $use_icc = "yes" -a $CC = "icc"; then
  AC_CHECK_HEADERS(mathimf.h)
fi

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_LONG_LONG_INT
AC_TYPE_UNSIGNED_LONG_LONG_INT
AC_STRUCT_TM
AC_TYPE_UID_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MMAP
AC_FUNC_STAT
AC_FUNC_STRFTIME
AC_CHECK_FUNCS([atexit getenv memcpy memmove memset mkdir munmap strstr \
sysconf setlinebuf dlopen dlsym dlclose])
AC_CHECK_FUNC([isnan], AC_DEFINE_UNQUOTED([HAVE_ISNAN2], 1,
		[Second isnan check]))

# Check support for large files
AC_SYS_LARGEFILE
AC_FUNC_FSEEKO

datadir2=$(URBI_RESOLVE_DIR([$datadir]))

# Provide special options for the INTEL MKL
AC_ARG_WITH(mkl-libdir,
        [AS_HELP_STRING([--with-mkl-libdir=<MKL library path>],
	[Alternative directory for the MKL library])],
    mkl_fftw_libdir=$withval
    mkl_lapacke_libdir=$withval,
    mkl_fftw_libdir=""
    mkl_lapacke_libdir=""
    )
AC_ARG_WITH(mkl-incdir,
        [AS_HELP_STRING([--with-mkl-incdir=<MKL include dir>],
	[Alternative directory for MKL include files])],
    fftw_incdir=$withval
    lapacke_incdir=$withval,
    fftw_incdir=""
    lapacke_incdir=""
    )

# Provide special options for ATLAS
AC_ARG_WITH(atlas-libdir,
        [AS_HELP_STRING([--with-atlas-libdir=<ATLAS library path>],
	[Alternative directory for the ATLAS library])],
    atlas_libdir=$withval,
    atlas_libdir=""
    )
AC_ARG_WITH(atlas-incdir,
        [AS_HELP_STRING([--with-atlas-incdir=<ATLAS include dir>],
	[Alternative directory for ATLAS include files])],
    atlas_incdir=$withval,
    atlas_incdir=""
    )

# Provide special options for FFTW
AC_ARG_WITH(fftw-libdir,
        [AS_HELP_STRING([--with-fftw-libdir=<FFTW library path>],
	[Alternative directory for the FFTW library])],
    fftw_libdir=$withval,
    fftw_libdir=""
    )
AC_ARG_WITH(fftw-incdir,
        [AS_HELP_STRING([--with-fftw-incdir=<FFTW include dir>],
	[Alternative directory for FFTW include files])],
    fftw_incdir=$withval,
    fftw_incdir=""
    )

# Provide a special option for CDSClient
AC_ARG_WITH(cdsclient-dir,
	[AS_HELP_STRING([--with-cdsclient-dir=<aclient execution path>],
	[Alternative default to the aclient execution path])],
    cdsclient_path="$withval",
    cdsclient_path=""
    )

# Provide special options for PLPlot
AC_MSG_CHECKING([whether we should use PLPlot])
use_plplot=yes
AC_ARG_ENABLE(plplot,
	[AS_HELP_STRING([--disable-plplot], [Disable PLPlot (on by default)])],
      if test "$enableval" = "no"; then
        use_plplot="no"
        AC_MSG_RESULT([no])
      else
        use_plplot="yes"
        AC_MSG_RESULT([yes])
      fi ,
      use_plplot="yes"
      AC_MSG_RESULT([yes])
    )
AC_ARG_WITH(plplot-libdir,
    [AS_HELP_STRING([--with-plplot-libdir=<PLPlot library dir>],
	[Alternative directory for the PLPlot library])],
      plplot_libdir=$withval,
      plplot_libdir=""
    )
AC_ARG_WITH(plplot-incdir,
        [AS_HELP_STRING([--with-plplot-incdir=<PLPlot include dir>],
	[Alternative directory for PLPlot include files])],
    plplot_incdir=$withval,
    plplot_incdir=""
    )

# Provide a special option for the default XSLT URL
AC_ARG_WITH(xsl_url,
	[AS_HELP_STRING([--with-xsl_url=<default URL for XSLT filter>],
	[Alternative default URL for the XSLT filter])],
    [xsl_url=$withval],
    [xsl_url="file://$datadir2/$PACKAGE_NAME/$PACKAGE_NAME.xsl"]
    )

AC_DEFINE_UNQUOTED([XSL_URL], "$xsl_url",[Default URL of the XSLT filter])

# Set flags for multithreading
n_pthreads=16
AC_ARG_ENABLE(threads,
	[AS_HELP_STRING([--enable-threads@<:@=<max_number_of_threads>@:>@],
	[Enable multhreading (on with up to 16 threads by default)])],
    if test "$enableval" = "no"; then
      use_pthreads="no"
    else
      use_pthreads="yes"
      if test "$enableval" != "yes"; then
        n_pthreads=$enableval
      fi
    fi,
    use_pthreads=yes
    )

# Provide special option for profiling
AC_MSG_CHECKING([for profiler mode])
AC_ARG_ENABLE(profiling,
	[AS_HELP_STRING([--enable-profiling],
	[Enable special mode for profiling (off by default)])],
        use_prof="yes"
	AC_MSG_RESULT([yes]),
        use_prof="no"
	AC_MSG_RESULT([no]))

# Enable linking options for making the executable as portable as possible.
AC_MSG_CHECKING([best linking option])
AC_ARG_ENABLE(best-link,
	[AS_HELP_STRING([--enable-best-link],
	[Choose the right combination of static and dynamic linking to make \
the executable as portable as possible (off by default)])],
        use_best="yes"
	AC_MSG_RESULT([yes]),
        use_best="no"
	AC_MSG_RESULT([no]))

################# Actions to complete in case of multhreading ################
AC_DEFINE_UNQUOTED(THREADS_NMAX, $n_pthreads,[Maximum number of POSIX threads])
if test "$use_pthreads" = "yes"; then
  AC_MSG_CHECKING([for multithreading])
  AC_MSG_RESULT([maximum of $n_pthreads thread(s)])
  AC_DEFINE(USE_THREADS, 1, [Triggers multhreading])
# CC, CFLAGS and LIBS are system and compiler-dependent
  ACX_PTHREAD
  CC="$PTHREAD_CC"
  [AM_CFLAGS="$AM_CFLAGS $PTHREAD_CFLAGS -D_REENTRANT"]
  LIBS="$PTHREAD_LIBS $LIBS"
fi
AM_CONDITIONAL(USE_THREADS, test $use_pthreads = "yes")

############ handle the INTEL MKL library (linear algebra+FFT) ###########
if test "$use_mkl" = "yes"; then
  ACX_MKL_LAPACKE($mkl_lapacke_libdir,$mkl_lapacke_incdir,$use_pthreads,
	[use_mkl_lapacke=yes],[use_mkl_lapacke=no])
  if test "$use_mkl_lapacke" = "yes"; then
    LIBS="$MKL_LAPACKE_LIBS $LIBS"
   [AM_CFLAGS="$AM_CFLAGS $MKL_LAPACKE_CFLAGS"]
    if test "$MKL_LAPACKE_WARN" != ""; then
      AC_MSG_WARN([$MKL_LAPACKE_WARN])
    fi
  else
    AC_MSG_ERROR([$MKL_LAPACKE_ERROR Exiting.])
  fi
  ACX_MKL_FFTW($mkl_fftw_libdir,$mkl_fftw_incdir,$use_pthreads,yes,
	[use_mkl_fftw=yes],[use_mkl_fftw=no])
  if test "$use_mkl_fftw" = "yes"; then
#    LIBS="$FFTW_LIBS $LIBS"
#   [AM_CFLAGS="$AM_CFLAGS $MKL_FFTW_CFLAGS"]
    if test "$MKL_FFTW_WARN" != ""; then
      AC_MSG_WARN([$MKL_FFTW_WARN])
    fi
  else
    AC_MSG_ERROR([$MKL_FFTW_ERROR Exiting.])
  fi
else
################ handle the FFTW library (Fourier transforms) ################
  ACX_FFTW($fftw_libdir,$fftw_incdir,$use_pthreads,yes,
	[use_fftw=yes],[use_fftw=no])
  if test "$use_fftw" = "yes"; then
    LIBS="$FFTW_LIBS $LIBS"
    if test "$FFTW_WARN" != ""; then
      AC_MSG_WARN([$FFTW_WARN])
    fi
  else
    AC_MSG_ERROR([$FFTW_ERROR Exiting.])
  fi

################## handle the ATLAS library(linear algebra) ##################
  ACX_ATLAS($atlas_libdir,$atlas_incdir,$use_pthreads,
	[use_atlas=yes],[use_atlas=no])
  if test "$use_atlas" = "yes"; then
    LIBS="$ATLAS_LIBS $LIBS"
    if test "$ATLAS_WARN" != ""; then
      AC_MSG_WARN([$ATLAS_WARN])
    fi
  else
    AC_MSG_ERROR([$ATLAS_ERROR Exiting.])
  fi
fi

################# handle the CDSclient package (CDS queries) #################
ACX_CDSCLIENT($cdsclient_path,[use_cdsclient=yes],[use_cdsclient=no])
if test x$cdsclient_path = x; then
  cdsclient="aclient_cgi"
else
  cdsclient="$cdsclient_path/aclient_cgi"
fi
if test "$use_cdsclient" = "no"; then
  AC_MSG_ERROR([$cdsclient not found. exiting.]
["aclient_cgi" from the CDSclient package is required to use SCAMP:]
[1. Install the CDSclient package: http://cdsweb.u-strasbg.fr/doc/cdsclient.html]
[2. Use the --with-cdsclient=<aclient_path> configure option if aclient_cgi is
   not in your execution PATH.])
fi
AC_DEFINE_UNQUOTED(CDSCLIENT, "$cdsclient", [CDSclient executable])

################# handle the PLPlot library (graphic plots) ##################
AC_MSG_CHECKING([PLPlot configure options])
if test "$use_plplot" = "yes"; then
  AC_MSG_RESULT([yes])
  ACX_PLPLOT($plplot_libdir,$plplot_incdir,,[use_plplot=no])
  if test "$use_plplot" = "yes"; then
    [AM_CPPFLAGS="$AM_CPPFLAGS $PLPLOT_CFLAGS"]
    [LIBS="$PLPLOT_LIBS $LIBS"]
  else
    AC_MSG_WARN(PLPlot library not found/incomplete! Plotting deactivated.)
  fi
else
  AC_MSG_RESULT([no])
  AC_MSG_WARN(Plotting deactivated.)
fi
AM_CONDITIONAL(USE_PLPLOT, test $use_plplot = "yes")

# Link with gprof option
if test "$use_prof" = "yes"; then
  if test "$use_icc" = "yes"; then
    AM_CFLAGS="$AM_CFLAGS -pq"
  else
    AM_CFLAGS="$AM_CFLAGS -pg"
  fi
  use_best="no"
fi

# "Best" linking option
if test "$use_best" = "yes"; then
  AM_LDFLAGS="-shared-libgcc -static-libtool-libs -avoid-version $AM_LDFLAGS"
fi

AC_SUBST(AM_CFLAGS)
AC_SUBST(AM_CPPFLAGS)
AC_SUBST(AM_LDFLAGS)

# Display compiler and linker flags
AC_MSG_RESULT([***************************************************************])
AC_MSG_RESULT([Compile cmdline:  $CC $AM_CPPFLAGS $CPPFLAGS $AM_CFLAGS $CFLAGS])
AC_MSG_RESULT([Link    cmdline:  $CC $AM_LDFLAGS $LDFLAGS $LIBS])
AC_MSG_RESULT([Default XSLT URL: $xsl_url])
AC_MSG_RESULT([***************************************************************])

AC_CONFIG_FILES([
Makefile
man/Makefile
src/Makefile
src/fits/Makefile
src/wcs/Makefile
scamp.spec
man/scamp.1])
AC_OUTPUT
