.. File Using.rst

.. _using_SCAMP:

Using SCAMP
===========

|SCAMP| is run from the shell with the following syntax:

.. code-block:: console

 $ scamp Catalog1 [Catalog2 ...] -c configuration-file [-Parameter1 Value1 -Parameter2 Value2 ...]

The parts enclosed within brackets are optional. The file names of input
catalogues can be directly provided in the command line, or in lists that are
ASCII files with each catalogue name preceded by ``@`` (one per line). One
should use lists instead of the catalogue file names if the number of input
catalogues is too large to be handled directly by the shell. Any
`-Parameter Value` statement in the command-line overrides the corresponding
definition in the configuration file or any default value (see below).

Input files
-----------

Catalogues
~~~~~~~~~~

Catalogue files read by |SCAMP| must be in |SExtractor|_'s "FITS_LDAC" binary
format. It is strongly advised to use |SExtractor| version 2.4.4 or later.
The catalogues *must* contain all the following measurements in order to
be processable by |SCAMP|:

* Centroid coordinates. They must be specified with the ``CENTROID_KEYS``
  configuration parameter (default: ``XWIN_IMAGE`` and ``YWIN_IMAGE``).

* Centroid error ellipse, as defined by the ``CENTROIDERR_KEYS``
  configuration parameter (default: ``ERRAWIN_IMAGE``, ``ERRBWIN_IMAGE`` and
  ``ERRTHETAWIN_IMAGE``).

* Factors controlling astrometric distortion. These are set with the
  ``DISTORT_KEYS`` configuration parameter (default: ``XWIN_IMAGE`` and
  ``YWIN_IMAGE``).

* Flux measurements, defined by the ``PHOTFLUX_KEY`` configuration parameter
  (default: ``FLUX_AUTO``).

-  Flux uncertainties, defined by the ``PHOTFLUXERR_KEY`` configuration
   parameter (default: ``FLUXERR_AUTO``).

In addition, it is strongly advised (but not mandatory) to include the following
optional |SExtractor| measurements:

* ``FLAGS``, ``FLAGS_WEIGHT`` and/or ``IMAFLAGS_ISO`` flags for filtering out
  blended and corrupted detections.

* the ``FLUX_RADIUS`` half-light radius estimation for filtering out both small
  glitches and extended objects.

* ``ELONGATION`` to help filtering out objects such as trails and diffraction
  spikes

* a measurement of the object "spread" compared to that of the PSF model:    
  ``SPREAD_MODEL`` and its uncertainty ``SPREADERR_MODEL``. These measurements
  are not used by |SCAMP| for selection, but they get propagated to the output
  catalogues.

:file:`.ahead` header files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The binary catalogues in "FITS_LDAC" format read by |SCAMP| contain a
copy of the original FITS image headers. These headers provide
fundamental information such as frame dimensions, World Coordinate
System (|WCS|) data and many other FITS keywords which |SCAMP| uses to
derive a full astrometric and photometric calibration. It is often
needed to change or add keywords in some headers. Editing FITS files is
not convenient, so |SCAMP| provides read (and write) support for
"external" header files. External headers may either be real FITS header
cards (no carriage-return), or ASCII files containing lines in FITS-like
format, with the final line starting with the ``END     `` keyword.
Multiple extensions must be separated by an ``END     `` line. External
"headers" need not contain all the FITS keywords normally
required. The keywords present in external headers are only there to
override their counterparts in the original image headers or to add new
ones.

Hence for every input (say, :file:`xxxx.cat`) FITS catalogue, |SCAMP| looks for
a :file:`xxxx.ahead` header file, loads it if present, and overrides or adds to
image header keywords those found there. ``.ahead`` is the default suffix;
it can be changed using the ``AHEADER_SUFFIX`` configuration parameter.

It is often useful to add/modify the same FITS keyword(s) in *all*
input catalogues. |SCAMP| offers the possibility to put these keywords in
one single external header file, which is read before all other
:file:`.ahead` files (but after reading the catalogue headers). The name of this
file is ``scamp.ahead`` by default; it can be changed using the
``AHEADER_GLOBAL`` configuration parameter. **show an example of a typical
.ahead file**

Output files
------------

:file:`.head` header files
~~~~~~~~~~~~~~~~~~~~~~~~~~

|SCAMP| itself generates FITS header keywords, containing updated
astrometric and photometric information. These keywords are written in
ASCII to external header files, with the :file:`.head` filename extension by
default (the suffix can be changed with the ``HEADER_SUFFIX`` configuration
parameter). In combination with the original image files, these .head
headers are ready to be used by the |SWarp|_ image stacking tool
:cite:`2002ASPC..281..228B`.

The astrometric engine at the heart of |SCAMP| and |SWarp| is based on
Mark Calabretta’s |WCSLIB|_ library
:cite:`2002A&A...395.1061G,2002A&A...395.1077C`, to which we added
support for the |TPV|_ description of polynomial distortions [#about_tpv]_. All 
celestial coordinate computations are performed in the equatorial system,
although galactic and ecliptic coordinates are supported in input and output.

Output catalogues
~~~~~~~~~~~~~~~~~

|SCAMP| can save three kinds of catalogues: local copies of the reference
catalogues downloaded from the Vizier server (see §[chap:astref]), a
“merged”, calibrated version of input catalogues (§[chap:mergedcat]),
and a “full” calibrated version of input a reference catalogues (§[chap:fullcat]).

Diagnostic files
~~~~~~~~~~~~~~~~

Two types of files can be generated by |SCAMP|, providing diagnostics
about the calibrations:

* *Check-plots* are graphic charts generated by |SCAMP|, showing scatter
  plots or calibration maps. The ``CHECKPLOT_TYPE`` and ``CHECKPLOT_NAME``
  configuration parameters allow the user to provide a list of
  check-plot types and file names, respectively. A variety of raster
  and vector file formats, from JPEG to Postscript, can be set with
  ``CHECKPLOT_DEV``. PNG is the default. See the
  :ref:`CHECKPLOT <param_checkplot>` section of
  the :ref:`configuration parameter list <param_list>` below for details.

* An |XML|_ file providing a processing summary and various statistics in
  |VOTable|_ format is written if the ``WRITE_XML`` switch is set to ``Y`` (the
  default). The ``XML_NAME`` parameter can be used to change the default
  file name :file:`scamp.xml`. The |XML| file can be displayed with any recent
  web browser; the |XSLT|_ stylesheet installed together with |SCAMP| will
  automatically translate it into a dynamic, user-friendly web-page.
  For more advanced usages (e.g., access from a remote web server),
  alternative |XSLT| translation URLs may be specified using the ``XSL_URL``
  configuration parameter.

The Configuration file
----------------------

Each time it is run, |SCAMP| looks for a configuration file. If no
configuration file is specified in the command-line, it is assumed to be
called :file:`scamp.conf` and to reside in the current directory. If no
configuration file is found, |SCAMP| uses its own internal default
configuration.

Creating a configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

|SCAMP| can generate an ASCII dump of its internal default configuration,
using the ``-d`` option. By redirecting the standard output of |SCAMP| to a
file, one creates a configuration file that can easily be modified
afterward:

.. code-block:: console

$ scamp -d >scamp.conf

A more extensive dump with less commonly used parameters can be
generated by using the ``-dd`` option.

Format of the configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The format is ASCII. There must be only one parameter set per line,
following the form:

     *Config-parameter      Value(s)*

Extra spaces or linefeeds are ignored. Comments must begin with a ``#``
and end with a linefeed. Values can be of different types: strings (can
be enclosed between double quotes), floats, integers, keywords or
Boolean (`Y/y` or `N/n`). Some parameters accept zero or several values,
which must then be separated by commas. Values separated by commas,
spaces, tabs or linefeeds may also be read from an ASCII file if what is
given is a filename preceded with ``@`` (e.g. `@values.txt`). Integers can be
given as decimals, in octal form (preceded by digit O), or in
hexadecimal (preceded by `0x`). The hexadecimal format is particularly
convenient for writing multiplexed bit values such as binary masks.
Environment variables, written as ``$HOME`` or ``${HOME}`` are expanded.

.. _param_list:

Configuration parameter list
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here is a list of all the parameters known to |SCAMP|. Please refer to
next section for a detailed description of their meaning. Some
“advanced” parameters (indicated with an asterisk) are also listed. They
must be used with caution, and may be rescoped or removed without notice
in future versions.

.. [#about_tpv] The |TPV| description was originally proposed by E. Greisen
and M. Calabretta in a 2000 draft, but abandoned in later 
versions :cite:`2004ASPC..314..551C`. Following adoption in |SCAMP| and in
large data processing centers it eventually `became a registered FITS convention
in 2012 <http://fits.gsfc.nasa.gov/registry/tpvwcs.html>`_, and is now included
in recent versions of the |WCSLIB|.


.. include:: keys.rst

