.. File Installing.rst

Installing the software
=======================

Hardware requirements
---------------------

|SCAMP| runs in (ANSI) text-mode from a shell. A window system is
not necessary for basic operation.

The amount of memory required depends mostly on the size of the input
catalogues and on the number of exposures and
astrometric “contexts” (stable instruments) involved in the astrometric
solution. Each detection in the input catalogues amounts to about 140 bytes,
plus a few tens of kbytes for every FITS table. To this one
should add the memory space used by the normal equation matrix, which is
:math:`8 \times N_{\rm T}^2` bytes, with, in the default |SCAMP|
configuration,

.. math::

   N_{\rm T} = N_{\rm ast} \times N_{\rm arr} \times N_{\rm P_{arr}}
       + (N_{\rm exp} - N_{\rm ast}) \times N_{\rm P_{foc}},

where :math:`N_{\rm arr}` is the number of focal plane arrays
(extensions) in each exposure, :math:`N_{\rm P_{arr}}` the number of
polynomial terms describing the static distortion pattern of each array
(20 for a :math:`3^{\rm rd}` in :math:`x` and :math:`y`),
:math:`N_{\rm P_{\rm foc}}` the number of polynomial terms for the
exposure\-dependent focal plane distortion pattern, and :math:`N_{\rm exp}`
the number of exposures. Actually one should probably double the memory space
used by the normal equation matrix to account for buffers in the ATLAS library.
It is not uncommon to see memory usage amounting to gigabytes when
many large mosaic exposures are involved. For instance, computing an astrometric
solution with :math:`N_{\rm P_{arr}}=60` and :math:`N_{\rm P_{\rm foc}}=6` for
a set of 500 exposures of a 60-CCD camera, each with 10,000 detections, spread
over three runs and five bands, may consume as much as 8GB of memory.

Although multiple CPU cores are not required for running |SCAMP|, they can
dramatically reduce execution time, especially when the solution is computed
over a large number of exposures.

Obtaining |SCAMP|
-----------------

For Linux users, the simplest way to have |SCAMP| up and running is to install
the standard binary package the comes with your Linux distribution. Run, e.g.,
``apt-get scamp`` (on Debian) or ``dnf scamp`` (Fedora) and
|SCAMP|, as well as all its dependencies, will automatically be installed. If
you decided to install the package this way you may skip the following and move
straight to the :ref:`next section <using_SCAMP>`.

However if |SCAMP| is not available in your distribution, or to obtain the most
recent version, the |SCAMP| source package can be downloaded from `the official
GitHub repository <https://github.com/astromatic/scamp>`_ . One may choose
`one of the stable releases <https://github.com/astromatic/scamp/releases>`_,
or for the fearless, `a copy of the current master development branch
<https://github.com/astromatic/scamp/archive/master.zip>`_.

Software requirements
---------------------

|SCAMP| has been developed on `GNU/Linux <http://en.wikipedia.org/wiki/Linux>`_
machines and should compile on any
`POSIX <http://en.wikipedia.org/wiki/POSIX>`_-compliant system (this includes
|OSX|_ and `Cygwin <http://www.cygwin.com>`_ on |Windows|_, at the price of
some difficulties with the configuration), provided
that the *development* packages of the following libraries have been installed:

* |ATLAS|_ V3.6 and above [#atlas_install]_,
* |FFTw|_ V3.0 and above [#fftw_install]_, 
* |PLPlot|_ V5.9 and above.

On Fedora/Redhat distributions for instance, the development packages above are
available as ``atlas-devel``, ``fftw-devel`` and ``plplot-devel``.
|PLPlot| is only required for producing diagnostic plots. Note
that |ATLAS| and |FFTw| are not necessary if |SCAMP| is linked with
|Intel|'s |MKL|_ library.

Installation
------------

To install from the |GitHub| source package, you must first uncompress the
archive:

.. code-block:: console

  $ unzip scamp-<version>.zip

A new directory called :file:`scamp-<version>` should now appear at the current
location on your disk. Enter the directory and generate the files required by
the `autotools <http://en.wikipedia.org/wiki/GNU_Build_System>`_, which the
package relies on:

.. code-block:: console

  $ cd scamp-<version>
  $ sh autogen.sh

A :program:`configure` script is created. This script has many options, which
may be listed with the ``--help`` option:

.. code-block:: console

  $ ./configure --help

No options are required for compiling with the default GNU C compiler
(:program:`gcc`) if all the required libraries are installed at their default
locations:

.. code-block:: console

  $ ./configure

Compared to :program:`gcc` and the librairies above, the combination of the
|Intel| compiler (:program:`icc`) and the |MKL|_ libraries can give the |SCAMP|
executable a strong boost in performance, thanks to better vectorized code.
If :program:`icc` and the |MKL| are installed on your system [#geticc]_ , you
can take advantage of them using

.. code-block:: console

  $ ./configure --enable-mkl

Additionally, if the |SCAMP| binary is to be run on a different machine
that does not have :program:`icc` and the |MKL| installed (e.g., a cluster
computing node), you must configure a partially statically linked executable
using

.. code-block:: console

  $ ./configure --enable-mkl --enable-auto-flags --enable-best-link

In all cases, |SCAMP| can now be compiled with

.. code-block:: console

  $ make -j

An :file:`src/scamp` executable is created. For system-wide installation, run
the usual

.. code-block:: console

  $ sudo make install

You may now check that the software is properly installed by simply
typing in your shell:

.. code-block:: console

  $ scamp

which will return the version number and other basic information (note that
some shells require the :program:`rehash` command to be run before making a
freshly installed executable accessible in the execution path).

.. [#mac_install] Mac OS X |.dmg|_ packages should be available soon.
.. [#atlas_install] Use the ``--with-atlas`` and/or
   ``--with-atlas-incdir`` options of the |SCAMP| :command:`configure`
   script to specify the |ATLAS| library and include paths if |ATLAS| files are 
   installed at unusual locations.
.. [#fftw_install] Make sure that |FFTW| has been compiled with
   :command:`configure` options ``--enable-threads --enable-float``.
.. [#geticc] The Linux versions of the |Intel| compiler and |MKL| are
   `available for free to academic researchers, students, educators and open
   source contributors <http://software.intel.com/qualify-for-free-software>`_.

.. include:: keys.rst

