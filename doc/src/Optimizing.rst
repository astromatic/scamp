.. File Optimizing.rst

Optimizing scamp with ICC
=========================

Using Intel Profile Guided Optimization
---------------------------------------
First compile with -prof-gen activated:
.. code-block:: console

  $ ./configure --enable-mkl CFLAGS="-prof-gen=srcpos prof-dir=/prof/file/dir
  $ make
   
Start the program. This will generate some files:
.. code-block:: console

  $ ./src/scamp myfiles*.cat

Compile again with -prof-use
.. code-block:: console

  $ ./configure --enable-mkl CFLAGS="-prof-use prof-dir=/prof/file/dir


Optimization example
--------------------
The best speed is presently obtained with these flags:
.. code-block:: console

  $ ./configure --enable-mkl CFLAGS=" -w3 -diag-disable:remark -g -O3 -ansi-alias -xCORE-AVX2 -ipo"

