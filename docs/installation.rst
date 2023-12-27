Installation
============

Requirements
------------

Before you install REDUCEME, make sure that PGPLOT and CFITSIO are
already installed in your system. If this is not the case, you must download
`PGPLOT <http://www.astro.caltech.edu/~tjp/pgplot/>`_, and
`CFITSIO <http://heasarc.gsfc.nasa.gov/fitsio/>`_.

Some details about how I do typically install PGPLOT and CFITSIO under
``Linux`` and ``Mac OS X`` are given `here for PGPLOT
<https://guaix.fis.ucm.es/~ncl/howto/howto-pgplot>`_, and
`here for CFITSIO
<https://guaix.fis.ucm.es/~ncl/howto/howto-cfitsio>`_.

REDUCEME Installation
---------------------

To install REDUCEME you need to perform the following steps:


1.- Download the latest distribution from github:

::

    $ git clone https://github.com/nicocardiel/reduceme

2.- Enter into the directory ``reduceme`` and prepare the code to be compiled

::

   $ autoreconf -s -i -f
   $ ./configure --program-prefix=R5-

.. note:: Mac users can easily indicate a different Fortran compiler using
      ``./configure F77=gfortran-mp-13 CC=gcc-mp-13 --program-prefix=R5-``.

.. note:: If you find problems at this step detecting PGPLOT, you can help
   ``configure`` by setting the expected location. For example:

   ::

      $ ./configure  --program-prefix=R5- LDFLAGS="-L/opt/local/lib"

.. warning:: I strongly suggest to use ``--program-prefix=R5-`` in order to add
   a prefix to all the REDUCEME programs. Although it is not strictly
   necessary, it avoids potential name collisions with other software packages.

.. warning:: Since Fortran 77 statically declares the dimensions of the arrays 
   at compilation time, you may need to declare the maximum size of the
   expected arrays while running ``configure``:
   
   ::

      $ ./configure --program-prefix=R5- NCMAX=4096 NSMAX=4096

3.- Compile the code:

::

   $ make


4.- Install the package in the system (you may need root privileges):

::

   $ sudo make install

