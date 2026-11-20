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

Optionally, if the `GNU readline
<https://tiswww.case.edu/php/chet/readline/rltop.html>`_ library is available,
REDUCEME uses it when asking for file names: pressing ``TAB`` completes the
file name up to the longest unambiguous prefix, pressing it again lists the
possible completions, and previous entries can be recalled with the arrow
keys. If readline is not found, REDUCEME is compiled without this feature and
works as usual. Note that the development files (headers) are required:

- Debian/Ubuntu: ``sudo apt install libreadline-dev``
- Fedora/RHEL: ``sudo dnf install readline-devel``
- macOS with Homebrew: ``brew install readline``
- macOS with MacPorts: ``sudo port install readline``

On macOS without GNU readline, ``configure`` may detect the readline
compatibility layer of the system library ``libedit``, which also provides
file name completion with minor differences in behaviour.

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

.. note:: ``configure`` reports at the end whether file name completion with
   readline has been enabled. On macOS, the Homebrew and MacPorts
   installations of GNU readline are detected automatically. Otherwise, the
   location of the library can be indicated with ``--with-readline=DIR`` (the
   files are searched in ``DIR/include`` and ``DIR/lib``), e.g.:

   ::

      $ ./configure --program-prefix=R5- --with-readline=$HOME/local

   Use ``--without-readline`` to disable this feature. Readline can also
   be configured through the ``~/.inputrc`` file; for example, the line
   ``set show-all-if-ambiguous on`` lists the possible completions with the
   first ``TAB``.

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

