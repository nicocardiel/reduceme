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
   $ ./configure --program-prefix=R6-

.. note:: Mac users can easily indicate a different Fortran compiler using
      ``./configure FC=gfortran-mp-13 CC=gcc-mp-13 --program-prefix=R6-``.

.. note:: The options of the Fortran compiler are set with ``FCFLAGS``,
   which replaces the default ``-g -O2``. If you set it, keep the
   optimisation option, e.g. ``FCFLAGS="-g -O2 -Wall"``.

.. note:: If you find problems at this step detecting PGPLOT, you can help
   ``configure`` by setting the expected location. For example:

   ::

      $ ./configure  --program-prefix=R6- LDFLAGS="-L/opt/local/lib"

.. note:: ``configure`` reports at the end whether file name completion with
   readline has been enabled. On macOS, the Homebrew and MacPorts
   installations of GNU readline are detected automatically. Otherwise, the
   location of the library can be indicated with ``--with-readline=DIR`` (the
   files are searched in ``DIR/include`` and ``DIR/lib``), e.g.:

   ::

      $ ./configure --program-prefix=R6- --with-readline=$HOME/local

   Use ``--without-readline`` to disable this feature. Readline can also
   be configured through the ``~/.inputrc`` file; for example, the line
   ``set show-all-if-ambiguous on`` lists the possible completions with the
   first ``TAB``.

.. warning:: I strongly suggest to use ``--program-prefix=R6-`` in order to add
   a prefix to all the REDUCEME programs. Although it is not strictly
   necessary, it avoids potential name collisions with other software packages.
   It also allows REDUCEME 6 to be installed next to a previous installation
   of REDUCEME 5 (usually installed with the prefix ``R5-``) without
   overwriting its programs.

.. note:: The images are allocated at run time with their actual size, so the
   image dimensions no longer have to be fixed while running ``configure``
   (the ``NCMAX`` and ``NSMAX`` variables of previous versions have been
   removed). The maximum dimensions are 65536 channels and 65536 scans.

3.- Compile the code:

::

   $ make


4.- Install the package in the system (you may need root privileges):

::

   $ sudo make install

Environment
-----------

Before running the programs, define the environment variables:

- ``PGPLOT_DIR``: the directory of the PGPLOT library (every program checks
  that it is defined).
- ``reduceme_dir``: the root directory of the REDUCEME source code, i.e. the
  directory ``reduceme`` created by ``git clone``. The programs read from it
  the help texts (``help/programs/``) and the data files in ``files/`` (arc
  line lists, extinction curves, standard star tables, index definitions).
  ``make install`` does not copy these files, so the variable must point to
  the source directory. If REDUCEME 5 is also installed, point it to the
  REDUCEME 6 source directory when using the ``R6-`` programs.

For example (bash, macOS with Homebrew), if the repository was cloned in
your home directory (``cd $HOME; git clone ...``):

::

   $ export PGPLOT_DIR=/opt/homebrew/opt/pgplot/lib
   $ export reduceme_dir=$HOME/reduceme

In general, run ``pwd`` inside the cloned ``reduceme`` directory and use the
result as the value of ``reduceme_dir``. Add these lines (without the ``$``
prompt) to your shell start-up file (e.g. ``~/.bashrc`` or ``~/.zshrc``) to
define them in every session.

Installing REDUCEME 5.1 (``R5-`` programs)
------------------------------------------

REDUCEME 5.1 is available as the git tag ``5.1``. To install it alongside
REDUCEME 6, clone it in a separate directory (here ``reduceme5``) and use the
prefix ``R5-``:

::

   $ git clone --branch 5.1 https://github.com/nicocardiel/reduceme reduceme5
   $ cd reduceme5
   $ ./autogen.sh
   $ ./configure --program-prefix=R5-
   $ make
   $ sudo make install

The installed files do not overlap with those of REDUCEME 6: only the
programs and the scripts ``RUN`` and ``exheads.sh`` are installed, all of
them with their prefix. The prefix must be exactly ``R5-``, since the script
``RUN`` of version 5.1 calls ``R5-exheads.sh`` by that name. Note the
differences in the ``configure`` options of version 5.1:

- The Fortran compiler is selected with ``F77`` and its options with
  ``FFLAGS`` (e.g. ``./configure F77=gfortran --program-prefix=R5-``).
- The maximum image dimensions are fixed at compilation time (2400x2400 by
  default); larger images require, for example,
  ``./configure --program-prefix=R5- NCMAX=4096 NSMAX=4096``.
- There are no regression tests (``make check``).

When using the ``R5-`` programs, the environment variable ``reduceme_dir``
must point to the REDUCEME 5.1 source directory (e.g.
``export reduceme_dir=$HOME/reduceme5``), since they read the help texts and
data files from it.

