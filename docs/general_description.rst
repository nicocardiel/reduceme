General description
===================

REDUCEME is an astronomical data reduction package, specially devoted to the
analysis of long-slit spectroscopic data. This software was created by N.
Cardiel as part of his thesis work, developed under the supervision of J.
Gorgas, at the Departamento de Astrofísica of the Universidad Complutense de
Madrid.

The initial help of S. Pedraz, J. Cenarro and O. Alonso was very important in
finding bugs and suggesting modifications, which have improved most of the
programs. Since then, more people have also contributed with comments and
suggestions (A. Gil de Paz, E. García-Dabó, P. Sánchez-Blázquez, E.
Mármol-Queraltó, and E. Toloba). Many thanks to all of them.  

*I am making this software freely available because I think some of the
programs may be useful for other people, but I am far from providing strong
support (if any!) since I am currently very busy with my teaching duties and
with additional research activities.  So, please, if you are planning to use
it, be aware of that.*

This package does not intend to be a general and complete spectroscopic data
reduction package, and it has been created strongly biased to solve the needs
we have encountered during the reduction of our long-slit spectroscopic data.
Unfortunately, and due to historical reasons, this software does not use FITS
as the data format, but the unformatted Fortran raw format (which, in old
times, used to be a very fast way to access files with slow computers). This
means that any potential user should transform the FITS files to REDUCEME
format. The reverse operation (from REDUCEME to FITS format) is also available.

This package consists in a set of programs written in Fortran 77, and also
includes some shell scripts (using the C shell syntax) to perform repetitive
tasks. This document contains an user description of the programs, giving some
guidelines to the reader who wants to extend their capabilities, with the
inclusion of own external programs.

Graphics (line plots and images) in this package are done with the help of the
excelent library `PGPLOT <http://www.astro.caltech.edu/~tjp/pgplot/>`_. A
subset of subroutines, called `button <http://button.readthedocs.io/>`_,
have been specially written to enable the user to communicate interactively
with the image display employing graphic buttons.

Error handling
--------------

One of the most interesting advantages of using REDUCEME programs is that for
each image an associated error image can also be processed throughout the
reduction process, allowing for a careful control of the error propagation. A
more detailed description of this technique, and its application to the
measurement of line-strength indices, is given in these papers:

* `Cardiel et al. (1998)
  <http://cdsads.u-strasbg.fr/abs/1998A%26AS..127..597C>`_: 
  Reliable random error estimation in the measurement of line-strength indices.
* `Cardiel et al. (2002) 
  <http://cdsads.u-strasbg.fr/abs/2002SPIE.4847..297C>`_: 
  Proper handling of random errors and distortions in astronomical data 
  analysis
* `Cardiel et al. (2003a)
  <http://cdsads.u-strasbg.fr/abs/2003RMxAC..16...73C>`_:
  A New Approach in Data Reduction: Proper Handling of Random Errors and Image 
  Distortions
* `Cardiel et al. (2003b)
  <http://cdsads.u-strasbg.fr/abs/2003A%26A...409..511C>`_:
  Using spectroscopic data to disentangle stellar population properties

.. note:: You can also find a general description of REDUCEME in Chapter 3 of
   my :download:`PhD Thesis<aux/thesis_ncardiel.pdf>` (see also appendix B and
   C).

.. warning:: Since there is no explicit publication describing this package,
   you can acknowledge the use of REDUCEME by citing Cardiel (1999): Cardiel,
   N., 1999, Ph.D. thesis, Universidad Complutense de Madrid

In addition, you can find :download:`here<aux/logo_reduceme.tex>` how to create
the REDUCEME logo in LaTeX.




