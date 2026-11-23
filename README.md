# REDUCEME

REDUCEME is an astronomical data reduction package, specially devoted to the
analysis of long-slit spectroscopic data. This software was created by N.
Cardiel as part of his thesis work, developed under the supervision of J.
Gorgas, at the Departamento de Astrofísica of the Universidad Complutense de
Madrid.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.

- Maintainer: Nicolás Cardiel, cardiel@ucm.es
- Webpage (documentation): http://reduceme.readthedocs.io/en/latest/
- Webpage (source): https://github.com/nicocardiel/reduceme

## From REDUCEME 5 (R5) to REDUCEME 6 (R6)

REDUCEME 6 is a modernisation of the code of REDUCEME 5.1. The programs ask
the same questions in the same order and read and write the same files, so
the way of working does not change. The main differences are:

- **No image size limits fixed at compilation time.** In REDUCEME 5 the
  images were static arrays whose maximum size (by default 2400x2400) had to
  be chosen when running `configure` (`NSMAX=... NCMAX=...`). REDUCEME 6
  allocates each image with its actual size when it is read; the maximum
  dimensions are 65536 scans and 65536 channels, and the `NSMAX`/`NCMAX`
  options of `configure` no longer exist.
- **Modern Fortran.** The sources were converted from Fortran 77 (fixed form)
  to free-form Fortran (`.f90`), the libraries are Fortran modules (the
  compiler now checks the arguments of every call), and the compiler
  extensions were replaced by standard Fortran 2018. The Fortran compiler is
  selected with `FC` and its options with `FCFLAGS` (instead of `F77` and
  `FFLAGS`).
- **Regression tests.** `make check` runs the programs on synthetic data and
  compares the results with those of REDUCEME 5.1.
- **Bug fixes.** Several bugs of version 5.1 found during the modernisation
  have been fixed; they are listed in `MODERNIZATION.md`.
- **File format.** The REDUCEME file format has not changed: files created
  with REDUCEME 5 can be used with REDUCEME 6 and vice versa.
- **Program names.** The recommended prefix of the installed programs is now
  `R6-` (it was `R5-`), so that both versions can be installed at the same
  time without overwriting each other (`R5-istat` and `R6-istat`, for
  example). Scripts written for REDUCEME 5 only need the new prefix in the
  program names. `R6-RUN` and `R6-exheads.sh` call the programs with their
  own prefix.

## Installation

### Requirements

- A Fortran compiler (gfortran) and a C compiler.
- [PGPLOT](http://www.astro.caltech.edu/~tjp/pgplot/) and
  [CFITSIO](http://heasarc.gsfc.nasa.gov/fitsio/).
- GNU autotools (`autoconf`, `automake`) to generate `configure`.
- Optionally, [GNU readline](https://tiswww.case.edu/php/chet/readline/rltop.html)
  (with its development files) for file name completion with `TAB`; if it is
  not found, REDUCEME is compiled without this feature.
- Python 3, only to run the regression tests (`make check`).

### Steps

```
git clone https://github.com/nicocardiel/reduceme
cd reduceme
./autogen.sh
./configure --program-prefix=R6-
make
make check            # optional: regression tests
sudo make install
```

`configure` accepts the usual variables, for example:

- `FC=gfortran-mp-13 CC=gcc-mp-13` to choose the compilers;
- `FCFLAGS="-g -O2 -Wall"` to change the compiler options (keep `-O2`);
- `LDFLAGS="-L..."` when the libraries are not in a standard location;
- `--with-readline=DIR` or `--without-readline`.

On macOS with Homebrew, for example:

```
./configure FC=gfortran \
  'LDFLAGS=-L/opt/homebrew/opt/libx11/lib -L/opt/homebrew/opt/pgplot/lib -L/opt/homebrew/opt/cfitsio/lib' \
  --program-prefix=R6-
```

### Environment

Before running the programs, define the environment variables:

- `PGPLOT_DIR`: the directory of the PGPLOT library (every program checks
  that it is defined).
- `reduceme_dir`: the root directory of the REDUCEME source code, i.e. the
  directory `reduceme` created by `git clone`. The programs read from it the
  help texts (`help/programs/`) and the data files in `files/` (arc line
  lists, extinction curves, standard star tables, index definitions).
  `make install` does not copy these files, so the variable must point to
  the source directory. If REDUCEME 5 is also installed, point it to the
  REDUCEME 6 source directory when using the `R6-` programs.

For example (bash, macOS with Homebrew), if the repository was cloned in
your home directory (`cd $HOME; git clone ...`):

```
export PGPLOT_DIR=/opt/homebrew/opt/pgplot/lib
export reduceme_dir=$HOME/reduceme
```

In general, run `pwd` inside the cloned `reduceme` directory and use the
result as the value of `reduceme_dir`. Add these lines to your shell start-up
file (e.g. `~/.bashrc` or `~/.zshrc`) to define them in every session.

See the [documentation](http://reduceme.readthedocs.io/en/latest/) for more
details.
