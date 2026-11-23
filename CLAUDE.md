# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

REDUCEME is a long-slit spectroscopic data reduction package written in
Fortran (free-form `.f90` sources, still mostly Fortran 77 style, compiled
with gfortran), with a few C files. Each
program in `src/` is an interactive terminal application that prompts the user
for parameters and uses PGPLOT for graphics; some programs also use CFITSIO.
Documentation (Sphinx, Read the Docs theme) lives in `docs/`.

The code is being modernised in phases (Fortran 77 -> modern Fortran):
read `MODERNIZATION.md` for the plan, the current status, pending work and
known bugs before working on it, and update it at the end of each phase.

## Build

Autotools project. Generated files (`configure`, `Makefile.in`, `config.h.in`,
...) are not tracked and must be regenerated after editing `configure.ac` or any
`Makefile.am`:

```
./autogen.sh                      # autoreconf -s -i -f
./configure --program-prefix=R6-  # installed programs get the R6- prefix
make                              # builds libraries and all programs in src/
make install
```

- Image sizes are not fixed at build time (images are allocatable); the
  maximum dimensions `NSMAX`/`NCMAX` are 65536 (`futils/redlib_m.f90.in`).
- On this Mac (Homebrew) the package has been configured with:
  `./configure FC=gfortran 'LDFLAGS=-L/opt/homebrew/opt/libx11/lib -L/opt/homebrew/opt/pgplot/lib -L/opt/homebrew/opt/cfitsio/lib' --program-prefix=R6-`
- The Fortran compiler is selected with `FC` and its options with `FCFLAGS`,
  which replaces the default `-g -O2`: keep `-O2` (e.g.
  `FCFLAGS="-g -O2 -Wall"`), since without optimisation the floating-point
  rounding changes and some regression tests fail.
- Optional GNU readline (file name completion): `--with-readline[=DIR]`,
  `--without-readline`; Homebrew and MacPorts locations are detected
  automatically on macOS. `configure` reports whether it was enabled.
- Build a single program (libraries must already be built): `make -C src findarc`.
- Check warnings of one file: `gfortran -Wall -I src -fsyntax-only src/findarc.f90`
  (`-Wcharacter-truncation` is relevant; the many `1/255` warnings from
  `X(1:1)=READC(...)` single-character answers are expected).
- Standard conformance: `make check-std` compiles every Fortran file with
  `-std=f2018 -Werror=intrinsics-std` (gfortran) and fails if any compiler
  extension is used; obsolescent features (COMMON, `CHARACTER*n`...) are
  allowed.
- Documentation: `cd docs && make html` (requirements in `docs/requirements.txt`).

## Tests

Regression tests live in `tests/` (`test/` and `specials/` are unrelated,
unbuilt Fortran sources). They run programs on synthetic data and compare
stdout and output files with reference results from version 5.1 (git tag
`5.1`); see `tests/README`.

```
make check                                   # whole suite (needs python3)
python3 tests/run_tests.py irevx istat       # selected cases
python3 tests/run_tests.py --keep irevx      # keep the temporary directory
```

- A case (`tests/cases/<name>.case`) lists the program and the answers fed to
  its prompts on stdin; inputs come from `tests/mkdata.py`.
- References are created with `--update`, only with binaries built from the
  reference version, never from modified code.
- Programs that use the graphic cursor or buttons cannot be automated; they
  are listed in `tests/MANUAL_CHECKS`.
- Programs read from stdin, so they can also be driven with a pipe or a
  pseudo-terminal (readline is only used when stdin/stdout are terminals).

## Architecture

Build order (`Makefile.am`): `files futils button fspec red src scripts
tests`; each library depends only on the previous ones (futils <- button,
futils <- fspec <- red), and the programs in `src/` use all of them.

Each library is a module (`futils_m`, `button_m`, `fspec_m`, `red_m`): the
file `<lib>/<lib>_m.f90` includes, after `CONTAINS`, the files with the
procedures, which keep one file per routine. To add a routine, create its
file, add an `INCLUDE` to the module file and add the file to
`<LIB>_ROUTINES` in `<lib>/Makefile.am`. Callers get explicit interfaces
through `USE` (never declare the type of a library function nor use
`EXTERNAL` for it). `make` does not track module dependencies of the
programs: after changing a procedure interface, run `make clean`.

- `futils/` (`libfutils.a`, module `futils_m`): terminal input/output
  utilities. `READC`, `READI`, `READF`, `READ2I`, `READILIM` read values with a
  default (`'@'` = no default) and, for `READC`, a string of valid characters
  (`'@'` = any). `READC_FILE` (same interface as `READC`) reads file names with
  readline completion (`freadline.c`, C/Fortran interop via `BIND(C)`) and
  lists files when the answer contains `*` or `?`. Typing `endprogram` at any
  prompt asks to stop the program (`ENDPROGRAM`). It also compiles the module
  `redlib_m` (see below).
- `button/` (`libbutton.a`, `button_m`): PGPLOT graphic buttons (`BUTTON`,
  `IFBUTTON`, ...); their state is kept in private module variables.
- `fspec/` (`libfspec.a`, `fspec_m`): spectral analysis routines (FFT,
  fitting, indices, bands...).
- `red/` (`libred.a`, `red_m`): REDUCEME image I/O and common routines.
  - `INFILEX`/`OUTFILEX` prompt for, check and open files and read/write the
    REDUCEME header. MODE 1 = REDUCEME unformatted format, 2 = unformatted
    without header, 3 = ASCII, 4 = FITS records (`INFILEX` only); MODE+10
    opens the given name directly without prompting; MODE 21 (`INFILEX`)
    requires NSCAN/NCHAN/STWV/DISP to match the values passed in. Both return
    `CHARACTER*255`: variables receiving them must be `CHARACTER*255` to avoid
    truncating long paths.
  - `WELCOME(date)` prints the banner; `SHOWHLP` prints help text from
    `$reduceme_dir/help/programs/<program>.hlp`.
- `src/`: one executable per `.f90` file. A new program must be added to
  `bin_PROGRAMS` and given `<name>_SOURCES` in `src/Makefile.am`; programs
  using CFITSIO set `<name>_LDADD = $(ADDCFITSIO)`.
- `files/`: data files (arc line lists `*.ldat`, extinction curves, standard
  star tables) found through the `reduceme_dir` environment variable.
- `scripts/RUN`: csh script that runs a REDUCEME script over several files.
  Special run modes are switched on by the presence of `.running_RUN`,
  `.running_MANUAL` or `.running_HLPHTML` in the current directory (`LRUNX`);
  the `READ*` routines then echo the answers (coloured, LaTeX or HTML).
  `RUN` and `exheads.sh` call other installed programs with the prefix
  taken from their own installed name (`R6-RUN` -> `R6-exheads.sh`); never
  write a fixed prefix in them.

### Global data and include files

- Module `redlib_m` (`futils/redlib_m.f90`, generated by `configure` from
  `futils/redlib_m.f90.in`; edit the `.in` file): global parameters (`NSMAX`,
  `NCMAX`, `NMAXFFT`, `CLAVE_RED`, and `CREDUCEVERSION` from `AC_INIT` in
  `configure.ac`)
  and the global image header variables (`AIRMASS`, `TIMEXPOS`, `OBJECT`,
  `FITSFILE`, `COMMENT`, `THISPROGRAM`).
- `red/redlib.inc` declares `NSCAN`, `NCHAN`, `STWV`, `DISP`, which are
  local variables of every unit that includes it (together with
  `USE REDLIB_M`); `configure` links it into `fspec/` and `src/` (edit the
  original in `red/`). `CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)` only
  silences unused-variable warnings.
- `red/pgplot.inc` reproduces PGPLOT internal COMMON blocks (used by
  `PGIDEN_RED`).

### REDUCEME file format

Fortran unformatted sequential file: record `'abcdefghijkl'` (`CLAVE_RED`),
then `NSCAN,NCHAN`, `STWV,DISP`, `AIRMASS`, `TIMEXPOS`, and length-prefixed
`OBJECT`, `FITSFILE`, `COMMENT`, followed by one record per scan. Error frames
append `' @ERROR@'` to `OBJECT`; `GUESSEF` proposes the error file name by
inserting `e` before the extension (`spec.u` -> `spece.u`).

## Conventions

- Free-form source (`.f90`, including the `.inc` files), converted from fixed
  form keeping its layout: code starts in column 9, a continued line ends
  with `&` and the continuation line starts with `&` in column 6, comments
  use `!`. Some files are Latin-1 encoded; preserve the encoding when
  editing them.
- Each file starts with a header whose second line is
  `! Version DD-Month-YYYY ... file: name.f90`, with the `file:` field
  right-aligned at column 79. Programs also pass a version date to
  `CALL WELCOME('DD-Month-YYYY')`. When a file is modified, update both dates
  (English month names, both dates identical).
- No compiler extensions (checked by `make check-std`): use
  `WRITE(*,100,ADVANCE='NO')` with a format without `$` for prompts,
  `GET_ENVIRONMENT_VARIABLE` instead of `GETENV`, `SYSTEMFUNCTION`
  (futils, based on `EXECUTE_COMMAND_LINE`) to run shell commands, `RANRED`
  (fspec) for random numbers, and `INTEGER`/`DOUBLE PRECISION` instead of
  `INTEGER*4`/`REAL*8`.
- Image-size arrays (phase 4 of MODERNIZATION.md): declare them
  `REAL, ALLOCATABLE :: A(:,:)`, allocate them with the actual dimensions
  `A(NCHAN,NSCAN)` after reading the header and set them to zero (static
  arrays were zero in practice and some programs rely on it). Library
  routines receive images as assumed-shape arguments (`A(:,:)`). Images
  that a program shared with its internal procedures through COMMON go
  to a data module defined in the program file (`<PROGRAM>_DATA_M`).
  `NSMAX`/`NCMAX` (65536) are only limits: 1-D work arrays (one spectrum,
  one scan profile) may be dimensioned with them, but arrays with another
  large dimension must be allocatable, and they must not be used as a
  number of points (e.g. to draw a curve).
- Dummy character arguments are declared `CHARACTER*(*)`; keep it that way so
  callers can pass longer strings.
- COMMON blocks must be declared with identical sizes in every program unit
  that uses them (correspondence is by position, not by name).
