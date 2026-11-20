# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

REDUCEME is a long-slit spectroscopic data reduction package written in
Fortran 77 (fixed form, compiled with gfortran), with a few C files. Each
program in `src/` is an interactive terminal application that prompts the user
for parameters and uses PGPLOT for graphics; some programs also use CFITSIO.
Documentation (Sphinx, Read the Docs theme) lives in `docs/`.

## Build

Autotools project. Generated files (`configure`, `Makefile.in`, `config.h.in`,
...) are not tracked and must be regenerated after editing `configure.ac` or any
`Makefile.am`:

```
./autogen.sh                      # autoreconf -s -i -f
./configure --program-prefix=R5-  # installed programs get the R5- prefix
make                              # builds libraries and all programs in src/
make install
```

- Array sizes are fixed at configure time: `./configure NSMAX=4096 NCMAX=4096`
  (spatial / spectral maximum dimensions, substituted into `red/redlib.inc`).
- On this Mac (Homebrew) the package has been configured with:
  `./configure FC=gfortran FCFLAGS=-Wall 'LDFLAGS=-L/opt/homebrew/opt/libx11/lib -L/opt/homebrew/opt/pgplot/lib -L/opt/homebrew/opt/cfitsio/lib' --program-prefix=R5-`
- Optional GNU readline (file name completion): `--with-readline[=DIR]`,
  `--without-readline`; Homebrew and MacPorts locations are detected
  automatically on macOS. `configure` reports whether it was enabled.
- Build a single program (libraries must already be built): `make -C src findarc`.
- Check warnings of one file: `gfortran -Wall -I src -fsyntax-only src/findarc.f`
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

Build order (`Makefile.am`): `files button futils fspec red src scripts`.

- `futils/` (`libfutils.a`): terminal input/output utilities. `READC`, `READI`,
  `READF`, `READ2I`, `READILIM` read values with a default (`'@'` = no
  default) and, for `READC`, a string of valid characters (`'@'` = any).
  `READC_FILE` (same interface as `READC`) reads file names with readline
  completion (`freadline.c`, C/Fortran interop via `BIND(C)`) and lists files
  when the answer contains `*` or `?`. Typing `endprogram` at any prompt asks
  to stop the program (`ENDPROGRAM`).
- `red/` (`libred.a`): REDUCEME image I/O and common routines.
  - `INFILEX`/`OUTFILEX` (declared in `iofile.inc`) prompt for, check and open
    files and read/write the REDUCEME header. MODE 1 = REDUCEME unformatted
    format, 2 = unformatted without header, 3 = ASCII, 4 = FITS records
    (`INFILEX` only); MODE+10 opens the given name directly without prompting;
    MODE 21 (`INFILEX`) requires NSCAN/NCHAN/STWV/DISP to match the values
    passed in. Both return `CHARACTER*255`: variables receiving them must be
    `CHARACTER*255` to avoid truncating long paths.
  - `redlib.inc` (generated from `redlib.inc.in`) holds the global parameters
    (`NSMAX`, `NCMAX`, `CLAVE_RED`, `CREDUCEVERSION`) and COMMON variables
    describing the current image (`NSCAN`, `NCHAN`, `STWV`, `DISP`, `AIRMASS`,
    `TIMEXPOS`, `OBJECT`, `FITSFILE`, `COMMENT`, `THISPROGRAM`).
  - `WELCOME(date)` prints the banner; `SHOWHLP` prints help text from
    `$reduceme_dir/help/programs/<program>.hlp`.
- `button/` (`libbutton.a`): PGPLOT graphic buttons (`BUTTON`, `IFBUTTON`, ...).
- `fspec/` (`libfspec.a`): spectral analysis routines (FFT, indices, bands...).
- `src/`: one executable per `.f` file. A new program must be added to
  `bin_PROGRAMS` and given `<name>_SOURCES` in `src/Makefile.am`; programs
  using CFITSIO set `<name>_LDADD = $(ADDCFITSIO)`.
- `files/`: data files (arc line lists `*.ldat`, extinction curves, standard
  star tables) found through the `reduceme_dir` environment variable.
- `scripts/RUN`: csh script that runs a REDUCEME script over several files.
  Special run modes are switched on by the presence of `.running_RUN`,
  `.running_MANUAL` or `.running_HLPHTML` in the current directory (`LRUNX`);
  the `READ*` routines then echo the answers (coloured, LaTeX or HTML).

### Include files

`configure` (`AC_CONFIG_LINKS`) links `futils/futils.inc` into `red/`,
`fspec/` and `src/`, `red/redlib.inc` into `fspec/` and `src/`, and
`red/iofile.inc` into `src/`. Edit the originals (`futils/futils.inc`,
`red/redlib.inc.in`, `red/iofile.inc`), never the links. The version string comes from `AC_INIT` in `configure.ac`
(substituted into `CREDUCEVERSION`, `CHARACTER*12`).

### REDUCEME file format

Fortran unformatted sequential file: record `'abcdefghijkl'` (`CLAVE_RED`),
then `NSCAN,NCHAN`, `STWV,DISP`, `AIRMASS`, `TIMEXPOS`, and length-prefixed
`OBJECT`, `FITSFILE`, `COMMENT`, followed by one record per scan. Error frames
append `' @ERROR@'` to `OBJECT`; `GUESSEF` proposes the error file name by
inserting `e` before the extension (`spec.u` -> `spece.u`).

## Conventions

- Fixed-form source: statements must end before column 72 (use continuation
  lines with `+` in column 6). Some files are Latin-1 encoded; preserve the
  encoding when editing them.
- Each file starts with a header whose second line is
  `C Version DD-Month-YYYY ... file: name.f`, with the `file:` field
  right-aligned at column 79. Programs also pass a version date to
  `CALL WELCOME('DD-Month-YYYY')`. When a file is modified, update both dates
  (English month names, both dates identical).
- No compiler extensions (checked by `make check-std`): use
  `WRITE(*,100,ADVANCE='NO')` with a format without `$` for prompts,
  `GET_ENVIRONMENT_VARIABLE` instead of `GETENV`, `SYSTEMFUNCTION`
  (futils, based on `EXECUTE_COMMAND_LINE`) to run shell commands, `RANRED`
  (fspec) for random numbers, and `INTEGER`/`DOUBLE PRECISION` instead of
  `INTEGER*4`/`REAL*8`.
- Dummy character arguments are declared `CHARACTER*(*)`; keep it that way so
  callers can pass longer strings.
- COMMON blocks must be declared with identical sizes in every program unit
  that uses them (correspondence is by position, not by name).
