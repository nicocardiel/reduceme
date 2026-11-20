# Modernisation of REDUCEME: status and plan

This file records the plan to move REDUCEME from Fortran 77 to modern
Fortran, what has been done so far, and the pending work, so that the work
can be resumed in any session or machine. Keep it up to date at the end of
each phase.

## Goals and approach

Goals:

1. Remove the image size limits: images are static arrays `A(NCMAX,NSMAX)`
   whose size is fixed at configure time (default 2400x2400).
2. Maintainability: modules instead of `COMMON`/`INCLUDE`, explicit
   interfaces, free-form source.
3. Long-term portability: no compiler extensions.

Approach: gradual, in phases. At the end of each phase the package must
build and behave exactly as before:

- `make check`: regression tests (`tests/`), whose references were
  generated with REDUCEME 5.1 (git tag `5.1`, commit `1562f4d`).
  References are never regenerated from modified code.
- `make check-std` (from phase 1 on): no compiler extensions
  (`gfortran -std=f2018 -Werror=intrinsics-std`).
- One branch per phase, small commits (one pattern or group of programs
  per commit). Before starting each phase, its scope is reviewed and agreed
  with the maintainer.
- Version dates (file header and `WELCOME` call) are updated in every
  modified file (see CLAUDE.md).

## Status

| Phase | Description | Status |
|---|---|---|
| 0 | Regression tests (`make check`) | Done, merged in `master` (`711f48e`) |
| 1 | Remove compiler extensions (`make check-std`) | Done, merged in `master` (`da89e10`..`0475390`) |
| 2 | Free-form source (`.f` -> `.f90`) | Done, merged in `master` |
| 3 | Modules in the libraries | Done on branch `phase3-modules` |
| 4 | Dynamic arrays (remove `NSMAX`/`NCMAX`) | Not started |
| 5 | Structured control flow, internal COMMON blocks | Ongoing/opportunistic, not started |

## Baseline (measured before phase 1)

- 227 `.f` files, ~65,700 code lines (`src/`: 104 programs, ~55,500).
- 1213 `GOTO`, 947 `COMMON` statements (235 distinct blocks), 700
  references to `NSMAX`/`NCMAX`, ~4000 PGPLOT calls, 139 CFITSIO calls.
- 85 of the 104 programs do not use the graphic cursor or buttons and can
  be driven through stdin; the other 19 are listed in `tests/MANUAL_CHECKS`.

## Done

### Phase 0: regression tests

- `tests/mkdata.py` writes synthetic REDUCEME files directly (gfortran
  sequential unformatted layout), independently of the code under test.
- `tests/run_tests.py` runs each case (`tests/cases/*.case`: program and
  answers fed to stdin, optional command-line arguments, multi-program
  cases) with PGPLOT device `/NULL`, and compares stdout and output files:
  numbers in text with relative tolerance 1e-6, REDUCEME files byte by
  byte, FITS files ignoring date cards. gfortran runtime error locations
  and backtraces are normalised.
- 32 cases (27 in phase 0, 5 added at the start of phase 1).
- The references were generated on macOS arm64 with gfortran 16.2
  (Homebrew). On another architecture or compiler version, floating-point
  results may differ in the last bit and the byte-by-byte comparison of
  REDUCEME files may fail although the code is correct: in that case,
  regenerate the references there from tag `5.1` (see `tests/README`) or
  add a numerical tolerance to the comparison of REDUCEME files.

### Phase 1: no compiler extensions

- `RANRED` (fspec) reimplemented in standard Fortran with the same
  algorithm as the GNU `RAND`/`SRAND` (Park-Miller, multiplier 16807;
  value `((x-1)/2**9)/2**22`; seed 0 -> 123459876): identical sequences,
  checked bit by bit for 1.4 million values. Only the clock seed
  (`NSEED<0`) now comes from `SYSTEM_CLOCK` instead of `TIME`.
- `GETENV` -> `GET_ENVIRONMENT_VARIABLE`; `SYSTEM` -> `EXECUTE_COMMAND_LINE`
  (in `SYSTEMFUNCTION`, with `CMDSTAT`); `IARGC`/`GETARG` ->
  `COMMAND_ARGUMENT_COUNT`/`GET_COMMAND_ARGUMENT`.
- `INTEGER*4`/`REAL*8` -> `INTEGER`/`DOUBLE PRECISION`; `DATA` for COMMON
  variables (cleanest, fitshead) -> array assignments.
- `$` edit descriptor -> `ADVANCE='NO'` in 3388 `WRITE` statements (228
  `FORMAT` labels plus inline formats); converted by a script working on
  whole statements and program units, and checked statement by statement.
- Remaining `-std=f2018` diagnostics are warnings about obsolescent
  features only (COMMON 3698, `CHARACTER*n` 2384, labelled DO 30, DATA
  after executable statements 3): they belong to later phases.

### Phase 2: free-form source

- Two steps, so that git keeps the history of every file: first the files
  were only renamed (`git mv x.f x.f90`, compiled with `-ffixed-form`),
  then their content was converted.
- Conversion with an ad hoc script (not findent) keeping a 1:1
  correspondence of lines: `C`/`c`/`*` comments -> `!` (`Comment` ->
  `!Comment`); continuation marker in column 6 -> `&` plus a trailing `&`
  on the previous line; lines ending inside a character constant and
  continued FORMAT statements padded with blanks to column 72 (their
  content includes those blanks in fixed form); `STOP'x'` -> `STOP 'x'`.
  The include files (`*.inc`, `red/redlib.inc.in`) were converted too;
  `test/` and `specials/` (not built) were left unchanged.
- Verification: each file compiled (`-O2`, same file name, so that the
  embedded line numbers and file name coincide) before and after the
  conversion: the 227 object files are identical byte by byte. Same
  `-Wall` diagnostics; `make check` and `make check-std` pass.
- Build system: `AC_PROG_F77` -> `AC_PROG_FC`; the compiler is selected
  with `FC` and its options with `FCFLAGS`, which replaces the default
  `-g -O2`. Keep `-O2` when setting `FCFLAGS` (e.g. `"-g -O2 -Wall"`):
  without optimisation the floating-point rounding changes and some
  regression tests fail (documented in CLAUDE.md, tests/README and
  docs/installation.rst).
- creahelp now accepts the free-form header markers (`!Comment`,
  `! Program:`...) and the `.f90` extension; its output for the program
  headers is identical to that obtained from the fixed-form sources.
- Version date and file name (`file: name.f90`) updated in every header.

### Phase 3: modules in the libraries

- Preparation: dependency cycle red <-> fspec removed (AVOID_WARNINGS moved
  to futils; FFTCORRZOOM and FFTPREP, plotting routines, moved to red);
  build order futils, button, fspec, red. Duplicated procedure names
  (tolerated by the linker in static libraries) renamed: the integer
  `ORDENA1F` of ordena1i.f90 -> `ORDENA1I`, `FUNKGC` of gauscfit_movel ->
  `FUNKGC_MOVEL`, `FPOLY` of polfitsig -> `FPOLY_SIG` (the link map
  confirmed which version each program used), and the program versions of
  `FMEAN2` (basicred, istat), `SEXAG` (fitshead) and `TRUELEN` (rvel) ->
  `*_LOCAL`. Only change of behaviour: in rvel, library routines now use the
  library `TRUELEN` (it differs only for strings ending in control
  characters).
- Module `REDLIB_M` (futils/redlib_m.f90.in, generated by configure):
  global parameters and image header variables (formerly COMMON
  /BLKRED01-06/). `red/redlib.inc` only declares the local variables
  NSCAN, NCHAN, STWV, DISP.
- One module per library (`futils_m`, `button_m`, `fspec_m`, `red_m`) that
  includes the procedure files after CONTAINS. Removed: 414 declarations
  of library functions, 11 EXTERNAL statements, futils.inc, iofile.inc,
  readcbutt.inc, button.inc (its COMMON blocks are now private variables of
  `button_m`). `READC`/`READCBUTT` return `CHARACTER*255` (assumed-length
  character functions cannot be module procedures).
- Explicit interfaces revealed no argument mismatch between the procedures
  of the package (checked beforehand compiling every program together with
  all the libraries).
- Verification: make check, make check-std (now compiling the modules in
  order and then the programs), same -Wall diagnostics, interactive check.
  Object code cannot be compared (module procedure names are mangled).

## Next phases (agreed scope; review before starting each one)

### Phase 4: dynamic arrays

- Library routines first: the 14 files (11 in `fspec/`, 3 in `red/`) with
  arrays dimensioned by `NCMAX` (cauchyfit, broaden, chrebin, gauscfit,
  gauss2afit, gaussfit, mideind, gauss2bfit, pseudofit, rvrebin, splfit;
  fftcorrzoom, infilex, zscale): dimension by argument or `ALLOCATABLE`; arrays
  shared with callers through COMMON go to a module.
- Then programs one by one (extend their tests first):
  `REAL A(NCMAX,NSMAX)` -> `REAL, ALLOCATABLE :: A(:,:)` allocated as
  `A(NCHAN,NSCAN)` after reading the header (keep the (channel, scan)
  index order); limits and loops use the actual sizes; large arrays in
  COMMON go to a module; no large automatic arrays (stack).
- Watch: programs reading several images (`INFILEX` mode 21), FFT
  programs (`NMAXFFT`), cleanest and plots (size, many COMMON blocks).
- Finally remove `NSMAX`/`NCMAX` from `configure.ac`, `red/redlib.inc.in`,
  `docs/installation.rst` and CLAUDE.md, and test an image larger than the
  old limit (e.g. 4096x4096) without recompiling.

### Phase 5: opportunistic

`GOTO` -> `DO`/`EXIT`/`CYCLE`/`SELECT CASE`, internal COMMON blocks ->
modules or arguments, only when a program is modified for another reason
and has tests covering it. No mass restructuring.

### Out of scope

- PGPLOT replacement (e.g. the PGPLOT-compatible library giza; compatibility
  to be verified): it is the main long-term portability risk but is not
  solved by modernising the Fortran code.
- The REDUCEME file format must stay bit-compatible.
- Rewriting in another language.

## Other pending work

- File name completion (`READC_FILE`) in the programs that still read file
  names with `READ(*,'(A)')`: genimage (`FILELINES`), mcolumns
  (`NEWFILEIN`), movel (`ERRFILE1`, `OUTFILE`), multpol (`POLFILE`,
  `OUTFILE`), offindex (`INFILE`), replacefits (`FILENAME`, `FILECOPY`),
  writefits and writefits_log (`OUTFILE`, `TEMPFILE`). Review case by
  case: an empty answer may have a meaning (e.g. movel uses it to guess
  the error file), and some programs echo the answer themselves in the
  LRUN/manual/HTML modes (`READC_FILE` already does it).
- Truncation warnings not related to file names (left on purpose):
  ~400 `1/255` from `X(1:1)=READC(...)` (single-character answers,
  harmless), labels in fitdata and creahelp, `OBJECT0` in midez.

## Known bugs (present in REDUCEME 5.1, not fixed yet)

- fitshead, outfitshead: gfortran runtime error "Sequential READ or WRITE
  not allowed after EOF marker" when appending to the output file (they
  write after reading it to the end without `BACKSPACE`). The test case
  `fitshead_known_error` documents the current behaviour; update its
  reference when this is fixed.
- writefits: the year is written as `**` in the `HISTORY Date:` card (the
  year does not fit in its two-digit field).
- `red/my_pgend.f`, `futils/endprogram.f`: shell commands written as
  `'\\rm ...\0'`; since gfortran does not process backslashes, the
  temporary files of the GIF captures are probably not removed (not
  verified).
- `READ*` routines (futils): when an answer is rejected and there is no
  default value, the prompt `? ` is printed twice (`? ? `).
- End of file on stdin makes all the `READ*` routines stop with a Fortran
  runtime error.
- `READC` accepts answers longer than one character when all their
  characters are valid (e.g. `ny` for a y/n question is taken as `n`).
- Non-conforming calls to PGPLOT (they work in practice; found in phase 3
  by compiling all the sources together): `CALL PGBAND(...)` although
  PGBAND is an INTEGER function (several programs); `PGSCR` called with
  integer constants where REAL is expected (e.g. `CALL PGSCR(3,0,0.5,0)`);
  `PGPOINT` called with scalar expressions where arrays are expected.
