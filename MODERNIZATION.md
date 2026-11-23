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
| 4 | Dynamic arrays (remove `NSMAX`/`NCMAX`) | Done (subphases 4.1-4.5; 4.5 on branch `phase4-5-limits`, visual checks pending) |
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

Agreed design ("two levels"): arrays whose size grows with the image
(2-D `(NCMAX,NSMAX)`, square `(NSMAX,NSMAX)`, flattened `NSMAX*NCMAX`)
become allocatable with the actual image size; 1-D arrays stay static.
`NSMAX`/`NCMAX` keep their configure value during the transition (static
2-D arrays still exist in unconverted programs) and become a fixed 1-D
limit of 65536 in subphase 4.5, after checking that no static 2-D array
remains; then the configure options are removed. Subphases, each in its
own branch and merged separately: 4.1 libraries, 4.2 simple programs
(2-D arrays only local to the main unit), 4.3 complex programs without
cursor, 4.4 interactive programs (visual checks by the maintainer), 4.5
limits, configure and documentation. Extend the tests of each program
before converting it.

Progress and notes:

- 4.1 (done): `ZSCALE` takes an assumed-shape image (valid for static and
  allocatable callers) and allocates its auxiliary array with the size of
  the region (checked bit by bit against the previous version).
- 4.2 (done): the 23 simple programs (2-D arrays only local to the main
  unit) converted. Pattern: all the image-size arrays of a program are
  allocated together, with the actual dimensions, right after reading
  the image header (or after asking the output dimensions), and set to
  zero (`.FALSE.` for logical masks), as the previous static arrays were
  in practice; a previous allocation is released first. Special cases:
  outputs that grow one scan/channel at a time (adnsc, adnch) start with
  the input size and double when full (`MOVE_ALLOC`); outputs larger than
  the input (resample, resample_flux) are enlarged keeping the contents;
  rotate allocates the result with swapped dimensions; genimage allocates
  its auxiliary image at least as large as the output one. 17 test cases
  added (49 in total), with references generated from a build of tag 5.1
  in a git worktree, covering every converted program and the new code
  paths (growth, resampling to more channels, rotations, genimage
  options).
- 4.3 (done): the 23 programs without cursor whose image arrays are not
  only local to the main unit, or need special handling (index, snratio
  and `MIDEIND` are left for 4.4, see below). Group A (12 programs with
  local arrays): basicred, corrext, fit1dpol, gain, interp, istat,
  leefits_new, miderot, nortwi, replacefits, writefits, writefits_log;
  flattened `PIXEL` arrays (basicred, istat) allocated with
  `NCHAN*NSCAN`, and the CFITSIO calls `FTP2DE`/`FTP2DJ` receive the
  actual first dimension instead of `NCMAX`. Group B (11 programs):
  images shared through COMMON with the internal procedures (autocos,
  autocos2, fit2dpol, leefits, portable, prfcal, sdistor) moved to a data
  module defined in the program file (`<PROGRAM>_DATA_M`, before the main
  program); arrays whose second index is not a scan (midez templates,
  prfcal spectra, autocos files) allocated with their actual number of
  columns; sdistor's `PGGRAY` receives the actual dimensions. 20 test
  cases added (68 in total), with new synthetic inputs (exposures with
  cosmic rays, an arc with curved lines, template stars and a galaxy with
  known velocities, an S-distorted point source); references from tag
  5.1. Besides `make check`, the suite was run with a build using
  `-fcheck=bounds`: no out-of-bounds access in the converted programs
  (fit2dpol differs in the 6th digit in that build also before the
  conversion: different rounding); it revealed the `PGIDEN_RED` and
  `FFTPREP` problems listed in the known bugs.
- 4.4 (done): the 16 interactive programs and
  `MIDEIND`; no static image-size array remains in the package. Group A:
  index, snratio and `MIDEIND`, whose shared image moves to the data
  module `INDEX_DATA_M` (`fspec/index_data_m.f90`, apart from `FSPEC_M`
  so that `A` and `ERR` are only visible where it is used); snratio's
  `PLOTNS` becomes `(NSCAN,NSCAN)`. Group B (local arrays): adnhand,
  interlines, midelines, multfit, plotsp3d, skysubm, vaucoul. Group C
  (COMMON -> `<PROGRAM>_DATA_M`): fillimage, interpmove, leeposs, testwc,
  fit2dspl. Group D: cleanest (static local arrays that keep their values
  between calls, such as the UNDO copies, become `ALLOCATABLE, SAVE` and
  are allocated in the first call) and plots (buffers enlarged by
  `GROWBUFF` when a file needs room for more scans). All the PGPLOT calls
  with image arrays receive the actual dimensions. 7 test cases added (74
  in total) for the parts that do not need the cursor (index, snratio,
  adnhand up to its known error, skysubm, fit2dspl). The suite was also
  run with `-fcheck=bounds`, and every program was started with that
  build and the `/NULL` device. Fixed on the way: snratio output changed
  by phase 1 (restored), skysubm sky errors (5.1 bug, reference
  regenerated), pixel value under the cursor outside the image (fillimage,
  snratio). Visual checks by the maintainer, comparing with 5.1 and with
  the `-fcheck=bounds` build: all the interactive programs except leeposs
  (needs a Digital Sky Survey FITS image; not checked by hand) and index
  (no cursor); they found the skysubm option `-3` bug (fixed).
- 4.5 (done, visual checks pending): `NSMAX` and `NCMAX` fixed to 65536
  (the same as `NMAXFFT`) in `futils/redlib_m.f90.in`, and the configure
  variables removed; they are now only the maximum image dimensions
  checked by `INFILEX` and the size of the 1-D work arrays. A build with
  65536 revealed what had to be changed first:
  - `CHREBIN`/`RVREBIN` declared their spectra `(NCMAX)` and received the
    8192-point FFT arrays of midez and nortwi (compilation error; also
    non-conforming with allocatable images): now `(NCHAN)`.
  - FFT work arrays of 8192 points (local `NMAX` in `FFT2POWER`,
    `FFTCORREL`, `FFTPREP`, corrfft, miderot, midez, nortwi, tfourier,
    wcnoarc), smaller than the image limit (midez, miderot and nortwi
    copy whole spectra into them; miderot and nortwi stopped at start-up
    with `NCMAX>8192`): now `NMAXFFT`.
  - Curves drawn with `NCMAX` points (fftcorrzoom, findmax_wl, sdistor):
    now 1000; fftcorrzoom limits the points around the maximum by the
    size of the correlation function instead of `(NCMAX-1)/2`.
  - leeposs (own `NSMAX=NCMAX=2400`) and fitshead (`NCMAX_LOCAL=2048`) use
    the global limits.
  - Static 2-D arrays with an `NSMAX`/`NCMAX` dimension and a large fixed
    one become allocatable with the actual size: skysubm `SKYRAN`,
    midelines `POL_SIMUL`, fitcdis and fitcdis_new `XPOS`/`XXPOS`, testwc
    `DWV0_MEAS` (to `TESTWC_DATA_M`) and the response curves `FLUX` of
    index, snratio and `MIDEIND` (COMMON `/BLKINDEX1/`, linked into every
    program, to `INDEX_DATA_M`). Static data with 65536 (zero-filled):
    at most 52 MB (movel), median 10 MB; without these changes, up to
    290 MB.
  - Uninitialised variables found by the different memory layout (index,
    snratio) and a short path buffer (corrext); see the bugs fixed.
  - The test runner replaces the `READILIM` upper limits of 1000 or more
    (`(1,...,2400)` in the references, `(1,...,65536)` now, or derived
    values such as `(NSMAX-2)/2` in midez) by a placeholder.

  Verification: `make check` (74 cases) and `make check-std` pass; the
  suite with `-fcheck=bounds` (plus the 4 local patches for the known
  bugs `PGIDEN_RED`, `FFTPREP`, `SPLFIT` and snratio `SNSIM(0)`) differs
  only in fit2dpol, as in 4.3. Nine programs (irevx, istat, imath,
  rotate, binning, isubset, adnsc, resample, writefits_log) run on
  4096x4096 images give the same stdout and output files as 5.1 built
  with `NSMAX=NCMAX=4096`, also with the `-fcheck=bounds` build; irevx
  reverses a 10x20000 image correctly, and an image with 70000 channels
  is rejected by `INFILEX`. Pending: visual checks of the changed plots
  and interactive paths not covered by the tests (sdistor and findmax_wl
  fitted curves, fftcorrzoom in corrfft, midelines errors, testwc).
- rotate checks `NSCAN.GT.NSMAX` and `NCHAN.GT.NCMAX` before allowing
  rotations, which can never be true (INFILEX rejects those images); the
  intended check was probably with the dimensions swapped. Left as it
  is; it becomes irrelevant with allocatable arrays.
- `MIDEIND` (fspec) shares the image with index and snratio through
  COMMON /BLKINDEX0/ (and other blocks): convert it together with those
  two programs, replacing the COMMON blocks by a data module. Done in 4.4
  (`/BLKINDEX1/` in 4.5).
- Before 4.5, review the uses of `NCMAX`/`NSMAX` as constants rather than
  limits, whose meaning would change with the new value: e.g.
  `fftcorrzoom` draws a curve with `NCMAX` points (`XPOL(NCMAX)`,
  `PGLINE(NCMAX,...)`), loops `DO I=1,NCMAX`, `READILIM(...,NCMAX)`.
  Done in 4.5.
- leeposs defines its own `NSMAX`/`NCMAX` (2400) instead of using
  `REDLIB_M`; its 1-D arrays and size checks depend on them. Done in 4.5.
- The 32 PGPLOT calls that pass `NCMAX,NSMAX` as declared dimensions
  (`PGGRAY`, `PGIMAG`...) must pass the actual dimensions when the array
  becomes allocatable; they are not covered by the regression tests
  (sdistor's `PGGRAY`, converted in 4.3, is the only one done). Done in
  4.4.

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

## Bugs of REDUCEME 5.1 fixed during the modernisation

- skysubm: with the errors computed from the approximated formulae, the
  sky errors were stored in `ESKY(NI,I)` with `I` left over from a
  previous loop (`NSCAN+1`) instead of `ESKY(NI,J)`, so the output sky
  error frame was all zeros (and, with allocatable arrays, it would have
  written outside the array). Fixed in phase 4.4; reference `skye.u`
  regenerated.
- skysubm: option `-3` (change the limits of the previous plot) replotted
  the channel `NI`, which holds the last option entered; after `-1`, `-2`
  or `-4` it read the image at a negative channel (garbage plot in 5.1,
  out of the allocatable array now). It now replots the last fitted
  channel. Found in the visual checks of phase 4.4.
- index: `RVELERR0` was undefined with a single radial velocity and no
  errors, or with a table file of velocities; it is printed in the output
  when it is not zero. Now set to zero. Found in phase 4.5.
- snratio: `LINDOK` (indices that can be measured) was undefined for the
  index types without a branch in the selection loop (slopes, `ITI<0`),
  so those indices could appear in the list. Now `.FALSE.`. Found in
  phase 4.5.
- corrext: `REDUCEMEDIR` was `CHARACTER*75`, so the default extinction
  file was not found when `reduceme_dir` was longer. Now `CHARACTER*255`.
  Found in phase 4.5.
- sdistor: the Cauchy fit was computed with `NCMAX` points spread over
  the plot range, but only the first 1000 were drawn (part of the
  curve). Now 1000 points over the whole range. Fixed in phase 4.5.

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
- macOS with a shared PGPLOT library (`libpgplot.dylib`): `red/pgplot.inc`
  reproduces the PGPLOT COMMON blocks, but the executable gets its own
  copy of `/PGPLT1/` and `/PGPLT2/` (ld64 default `-commons ignore_dylibs`)
  instead of sharing the one in the library, so `PGIDEN_RED` always sees
  `PGID=0` and reads `PGNY(0)`, `PGXSZ(0)`... out of bounds (wrong size or
  position of the identification label). Found with `-fcheck=bounds` in
  phase 4.3. Possible fixes: link with `-Wl,-commons,use_dylibs` on macOS,
  or obtain the values through the PGPLOT query routines (`PGQVSZ`...)
  instead of the COMMON blocks.
- `FFTPREP` (red) declares `SP(N)` with `N` the FFT size, but accesses
  `SP(NC1:NC2)` (before moving the region to the start), beyond `N`
  whenever `NC2 > N`.
  The actual arguments (`NMAX=NMAXFFT` arrays in miderot, nortwi...) are
  large enough, so no memory is overwritten, but the declaration is not
  conforming (found with `-fcheck=bounds`).
- autocos2: when the files are opened again to remove the cosmic rays,
  it prints `FILEN(I)` instead of `FILEN(II)` (`I` is left over from a
  previous loop and points to an unused element), so garbage is written
  after `* Opening:`. The test case `autocos2` does not compare stdout
  because of this.
- snratio: in the table of the symmetric binning the scan number is
  written with `(I4,3X,$)`; with `$` gfortran drops the trailing `3X`, so
  the columns are not aligned with the header. Phase 1 changed it to
  `ADVANCE='NO'`, which keeps the blanks (output different from 5.1); the
  5.1 output was restored in phase 4.4 by removing the `3X`.
- snratio: `DO WHILE(((IL1-1).GE.1).AND.(SNSIM(IL1-1).GT.SNSIM(IL1)))`
  evaluates `SNSIM(0)` when `IL1=1` (Fortran does not guarantee
  short-circuit evaluation); harmless in practice, found with
  `-fcheck=bounds`.
- adnhand: `GIVESN` writes REAL numbers with list-directed output
  (`WRITE(CMEAN,*) 0.0`) into `CHARACTER*15` variables; gfortran needs 17
  characters, so the program always stops with "End of record" after
  asking the channel region. The test case `adnhand_known_error`
  documents the current behaviour.
- `FFTCORRZOOM` (red): the polynomial fit around the maximum of the
  correlation uses `X(IMAX-NSIDE:IMAX+NSIDE)` without checking that it
  stays inside `1..N` (only if the maximum is near the edges).
- `SPLFIT` (fspec): `DO WHILE((X(J).LE.XD(I)).AND.(J.LE.N))` evaluates
  `X(N+1)` (no short-circuit evaluation in Fortran); the actual arguments
  are larger static arrays, so it is harmless in practice. Found with
  `-fcheck=bounds` (fit2dspl).
