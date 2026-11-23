#!/usr/bin/env python3
"""Regression tests for the REDUCEME programs.

Each test case (``cases/<name>.case``) runs one program in a temporary
directory containing the synthetic data created by ``mkdata.py``, feeding
the program's prompts from standard input. The standard output (without the
version banner) and the requested output files are compared with the
reference results stored in ``references/<name>/``.

Case file format::

    # comments and blank lines are ignored before the stdin section
    program = irevx
    outputs = out.u            # optional, space-separated
    timeout = 30               # optional, seconds (per program run)
    rtol = 1e-6                # optional, relative tolerance for numbers
    atol = 0                   # optional, absolute tolerance for numbers
    args = file.u              # optional, command-line arguments
    compare_stdout = no        # optional (default yes); "no" for programs
                               # whose output is random (only the return
                               # code and the output files are compared)
    fits = no                  # optional (default yes); "no" excludes the
                               # case from the FITS mode
    --- stdin
    img_a.u
    out.u

A case may also run several programs in sequence in the same directory
(e.g. a format conversion and its inverse): instead of ``program`` and
``--- stdin``, use one ``--- run PROGRAM [ARGS]`` section per run, each
followed by the lines fed to that program.

Comparisons:

- stdout and ASCII files: line by line; numbers are compared with a relative
  tolerance (``--rtol``), the remaining text must be identical. FITS header
  cards with the current date copied into text files (``DATE...= ``,
  ``HISTORY Date:``) are not compared.
- FITS files (``.fits``, ``.fit``): header cards compared ignoring the
  cards that contain the current date and the ``HISTORY REDUCEMEv...``
  cards (REDUCEME version); data compared byte by byte.
- other binary files (REDUCEME format): byte by byte.

In the standard output, the source directory (``reduceme_dir``) and the
temporary working directory are replaced by ``$reduceme_dir`` and
``$workdir``, and the version banner is removed. In gfortran runtime error
messages, the source line number is replaced by ``#``, the source file
extension is written as ``.f`` (sources were renamed from ``.f`` to ``.f90``)
and the backtrace lines are removed, so that a case documenting a known
runtime error does not depend on line numbers, file names or the call stack.

The upper limits shown by ``READILIM`` prompts, ``(N1,...,N2)``, that are
1000 or larger are replaced by ``(N1,...,<limit>)`` in the standard output
and in the reference: they derive from the maximum image dimensions
(``NSMAX``, ``NCMAX``), which were 2400 when the references were created
from version 5.1 and are 65536 now (e.g. ``(1,...,2400)``, or
``(1,...,1199)`` for ``(NSMAX-2)/2`` in midez). The synthetic images are
much smaller, so no other limit reaches that value.

FITS mode (``--fits``): the cases whose answers or outputs contain REDUCEME
file names ``NAME.u`` are run again with ``NAME.fits`` instead (the inputs
created by ``mkdata.py`` are converted into FITS files keeping all their
header values), and compared with the same references: the FITS files
written instead of REDUCEME reference files must have the same dimensions,
REDUCEME header values and data (bit by bit); in the standard output and in
the text files the FITS names are written back as the original names.

Usage
-----
    python3 run_tests.py [--builddir DIR] [--update | --fits] [case ...]

``--update`` (re)creates the reference results instead of comparing; it must
only be used with binaries built from the reference version.

Only the Python standard library is used. Exit status: 0 if all the cases
pass, 1 otherwise.
"""

import argparse
import array
import difflib
import os
import re
import shutil
import struct
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mkdata  # noqa: E402

NUMBER = re.compile(r'[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[EeDd][-+]?\d+)?')
# version banner printed by WELCOME (red/welcome.f): version and date change
BANNER = re.compile(r'REDUCEMEv\S+\s+Welcome to .*Version:')
# gfortran runtime error location and backtrace lines
RUNTIME_LINE = re.compile(r'At line \d+ of file (\S+?)\.f(?:90)?\b')
BACKTRACE = re.compile(r'^#\d+\s+0x[0-9a-fA-F]+.*$')
# READILIM upper limits derived from the maximum image dimensions
SIZE_LIMIT = re.compile(r'\((-?\d+),\.\.\.,\d{4,}\)')


def parse_case(path):
    """Read a case file.

    Parameters
    ----------
    path : str
        Case file name.

    Returns
    -------
    dict
        Keys 'steps' (list of (program, args, stdin) tuples), 'outputs' (list),
        'timeout' (float) and 'rtol', 'atol' (float or None).
    """
    case = {'outputs': [], 'timeout': 30.0, 'steps': [], 'rtol': None,
            'atol': None, 'args': [], 'compare_stdout': True, 'fits': True}
    with open(path) as fh:
        lines = fh.read().split('\n')
    program = None
    i = 0
    while i < len(lines):
        line = lines[i]
        if line.startswith('--- '):
            break
        line = line.split('#', 1)[0].strip()
        i += 1
        if not line:
            continue
        key, _, value = line.partition('=')
        key, value = key.strip(), value.strip()
        if key == 'program':
            program = value
        elif key == 'outputs':
            case['outputs'] = value.split()
        elif key in ('timeout', 'rtol', 'atol'):
            case[key] = float(value)
        elif key == 'args':
            case['args'] = value.split()
        elif key in ('compare_stdout', 'fits'):
            case[key] = value.lower() not in ('no', 'false', '0')
        else:
            raise ValueError('%s: unknown key %r' % (path, key))
    # stdin sections
    while i < len(lines):
        head = lines[i].strip()
        body = []
        i += 1
        while i < len(lines) and not lines[i].startswith('--- '):
            body.append(lines[i])
            i += 1
        if head == '--- stdin':
            if program is None:
                raise ValueError('%s: "--- stdin" requires "program"' % path)
            prog, args = program, case['args']
        elif head.startswith('--- run '):
            words = head[len('--- run '):].split()
            prog, args = words[0], words[1:]
        else:
            raise ValueError('%s: unknown section %r' % (path, head))
        case['steps'].append((prog, args, '\n'.join(body)))
    if not case['steps']:
        raise ValueError('%s: no "--- stdin" or "--- run" section' % path)
    return case


def normalise_limits(text):
    """Replace the READILIM upper limits of 1000 or more by a placeholder.

    Parameters
    ----------
    text : str
        Program output or reference output.

    Returns
    -------
    str
        Text with ``(N1,...,N2)``, N2 >= 1000, written as ``(N1,...,<limit>)``.
    """
    return SIZE_LIMIT.sub(r'(\1,...,<limit>)', text)


def normalise_stdout(text, replacements=()):
    """Remove run-dependent content from the program output.

    Parameters
    ----------
    text : str
        Program output.
    replacements : sequence of (str, str)
        Literal substitutions (e.g. absolute paths by placeholders), applied
        in order.

    Returns
    -------
    str
        Output without the version banner, with the substitutions done and
        with the size-dependent READILIM limits normalised.
    """
    for old, new in replacements:
        text = text.replace(old, new)
    text = normalise_limits(text)
    out = []
    for line in text.split('\n'):
        if BANNER.search(line):
            continue
        line = RUNTIME_LINE.sub(r'At line # of file \1.f', line)
        if BACKTRACE.match(line):
            continue
        out.append(line.rstrip())
    return '\n'.join(out).rstrip() + '\n'


def is_binary(data):
    """True if data does not look like ASCII text."""
    return b'\0' in data or any(b > 127 for b in data[:4096])


def compare_text(ref, new, rtol, atol=0.0):
    """Compare two texts line by line with numeric tolerance.

    Returns
    -------
    str or None
        Description of the first difference, or None if equivalent.
    """
    rl, nl = ref.split('\n'), new.split('\n')
    if len(rl) != len(nl):
        diff = '\n'.join(list(difflib.unified_diff(rl, nl, 'reference', 'new',
                                                   lineterm='', n=1))[:40])
        return 'different number of lines (%d vs %d)\n%s' % (len(rl), len(nl),
                                                             diff)
    for k, (a, b) in enumerate(zip(rl, nl), 1):
        if a == b:
            continue
        # FITS header cards with the date copied into text files (outfitshead)
        if FITS_DATE_CARD_TEXT.match(a) and FITS_DATE_CARD_TEXT.match(b):
            continue
        if NUMBER.sub('#', a) != NUMBER.sub('#', b):
            return 'line %d differs:\n  ref: %s\n  new: %s' % (k, a, b)
        for x, y in zip(NUMBER.findall(a), NUMBER.findall(b)):
            fx = float(x.replace('D', 'E').replace('d', 'e'))
            fy = float(y.replace('D', 'E').replace('d', 'e'))
            if abs(fx - fy) > max(atol, rtol * max(abs(fx), abs(fy))):
                return 'line %d: %s != %s\n  ref: %s\n  new: %s' % (k, x, y,
                                                                    a, b)
    return None


# FITS header cards that depend on the date of execution
# (and the HISTORY card with the REDUCEME version written by SCRATCH2FITS)
FITS_DATE_CARD = re.compile(r'^(DATE|HISTORY Date:|HISTORY REDUCEMEv)')
# the same cards in text files: stricter (keyword padded to 8 characters and
# followed by '= '), so that other lines starting with DATE are compared
FITS_DATE_CARD_TEXT = re.compile(r'^(DATE[-\w ]{0,4}= |HISTORY Date:)')


def split_fits(data):
    """Split a FITS file into (list of header cards, rest of the file).

    Only the primary header is parsed; everything after it is returned as
    raw bytes.
    """
    cards, pos = [], 0
    while pos + 80 <= len(data):
        card = data[pos:pos + 80].decode('ascii', 'replace')
        cards.append(card)
        pos += 80
        if card.startswith('END') and card[3:].strip() == '':
            break
    pos = ((pos + 2879) // 2880) * 2880   # header padded to 2880-byte blocks
    return cards, data[pos:]


def compare_fits(ref, new):
    """Compare two FITS files ignoring the date-dependent header cards."""
    rcards, rdata = split_fits(ref)
    ncards, ndata = split_fits(new)
    rcards = [c for c in rcards if not FITS_DATE_CARD.match(c)]
    ncards = [c for c in ncards if not FITS_DATE_CARD.match(c)]
    if rcards != ncards:
        for a, b in zip(rcards, ncards):
            if a != b:
                return 'FITS header differs:\n  ref: %s\n  new: %s' % (
                    a.rstrip(), b.rstrip())
        return 'FITS header differs in number of cards (%d vs %d)' % (
            len(rcards), len(ncards))
    if rdata != ndata:
        return 'FITS data differ'
    return None


def compare_file(ref_path, new_path, rtol, atol):
    """Compare a reference file with a newly generated one."""
    if not os.path.exists(new_path):
        return 'output file not created'
    with open(ref_path, 'rb') as fh:
        ref = fh.read()
    with open(new_path, 'rb') as fh:
        new = fh.read()
    if new_path.lower().endswith(('.fits', '.fit')):
        return compare_fits(ref, new)
    if is_binary(ref) or is_binary(new):
        if ref == new:
            return None
        if len(ref) != len(new):
            return 'binary files differ in size (%d vs %d bytes)' % (len(ref),
                                                                     len(new))
        first = next(i for i in range(len(ref)) if ref[i] != new[i])
        return 'binary files differ (first difference at byte %d)' % first
    return compare_text(ref.decode('latin-1'), new.decode('latin-1'), rtol,
                        atol)


# ---------------------------------------------------------------------------
# FITS mode (--fits): the cases are rerun with FITS files instead of REDUCEME
# files and the results are compared with the same references.

# REDUCEME file names in the answers: NAME.u (not part of a path nor followed
# by another extension). Files with other extensions (e.g. .img) are not
# replaced: some programs derive other file names from them (gal.img.log)
# or print them in fixed-width columns.
REDUCEME_NAME = re.compile(r'(?<![\w.\-/])([\w\-]+)\.u(?![\w.])')


def read_records(data):
    """Split a Fortran sequential unformatted file into its records."""
    records, pos = [], 0
    while pos + 4 <= len(data):
        n = struct.unpack('=i', data[pos:pos + 4])[0]
        records.append(data[pos + 4:pos + 4 + n])
        pos += n + 8
    return records


def read_reduceme(data):
    """Parse a file in REDUCEME format.

    Parameters
    ----------
    data : bytes
        File contents.

    Returns
    -------
    dict or None
        None if the file is not in REDUCEME format. Otherwise: 'nscan',
        'nchan' (int); 'stwv', 'disp', 'airmass', 'timexpos' (the REAL values
        as 4 bytes in native order, for exact comparisons); 'object',
        'fitsfile', 'comment' (str); 'scans' (list of bytes, one record per
        scan).
    """
    if data[4:16] != mkdata.CLAVE_RED:
        return None
    rec = read_records(data)
    nscan, nchan = struct.unpack('=ii', rec[1])
    strings, k = [], 5
    for _ in range(3):
        n = struct.unpack('=i', rec[k])[0]
        k += 1
        if n > 0:
            strings.append(rec[k].decode('latin-1'))
            k += 1
        else:
            strings.append('')
    return {'nscan': nscan, 'nchan': nchan, 'stwv': rec[2][:4],
            'disp': rec[2][4:8], 'airmass': rec[3], 'timexpos': rec[4],
            'object': strings[0], 'fitsfile': strings[1],
            'comment': strings[2], 'scans': rec[k:k + nscan]}


def reduceme_to_fits(src, dst):
    """Convert a REDUCEME file into a FITS file as written by SCRATCH2FITS.

    The header keeps all the REDUCEME header values: CRPIX1=1, CRVAL1=STWV,
    CDELT1=DISP, OBJECT, EXPTIME, AIRMASS, FITSFILE (blank when it is blank
    in the REDUCEME file, so that FITS2SCRATCH does not replace it by the
    file name), REDCOMM and ERRFRAME for error frames.

    Returns
    -------
    bool
        False if src is not a REDUCEME file (nothing is written).
    """
    with open(src, 'rb') as fh:
        head = read_reduceme(fh.read())
    if head is None:
        return False

    def real(raw):
        return struct.unpack('=f', raw)[0]

    def num(value):
        return '%20s' % ('%.9g' % value)   # 9 digits: exact REAL round trip

    def text(value):
        if len(value) > 68:
            raise ValueError('%s: string too long for a FITS card' % src)
        return "'%-8s'" % value.replace("'", "''")

    obj = head['object']
    error_frame = obj.endswith(mkdata.ERROR_TAG)
    if error_frame:
        obj = obj[:-len(mkdata.ERROR_TAG)]
    cards = [('CTYPE1', text('WAVE')), ('CRPIX1', num(1.0)),
             ('CRVAL1', num(real(head['stwv']))),
             ('CDELT1', num(real(head['disp']))), ('OBJECT', text(obj)),
             ('EXPTIME', num(real(head['timexpos']))),
             ('AIRMASS', num(real(head['airmass']))),
             ('FITSFILE', text(head['fitsfile']))]
    if head['comment']:
        cards.append(('REDCOMM', text(head['comment'])))
    if error_frame:
        cards.append(('ERRFRAME', '%20s' % 'T'))
    nchan = head['nchan']
    data = [list(struct.unpack('=%df' % nchan, scan))
            for scan in head['scans']]
    mkdata.write_fits(dst, data, cards=cards,
                      naxis=1 if head['nscan'] == 1 else 2)
    return True


def fits_keywords(cards):
    """Return the values of the FITS header keywords as strings.

    String values are unquoted, without trailing blanks, and long strings
    (CONTINUE convention) are joined.
    """
    values, last = {}, None
    string = re.compile(r"'((?:[^']|'')*)'")
    for card in cards:
        key = card[:8].strip()
        if key == 'CONTINUE' and last is not None and \
                values[last].endswith('&'):
            m = string.match(card[10:].strip())
            values[last] = values[last][:-1] + (
                m.group(1).replace("''", "'").rstrip() if m else '')
            continue
        last = None
        if card[8:10] != '= ':
            continue
        value = card[10:].strip()
        m = string.match(value)
        if m:
            values[key] = m.group(1).replace("''", "'").rstrip()
            last = key
        else:
            values[key] = value.split('/')[0].strip()
    return values


def compare_fits_reduceme(ref, new):
    """Compare a FITS file with a reference file in REDUCEME format.

    The image dimensions, the REDUCEME header values (from CRVAL1, CRPIX1,
    CDELT1, EXPTIME, AIRMASS, OBJECT, ERRFRAME, FITSFILE, REDCOMM) and the
    data (bit by bit) must be identical.
    """
    head = read_reduceme(ref)
    cards, body = split_fits(new)
    kw = fits_keywords(cards)
    nscan, nchan = head['nscan'], head['nchan']
    found = [kw.get('BITPIX'), kw.get('NAXIS'), kw.get('NAXIS1'),
             kw.get('NAXIS2', '1')]
    expected = ['-32', '1' if nscan == 1 else '2', str(nchan), str(nscan)]
    if found != expected:
        return 'FITS BITPIX, NAXIS, NAXIS1, NAXIS2 = %s (expected %s)' % (
            found, expected)

    def real(key, default=None):
        return struct.pack('=f', float(kw.get(key, default)))

    try:
        stwv = struct.pack('=f', float(kw['CRVAL1']) + (
            1.0 - float(kw.get('CRPIX1', '1'))) * float(kw['CDELT1']))
        reals = [('STWV (CRVAL1)', stwv, head['stwv']),
                 ('DISP (CDELT1)', real('CDELT1'), head['disp']),
                 ('EXPTIME', real('EXPTIME'), head['timexpos']),
                 ('AIRMASS', real('AIRMASS'), head['airmass'])]
    except (KeyError, TypeError, ValueError) as exc:
        return 'FITS header: missing or invalid keyword (%s)' % exc
    problems = ['%s differs' % k for k, a, b in reals if a != b]
    obj = head['object']
    error_frame = obj.endswith(mkdata.ERROR_TAG)
    if error_frame:
        obj = obj[:-len(mkdata.ERROR_TAG)]
    for key, ref_value in (('OBJECT', obj), ('FITSFILE', head['fitsfile']),
                           ('REDCOMM', head['comment'])):
        if kw.get(key, '') != ref_value.rstrip():
            problems.append('%s = %r (expected %r)' % (key, kw.get(key, ''),
                                                       ref_value.rstrip()))
    if (kw.get('ERRFRAME', 'F') == 'T') != error_frame:
        problems.append('ERRFRAME = %s' % kw.get('ERRFRAME', 'F'))
    values = array.array('f')
    values.frombytes(body[:4 * nscan * nchan])
    if sys.byteorder == 'little':
        values.byteswap()                     # FITS data are big-endian
    if values.tobytes() != b''.join(head['scans']):
        problems.append('FITS data differ')
    return '; '.join(problems) or None


def fits_case_names(workdir, case):
    """Prepare a case for the FITS mode.

    Converts the REDUCEME files NAME.u created by mkdata.py into NAME.fits
    and maps every name NAME.u found in the answers and in the outputs of
    the case to NAME.fits.

    Returns
    -------
    dict
        Original name -> FITS name.
    """
    names = {}
    for fname in sorted(os.listdir(workdir)):
        m = REDUCEME_NAME.fullmatch(fname)
        if m and reduceme_to_fits(os.path.join(workdir, fname),
                                  os.path.join(workdir,
                                               m.group(1) + '.fits')):
            names[fname] = m.group(1) + '.fits'
    texts = [' '.join(args) + '\n' + stdin for _, args, stdin in case['steps']]
    for m in REDUCEME_NAME.finditer('\n'.join(texts + case['outputs'])):
        names.setdefault(m.group(0), m.group(1) + '.fits')
    fits_names = list(names.values())
    if len(set(fits_names)) != len(fits_names):
        raise RuntimeError('FITS mode: two files map to the same FITS name')
    return names


def fits_eligible(case):
    """True if the case can be run in FITS mode (it uses REDUCEME files)."""
    return case['fits'] and any(
        REDUCEME_NAME.search(' '.join(args) + '\n' + stdin)
        for _, args, stdin in case['steps'])


def rename_files(text, names):
    """Replace the file names in text (keys of names) by their values."""
    if not names:
        return text
    keys = sorted(names, key=len, reverse=True)
    pattern = re.compile(r'(?<![\w\-])(%s)(?![\w])'
                         % '|'.join(re.escape(k) for k in keys))
    return pattern.sub(lambda m: names[m.group(1)], text)


def compare_file_fits(ref_path, new_path, names, rtol, atol):
    """Compare an output of a case run in FITS mode with its reference.

    A REDUCEME reference is compared with the FITS file written instead of
    it; text files are compared after writing the FITS names back as the
    original names; other files as in compare_file.
    """
    if not os.path.exists(new_path):
        return 'output file not created'
    with open(ref_path, 'rb') as fh:
        ref = fh.read()
    with open(new_path, 'rb') as fh:
        new = fh.read()
    if read_reduceme(ref) is not None and new_path.endswith('.fits'):
        return compare_fits_reduceme(ref, new)
    if not (new_path.lower().endswith(('.fits', '.fit')) or
            is_binary(ref) or is_binary(new)):
        back = {v: k for k, v in names.items()}
        return compare_text(ref.decode('latin-1'),
                            rename_files(new.decode('latin-1'), back),
                            rtol, atol)
    return compare_file(ref_path, new_path, rtol, atol)


def run_case(name, case, builddir, srcdir, keep, fits=False):
    """Run a case in a temporary directory.

    Parameters
    ----------
    fits : bool
        FITS mode: the REDUCEME files of the answers are replaced by FITS
        files (see fits_case_names); the FITS names in the standard output
        are written back as the original names.

    Returns
    -------
    tuple
        (workdir, normalised stdout, return code, names). names maps the
        original file names to the FITS names (empty without fits). The
        caller must remove workdir.
    """
    workdir = tempfile.mkdtemp(prefix='reduceme_%s_' % name)
    mkdata.main(workdir)
    steps, names = case['steps'], {}
    if fits:
        names = fits_case_names(workdir, case)
        steps = [(program, [rename_files(a, names) for a in args],
                  rename_files(stdin, names))
                 for program, args, stdin in steps]
    env = dict(os.environ)
    env['PGPLOT_DEV'] = '/NULL'
    env.setdefault('PGPLOT_DIR', '/nonexistent')  # WELCOME only checks it is set
    env['reduceme_dir'] = srcdir
    outputs, codes = [], []
    for k, (program, args, stdin) in enumerate(steps, 1):
        exe = os.path.join(builddir, 'src', program)
        try:
            proc = subprocess.run([exe] + args, input=stdin.encode(),
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.STDOUT, cwd=workdir,
                                  env=env, timeout=case['timeout'])
            out, code = proc.stdout.decode('latin-1'), proc.returncode
        except subprocess.TimeoutExpired as exc:
            out = (exc.stdout or b'').decode('latin-1') + '\n[TIMEOUT]'
            code = 'timeout'
        if len(case['steps']) > 1:
            out = '=== run %d: %s\n%s' % (k, program, out)
        outputs.append(out)
        codes.append(str(code))
    out, code = '\n'.join(outputs), ' '.join(codes)
    if fits:
        out = rename_files(out, {v: k for k, v in names.items()})
    if keep:
        print('  workdir kept: %s' % workdir)
    replacements = ((os.path.realpath(workdir) + os.sep, '$workdir/'),
                    (workdir + os.sep, '$workdir/'),
                    (srcdir, '$reduceme_dir'))
    return workdir, normalise_stdout(out, replacements), code, names


def main(argv=None):
    """Run the regression tests (or update the references)."""
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('cases', nargs='*', help='case names (default: all)')
    parser.add_argument('--builddir', default=os.environ.get(
        'top_builddir', os.path.join(HERE, '..')),
        help='top build directory (contains src/)')
    parser.add_argument('--srcdir', default=os.environ.get(
        'top_srcdir', os.path.join(HERE, '..')),
        help='top source directory (used as reduceme_dir)')
    parser.add_argument('--update', action='store_true',
                        help='create/overwrite the reference results')
    parser.add_argument('--rtol', type=float, default=1e-6,
                        help='relative tolerance for numbers in text output')
    parser.add_argument('--atol', type=float, default=0.0,
                        help='absolute tolerance for numbers in text output')
    parser.add_argument('--keep', action='store_true',
                        help='keep the temporary working directories')
    parser.add_argument('--fits', action='store_true',
                        help='FITS mode: rerun the cases that use REDUCEME '
                        'files with FITS files instead')
    args = parser.parse_args(argv)
    if args.fits and args.update:
        parser.error('--update cannot be used with --fits')
    builddir = os.path.abspath(args.builddir)
    srcdir = os.path.abspath(args.srcdir)
    casedir = os.path.join(HERE, 'cases')
    refroot = os.path.join(HERE, 'references')
    names = args.cases or sorted(f[:-5] for f in os.listdir(casedir)
                                 if f.endswith('.case'))
    failures, nrun = [], 0
    for name in names:
        case = parse_case(os.path.join(casedir, name + '.case'))
        if args.fits and not fits_eligible(case):
            continue
        nrun += 1
        workdir, out, code, fnames = run_case(name, case, builddir, srcdir,
                                              args.keep, args.fits)
        refdir = os.path.join(refroot, name)
        rtol = args.rtol if case['rtol'] is None else case['rtol']
        atol = args.atol if case['atol'] is None else case['atol']
        try:
            if args.update:
                if os.path.isdir(refdir):
                    shutil.rmtree(refdir)
                os.makedirs(refdir)
                with open(os.path.join(refdir, 'stdout.txt'), 'w',
                          encoding='latin-1') as fh:
                    fh.write(out)
                with open(os.path.join(refdir, 'returncode.txt'), 'w') as fh:
                    fh.write('%s\n' % code)
                for fname in case['outputs']:
                    src = os.path.join(workdir, fname)
                    if not os.path.exists(src):
                        raise RuntimeError('%s: output %s not created'
                                           % (name, fname))
                    shutil.copyfile(src, os.path.join(refdir, fname))
                print('UPDATED %s' % name)
                continue
            problems = []
            with open(os.path.join(refdir, 'returncode.txt')) as fh:
                ref_code = fh.read().strip()
            if ref_code != str(code):
                problems.append('return code %s (reference %s)'
                                % (code, ref_code))
            if case['compare_stdout']:
                with open(os.path.join(refdir, 'stdout.txt'),
                          encoding='latin-1') as fh:
                    msg = compare_text(normalise_limits(fh.read()), out,
                                       rtol, atol)
                if msg:
                    problems.append('stdout: ' + msg)
            for fname in case['outputs']:
                if args.fits:
                    msg = compare_file_fits(
                        os.path.join(refdir, fname),
                        os.path.join(workdir, fnames.get(fname, fname)),
                        fnames, rtol, atol)
                else:
                    msg = compare_file(os.path.join(refdir, fname),
                                       os.path.join(workdir, fname), rtol,
                                       atol)
                if msg:
                    problems.append('%s: %s' % (fname, msg))
            if problems:
                failures.append(name)
                print('FAIL %s' % name)
                for p in problems:
                    print('  ' + p.replace('\n', '\n  '))
            else:
                print('PASS %s' % name)
        finally:
            if not args.keep:
                shutil.rmtree(workdir, ignore_errors=True)
    if not args.update:
        print('\n%d cases%s, %d failed%s' % (
            nrun, ' (FITS mode)' if args.fits else '', len(failures),
            (': ' + ' '.join(failures)) if failures else ''))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
