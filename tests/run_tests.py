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

Usage
-----
    python3 run_tests.py [--builddir DIR] [--update] [case ...]

``--update`` (re)creates the reference results instead of comparing; it must
only be used with binaries built from the reference version.

Only the Python standard library is used. Exit status: 0 if all the cases
pass, 1 otherwise.
"""

import argparse
import difflib
import os
import re
import shutil
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
            'atol': None, 'args': [], 'compare_stdout': True}
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
        elif key == 'compare_stdout':
            case['compare_stdout'] = value.lower() not in ('no', 'false', '0')
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


def run_case(name, case, builddir, srcdir, keep):
    """Run a case in a temporary directory.

    Returns
    -------
    tuple
        (workdir, normalised stdout, return code). The caller must remove
        workdir.
    """
    workdir = tempfile.mkdtemp(prefix='reduceme_%s_' % name)
    mkdata.main(workdir)
    env = dict(os.environ)
    env['PGPLOT_DEV'] = '/NULL'
    env.setdefault('PGPLOT_DIR', '/nonexistent')  # WELCOME only checks it is set
    env['reduceme_dir'] = srcdir
    outputs, codes = [], []
    for k, (program, args, stdin) in enumerate(case['steps'], 1):
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
    if keep:
        print('  workdir kept: %s' % workdir)
    replacements = ((os.path.realpath(workdir) + os.sep, '$workdir/'),
                    (workdir + os.sep, '$workdir/'),
                    (srcdir, '$reduceme_dir'))
    return workdir, normalise_stdout(out, replacements), code


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
    args = parser.parse_args(argv)
    builddir = os.path.abspath(args.builddir)
    srcdir = os.path.abspath(args.srcdir)
    casedir = os.path.join(HERE, 'cases')
    refroot = os.path.join(HERE, 'references')
    names = args.cases or sorted(f[:-5] for f in os.listdir(casedir)
                                 if f.endswith('.case'))
    failures = []
    for name in names:
        case = parse_case(os.path.join(casedir, name + '.case'))
        workdir, out, code = run_case(name, case, builddir, srcdir, args.keep)
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
                msg = compare_file(os.path.join(refdir, fname),
                                   os.path.join(workdir, fname), rtol, atol)
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
        print('\n%d cases, %d failed%s' % (len(names), len(failures),
              (': ' + ' '.join(failures)) if failures else ''))
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
