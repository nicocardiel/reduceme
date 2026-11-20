#!/usr/bin/env python3
"""Generate the synthetic input data used by the REDUCEME regression tests.

The REDUCEME files are written directly (without using REDUCEME itself), so
that the test inputs do not depend on the code being tested. The layout
reproduces what OUTFILEX (red/outfilex.f) writes with gfortran: Fortran
sequential unformatted records, each one surrounded by a 4-byte record
length marker (native byte order).

Usage
-----
    python3 mkdata.py OUTDIR

Only the Python standard library is used.
"""

import math
import os
import struct
import sys

CLAVE_RED = b'abcdefghijkl'
ERROR_TAG = ' @ERROR@'


def _record(payload):
    """Return a Fortran sequential unformatted record (gfortran layout).

    Parameters
    ----------
    payload : bytes
        Record contents.

    Returns
    -------
    bytes
        Payload surrounded by its length as 4-byte integers.
    """
    marker = struct.pack('=i', len(payload))
    return marker + payload + marker


def write_reduceme(path, data, stwv, disp, airmass=1.0, timexpos=0.0,
                   obj='', fitsfile='', comment='', error_frame=False):
    """Write an image in REDUCEME format.

    Parameters
    ----------
    path : str
        Output file name.
    data : list of list of float
        Image values, ``data[i][j]`` for scan ``i`` (spatial direction,
        0-based) and channel ``j`` (wavelength direction, 0-based). This
        matches the Fortran array ``A(NCHAN,NSCAN)`` written one scan per
        record.
    stwv : float
        Central wavelength of the first pixel (Angstrom).
    disp : float
        Dispersion (Angstrom/pixel).
    airmass, timexpos : float
        Airmass and exposure time (s).
    obj, fitsfile, comment : str
        Header strings.
    error_frame : bool
        If True, append ' @ERROR@' to the object name, as OUTFILEX does for
        error frames.
    """
    nscan = len(data)
    nchan = len(data[0])
    obj_out = obj + ERROR_TAG if error_frame else obj
    out = bytearray()
    out += _record(CLAVE_RED)
    out += _record(struct.pack('=ii', nscan, nchan))
    out += _record(struct.pack('=ff', stwv, disp))
    out += _record(struct.pack('=f', airmass))
    out += _record(struct.pack('=f', timexpos))
    for text in (obj_out, fitsfile, comment):
        encoded = text.encode('ascii')
        out += _record(struct.pack('=i', len(encoded)))
        if encoded:
            out += _record(encoded)
    for row in data:
        if len(row) != nchan:
            raise ValueError('all scans must have the same number of channels')
        out += _record(struct.pack('=%df' % nchan, *row))
    with open(path, 'wb') as fh:
        fh.write(out)


def gaussian(x, x0, sigma, amp):
    """Gaussian profile evaluated at x."""
    return amp * math.exp(-0.5 * ((x - x0) / sigma) ** 2)


def image_a(nscan, nchan):
    """Continuum with a slope, an absorption and an emission line.

    Channels and scans are 1-based in the formulas (Fortran convention).
    """
    return [[100.0 + 0.5 * j + 5.0 * i
             - gaussian(j, 15.0, 2.0, 30.0)
             + gaussian(j, 30.0, 1.5, 40.0)
             for j in range(1, nchan + 1)]
            for i in range(1, nscan + 1)]


def image_b(nscan, nchan):
    """Smooth image with a sinusoidal modulation."""
    return [[50.0 + 0.2 * j - 2.0 * i + 3.0 * math.sin(j / 4.0)
             for j in range(1, nchan + 1)]
            for i in range(1, nscan + 1)]


def errors_of(data):
    """Poisson-like error frame: sqrt(|value|)."""
    return [[math.sqrt(abs(v)) for v in row] for row in data]


def main(outdir):
    """Create all the test input files in outdir."""
    os.makedirs(outdir, exist_ok=True)
    nscan, nchan = 6, 40
    stwv, disp = 5000.0, 1.5
    a = image_a(nscan, nchan)
    b = image_b(nscan, nchan)
    common = dict(stwv=stwv, disp=disp, airmass=1.25, timexpos=900.0)
    write_reduceme(os.path.join(outdir, 'img_a.u'), a, obj='Synthetic A',
                   **common)
    write_reduceme(os.path.join(outdir, 'img_ae.u'), errors_of(a),
                   obj='Synthetic A', error_frame=True, **common)
    write_reduceme(os.path.join(outdir, 'img_b.u'), b, obj='Synthetic B',
                   **common)
    write_reduceme(os.path.join(outdir, 'img_be.u'), errors_of(b),
                   obj='Synthetic B', error_frame=True, **common)
    # single spectrum (NSCAN=1) and its error frame
    write_reduceme(os.path.join(outdir, 'spec.u'), a[:1], obj='Spectrum',
                   **common)
    write_reduceme(os.path.join(outdir, 'spece.u'), errors_of(a[:1]),
                   obj='Spectrum', error_frame=True, **common)
    # single spatial cross section (NCHAN=1)
    write_reduceme(os.path.join(outdir, 'xsec.u'),
                   [[row[19]] for row in a], obj='Cross section', **common)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('usage: mkdata.py OUTDIR')
    main(sys.argv[1])
