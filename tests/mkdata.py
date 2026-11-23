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
C_KMS = 299792.46          # speed of light (km/s), as in midez and movel


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


def write_fits16(path, data, bscale, bzero, cards=()):
    """Write an image as a FITS file with BITPIX=16 (as read by fitshead).

    Parameters
    ----------
    path : str
        Output file name.
    data : list of list of float
        Image values, ``data[i][j]`` for row ``i`` (NAXIS2) and column ``j``
        (NAXIS1). They are stored as ``round((value - bzero) / bscale)``.
    bscale, bzero : float
        Scaling written in the header (physical = stored * bscale + bzero).
    cards : sequence of (str, str)
        Additional keywords and their values already formatted as FITS
        fixed-format fields (strings with their quotes, starting in column
        11; numbers right-justified to column 30).
    """
    naxis2 = len(data)
    naxis1 = len(data[0])

    def card(key, value):
        return '%-8s= %s' % (key, value)

    def num(value):
        return '%20s' % value

    header = [card('SIMPLE', num('T')), card('BITPIX', num(16)),
              card('NAXIS', num(2)), card('NAXIS1', num(naxis1)),
              card('NAXIS2', num(naxis2)), card('BSCALE', num(repr(bscale))),
              card('BZERO', num(repr(bzero)))]
    header += [card(k, v) for k, v in cards]
    header.append('END')
    text = ''.join('%-80s' % c for c in header)
    text += ' ' * (-len(text) % 2880)                # header blocks of 2880
    values = [int(round((v - bzero) / bscale)) for row in data for v in row]
    body = struct.pack('>%dh' % len(values), *values)  # big-endian int16
    body += b'\0' * (-len(body) % 2880)
    with open(path, 'wb') as fh:
        fh.write(text.encode('ascii') + body)


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
    # four flat-field like images with different levels and a list file
    # (extension .img, so that they do not appear in the listings of *.u)
    names = []
    for k in range(1, 5):
        flat = [[(1000.0 * k) * (1.0 + 0.01 * math.sin(j / 3.0 + i))
                 for j in range(1, nchan + 1)] for i in range(1, nscan + 1)]
        name = 'flat%d.img' % k
        write_reduceme(os.path.join(outdir, name), flat, obj='Flat %d' % k,
                       **common)
        names.append(name)
    with open(os.path.join(outdir, 'flats.lis'), 'w') as fh:
        fh.write('\n'.join(names) + '\n')
    # five exposures of the same object with small differences between
    # them and some cosmic rays (one of them with a tail), and a list file
    names = []
    cosmics = {2: [(3, 12, 3000.0), (3, 13, 400.0), (4, 12, 300.0)],
               4: [(5, 30, 2000.0)]}             # frame: (scan, channel, counts)
    for k in range(1, 6):
        img = [[b[i - 1][j - 1] *
                (1.0 + 0.03 * math.sin(1.7 * k + 0.9 * j + 2.3 * i))
                for j in range(1, nchan + 1)] for i in range(1, nscan + 1)]
        for i, j, counts in cosmics.get(k, []):
            img[i - 1][j - 1] += counts
        name = 'cr%d.img' % k
        write_reduceme(os.path.join(outdir, name), img, obj='Exposure %d' % k,
                       **common)
        names.append(name)
    with open(os.path.join(outdir, 'crs.lis'), 'w') as fh:
        fh.write('\n'.join(names) + '\n')
    # arc with three emission lines curved along the slit (C-distortion)
    # and a file with their approximate positions (channels, 1-based)
    lines = (8.0, 20.0, 32.0)
    arc = [[10.0 + sum(gaussian(j, x0 + 0.1 * (i - 3.5) ** 2, 1.2, 500.0)
                       for x0 in lines)
            for j in range(1, nchan + 1)] for i in range(1, nscan + 1)]
    write_reduceme(os.path.join(outdir, 'arc.img'), arc, obj='Arc', **common)
    with open(os.path.join(outdir, 'lincdis.dat'), 'w') as fh:
        fh.write(''.join('%d\n' % x0 for x0 in lines))
    # wavelength calibration polynomial of the arc, as written by fitlin
    # (term number and coefficient; wavelength as a function of channel)
    with open(os.path.join(outdir, 'arc.pol'), 'w') as fh:
        fh.write('1 4998.0\n2 1.5\n3 0.002\n')
    # spectra for radial velocity and velocity dispersion: two template
    # stars (separately and in a single frame) and a galaxy frame made from
    # the first star, shifted and broadened (velocities in km/s)
    nch, wv0, dwv = 300, 4800.0, 1.0
    wlines = (4825.0, 4861.0, 4890.0, 4920.0, 4957.0, 4990.0, 5015.0,
              5040.0, 5070.0)
    depths = {1: (0.30, 0.50, 0.20, 0.25, 0.40, 0.15, 0.35, 0.20, 0.30),
              2: (0.20, 0.60, 0.15, 0.30, 0.30, 0.25, 0.25, 0.30, 0.20)}

    def star(k, vel=0.0, sigma=0.0):
        """Template k redshifted by vel and broadened by sigma (km/s)."""
        spec = []
        for j in range(1, nch + 1):
            wv = wv0 + dwv * (j - 1)
            flux = 1000.0 * (1.0 + 0.0005 * (j - nch / 2))
            for wl, dep in zip(wlines, depths[k]):
                wlz = wl * (1.0 + vel / C_KMS)
                sig = math.hypot(1.5, wl * sigma / C_KMS)
                flux -= 1000.0 * dep * 1.5 / sig * \
                    math.exp(-0.5 * ((wv - wlz) / sig) ** 2)
            spec.append(flux)
        return spec

    sp = dict(stwv=wv0, disp=dwv, airmass=1.0, timexpos=600.0)
    for k in (1, 2):
        write_reduceme(os.path.join(outdir, 'templ%d.img' % k), [star(k)],
                       obj='Star %d' % k, **sp)
    write_reduceme(os.path.join(outdir, 'templs.img'), [star(1), star(2)],
                   obj='Stars', **sp)
    gal = [[0.5 * f for f in star(1, 300.0 + 30.0 * (i - 3), 150.0)]
           for i in range(1, 6)]
    write_reduceme(os.path.join(outdir, 'gal.img'), gal, obj='Galaxy', **sp)
    write_reduceme(os.path.join(outdir, 'gale.img'), errors_of(gal),
                   obj='Galaxy', error_frame=True, **sp)
    write_reduceme(os.path.join(outdir, 'gal1.img'), gal[2:3],
                   obj='Galaxy centre', **sp)
    # point source whose spatial profile (along the slit) is displaced
    # with wavelength (S-distortion), and its error frame
    star = [[20.0 + gaussian(i, 12.0 + 1.5 * math.sin((j - 20) / 8.0),
                             1.8, 1000.0)
             for j in range(1, nchan + 1)] for i in range(1, 26)]
    write_reduceme(os.path.join(outdir, 'star.img'), star, obj='Star',
                   **common)
    write_reduceme(os.path.join(outdir, 'stare.img'), errors_of(star),
                   obj='Star', error_frame=True, **common)
    # scan regions of each galaxy spectrum, as written by adnsc (read by movel)
    with open(os.path.join(outdir, 'gal.img.log'), 'w') as fh:
        fh.write(''.join('%d,%d\n' % (i, i) for i in range(1, 6)))
    # 16-bit integer FITS image (fitshead only reads BITPIX=16)
    write_fits16(os.path.join(outdir, 'img16.fits'), a, 0.01, 100.0,
                 cards=(('BUNIT', "'ADU     '"), ('BLANK', '%20d' % -32768),
                        ('OBJECT', "'Synthetic A'"),
                        ('EXPTIME', '%20s' % '900.0'),
                        ('AIRMASS', '%20s' % '1.25'),
                        ('RA', "'12:34:56.78'"), ('DEC', "'+12:34:56.7'"),
                        ('GAIN', '%20s' % '1.5'),
                        ('READNOIS', '%20s' % '5.0')))


if __name__ == '__main__':
    if len(sys.argv) != 2:
        sys.exit('usage: mkdata.py OUTDIR')
    main(sys.argv[1])
