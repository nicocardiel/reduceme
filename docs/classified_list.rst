Classified list of programs
===========================

Input/Output
------------

Since version 6.1 the programs read and write FITS files directly (see
:ref:`fits_files`), so the conversions between FITS and REDUCEME format
(``R6-leefits``, ``R6-writefits``...) are no longer necessary.

``R6-interchange_32_64``: Transforms REDUCEME images from 32 bits to 64 bits
architecture (and viceversa).

``R6-leefits``: Reads a FITS file and creates a new file with REDUCEME format.

``R6-leefits_new``: Reads a FITS file and creates a new file with REDUCEME
format.

``R6-portable``: Transforms REDUCEME images into a more portable ASCII format.

``R6-readfits``: Reads a FITS file and creates a new file with REDUCEME format
(this program is obsolete; use leefits).

``R6-replacefits``: Put data from a REDUCEME file into a previously existing
FITS file. Note: this program requires the FITSIO subroutine package

``R6-sunalpha``: Transforms REDUCEME images from sun architecture to alpha
architecture (and viceversa).

``R6-writeascii``: Reads a file with REDUCEME format and creates a new file
with ASCII format.

``R6-writefits``: Reads a file with REDUCEME format and creates a new file with
FITS format. Note: this program requires the FITSIO subroutine package

``R6-writefits_log``: Reads a file with REDUCEME format and creates a new file
with FITS format and CTYPE1="WAVE-LOG". Note: this program requires the FITSIO
subroutine package

Examination and statistical analysis
------------------------------------

``R6-exhead``: Shows and changes image header information. Image data are not
modified.

``R6-exheads``: Shows image header information without prompting.

``R6-histogram``: Performs a histogram of data tabulated in a single column
ASCII file.

``R6-istat``: Provides some statistics of an image.

Arithmetic and manipulations
----------------------------

``R6-addnf``: Add several images, taking into account offsets in the spatial
direction.

``R6-adnch``: Add channels of an image, generating a new image with the same
NSCAN and variable NCHAN.

``R6-adnhand``: Create a new spectrum by adding spectra from different images
(with interactive monitoring of the S/N ratio).

``R6-adnsc``: Add scans of an image, generating a new image with the same NCHAN
and variable NSCAN.

``R6-adnscgrad``: Add scans of an image, generating a new image with the same
NCHAN and variable NSCAN.

``R6-basicred``: Determine the BIAS level and output the useful region of a
frame after dividing it by the required flatfields. This program also generates
the error frames (gain and readout noise of the employed detector must be
known).

``R6-binning``: Perform a constant binning in the spatial and wavelength
direction.

``R6-broadima``: Broaden selected spectra of an image.

``R6-broadres``: Broadens selected spectra of an image to produce spectra at a
fixed spectral resolution. It needs a table with velocity dispersions for each
spectrum (scan).

``R6-broadsp``: Broaden a single spectrum by convolving with a gaussian. If the
input file is an image, the program extracts a single spectrum prior
broadening.

``R6-corrfft``: Cross-correlates spectra using FFT.

``R6-fillimage``: Fill an image region by extrapolating the signal in a
contigous region. The program detects the transition between the extrapolated
region and the fitted region as a sudden signal change.

``R6-fit1dpol``: Fits 1-D polynomials.

``R6-fit2dpol``: Fits 2-D polynomials.

``R6-fit2dspl``: Fits 2-D splines and local polynomials.

``R6-gain``: Measure the gain using several flatfield images.

``R6-genimage``: Creates an artificial image using, if necessary, regions from
other images and/or tabulated data from ASCII files.

``R6-gluesc``: Takes different spectra and creates a new frame (in which each
individual spectrum correspond to the former single spectra).

``R6-growx``: Expands a single spectrum into an image.

``R6-growy``: Expands a single spatial cross section into an image.

``R6-ifilter``: Applies a filter to an image.

``R6-imath``: Performs arithmetic with images, spectra, spatial cross sections
and constants.

``R6-interp``: Interpolation/extrapolation of data in an image by using
polynomials.

``R6-interpmove``: Interpolation/extrapolation of data in an image by using
polynomials around a previously fitted polynomial.

``R6-irevx``: Reverses an image (or spectrum) in the X-direction.

``R6-irevy``: Reverse an image (or spatial direction) in the Y-direction.

``R6-isubset``: Produces a subset of an image.

``R6-movel.f``: Computation of radial velocity and velocity dispersion using
the MOVEL and OPTEMA algorithms

``R6-multfit``: perform simultaneous fits.

``R6-resample``: Transforms an image with an initial STWV and DISP into another
image with a different STWV and DISP.

``R6-resample_flux``: Transforms an image with an initial wavelength
calibration (linear or log-linear) into another image with a new and linear
wavelength calibration.

``R6-rotate``: Rotates an image.

``R6-tfourier``: Filters a spectrum using FFT.

Distortion
----------

``R6-cdisc``: Correct C-distortion using the polynomial fits of fitcdis.

``R6-fitcdis``: Calculates the C-distortion of an image.

``R6-fitcdis``: Calculates the C-distortion of an image.

``R6-nortwi``: Wavelength normalization of twilight flatfields calculating
C-distortion with cross correlation.

``R6-sdistor``: Corrects S-distortion.

Graphic display
---------------

``R6-plots``: General program to produce line plots and display images.

``R6-plotsp3d``: Plots successive spectra with a 3-D perspective.

``R6-stplot``: Plots the spectrophotometric spectra.

Error handling
--------------

``R6-esterror``: Creates an error file from an initial image, taking into
account the r.m.s. measured in each spectrum in a given wavelength range.

``R6-generror``: Creates an error file from an initial image, taking into
account the readout noise and gain of the detector.

``R6-randsc``: Creates a fake image through bootstraping from an original image
and its associated error frame. In each pixel, the original signal is randomly
modified using the correspoding error.

``R6-snratio``: Determines the S/N ratio as a function of binning in the
spatial direction.

Cosmic Rays
-----------

``R6-autocos2``: Automatic removal of cosmic rays in many similar images
simultaneously (maximum number of images=50, maximum number of cosmic
rays=1000).

``R6-autocos``: Automatic removal of cosmic rays in many similar images
simultaneously (maximum number of images=50, maximum number of cosmic
rays=1000). The program detects cosmic rays by comparing the number of counts
in the same pixel in all the available frames.

``R6-cleanest``: Automatic and interactive removal of cosmic rays.

Wavelengths
-----------

``R6-air2vacuum``: Transforms wavelengths from air to vacuum using the
equations from Greisen et al. 2002 (Representation of spectral coordinates in
FITS).

``R6-calambda``: determine the wavelength as a function of the channel number,
using the wavelength calibration polynomial.

``R6-findarc``: Interactive arc line identification.

``R6-findmax``: Automatic detection of line peaks in a spectrum.

``R6-findmax``: Automatic detection of line peaks in a spectrum and estimation
of the peaks through the fit to a gaussian

``R6-fitlin``: Calculates the wavelength calibration polynomials.

``R6-predchan``: determine the final position (channel and wavelength) of a
pixel in the wavelength axis corresponding to a given channel position (known
before the C-distortion correction and the wavelength calibration).

``R6-rebincw``: Performs the wavelength calibration of an image by using the
polynomial fits obtained with fitlin.

``R6-rvshift``: Applies a radial velocity shift to an image.

``R6-rvshiftrot``: Applies radial velocity shifts to the different scans of an
image. This program can be used to correct 2D spectroscopic images from
rotation curves. It is similar to rvshift, but using different radial
velocities for each scan (these velocities are read from an ASCII file).

``R6-rvz``: Computes the corresponding radial velocity for a given redshift.

``R6-shiftch``: Applies a constant shift (channels) to an image.

``R6-shpol``: Determines the new coefficients of a polynomial after changing
the x-axis origin.

``R6-testwc``: Test the wavelength calibration by measuring the positions of
sky lines.

``R6-wcnoarc``: Performs the wavelength calibration of an image by using the
polynomial fits obtained with fitlin for other frames.

Flux calibration
----------------

``R6-fcalspl``: Generates a flux calibration spectrum using splines.

``R6-fcalsplnew``: Generates a flux calibration spectrum using splines,
allowing an absolute flux calibration.

``R6-prfcal``: Creates a mean flux calibration curve from individual curves,
and generates an image containing the mean curve (as the first spectrum) and
the individual ones (in successive spectra).

Extinction correction
---------------------

``R6-airmass``: Calculate the airmass for fixed observing conditions.

``R6-corrext``: Corrects spectra from atmospheric and interstellar extinction.

Sky subtraction
---------------

``R6-skysubm``: Calculates and subtracts a sky image by fitting polynomials to
each channel.

Measurement
-----------

``R6-index``: Measures indices in a fully calibrated image. The program
requires an external file containing the index definitions. For this purpose,
the program looks first for a file called 'myindex.dat' in the current
directory. If this file does not exist, the program then looks for a file
called 'index.dat' (located in the subdirectory 'files' of the distribution
package). If this last file is also missing, the program stops.

``R6-midelines``: Measures lines (EW) in a spectrum selecting regions manually.

``R6-miderot``: Determines the rotation curve of a galaxy by using cross
correlation.

``R6-midez``: Determines redshifts by using cross correlation.

Miscellany
----------

``R6-chancoor``: Transform r.a. and dec. coordinates from an initial equinox to
another equinox.

``R6-creahelp``: Creates the help file 'helpred.txt', and the auxiliary files
'allfiles.tex', and libraries.tex.

``R6-ecugal``: Transforms r.a. and dec. to galactic coordinates.

``R6-fitdata``: Fits polynomials using tabulated data.

``R6-fitshead``: Reads a FITS file and creates an output file with keyword
information.

``R6-fitstex``: Reads the output from fitshead and generates an output LaTeX
file with a table.

``R6-giveerrf.f``: Given a character string corresponding to a file name, this
program returns the expected error frame associated to that file.

``R6-interlines``: Interpolates lines in a spectrum by using data from a
similar spectrum.

``R6-leeposs``: Reads FITS images created with getimage from the Digital Sky
Survey, allowing the determination of accurate coordinates and plotting of
marks.

``R6-loglin``: Rebins data from a linear wavelength coverage to logarithmic
scale

``R6-mcolumns``: Manipulate columns of data in ASCII files.

``R6-offindex``: Computes index offsets due to differences in the
spectrophotometric system.

``R6-outfitshead``: Outputs the whole header of a FITS file into a text file.

``R6-rvel``: Determines the Earth`s velocity at a fixed time to correct radial
velocities.

``R6-rvrelat``: Transforms radial velocities computed from classical formulae
into radial velocities computed with relativistic expressions.

``R6-testrandom``: Program to test the random number generator.

``R6-vaucoul``: Fits a de Vaucouleurs profile to a galaxy frame.

