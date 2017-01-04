Auxiliary libraries
===================

Libray ``libred.a``
-------------------

``autoplot``
............

::

   SUBROUTINE AUTOPLOT(N,X,Y,N1,N2,
                       CX,CY,CT,
                       LLIMIT,X1,X2,Y1,Y2,FEXPAND,
                       JUST,LCLEAR,XOPT,YOPT,
                       IDATA,ICOLOR,
                       XV1,XV2,YV1,YV2)

   Input: N,X,Y,N1,N2,
          CX,CY,CT,
          LLIMIT,X1,X2,Y1,Y2,FEXPAND
          JUST,LCLEAR,XOPT,YOPT,
          IDATA,ICOLOR,
          XV1,XV2,YV1,YV2
   Output: X1,X2,Y1,Y2
   
   Plot X(N),Y(N) (with N in the range N1,...,N2) in the viewport region defined
   by XV1,XV2,YV1,YV2, using common routines from PGPLOT.
   
   INTEGER N -> No. of data points to be plotted
   REAL X(N) -> X-coordinates of the data to be plotted
   REAL Y(N) -> Y-coordinates of the data to be plotted
   INTEGER N1,N2 -> subset of X(N),Y(N) to be plotted (N in the range N1,...,N2)
   CHARACTER*(*) CX -> X-axis label \
   CHARACTER*(*) CY -> Y-axis label  |---> PGLABEL(CX,CY,CT)
   CHARACTER*(*) CT ->   plot title /
   LOGICAL LLIMIT -> if .TRUE., compute new window limits X1,X2,X3,X4
                     if .FALSE. employ input values
   REAL X1,X2,Y1,Y2 -> window limits
   REAL FEXPAND -> fraction to be employed to expand limits (only if
                   LLIMIT=.TRUE.)
   INTEGER JUST -> if JUST=1, the scales of the x and y world coordinates
                   will be equal
   LOGICAL LCLEAR -> if .TRUE. the viewport rectangle XV1,XV2,YV1,YV2 is cleared
                    prior plotting
   CHARACTER*(*) XOPT,YOPT -> controls the plotting of axes (see PGBOX):
      A: draw Axis (X axis is line Y=0, Y axis is line X=0)
      B: draw bottom (X) or left (Y) edge of frame
      C: draw top (X) or right (Y) edge of frame
      G: draw Grid of vertical (X) or horizontal (Y) lines
      I: invert the tick marks; i.e. draw them outside viewport
      L: label axis logarithmically
      N: write numeric lables in the conventional location
      P: extend major tick marks outside the box
      M: write numeric lables in the unconventional location
      T: draw major tick marks at the major coordinate interval
      S: draw minor tick marks (subticks)
      V: orient numeric labels vertically (only to Y)
      1: force decimal labelling, instead of automatic choice (PGNUMB)
      2: force exponential labelling, instead of automatic
   INTEGER IDATA -> controls the plotting of data:
      IDATA=-8,...,31: plot symbols with PGPOINT(N,X,Y,IDATA)
      IDATA=100: plot line with PGLINE(N,X,Y)
      IDATA=101: plot line with PGBIN(N,X,Y,.TRUE.)
      IDATA=102: do not plot data (only box, if required)
   INTEGER ICOLOR -> PGPLOT color for data
   REAL XV1,XV2,YV1,YV2 -> viewport limits

``avoid_warnings``
..................

::
   
   SUBROUTINE AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)

   Input: STWV,DISP,NSCAN,NCHAN
   Output: (nothing)
   
   Dummy function to avoid compilation warnings "unused variable".

``guessef``
...........

::

   SUBROUTINE GUESSEF(INFILE,OUTFILE)

   Input: INFILE
   Ouput: OUTFILE
   
   This subroutine determines the expected error file name OUTFILE from INFILE.
   An additional "e" character is located between the portion of the file name
   (INFILE) preceding the last period (if present) and the last period itself.
   
   CHARACTER*(*) INFILE
   CHARACTER*(*) OUTFILE

``infilex``
...........

::
   
   CHARACTER*75 FUNCTION INFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)

   Input: NF,FILENAME,MODE,LERR  (and NSCAN,NCHAN,STWV,DISP if MODE=21)
   Output: NSCAN,NCHAN,STWV,DISP
   Output (COMMON): AIRMASS,TIMEXPOS,OBJECT,FITSFILE,COMMENT
   
   This function opens an input file (if the file exist), reading the header
   keywords, and verifying whether the required file has REDUCEME format and
   that the image dimensions do not exceed the maximum expected values (defined
   in NCMAX,NSMAX --see the file redlib.inc--). The function returns the name
   of the file opened. This function does NOT read the data records (this
   action must be performed after a call to this function).
   
   INTEGER       NF -> logical unit number of the file to be opened
   CHARACTER*(*) FILENAME -> default file name ('@' means there is not default)
   INTEGER       NSCAN -> no. of scans (pixels in the spatial direction)
   INTEGER       NCHAN -> no. of channels (pixels in the wavelength direction)
   REAL          STWV -> central wavelength of the first pixel
   REAL          DISP -> dispersion (Angstroms/pixel)
   INTEGER       MODE -> indicates the expected file format:
                 MODE=1: unformatted (with full header), i.e. REDUCEME format
                 MODE=2: unformatted (without header)
                 MODE=3: formatted (ascii file)
                 MODE=4: formatted, with RECL=2880 to read FITS files
                 MODE=11,12,13,14: like MODE=1,2,3,4 but, if the file exist, it
                 is opened directly (without prompting)
                 MODE=21: like MODE=1 but the function verifies whether the
                 input values of NSCAN,NCHAN,STWV, and DISP are identical with
                 those in the header of the file.
   LOGICAL       LERR -> if .TRUE. the input file corresponds to an error frame;
                 if .FALSE. the input file does not corresponds to an error
                 frame (only when MODE=1, 11 or 21; otherwise it has no effect)
   
   Apart from NSCAN, NCHAN, STWV, and DISP, other global variables (declared
   through COMMON blocks in redlib.inc), are also (re)declared: AIRMASS,
   TIMEXPOS, OBJECT, FITSFILE and COMMENT.

``my_pgend``
............

::

   SUBROUTINE MY_PGEND

   If the file .running_HLPHTML exist, this subroutine allows to capture the
   last XServe image before calling PGEND.

``outfilex``
............

::
   
   CHARACTER*75 FUNCTION OUTFILEX(NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR)

   Input: NF,FILENAME,NSCAN,NCHAN,STWV,DISP,MODE,LERR
   Input (COMMON): AIRMASS,TIMEXPOS,OBJECT,FITSFILE,COMMENT
   Output: OUTFILEX
   
   This function opens an output file (if the file does NOT exist), writing
   the header keywords (usually declared in a previous call to INFILEX).
   The function returns the name of the file opened. This function does NOT
   write the data records (this action must be performed after a call to this
   function).
   
   INTEGER       NF -> logical unit number of the file to be opened
   CHARACTER*(*) FILENAME -> default file name ('@' means there is not default)
   INTEGER       NSCAN -> no. of scans (pixels in the spatial direction)
   INTEGER       NCHAN -> no. of channels (pixels in the wavelength direction)
   REAL          STWV -> central wavelength of the first pixel
   REAL          DISP -> dispersion (Angstroms/pixel)
   INTEGER       MODE -> indicates the expected file format:
                 MODE=1: unformatted (with full header), i.e. REDUCEME format
                 MODE=2: unformatted (without header)
                 MODE=3: formatted (ascii file)
                 MODE=11,12,13: like MODE=1,2,3 but, if the file exist, it
                 is opened directly (without prompting)
   LOGICAL       LERR -> if .TRUE. the input file corresponds to an error frame;
                 if .FALSE. the input file does not corresponds to an error
                 frame (only when MODE=1, 11 or 21; otherwise it has no effect)
   
   Apart from NSCAN, NCHAN, STWV, and DISP, other global variables (declared
   through COMMON blocks in redlib.inc), are also saved: AIRMASS,
   TIMEXPOS, OBJECT, FITSFILE and COMMENT.

``pgiden_red``
..............

::

   SUBROUTINE PGIDEN_RED
   
   Input (COMMON): THISPROGRAM,CREDUCEVERSION
   
   Write the current program name, the username, date, time and REDUCEME
   version at the bottom of the plot.

``pidegter``
............

::
   
   SUBROUTINE PIDEGTER(NTERM,IDN,LCOLOR)

   Output: NTERM,IDN,LCOLOR
   
   Open the graphic device(s), detecting whether color is available.
   
   INTEGER NTERM -> No. of opened graphic devices (maximum = MAX_ID_RED)
   INTEGER IDN(MAX_ID_RED) -> logical device number associated to each NTERM
   LOGICAL LCOLOR(MAX_ID_RED) -> .TRUE. if color is available

``showhlp``
...........

::

   SUBROUTINE SHOWHLP(CADENA)

   Input: CADENA
   Output: CADENA
   
   Show additional help in programs at running time. This routine is employed
   for maintenance purposes (i.e. the creation of the help WEB page). When
   running any REDUCEME program, this routine checks whether any of the
   following two files exist:
   .running_HLP     -> generates HELP info in the terminal
   .running_HLPHTML -> generates HELP info to create a WEB page
   
   CHARACTER*(*) CADENA -> character string to identify the piece of information
                           to be extracted from the help file, which must be
                           located in the directory $reduceme_dir/help/programs

``welcome``
...........

::

   SUBROUTINE WELCOME(CSTRING)

   Input: CSTRING
   Input (COMMON): THISPROGRAM,CREDUCEVERSION
   
   Write the welcome presentation of the programs. This subroutine also
   verifies whether the environment variables PGPLOT_DIR and reduceme_dir have
   been defined.
   
   CHARACTER*(*) CSTRING -> additional information to be shown as a centered
                 character string in the welcome presentation

Libray ``libfutils.a``
----------------------

``chlower``
...........

::

   SUBROUTINE CHLOWER(CADENA)
   
   Input: CADENA
   Output: CADENA
   
   Upper case characters in CADENA are transformed to lower case
   
   CHARACTER*(*) CADENA -> character string to be transformed

``chupper``
...........

::

   SUBROUTINE CHUPPER(CADENA)
   
   Input: CADENA
   Output: CADENA
   
   Lower case characters in CADENA are transformed to upper case
   
   CHARACTER*(*) CADENA -> character string to be transformed

``endprogram``
..............

::

   SUBROUTINE ENDPROGRAM(CADENA)
   
   Input: CADENA
   
   Multipurpose routine: stops the program, inserts CALL PGPAGE or
   captures XServe.
   
   CHARACTER*255 CADENA -> if CADENA='endprogram' CALL PGEND+STOP
                           if CADENA='newpagenew' CALL PGPAGE
                           if CADENA='capturegif' capture current X11 image
   
``findmm``
..........

::

   SUBROUTINE FINDMM(N,X,XMIN,XMAX)
   
   Input: N,X
   Output: XMIN,XMAX
   
   Return the maximum and minimum value of matrix X of N elements
   
   INTEGER N -> no. of elements of matrix X
   REAL    X(N) -> data matrix
   REAL    XMIN -> minimum value of X()
   REAL    XMAX -> maximum value of X()
   
``findmml``
...........

::

   SUBROUTINE FINDMML(N,N1,N2,X,XMIN,XMAX)
   
   Input: N,N1,N2,X
   Output: XMIN,XMAX
   
   Return the maximum and minimum value of matrix X of N elements (in the
   range from N1 to N2 exclusively)
   
   INTEGER N -> no. of elements of matrix X
   INTEGER N1 -> first element of X() to search minimum/maximum
   INTEGER N2 -> last element of X() to search minimum/maximum
   REAL    X(N) -> data matrix
   REAL    XMIN -> minimum value of X()
   REAL    XMAX -> maximum value of X()

``findmmlmask``
...............

::

   SUBROUTINE FINDMMLMASK(N,N1,N2,X,LMASK,XMIN,XMAX)
   
   Input: N,N1,N2,X,LMASK
   Output: XMIN,XMAX
   
   Return the maximum and minimum value of matrix X of N elements (in the
   range from N1 to N2 exclusively, making use of the boolean mask LMASK)
   
   INTEGER N -> no. of elements of matrix X
   INTEGER N1 -> first element of X() to search minimum/maximum
   INTEGER N2 -> last element of X() to search minimum/maximum
   REAL    X(N) -> data matrix
   LOGICAL LMASK(N) -> boolean mask
   REAL    XMIN -> minimum value of X()
   REAL    XMAX -> maximum value of X()

``lrunx``
.........

::

   SUBROUTINE LRUNX(LRUN,LMANUAL,LHTML)
   
   Output: LRUN,LMANUAL,LHTML
   
   Determine whether files .running_RUN, .running_MANUAL and .running_HLPHTML
   exist in current the directory.
   
   LOGICAL LRUN -> .TRUE. if file .running_RUN exists (.FALSE. otherwise)
   LOGICAL LMANUAL -> .TRUE. if file .running_MANUAL exists (.FALSE. otherwise)
   LOGICAL LHTML -> .TRUE. if file .running_HLPHTML exists (.FALSE. otherwise)

``read2i``
..........

::

   SUBROUTINE READ2I(CDEF,N1,N2)
   
   Input: CDEF
   Output: N1,N2
   
   Return 2 integers N1 and N2 entered by the user through the keyboard.
   
   CHARACTER*(*) CDEF -> character string with default values for N1 and N2
                 ('@' if there is no default)
   INTEGER N1 -> first integer
   INTEGER N2 -> second integer

``readc``
.........

::

   CHARACTER*(*) FUNCTION READC(CDEF,CVAL)
   
   Input: CDEF,CVAL
   Output: READC (function)
   
   Return a character string entered by the user through the keyboard.
   
   CHARACTER*(*) CDEF -> character with default value for READC
                 ('@' if there is no default)
   CHARACTER*(*) CVAL -> character string with valid characters
                 ('@' if all characters are valid)

``readf``
.........

::

   REAL FUNCTION READF(CDEF)
   
   Input: CDEF
   Output: READF (function)
   
   Return a float number entered by the user through the keyboard.
   
   CHARACTER*(*) CDEF -> character string with default value for READF
                 ('@' if there is no default)

``readi``
.........

::

   INTEGER FUNCTION READI(CDEF)
   
   Input: CDEF
   Output: READI (function)
   
   Return an integer number entered by the user through the keyboard.
   
   CHARACTER*(*) CDEF -> character string with default value for READI
                 ('@' if there is no default)

``readilim``
............

::

   INTEGER FUNCTION READILIM(CDEF,N1,N2)
   
   Input: CDEF,N1,N2
   Output: READILIM (function)
   
   Return an integer number entered by the user through the keyboard in the
   range from N1 to N2
   
   CHARACTER*(*) CDEF -> character string with default value for READILIM
                 ('@' if there is no default)
   INTEGER       N1 -> first limit for READILIM
   INTEGER       N2 -> second limit for READILIM

``rmblank``
...........

::

   SUBROUTINE RMBLANK(C1,C2,L)
   
   Input: C1
   Output: C2,L
   
   Remove blanks in character string C1, returning C2 with a true length L
   
   CHARACTER*(*) C1 -> input character string
   CHARACTER*(*) C2 -> output character string (C1 without blanks)
   INTEGER       L -> true len of C2

``showperc``
............

::

   SUBROUTINE SHOWPERC(N1,N2,ISTEP,I,NEXTINFO)
   
   Input: N1,N2,ISTEP,I,NEXTINFO
   Output: NEXTINFO
   
   Display the percentage of work performed in a loop, which has been defined
   in the range from N1 to N2, with an incremental step ISTEP, being I the
   current value of the loop control.
   
   INTEGER N1 -> first limit of the control variable of the loop
   INTEGER N2 -> last limit of the control variable of the loop
   INTEGER ISTEP -> the value by which the control variable is incremented
   INTEGER I -> current value of the control variable
   INTEGER NEXTINFO -> integer which stores the last fraction of 10 displayed;
                       this variable is initialized in the first call of this
                       routine

``truebeg``
...........

::

   INTEGER FUNCTION TRUEBEG(CADENA)

   Input: CADENA
   Output: TRUEBEG (function)
   
   Return the position of the first non-blank character in CADENA (ignoring
   also control characters with ASCII value < 32)
   
   CHARACTER*(*) CADENA -> input character string

``truelen``
...........

::

   INTEGER FUNCTION TRUELEN(CADENA)
   
   Input: CADENA
   Output: TRUELEN (function)
   
   Return the position of the last non-blank character in CADENA (ignoring also
   control characters with ASCII value < 32)
   
   CHARACTER*(*) CADENA -> input character string

Library ``libfspec.a``
----------------------

``bicubspl``
............

::

   SUBROUTINE BICUBSPL(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC)

   Input: X,Y,Z,NX,NY,NXDIM,NYDIM
   Output: AA,BB,CC
   
   This subroutine computes the spline coefficients required by the subroutine
   BICUBSPLX to compute a surface through bicubic spline interpolation.
   
   REAL X(NX) -> X-values to be fitted
   REAL Y(NY) -> Y-values to be fitted
   REAL Z(NX,NY) -> Z-values to be fitted
   INTEGER NX -> logical dimension of X and logical first dimension of Z
   INTEGER NY -> logical dimension of Y and logical second dimension of Z
   INTEGER NXDIM -> physical dimension of X and physical first dimension of Z
   INTEGER NYDIM -> pyshical dimension of Y and physical second dimension of Z
   REAL AA(NX,NY) -> spline coefficients of the one-dimensional cubic spline
                     fits to the rows of Z
   REAL BB(NX,NY) -> spline coefficients of the one-dimensional cubic spline
                     fits to the rows of Z
   REAL CC(NX,NY) -> spline coefficients of the one-dimensional cubic spline
                     fits to the rows of Z
   
``bicubsplx``
.............

::
   
   SUBROUTINE BICUBSPLX(X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC,X0,Y0,Z0)
   Input: X,Y,Z,NX,NY,NXDIM,NYDIM,AA,BB,CC,X0,Y0
   Output: Z0
   
   This subroutine computes the bicubic spline Z0 at X0,Y0, using the spline
   coefficients computed with BICUBSPL.
   
   REAL X(NX) -> X-values fitted with BICUBSPL
   REAL Y(NY) -> Y-values fitted with BICUBSPL
   REAL Z(NX,NY) -> Z-values fitted with BICUBSPL
   INTEGER NX -> logical dimension of X and logical first dimension of Z
   INTEGER NY -> logical dimension of Y and logical second dimension of Z
   INTEGER NXDIM -> physical dimension of X and physical first dimension of Z
   INTEGER NYDIM -> pyshical dimension of Y and physical second dimension of Z
   REAL AA(NX,NY) -> coefficients of the one-dimensional cubic spline fits
                    to the rows of Z computed with BICUBSPL
   REAL BB(NX,NY) -> coefficients of the one-dimensional cubic spline fits
                    to the rows of Z computed with BICUBSPL
   REAL CC(NX,NY) -> coefficients of the one-dimensional cubic spline fits
                    to the rows of Z computed with BICUBSPL
   REAL X0 -> X-value where the bicubic spline will be evaluated
   REAL Y0 -> Y-value where the bicubic spline will be evaluated
   REAL Z0 -> bicubic spline value at X0,Y0
   
``binsearch``
.............

::
   
   SUBROUTINE BINSEARCH(X,N,X0,N0)

   Input: X,N,X0,N0
   Output: N0
   
   Given the array X(N), and the test value X0, this subroutine returns an
   integer N0, such that X0 is between X(N0) and X(N0+1). As input N0 is
   employed to start the searching. If X0.LT.X(1) then N0=0 on output, whereas
   if X0.GT.X(N) then N0=N. If X0.EQ.X(K), N0=K on output.
   
   REAL    X(N) -> ordered input array (not necesarilly equally-spaced)
   INTEGER N -> no. of points in input array
   REAL    X0 -> argument to be searched for
   INTEGER N0 -> location of X0 in the input array
   
``broaden``
...........

::
   
   SUBROUTINE BROADEN(S1,S2,NCHAN,STWV,DISP,SIGMA,LERR)

   Input: S1,NCHAN,STWV,DISP,SIGMA,LERR
   Output: S2
   
   Broadens a single spectrum by convolving with a gaussian (variable along the
   wavelength scale).
   
   REAL    S1(NCMAX) -> input spectrum
   REAL    S2(NCMAX) -> output spectrum
   INTEGER NCHAN -> no. of channels
   REAL    STWV -> central wavelength of the first channel
   REAL    DISP -> dispersion (Angstrom/pixel)
   REAL    SIGMA(NCMAX) -> sigma value of the gaussian to be applied in each
                           channel
   LOGICAL LERR -> if .TRUE. S1 corresponds to an error spectrum
   
``cauchyfit``
.............

::
   
   SUBROUTINE CAUCHYFIT(X0,SIGMA,AMP,EEX0,EESIGMA,EEAMP,YRMSTOL)

   Input: YRMSTOL
   Input (COMMON): NP,XF,YF
   Output: X0,SIGMA,AMP,EEX0,EESIGMA,EEAMP
   
   Fit numerically a Cauchy function (using DOWNHILL):
   Y=AMP/[SIGMA^2+(X-X0)^2]
   
   REAL X0 -> center of the fitted Cauchy function
   REAL SIGMA -> sigma value of the fitted Cauchy function
   REAL AMP -> maximum of the fitted Cauchy function
   REAL EEX0,EESIGMA,EEAMP -> errors in X0,SIGMA,AMP (rms from DOWNHILL)
   REAL YRMSTOL -> stopping criterion for DOWNHILL
   
``cfftd``
.........

::
   
   SUBROUTINE CFFTD(N,XR,XI,IMODE)

   Input N,XR,XI,IMODE
   Output XR,YR
   
   If IMODE=1, this subroutine computes the FFT of the input complex vector
   XR + i XI, where XR and XI are both real variables. As output XR and XI are
   the real and imaginary part of the transform. If IMODE=-1 the subroutine
   evaluates the inverse FFT. See E. O. Brigham, The Fast Fourier Transform,
   pag.160.
   
   INTEGER N  -> number of points (must be a power of 2)
   DOUBLE PRECISION XR(N) -> real part of the input data
   DOUBLE PRECISION XI(N) -> imaginary part of the input data
   INTEGER IMODE -> +1 (direct FFT) or -1 (inverse FFT)
   
``cfft``
........

::
   
   SUBROUTINE CFFT(N,XR,XI,IMODE)

   Input N,XR,XI,IMODE
   Output XR,YR
   
   If IMODE=1, this subroutine computes the FFT of the input complex vector
   XR + i XI, where XR and XI are both real variables. As output XR and XI are
   the real and imaginary part of the transform. If IMODE=-1 the subroutine
   evaluates the inverse FFT. See E. O. Brigham, The Fast Fourier Transform,
   pag.160.
   
   INTEGER N  -> number of points (must be a power of 2)
   REAL XR(N) -> real part of the input data
   REAL XI(N) -> imaginary part of the input data
   INTEGER IMODE -> +1 (direct FFT) or -1 (inverse FFT)
   
``chequea_fileindex``
.....................
   
::
   
   SUBROUTINE CHEQUEA_FILEINDEX
   
   This subroutine verifies that there is an appropiate file containing
   the index definitions.
   
   
``chrebin``
...........

::
   
   SUBROUTINE CHREBIN(CHANSHIFT,NCHAN,S,SS)

   Input: CHANSHIFT,NCHAN,S
   Output: SS
   
   This subroutine applies a constant channel shift to a spectrum.
   
   REAL CHANSHIFT -> channel shift to be applied
   INTEGER NCHAN  -> number of channels
   REAL S(NCHAN)  -> initial spectrum
   REAL SS(NCHAN) -> shifted spectrum
   
``combpf``
..........

::
   
   DOUBLE PRECISION FUNCTION COMBPF(N,K)

   Input: N,K
   Output: COMBPF (function)
   
   Calculate the binomial coefficient N over K
   
   INTEGER N
   INTEGER K
   
``cubspl``
..........

::
   
   SUBROUTINE CUBSPL(X,Y,N,IMODE,S,A,B,C)

   Input: X,Y,N,IMODE,S
   Output: A,B,C,S
   
   This subroutine computes the coefficients of a cubic spline. See C.F. Gerald
   and P. O. Wheatley, in Applied Numerical Analysis, 4th edition, pag. 207.
   The subroutine returns the spline coefficients, where the spline defined
   in the interval between X(I),Y(I) and X(I+1),Y(I+1) is given by:
   
        Y = A(I)*(X-X(I))**3 + B(I)*(X-X(I))**2 + C(I)*(X-X(I)) + D(I)
   
   REAL X(N) -> X-values to be fitted
   REAL Y(N) -> Y-values to be fitted
   INTEGER N -> number of data points
   INTEGER IMODE -> End conditions mode: if S(I) represent the second derivative
                    at the point X(I),Y(I), the following four possibilites
                    are available:
                    1) IMODE=1: S(0)=0, S(N)=0. This is called natural cubic
                       spline. It is equivalent to assuming that the end cubics
                       aproach linearity at their extremities.
                    2) IMODE=2: S(0)=S(1), S(N)=S(N-1). This is equivalent to
                       assuming that the cubics approach parabolas at their
                       extremities.
                    3) IMODE=3: S(0) is a linear extrapolation from S(1) and
                       S(2), and S(N) is a linear extrapolation from S(N-2)
                       and S(N-1).
                    4) IMODE=4: Force the slopes at each end to assume certain
                       values.
   REAL S(N) -> if IMODE=4, in input S(1) and S(N) contain the first derivatives
                at X(1) and X(N). In output, this matrix contains the second
                derivatives
   REAL A(N) -> spline coefficients
   REAL B(N) -> spline coefficients
   REAL C(N) -> spline coefficients
   
``cubsplx``
...........

::
   
   SUBROUTINE CUBSPLX(X,Y,A,B,C,N,I0,X0,Y0)

   Input: X,Y,A,B,C,N,I0,X0
   Output: Y0
   
   The subroutine returns the cubic spline evaluated at X0, using the
   coefficients determined in a previous call to CUBSPL. The spline defined in
   the interval between X(I),Y(I) and X(I+1),Y(I+1) is given by:
   
        Y = A(I)*(X-X(I))**3 + B(I)*(X-X(I))**2 + C(I)*(X-X(I)) + D(I)
   
   If X0.LT.X(1), I=1 is employed (first computed spline)
   If X0.GT.X(N), I=N-1 is employed (last computed spline)
   
   REAL X(N) -> X-values fitted with CUBSPL
   REAL Y(N) -> Y-values fitted with CUBSPL
   REAL A(N) -> spline coefficients
   REAL B(N) -> spline coefficients
   REAL C(N) -> spline coefficients
   INTEGER N -> number of data points
   INTEGER I0 -> initial location to start the search of the place of X0 in
                 the X array
   REAL X0 -> X-value where the spline function will be evaluated
   REAL Y0 -> spline value at X0
   
``downhill``
............

::
   
   SUBROUTINE DOWNHILL(N,X0,DX0,YFUNK,A,B,G,YRMSTOL,XF,DXF,NEVAL)

   Input N,X0,DX0,YFUNK,A,B,G,YRMSTOL
   Output XF,NEVAL
   
   Minimization of the function YFUNK of N variables using the downhill
   simplex method, as explained by Nelder and Mead (1965, Computer Journal, 7,
   pags. 308-313). The routine returns when the stopping criterion is reached
   (the r.m.s. of the YFUNK values computed with all the vertices of the simplex
   is .LT. YRMSTOL), or when the number of function evaluations
   is too large (NEVAL.GT.NEVALMAX).
   
   INTEGER N -> number of variables
   REAL    X0(N) -> starting point (initial solution)
   REAL    DX0(N) -> N characteristic length scales, employed to derive N
                    aditional starting points which, together with X0, form
                    the (N+1) vertices of the simplex
   REAL    YFUNK -> function to be minimized
   REAL    A -> reflection coefficient (ALPHA); tipically, A=1.0
   REAL    B -> contraction coefficient (BETA); tipically, B=0.5
   REAL    G -> expansion coefficient (GAMMA); tipically, G=2.0
   REAL    YRMSTOL -> stopping criterion
   REAL    XF(N) -> final solution
   REAL    DXF(N) -> rms of XF evaluated from the final different points of the
                     simplex
   INTEGER NEVAL -> number of evaluations of YFUNK employed by DOWNHIL to reach
                    the solution; the routine returns NEVAL=-1 if something goes
                    wrong
   
``factorialpf``
...............

::
   
   DOUBLE PRECISION FUNCTION FACTORIALPF(N)

   Input: N
   Output: FACTORIALPF (function)
   
   Calculate N factorial
   
   INTEGER N
   
``fft2power``
.............

::
   
   SUBROUTINE FFT2POWER(N0,N)

   Input: N0
   Output: N
   
   Given an integer N0, this subroutine asks for a power of 2, such as
   N=2**K, being K integer and N.GE.N0. This value of N is employed by other
   subroutines to perform zero padding prior computing FFT.
   
``fftcorrel``
.............

::
   
   SUBROUTINE FFTCORREL(N,DATA1,DATA2,XCORR,FCORR)

   Input N,DATA1,DATA2
   Output XCORR,FCORR
   
   Compute the correlation function FCORR of two real data vectors DATA1 and
   DATA2, using FFT. We assume that both data sets have been properly filtered
   and dimensioned. We use the discrete correlation theorem, which says that
   the discrete correlation of two real functions f1 and f2 is one member of
   the discrete Fourier transform pair:
                         Corr(f1,f2) <==> F1 F2*
   where F1 and F2 are the discrete Fourier transforms of f1 and f2,
   respectively, and asterisk denotes complex conjugation.
   Note that DATA1 and DATA2 are not modified by this routine.
   
   INTEGER N  -> number of points (must be a power of 2)
   REAL DATA1(N) -> first data set to be correlated
   REAL DATA2(N) -> second data set to be correlated
   REAL XCORR(N) -> abcissas of the output correlation function
   REAL FCORR(N) -> output correlation function
   
``fftcorrzoom``
...............

::
   
   SUBROUTINE FFTCORRZOOM(X,Y,N,NL,NAME1,NAME2,NFIT,LPLOT,X0,Y0)

   Input: X,Y,N,NL,NAME1,NAME2,NFIT,LPLOT
   Output: X0,Y0
   
   This subroutine plots the correlation function given by X(N),Y(N), zooms in
   around the maximum and fit a second-order polynomial in order to obtain the
   exact location of the peak. The fitted value is returned through X0,Y0.
   
   REAL X(N) -> X-coordinates of the correlation function
   REAL Y(N) -> Y-coordinates of the correlation function
   INTEGER N -> dimension of X and N
   INTEGER NL -> number of pixels affected by zero padding
   CHARACTER*(*) NAME1 -> description of the first data set employed to compute
                          cross correlation
   CHARACTER*(*) NAME2 -> description of the second data set employed to compute
                          cross correlation
   INTEGER NFIT -> no. of points around peak to fit maximum (if NFIT=0, the
                   routine asks for this number; otherwise the routine returns
                   without prompting)
   LOGICAL LPLOT -> if .TRUE., plot zoomed region
   REAL X0 -> X-offset of the peak of the correlation function
   REAL Y0 -> peak value of the correlation function
   
``fftcosbell``
..............

::
   
   SUBROUTINE FFTCOSBELL(N,COSBELL,FL)

   Input: N, FL
   Output: COSBELL
   
   Compute a cosine bell expanding N pixels, being FL the fraction of pixels (at
   the beginning and at the end of the cosine bell) employed to perform the
   transition from zero to one. See Brault & White, A&A, 13, 169.
   
   INTEGER N -> dimension of COSBELL (pixels)
   REAL COSBELL(N) -> cosine bell
   REAL FL -> fraction of pixels at the borders of the cosine bell
   
``fftkfilter``
..............

::
   
   SUBROUTINE FFTKFILTER(N,KFILTER,K1,K2,K3,K4)

   Input: N,K1,K2,K3,K4
   Output: KFILTER
   
   Given K1, K2, K3 and K4, this subroutine computes the filter KFILTER, such
   as:
   
   KFILTER(I)=0., I=1,...,K1
   KFILTER(I) is a straight line from 0. to 1., I=K1,...,K2
   KFILTER(I)=1., I=K2,...,K3
   KFILTER(I) is a straight line from 1. to 0., I=K3,...,K4
   KFILTER(I)=0., I=K4,...,N/2+1
   
   KFILTER(I)=0., I=N,...,N-K1+2
   KFILTER(I) is a straight line from 0. to 1., I=N-K1+2,...,N-K2+2
   KFILTER(I)=1., I=N-K2+2,...,N-K3+2
   KFILTER(I) is a straight line from 1. to 0., I=N-K3+2,...,N-K4+2
   KFILTER(I)=0., I=N-K4+2,...,N/2+1
   
   N must be a power of 2, and 1.le.K1.le.K2.le.K3.le.K4.le.(N/2-1).
   This filter is employed by other routines to perform a frequency filter
   in the frequency domain of the FFT.
   
``fftprep``
...........

::
   
   SUBROUTINE FFTPREP(N,SP,NC1,NC2,LNORM,COSBELL,LFILT,KFILTER,LPLOT,CNAME)

   Input: N,SP,NC1,NC2,LNORM,COSBELL,LFILT,KFILTER,LPLOT,CNAME
   Output: SP
   
   This subroutine prepares spectrum SP for cross correlation. For this purpose
   the following steps are performed:
   - normalization of the spectrum, using only the data in the range
     [SP(NC1),SP(NC2)]
   - extraction of the useful spectral region: [SP(NC1),SP(NC2)] is converted
     into [SP(1),SP(NC2-NC1+1)]
   - multiplication of the spectrum by a cosine bell defined from 1 to NC2-NC1+1
   - zero padding from NC2-NC1+2 to N, where N=2**K, being K an integer
   - computation of the Fast Fourier Transform of SP(1),...,SP(N) and obtention
     of the power spectrum
   - filtering of the power spectrum by applying the user defined filter KFILTER
   - computation of the inverse FFT of the filtered power spectrum and obtention
     of the filtered spectrum SP(1),...,SP(N)
   - extraction of the initial spectral region: SP(1),...,SP(NC2-NC1+1)
   - multiplication by the cosine bell defined from 1 to NC2-NC1+1 (again!)
   - zero padding from NC2-NC1+2 to N (again!)
   
   INTEGER N -> number of points in output (must be a power of 2)
   REAL SP(N) -> spectrum to be prepared for cross correlation
   INTEGER NC1 -> index which indicates the first SP value to be employed
   INTEGER NC2 -> index which indicates the last SP value to be employed
   LOGICAL LNORM -> if .TRUE. SP is normalized prior any other manipulation
   REAL COSBELL(N) -> cosine bell (defined from point number 1 to NC2-NC1+1)
   LOGICAL LFILT -> if .TRUE. SP is filtered applying the filter KFILTER
   REAL KFILTER(N) -> filter (defined from point number 1 to N)
   LOGICAL LPLOT -> if .TRUE. plot intermediate plots
   CHARACTER*(*) CNAME -> string description for plots
   
``fintgauss`` and ``fintgausse``
................................

::
   
   REAL FUNCTION FINTGAUSS(X1,X2,N,X0,FACTOR1)  -> integral of Gaussian

   REAL FUNCTION FINTGAUSSE(X1,X2,N,X0,FACTOR1) -> integral of Gaussian^2

   Input: X1,X2,N,X0,FACTOR1
   Output: FINTGAUSS (function) or FINTGAUSSE (function)
   
   Calculate the integral of a Gaussian of width SIGMA, center at wavelength X0,
   using Simpson's rule, between wavelengths X1 and X2, with N intervals
   (N must be even).
   
   REAL X1,X2 -> integral limits (wavelengths)
   INTEGER N -> number of intervals to estimate integral
   REAL X0 -> Gaussian center (wavelength)
   REAL FACTOR1 -> = (c^2/(2*SIGMA^2))/(X0^2)
   
``fmean0``
..........

::
   
   REAL FUNCTION FMEAN0(N,X,SIGMA)

   Input: N,X
   Output: FMEAN0 (function),SIGMA
   
   Calculate the mean value of X(N) and its r.m.s.
   
   INTEGER N -> no. of elements
   REAL    X(N) -> input matrix
   REAL SIGMA -> r.m.s. around the mean value
   
``fmean1``
..........

::
   
   REAL FUNCTION FMEAN1(N,X)

   Input: N,X
   Output: FMEAN1 (function)
   
   Calculate the mean value of X(N)
   
   INTEGER N -> no. of elements
   REAL    X(N) -> input matrix
   
``fmean2``
..........

::
   
   REAL FUNCTION FMEAN2(N,X,TIMES)

   Input: N,X,TIMES
   Output: FMEAN2 (function)
   
   Calculate the mean value of X(N) rejecting points at TIMES sigma. The
   function can recover points rejected in previous iterations.
   
   INTEGER N -> no. of elements
   REAL    X(N) -> input matrix
   REAL    TIMES -> times sigma to reject points before calculating the mean
   
``fmedian1``
............

::
   
   REAL FUNCTION FMEDIAN1(N,X)

   Input: N,X
   Output: FMEDIAN (function), X (sorted)
   
   Calculate the median value of X(N). It is important to note that this
   subroutine rearranges the matrix X which is returned sorted.
   
   INTEGER N -> no. of elements
   REAL    X(N) -> input matrix
   
``fpercent``
............

::
   
   REAL FUNCTION FPERCENT(N,X,PERCENTILE)

   Input: N,X,PERCENTILE
   Output: FPERCENT (function)
   
   Calculate a fixed percentile of X(N)
   
   INTEGER N -> no. of elements
   REAL    X(N) -> input matrix
   REAL    PERCENTILE -> percentile to be computed
   
``fpoly``
.........

::
   
   REAL FUNCTION FPOLY(NDEG,COEFF,X)

   Input: NDEG,COEFF,X
   Output: FPOLY (function)
   
   Evaluate the polynomial of degree NDEG and coefficients COEFF at X.
   
   INTEGER NDEG -> polynomial degree
   REAL    COEFF(NDEG+1) -> polynomial coefficients
   REAL    X -> abscissa at which the polynomial is going to be evaluated
   
``gauscfit``
............

::
   
   SUBROUTINE GAUSCFIT(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,YRMSTOL)

   Input: YRMSTOL
   Input (COMMON): NP,XF,YF
   Output: X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAM,EEY0
   
   Fit numerically a gaussian + constant (using DOWNHILL):
   Y=Y0+AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
   
   REAL X0 -> center of the fitted gaussian
   REAL SIGMA -> sigma value of the fitted gaussian
   REAL AMP -> maximum of the fitted gaussian
   REAL Y0 -> constant
   REAL EEX0,EESIGMA,EEAMP,EEY0 -> errors in X0,SIGMA,AMP,Y0 (rms from DOWNHILL)
   REAL YRMSTOL -> stopping criterion for DOWNHILL
   
``gauscfit_movel``
..................

::
   
   SUBROUTINE GAUSCFIT_MOVEL(X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAMP,EEY0,YRMSTOL)

   Input: YRMSTOL
   Input (COMMON): NP,XF,YF
   Output: X0,SIGMA,AMP,Y0,EEX0,EESIGMA,EEAM,EEY0
   
   Fit numerically a gaussian + constant (using DOWNHILL):
   Y=Y0+AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
   
   REAL X0 -> center of the fitted gaussian
   REAL SIGMA -> sigma value of the fitted gaussian
   REAL AMP -> maximum of the fitted gaussian
   REAL Y0 -> constant
   REAL EEX0,EESIGMA,EEAMP,EEY0 -> errors in X0,SIGMA,AMP,Y0 (rms from DOWNHILL)
   REAL YRMSTOL -> stopping criterion for DOWNHILL
   
   This subroutine is identical to GAUSCFIT, but the physical dimensions
   of the data to be fitted is defined to be NMAXFFT. This change is
   required for the program movel.
   
``gauss2afit``
..............

::
   
   SUBROUTINE GAUSS2AFIT(NPFIT,XFIT,YFIT,EYFIT,DELTAX,X0,SIGMA,AMP,
                         EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,
                         YRMSTOL,NSIMUL)

   Input: NPFIT,XFIT,YFIT,EYFIT,DELTAX,YRMSTOL,NSIMUL
   Output: X0,SIGMA,AMP,EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP
   
   Numerical fit of 2 gaussians with the same width and area (using DOWNHILL),
   with a fixed separation given by DELTAX.
   Y=AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]+AMP*EXP[-((X-X0-DELTAX)^2/(2*SIGMA^2))]
   
   INTEGER NPFIT -> number of points to be fitted
   REAL XFIT,YFIT,EYFIT -> x, y and error
   REAL DELTAX -> separation between gaussians
   REAL X0 -> center of the fitted gaussian
   REAL SIGMA -> sigma value of the fitted gaussian
   REAL AMP -> maximum of the fitted gaussian
   REAL EX0,ESIGMA,EAMP -> errors in X0,SIGMA,AMP (due to EYFIT --simulations--)
   REAL EEX0,EESIGMA,EEAMP -> errors in X0,SIGMA,AMP (rms from DOWNHILL)
   REAL YRMSTOL -> stopping criterion for DOWNHILL
   INTEGER NSIMUL -> number of simulations to compute errors
   
``gauss2bfit``
..............

::
   
   SUBROUTINE GAUSS2BFIT(NPFIT,XFIT,YFIT,EYFIT,DELTAX,X0,SIGMA,AMP1,AMP2,
                         EX0,ESIGMA,EAMP1,EAMP2,EEX0,EESIGMA,EEAMP1,EEAMP2,
                         YRMSTOL,NSIMUL)

   Input: NPFIT,XFIT,YFIT,EYFIT,DELTAX,YRMSTOL,NSIMUL
   Output: X0,SIGMA,AMP1,AMP2,EX0,ESIGMA,EAMP1,EAMP2,EEX0,EESIGMA,EEAMP1,EEAMP2
   
   Numerical fit of 2 gaussians with the same width and different area
   (using DOWNHILL), with a fixed separation given by DELTAX.
   Y=AMP1*EXP[-((X-X0)^2/(2*SIGMA^2))]+AMP2*EXP[-((X-X0-DELTAX)^2/(2*SIGMA^2))]
   
   INTEGER NPFIT -> number of points to be fitted
   REAL XFIT,YFIT,EYFIT -> x, y and error
   REAL DELTAX -> separation between gaussians
   REAL X0 -> center of the fitted gaussian
   REAL SIGMA -> sigma value of the fitted gaussian
   REAL AMP1 -> maximum of the fitted gaussian #1
   REAL AMP2 -> maximum of the fitted gaussian #2
   REAL EX0,ESIGMA,EAMPn-> errors in X0,SIGMA,AMPn (due to EYFIT --simulations--)
   REAL EEX0,EESIGMA,EEAMPn-> errors in X0,SIGMA,AMPn (rms from DOWNHILL)
   REAL YRMSTOL -> stopping criterion for DOWNHILL
   INTEGER NSIMUL -> number of simulations to compute errors
   
``gaussfitamp``
...............

::
   
   SUBROUTINE GAUSSFITAMP(NPFIT,XFIT,YFIT,EYFIT,X0,SIGMA,EX0,ESIGMA,MINSIGMA,
                          AMP,EAMP,EEAMP,NSIMUL)

   Input: NPFIT,XFIT,YFIT,EYFIT,X0,SIGMA,EX0,ESIGMA,MINSIGMA,NSIMUL
   Output: AMP,EAMP,EEAMP
   
   Fit of AMP of a gaussian with X0 and SIGMA fixed:
   Y=AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
   
   INTEGER NPFIT -> number of points to be fitted
   REAL XFIT,YFIT,EYFIT -> x, y and error
   REAL X0, EX0 -> center of the gaussian and its error
   REAL SIGMA, ESIGMA -> sigma value of the gaussian and its error
   REAL MINSIGMA -> minimum SIGMA allowed in simulations (typically MINSIGMA
                    must be the spectral resolution)
   REAL AMP -> maximum of the fitted gaussian
   REAL EAMP -> error in SIGMA (due to EYFIT)
   REAL EEAMP -> error in SIGMA (due to EX0 and ESIGMA ---simulations---)
   INTEGER NSIMUL -> number of simulations to compute EEAMP
   
``gaussfit``
............

::
   
   SUBROUTINE GAUSSFIT(NPFIT,XFIT,YFIT,EYFIT,X0,SIGMA,AMP,
                       EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP,YRMSTOL,NSIMUL)

   Input: NPFIT,XFIT,YFIT,EYFIT,YRMSTOL,NSIMUL
   Output: X0,SIGMA,AMP,EX0,ESIGMA,EAMP,EEX0,EESIGMA,EEAMP
   
   Numerical fit of a gaussian (using DOWNHILL):
   Y=AMP*EXP[-((X-X0)^2/(2*SIGMA^2))]
   
   INTEGER NPFIT -> number of points to be fitted
   REAL XFIT,YFIT,EYFIT -> x, y and error
   REAL X0 -> center of the fitted gaussian
   REAL SIGMA -> sigma value of the fitted gaussian
   REAL AMP -> maximum of the fitted gaussian
   REAL EX0,ESIGMA,EAMP -> errors in X0,SIGMA,AMP (due to EYFIT --simulations--)
   REAL EEX0,EESIGMA,EEAMP -> errors in X0,SIGMA,AMP (rms from DOWNHILL)
   REAL YRMSTOL -> stopping criterion for DOWNHILL
   INTEGER NSIMUL -> number of simulations to compute errors
   
``integtab``
............

::
   
   REAL FUNCTION INTEGTAB(N,X,Y,X1,X2,IFLAG1,IFLAG2)

   Input: N,X,Y,X1,X2
   Output: INTEGTAB(function), IFLAG1,IFLAG2
   
   Performs the integration of a function given in a tabular form, between the
   limits X1 and X2. Note that the X matrix must be sorted in ascending order.
   
   INTEGER N -> input number of data in X and Y
   REAL    X(N) -> data matrix
   REAL    Y(N) -> data matrix
   REAL    X1 -> first limit of the integral
   REAL    X2 -> second limit of the integral
   INTEGER IFLAG1 -> = 0 : interpolation
                     = -1 : extrapolation towards lower X values
                     = +1 : extrapolation towards higher X values
                     = +9 : error (division by zero)
   INTEGER IFLAG2 -> = 0 : interpolation
                     = -1 : extrapolation towards lower X values
                     = +1 : extrapolation towards higher X values
                     = +9 : error (division by zero)
   
``lagrange``
............

::
   
   SUBROUTINE LAGRANGE(N,X,Y,X0,Y0)

   Input: N,X,Y,X0
   Output: Y0
   
   Calculate the Lagrangian polynomial and evaluate such polynomial at X=X0.
   We do not assume uniform spacing between the x-values, nor do we need the
   x-values arranged in a particular order. However, the x-values must all be
   distinct. We follow the algorithm described by B.P. Demidovich and I.A. Maron
   in Calculo Numerico Fundamental, Paraninfo 1988, pag. 593.
   
   INTEGER N -> no. of elements
   REAL    X(N) -> input x-matrix
   REAL    Y(N) -> input y-matrix
   REAL    X0 -> x-value where the Lagrangian polynomial is evaluated
   REAL    Y0 -> polynomial value at X=X0
   
``fitl``
........

::
   
   SUBROUTINE  FITL(X,Y,NX,sig,iw,xmn,xmx,A,B,SIGA,SIGB,STD)

   -----------------------------------------------------------------------C
   Linear Fit Routines                          J. Jesus Gonzalez G.
   -----------------------------------------------------------------------C
   -----------------------------------------------------------------------C
                        Fits the line y = b*x + a
       THIS ROUTINE ITERATES ELIMINATING HIGHLY DEVIANT POINTS
       INPUT:  X,   Y - Data arrays
                  SIG - Y-error of points (<=0 if a point is to be rejected)
                   IW - =0 if unweighted fit, weighted fit otherwise.
              xmn,xmx - Limits of fit.
   
       OUTPUT:
              B,    A - Slope and zero-ordinate.
           SIGB, SIGA - Stimated errors on B and A.
                  STD - Unbiased Std-Deviation from the fit.
   
       Remeber how to compute errors of predicted values:
       VAR(y(x)) = VAR(a) + VAR(b)*(x^2 - 2*xm*x), since the
       ab-error covariance is COV(ab)=-xm*VAR(b)
   
   -----------------------------------------------------------------------C

``lininterp``
.............

::
   
   REAL FUNCTION LININTERP(N,X,Y,X0,IFLAG,N1,N2)

   Input: N,X,Y,X0
   Output: LININTERP(function), IFLAG
   
   Performs a linear interpolation in the table X(N),Y(N) at x=X0. Note that the
   X matrix must be sorted in ascending order, although the
   
   INTEGER N -> input number of data in X and Y
   REAL    X(N) -> data matrix
   REAL    Y(N) -> data matrix
   REAL    X0 -> x-point at which the linear interpolation is evaluated
   INTEGER IFLAG -> = 0 : interpolation
                    = -1 : extrapolation towards lower X values
                    = +1 : extrapolation towards higher X values
                    = +2 : X0=X(N), which could produce some "border" effects
                    = +9 : error (division by zero)
   INTEGER N1 -> first data entry towards the left of X0
   INTEGER N2 -> first data entry towards the right of X0
   
``ludcmp``
..........

::
   
   SUBROUTINE LUDCMP(A,N,NDIM,ORDER,SCALEROW,IOK,IPAR)

   Input: A,N,NDIM
   Output: A,ORDER,SCALEROW,IOK,IPAR
   
   This subroutine computes the L and U triangular matrices equivalent to the
   A matrix, such that LU = A.  These matrices are returned in the space of A,
   in compact form. See C.F. Gerald and P. O. Wheatley, in Applied Numerical
   Analysis, 4th edition, pag. 106.
   
   REAL A(NDIM,NDIM) -> matrix of coefficients
   INTEGER N -> logical dimension of A
   INTEGER NDIM -> physical dimension of A in the calling program
   INTEGER ORDER(N) -> vector holding row order after pivoting
   REAL SCALEROW(N) -> vector holding scaling factors applied to each row
   INTEGER IOK -> returns 0 if everything works properly, +(the row number)
                  if all elements in a row are zero, or -(the row number) if
                  the pivot value is zero.
   INTEGER IPAR -> returns as +1 or -1 depending on whether the number of row
                   interchanges was even or odd, respectively
   
``lusolv``
..........

::
   
   SUBROUTINE LUSOLV(A,N,NDIM,ORDER,SCALEROW,B,X)

   Input: A,N,NDIM,ORDER,SCALEROW,B
   Output: X
   
   This subroutine solves the set of N linear equations A X = B, where the
   A matrix corresponds to the LU decomposition of the initial coefficient
   matrix. See C.F. Gerald and P. O. Wheatley, in Applied Numerical
   Analysis, 4th edition, pag. 110. The matrix A remains unchanged (also ORDER
   and SCALEROW), so subsequent calls to this subroutine, variying the B matrix,
   can be performed.
   
   REAL A(NDIM,NDIM) -> matrix of coefficients (LU in compact scheme)
   INTEGER N -> logical dimension of A
   INTEGER NDIM -> physical dimension of A in the calling program
   INTEGER ORDER(N) -> vector holding row order after pivoting in LUDCMP
   REAL SCALEROW(N) -> vector holding scaling factors applied to each row
   REAL B(N) -> right-hand side vector B
   REAL X(N) -> solution vector X
   
``mideind``
...........

::
   
   INTEGER FUNCTION MIDEIND(NS1,NS2,ITI,WV,FWV,CERR,RCVEL1,NCRES,FINDEX,EINDEX,
                            EJJGG,ESIMU,SN,FFPLAW,LONLY_SN)

   Input: NS1,NS2,ITI,WV,CERR,RCVEL1,NCRES,FFPLAW,LONLY_SN
   Output: FINDEX,EINDEX,EJJGG,ESIMU,SN
   Input/Output: (see COMMON blocks)
   
   Return MIDEIND=0 if atomic/molecular/D4000/B4000/generic index (and error)
   has been properly measured. MIDEIND=1 otherwise.
   
   INTEGER     NS1,NS2 -> scans to be coadded
   INTEGER     ITI -> index type:
                      ITI=   1: molecular index
                      ITI=   2: atomic index
                      ITI=   3: D4000 (Bruzual 1983)
                      ITI=   4: B4000 (own defintion)
                      ITI =  5: like B4000 but computing flux per angstrom
                      ITI=????: generic index:
                                ITI= C x 100 + L, where
                                     C: number of continuum regions
                                     L: number of absorption regions
                       Since Cmin=1, Cmax=99, Lmin=1, Lmax=99
                             ==> ITImin=101, ITImax=9999
   REAL        WV(NWVMAX) -> wavelength limits
   REAL        FWV(NWVMAX/4) -> constant factors to be used when computing
                             the index as multiplicative coefficients for
                             the absorption signal.
   CHARACTER*1 CERR -> 'y' if error spectrum is available, 'n' otherwise
   REAL        RCVEL1 -> 1 + v/c
   INTEGER     NCRES -> No. of response curve to be employed (1=averaged)
   REAL        FINDEX -> measured index
   REAL        EINDEX -> measured index error using own formulae
   REAL        EJJGG -> measured index error using JJGGs formulae
   REAL        ESIMU -> measured index error using numerical simulations
   REAL        SN -> averaged (S/N ratio)/Angstrom in the index region
   REAL        FFPLAW -> =0,...,10 indicates the fraction of light (in the
                      selected photometric band) to be used with a power law
                      of the form: F(lambda) = k * (lambda)**(alpha-2.0),
                      where alpha is defined through the variable ALPHAPLAW
                      (global variable ---COMMON---). If FFPLAW = -1 no
                      power law is employed.
   LOGICAL     LONLY_SN -> if .TRUE., the subroutine only computes the
                           mean S/N ratio and returns
   
``ordena1f1i``
..............

::
   
   SUBROUTINE ORDENA1F1I(N,A,B)

   Input: N,A,B
   Output: A,B
   
   Sorts the real array A(N) into ascending numerical order. The additional
   array B(N) is simultaneously changed in parallel with the array
   A(N). Note that the two input arrays are returned rearranged. We follow the
   Heapsort method, as described by D. Knuth, The Art of Computer Programming
   (pag.146, 5.2.3.).
   
   INTEGER N -> input number of data in A
   REAL    A(N) -> data matrix to be sorted
   INTEGER B(N) -> data matrix to be sorted in parallel with matrix A
   
``ordena1f2i``
..............

::
   
   SUBROUTINE ORDENA1F2I(N,A,B,C)

   Input: N,A,B,C
   Output: A,B,C
   
   Sorts the real array A(N) into ascending numerical order. The additional
   arrays B(N) and C(N) are simultaneously changed in parallel with the array
   A(N). Note that the three input arrays are returned rearranged. We follow the
   Heapsort method, as described by D. Knuth, The Art of Computer Programming
   (pag.146, 5.2.3.).
   
   INTEGER N -> input number of data in A
   REAL    A(N) -> data matrix to be sorted
   INTEGER B(N) -> data matrix to be sorted in parallel with matrix A
   INTEGER C(N) -> data matrix to be sorted in parallel with matrix A
   
``ordena1f``
............

::
   
   SUBROUTINE ORDENA1F(N,A)

   Input: N,A
   Output: A
   
   Sorts the real array A(N) into ascending numerical order. Note that the input
   array is returned rearranged. We follow the Heapsort method, as described by
   D. Knuth, The Art of Computer Programming (pag.146, 5.2.3.).
   
   INTEGER N -> input number of data in A
   REAL    A(N) -> data matrix to be sorted
   
``ordena1i``
............

::
   
   SUBROUTINE ORDENA1I(N,A)

   Input: N,A
   Output: A
   
   Sorts the integer array A(N) into ascending numerical order. Note that the
   input array is returned rearranged. We follow the Heapsort method, as
   described by C D. Knuth, The Art of Computer Programming (pag.146, 5.2.3.).
   
   INTEGER N -> input number of data in A
   INTEGER A(N) -> data matrix to be sorted
   
``ordena2f``
............

::
   
   SUBROUTINE ORDENA2F(N,A,B)

   Input: N,A,B
   Output: A,B
   
   Sorts the real array A(N) into ascending numerical order. The additional
   array B(N) is simultaneously changed in parallel with the array A(N).
   Note that both input arrays are returned rearranged. We follow the Heapsort
   method, as described by D. Knuth, The Art of Computer Programming (pag.146,
   5.2.3.).
   
   INTEGER N -> input number of data in A
   REAL    A(N) -> data matrix to be sorted
   REAL    B(N) -> data matrix to be sorted in parallel with matrix A
   
``ordena2i``
............

::
   
   SUBROUTINE ORDENA2I(N,A,B)

   Input: N,A,B
   Output: A,B
   
   Sorts the integer array A(N) into ascending numerical order. The additional
   array B(N) is simultaneously changed in parallel with the array A(N).
   Note that both input arrays are returned rearranged. We follow the Heapsort
   method, as described by D. Knuth, The Art of Computer Programming (pag.146,
   5.2.3.).
   
   INTEGER N -> input number of data in A
   INTEGER A(N) -> data matrix to be sorted
   INTEGER B(N) -> data matrix to be sorted in parallel with matrix A
   
``polfit``
..........

::
   
   SUBROUTINE  POLFIT(X,Y,SIGMAY,NPTS,NTERMS,MODE,A,CHISQR)

   >>> This subroutine is based on the subroutine from Data Reduction and
   Error Analysis for the Physical Sciences (Bevington, 1969) <<<
   
         LEAST-SQUARES FIT TO A POLYNOMIAL
         INPUT: X  -  ARRAY FOR INDEPENDENT VARIABLE
                Y  -  ARRAY FOR DEPENDENT VARIABLE
                SIGMAY  -  STANDARD DEVIATIONS FOR Y DATA POINTS
                NPTS  -  NUMBER OF PAIRS OF DATA POINTS
                NTERMS  - NUMBER OF COEFFICIENTS (DEGREE + 1)
                MODE  -  METHOD OF WEIGHTING (0 = NO WEIGHTING)
         OUTPUT:A  - ARRAY OF COEFFICIENTS
                CHISQR  -  REDUCED CHI SQUARE FOR FIT
   
         IT USES FUNCTION DETERM TO EVALUATE DETERMINANT OF MATRIX
         SUPPORTS NTERM UP TO 20
         FOR DETAILS SEE BEVINGTON(1969)
   
``pseudofit``
.............

::
   
   SUBROUTINE PSEUDOFIT(XF,YF,NF,NTERMS,YRMSTOL,WEIGHT,POWER,LUP,A)

   Input: XF,YF,NF,NTERMS,YRMSTOL,WEIGHT,POWER,LUP
   Output: A
   
   Calculate the polynomial fit to the upper/lower side of a set of data
   points.
   
   REAL XF(NF),YF(NF) -> data points to be fitted
   INTEGER NF -> number of data points
   INTEGER NTERMS -> number of coeffcients
   REAL YRMSTOL -> stopping criterion for DOWNHILL
   REAL WEIGHT -> weighting factor to enhance one side of the fit
   REAL POWER -> power to be used to compute distances
   LOGICAL LUP -> .TRUE.: fit upper side
                  .FALSE.: fit lower side
   REAL A(NTERMS) -> fitted coefficients
   
``ranred``
..........

::
   
   REAL FUNCTION RANRED(NSEED)

   Input: NSEED
   Output: RANRED (function)
   
   Return a random number in the range [0,1) using the intrinsic fortran
   function RAND(). If NSEED<0 a previous call to SRAND(TIME()) is also
   performed.
   
   INTEGER NSEED -> NSEED=0: RANRED returns the next random number in the
                             sequence.
                    NSEED<0: RANRED performs a previous call to the
                             intrinsic fortran function SRAND(TIME()), and
                             generates a random number in the new sequence.
                             In this case NSEED returns 0.
                    NSEED>0: RANRED performs a previous call to the
                             intrinsic fortran function SRAND(NSEED), and
                             generates a random number in the new sequence.
                             In this case NSEED returns 0.
   
``rebining``
............

::
   
   SUBROUTINE REBINING(X,Y,N,XX,YY,M,XINI,XINC)

   Input: X,Y,N,XINI,XINC
   Output: XX,YY,M
   
   Rebin a data table X(1:N),Y(1:N). If X(J)-X(J-1) .GT. XINC the routine
   performs a linear interpolation.
   
   REAL    X(N) -> ordered input array (not necesarilly equally-spaced)
   REAL    Y(N) -> input array
   INTEGER N -> no. of points in input array
   REAL    XX(M) -> equally-spaced output array
   REAL    YY(M) -> output array
   INTEGER M -> no. of points in output array
   REAL    XINI -> equal to XX(1)
   REAL    XINC -> equal to XX(J)-XX(J-1) for all J
   
``rinterp``
...........

::
   
   REAL FUNCTION RINTERP(LAMBDA)

   Input: LAMBDA
   Output: RINTERP (function)
   
   Calculate the extinction value A(lambda)/E(B-V) for a given wavelength
   LAMBDA. The tabulated data correspond to Savage & Mathis (1979, Ann.
   Rev. Astron. Astrophys., 17, 13).
   
   REAL LAMBDA -> input wavelength
   
``rvrebin``
...........

::

   
   SUBROUTINE RVREBIN(RADVEL,NCHAN,S,SS,STWV,DISP)

   Input: RADVEL,NCHAN,S,STWV,DISP
   Output: SS
   
   This subroutine applies a radial velocity shift to a spectrum.
   
   REAL RADVEL    -> radial velocity to be applied
   INTEGER NCHAN  -> number of channels
   REAL S(NCHAN)  -> initial spectrum
   REAL SS(NCHAN) -> shifted spectrum
   REAL STWV      -> central wavelength of the first pixel
   REAL DISP      -> dispersion (Angs./pixel) in the wavelength direction
   
``cubspl__``
............

::
   
   SUBROUTINE CUBSPL__(X,Y,N,IMODE,S,A,B,C)

   INTEGER N
   REAL X(N),Y(N)
   INTEGER IMODE
   REAL S(N)
   REAL A(N),B(N),C(N)
  
``cubsplx__``
.............

::

   SUBROUTINE CUBSPLX__(X,Y,A,B,C,N,I0,X0,Y0)
   
   Input: X,Y,A,B,C,N,I0,X0
   Output: Y0
   
   The subroutine returns the cubic spline evaluated at X0, using the
   coefficients determined in a previous call to CUBSPL__. The spline defined in
   the interval between X(I),Y(I) and X(I+1),Y(I+1) is given by:
   
        Y = A(I)*(X-X(I))**3 + B(I)*(X-X(I))**2 + C(I)*(X-X(I)) + D(I)
   
   If X0.LT.X(1), I=1 is employed (first computed spline)
   If X0.GT.X(N), I=N-1 is employed (last computed spline)
   
   REAL X(N) -> X-values fitted with CUBSPL__
   REAL Y(N) -> Y-values fitted with CUBSPL__
   REAL A(N) -> spline coefficients
   REAL B(N) -> spline coefficients
   REAL C(N) -> spline coefficients
   INTEGER N -> number of data points
   INTEGER I0 -> initial location to start the search of the place of X0 in
                 the X array
   REAL X0 -> X-value where the spline function will be evaluated
   REAL Y0 -> spline value at X0
   
``selbands``
............

::
   
   SUBROUTINE SELBANDS(CBAND,NPBAND,WV,RES)

   Input: CBAND
   Output: NPBAND,WV,RES
   
   Return the response curves of some common photometric bands.
   
   CHARACTER*1 CPBAND -> photometric band: U,B,V
   INTEGER     NPBAND -> no. of points which define the output table
   REAL        WV(NPBANDMAX) -> wavelengths
   REAL        RES(NPBANDMAX) -> response curve
   
``selindex``
............

::
   
   SUBROUTINE SELINDEX(NINDEX,WV,FWV,ITI,CLABEL)

   Input: NINDEX
   Output: WV,FWV,ITI,CLABEL
   
   Return the bandpass limits of atomic, and molecular indices (and the D4000).
   The subroutine looks first for a file called 'myindex.dat' in the current
   directory. If this file does not exist, the program then looks for a file
   called 'index.dat' (located in the subdirectory 'files' of the distribution
   package). If this last file is also missing, the program stops.
   
   INTEGER     NINDEX -> index number. If NINDEX=0 the routine returns ITI
               with the total number of defined indices.
   REAL        WV(NWVMAX) -> wavelength limits.
   REAL        FWV(NWVMAX/4) -> constant factors to be applied to the data in
                                the absorption bands.
   INTEGER     ITI -> index type:
               ITI =  -??: slope
               ITI =   1 : molecular
                   =   2 : atomic
                   =   3 : D4000
                   =   4 : B4000
                   =   5 : color
                   = ????: generic with
                           ITI= C x 100 + L, where C=No. continuum regions
                                                   L=No. absorption regions
                           Cmin=1, Cmax=99, Lmin=1, Lmax=99
                           (ITImin=101, ITImax=9999)
   CHARACTER*8 CLABEL -> character string with index identification
   
``sellines``
............

::
   
   SUBROUTINE SELLINES(NTYPE,NLINES,WV,CLABEL)

   Input: NTYPE
   Output: NLINES,WV,CLABEL
   
   Return the wavelength location of typical lines.
   
   INTEGER     NTYPE -> type of lines:
               NTYPE=0: Balmer serie
               NTYPE=1: typical emission lines
               NTYPE=2: typical sky lines
               NTYPE=3: typical absorption lines
   INTEGER     NLINES -> no. of returned lines
   REAL        WV(NLINMAX) -> wavelengths
   CHARACTER*8 CLABEL(NLINMAX) -> character string with line identification
   
``shindex``
...........

::
   
   SUBROUTINE SHINDEX(LINDOK,MODE)

   Input: LINDOK,MODE
   
   Show a list of available indices (with an 80 character width format) in the
   subroutine SELINDEX.
   
   LOGICAL LINDOK(NINDMAX) -> if .TRUE. the index is shown in the list
   INTEGER MODE -> if MODE=0 two extra options are displayed, namely
                   -1:EXIT and 0:ALL.
   
``splfit``
..........

::
   
   SUBROUTINE SPLFIT(N,X,Y,ND,XD,YRMSTOL,NOUT,XOUT,YOUT,XMIN,XMAX,SIGMA,LPLOTS)

   Input: N,X,Y,ND,XD,YRMSTOL,NOUT,XMIN,XMAX,SIGMA,LPLOTS
   Output: XOUT,YOUT
   
   Least-squares fit to splines, using ND knots located at XD().
   Input data are X(N), Y(N). XOUT(NOUT), YOUT(NOUT) are the output values which
   are computed in the range from XMIN to XMAX. The knot location determines the
   X(),Y() range employed in the fit (which is performed in the interval from
   XD(1) to XD(ND))
   
   INTEGER N -> initial number of points in input data
   REAL    X(N) -> sorted input data
   REAL    Y(N) -> input data
   INTEGER ND -> number of knots
   REAL    XD(ND) -> X location of each knot
   REAL    YRMSTOL ->  stopping criterion for DOWNHILL
   INTEGER NOUT -> number of points in output
   REAL    XOUT(NOUT) -> output data
   REAL    YOUT(NOUT) -> output data
   REAL    XMIN -> = XOUT(1)
   REAL    XMAX -> = XOUT(NOUT)
   REAL    SIGMA -> sigma of the fit
   LOGICAL LPLOTS -> if .TRUE. some plots are performed
   
``subprece``
............

::
   
   SUBROUTINE SUBPRECE(TII,RAI,DECI,TFF,RAF,DECF)

   Input: TII,RAI,DECI,TFF
   Output: RAF,DECF
   
   Transformation of coordinates given for an equinox to another equinox
   (precession effect).
   
   DOUBLE PRECISION TII -> initial equinox (year)
   DOUBLE PRECISION RAI -> initial right ascension (hours)
   DOUBLE PRECISION DECI -> initial declination (degrees)
   DOUBLE PRECISION TIF -> final equinox (year)
   DOUBLE PRECISION RAF -> final right ascension (hours)
   DOUBLE PRECISION DECF -> final declination (degrees)
   
``ulogreb``
...........

::
   
   SUBROUTINE ULOGREB(COPC,S,N,CRVAL,CRPIX,CDELT,SS,M,STWV,DISP)

   Input: COPC,S,N,CRVAL,CRPIX,CDELT,STWV,DISP
   Output: SS,M
   
   Transform a spectrum S(N) in logarithmic wavelength scale into a linear
   wavelength scale.
   
   CHARACTER*1 COPC -> type of wavelength calibration of input spectrum
                       COPC='1': CTYPE='WAVE'     (linear)
                       COPC='2': CTYPE='WAVE-LOG' (log10)
                       COPC='3': CTYPE='WAVE-LOG' (ln)
                       COPC='4': CTYPE='wavenumber'
   REAL    S(N) -> input spectrum
   INTEGER N -> no. of points in input spectrum
   REAL    CRVAL,CRPIX,CDELT -> wavelength calibration of input spectrum
                                if COPC='4', CRVAL and CRPIX are wavenumbers
   REAL    SS(M) -> output spectrum
   INTEGER M -> no. of points in output spectrum
   REAL    STWV, DISP -> linear wavelength calibration for output spectrum

Library ``libbutton.a``
-----------------------
   
``buttmis``
...........

::

   BUTTMIS

   This file contains subroutines and functions which are similar to other
   already defined for REDUCEME, namely TRUELEN, TRUEBEG, READI, READC, READF,
   LRUNX and RMBLANK. We duplicate them in order to create a library which must
   work independtly of REDUCEME.
   
   INTEGER FUNCTION TLENBUTT(CADENA)
   INTEGER FUNCTION TBEGBUTT(CADENA)
   INTEGER FUNCTION READIBUTT(CDEF)
   CHARACTER*(*) FUNCTION READCBUTT(CDEF,CVAL)
   REAL FUNCTION READFBUTT(CDEF)
   SUBROUTINE LRUNXBUTT(LRUN,LMANUAL,LHTML)
   SUBROUTINE RMBLANKBUTT(C1,C2,L)
   
``button``
..........

::
   
   SUBROUTINE BUTTON(N,TEXT,MODE)

   Input: N,TEXT,MODE
   Input (COMMON): global variables in button.inc
   
   Plot buttons and button text in different modes.
   
   INTEGER       N -> button number in the range of available buttons (which
                 runs from 1 to MAX_XBUTT x MAX_YBUTT)
   CHARACTER*(*) TEXT -> the text that will appear in the button
   INTEGER       MODE -> determine the button mode:
                 MODE=-2,-3,...: only text is plotted with PGPLOT color=-NMODE-1
                 (i.e. 1,2,3...)
                 MODE=-1 erase the button
                 MODE=0 whole button is plotted (text in black)
                 MODE=1 only text is plotted (white)
                 MODE=2 only text is plotted (black)
                 MODE=3 only text is plotted (gray, button disabled)
                 MODE=4 whole button with reversed colors (text in black)
                 MODE=5 whole button with reversed colors (text in white)
   
``buttqbr``
...........

::
   
   SUBROUTINE BUTTQBR(X1,X2,Y1,Y2)

   Output: X1,X2,Y1,Y2
   
   Return the button region limits.
   
   REAL X1 -> x-coordinate of the left hand edge of the button region viewport,
        in normalized device coordinates
   REAL X2 -> x-coordinate of the right hand edge of the button region viewport,
        in normalized device coordinates
   REAL Y1 -> y-coordinate of the bottom edge of the button region viewport,
        in normalized device coordinates
   REAL Y2 -> y-coordinate of the top edge of the button region viewport,
        in normalized device coordinates
   
``buttqcf``
...........

::
   
   SUBROUTINE BUTTQCF(FONT)

   Output: FONT
   
   Return the current character font type in buttons.
   
   INTEGER FONT -> the current font number (in range 1-4)
   
``buttqch``
...........

::
   
   SUBROUTINE BUTTQCH(SIZE)

   Output: SIZE
   
   Return the current character font size in buttons.
   
   REAL SIZE -> the current font size (dimensionless multiple of the default
        size)
   
``buttqex``
...........

::
   
   SUBROUTINE BUTTQEX(NBUT,LEXIST)

   Input: NBUT
   Output: LEXIST
   
   Return whether the asked button is active (currently available) or not.
   
   INTEGER NBUT -> button number
   LOGICAL LEXIST -> .TRUE. if the button is active, .FALSE. otherwise
   
``buttqit``
...........

::
   
   SUBROUTINE BUTTQIT(LOUTSIDE)

   Output: LOUTSIDE
   
   Return whether tick marks are drawn outside the viewport instead of inside.
   
   LOGICAL LOUTSIDE -> .TRUE. if ticks are drawn outside the viewport
   
``buttqpr``
...........

::
   
   SUBROUTINE BUTTQPR(X1,X2,Y1,Y2)

   Output: X1,X2,Y1,Y2
   
   Return the plot region limits.
   
   REAL X1 -> x-coordinate of the left hand edge of the plot region viewport,
        in normalized device coordinates
   REAL X2 -> x-coordinate of the right hand edge of the plot region viewport,
        in normalized device coordinates
   REAL Y1 -> y-coordinate of the bottom edge of the plot region viewport,
        in normalized device coordinates
   REAL Y2 -> y-coordinate of the top edge of the plot region viewport,
        in normalized device coordinates
   
``buttqxb``
...........

::
   
   SUBROUTINE BUTTQXB(NB)

   Output: NB
   
   Return MAX_XBUTT.
   
   INTEGER NB -> = MAX_XBUTT
   
``buttqyb``
...........

::
   
   SUBROUTINE BUTTQYB(NB)

   Output: NB
   
   Return MAX_YBUTT.
   
   INTEGER NB -> = MAX_YBUTT
   
``buttqytext``
..............

::
   
   SUBROUTINE BUTTQYTEXT(YTEXT)

   Output: YTEXT
   
   Return the current relative y-position of the text baseline in buttons
   (from 0 to 1)
   
   REAL YTEXT -> = YTEXT_BUTT
   
``buttsbr``
...........

::
   
   SUBROUTINE BUTTSBR(X1,X2,Y1,Y2)

   Input: X1,X2,Y1,Y2
   
   Set the button region limits.
   
   REAL X1 -> x-coordinate of the left hand edge of the button region viewport,
        in normalized device coordinates
   REAL X2 -> x-coordinate of the right hand edge of the button region viewport,
        in normalized device coordinates
   REAL Y1 -> y-coordinate of the bottom edge of the button region viewport,
        in normalized device coordinates
   REAL Y2 -> y-coordinate of the top edge of the button region viewport,
        in normalized device coordinates
   
``buttscf``
...........

::
   
   SUBROUTINE BUTTSCF(FONT)

   Input : FONT
   
   Set the character font type in buttons.
   
   INTEGER FONT -> the current font number (in range 1-4)
   
``buttsch``
...........

::
   
   SUBROUTINE BUTTSCH(SIZE)

   Input: SIZE
   
   Set the character height in buttons.
   
   REAL SIZE -> the current font size (dimensionless multiple of the default
        size)
   
``buttsex``
...........

::
   
   SUBROUTINE BUTTSEX(NBUT,LEXIST)

   Input: NBUT,LEXIST
   
   Set whether the asked button is active (currently available) or not.
   
   INTEGER NBUT -> button number
   LOGICAL LEXIST -> .TRUE. if the button is active, .FALSE. otherwise
   
``buttsit``
...........

::
   
   SUBROUTINE BUTTSIT(LOUTSIDE)

   Input: LOUTSIDE
   
   Set whether tick marks are drawn outside the viewport instead of inside.
   
   LOGICAL LOUTSIDE -> .TRUE. if ticks are drawn outside the viewport
   
``buttspr``
...........

::
   
   SUBROUTINE BUTTSPR(X1,X2,Y1,Y2)

   Input: X1,X2,Y1,Y2
   
   Set the plot region limits.
   
   REAL X1 -> x-coordinate of the left hand edge of the plot region viewport,
        in normalized device coordinates
   REAL X2 -> x-coordinate of the right hand edge of the plot region viewport,
        in normalized device coordinates
   REAL Y1 -> y-coordinate of the bottom edge of the plot region viewport,
        in normalized device coordinates
   REAL Y2 -> y-coordinate of the top edge of the plot region viewport,
        in normalized device coordinates
   
``buttsxb``
...........

::
   
   SUBROUTINE BUTTSXB(NB)

   Input: NB
   
   Set MAX_XBUTT.
   
   INTEGER NB -> = MAX_XBUTT
   
``buttsyb``
...........

::
   
   SUBROUTINE BUTTSYB(NB)

   Input: NB
   
   Set MAX_YBUTT.
   
   INTEGER NB -> = MAX_YBUTT
   
``buttsytext``
..............

::
   
   SUBROUTINE BUTTSYTEXT(YTEXT)

   Input: YTEXT
   
   Set the relative y-position of the text baseline in buttons (from 0 to 1)
   
   REAL YTEXT -> = YTEXT_BUTT
   
``ifbutton``
............

::
   
   SUBROUTINE IFBUTTON(XC,YC,NB)

   Input: XC,YC
   Output: NB
   
   Determine whether any button has been selected.
   
   REAL    XC -> world x-coordinate of the cursor
   REAL    YC -> world y-coordinate of the cursor
   INTEGER NB -> number of the selected button (if available). NB=0 if no
                 button has been selected.
   
``rpgband``
...........

::
   
   SUBROUTINE RPGBAND(MODE,POSN,XREF,YREF,XC,YC,CH)

   Input: MODE,POSN,XREF,YREF
   Output: XC,YC,CH
   
   This routine is similar to PGBAND, but it also allows the utilization of
   buttons in text mode.
   
   INTEGER     MODE -> display mode (see PGPLOT manual)
   INTEGER     POSN -> if POSN=1, the routine positions the cursor at the
               position specified by XREF,YREF
   REAL        XREF -> reference position
   REAL        YREF -> reference position
   REAL        XC -> the world x-coordinate of the cursor
   REAL        YC -> the world y-coordinate of the cursor
   CHARACTER*1 CH -> the character typed by the user
   
``rpgbegin``
............

::
   
   SUBROUTINE RPGBEGIN(NTERM,IDN,LCOLOR)

   Output: NTERM,IDN,LCOLOR
   Output (COMMON): all global variables in button.inc
   
   Open the graphic device(s) and assign the default values to the global
   variables:
         MAX_XBUTT=6
         MAX_YBUTT=2
         PGSCF_BUTT=2
         PGSCH_BUTT=1.
         YTEXT_BUTT=0.35
         X1VPORT=0.1
         X2VPORT=0.95
         Y1VPORT=0.1
         Y2VPORT=0.70
         X3VPORT=0.05
         X4VPORT=0.95
         Y3VPORT=0.80
         Y4VPORT=0.95
   
   INTEGER NTERM -> number of opened graphic devices to be employed
           simultaneously
   INTEGER IDN(8) -> identifier of the openned graphic devices
           (positive values returned by PGOPEN)
   LOGICAL LCOLOR(8) -> determines whether color is available or not
           in each opened graphic device
   
   
``rpgbegok``
............

::
   
   SUBROUTINE RPGBEGOK(TTERM,IMODE)

   Input: TTERM,IMODE
   Output (COMMON): all global variables in button.inc
   
   Open the graphic device TTERM and assign the default values to the global
   variables:
   MAX_XBUTT=6
   MAX_YBUTT=2
   PGSCF_BUTT=2
   PGSCH_BUTT=1.
   YTEXT_BUTT=0.35
   X1VPORT=0.1
   X2VPORT=0.95
   Y1VPORT=0.1
   Y2VPORT=0.70
   X3VPORT=0.05
   X4VPORT=0.95
   Y3VPORT=0.80
   Y4VPORT=0.95
   
   CHARACTER*(*) TTERM -> graphic device to be opened
   INTEGER       IMODE -> define the mode to operate with the buttons:
                          IMODE=0: graphics button
                          IMODE=1: text buttons (but plot graphic buttons)
                          IMODE=2: text buttons (without graphic buttons)
   
``rpgenv``
..........

::
   
   SUBROUTINE RPGENV(XMIN,XMAX,YMIN,YMAX,JUST,AXIS)

   Input: XMIN,XMAX,YMIN,YMAX,JUST,AXIS
   Input (COMMON): ITICKS_BUTT
   
   Perform the same functions than PGENV, although the plot surface is
   restricted to the rectangle defined by X1VPORT,X2VPORT,Y1VPORT,Y2VPORT.
   Other important difference with PGENV is that RPGENV does not clear the
   plot region of the new plot. A previous call to PGADVANCE, PGPAGE, PGERAS
   (RPGERAS, RPGERASB or RPGERASW) is required. The arguments of this routine
   are exactly the same than those in PGENV:
   
   REAL    XMIN -> the world x-coordinate at the bottom left corner of the
                   viewport
   REAL    XMAX -> the world x-coordinate at the top right corner of the
                   viewport
   REAL    YMIN -> the world y-coordinate at the bottom left corner of the
                   viewport
   REAL    YMAX -> the world y-coordinate at the top right corner of the
                   viewport
   INTEGER JUST -> if JUST=1, the scales of the x and y axes (in world
                   coordinates per inch) will be equal, otherwise they will be
                   scaled independently
   INTEGER AXIS -> controls the plotting of axes, tick marks, etc:
           AXIS = -2: draw no box, axes or labels
           AXIS = -1: draw box only
           AXIS =  0: draw box and label it with coordinates
           AXIS =  1: same as AXIS=0, but also draw the coordinate axes
           AXIS =  2: same as AXIS=1, but also draw grid lines
           AXIS = 10: draw box and label X-axis logarithmically
           AXIS = 20: draw box and label Y-axis logarithmically
           AXIS = 30: draw box and label both axes logarithmically
   
``rpgerasb``
............

::
   
   SUBROUTINE RPGERASB

   Input (COMMON) : X3VPORT,X4VPORT,Y3VPORT,Y4VPORT
   
   Clear the button region (preserving the plot region which does not overlap
   with the plot region).
   
``rpgeras``
...........

::
   
   SUBROUTINE RPGERAS

   Input (COMMON) : X1VPORT,X2VPORT,Y1VPORT,Y2VPORT
   
   Clear the plot region (preserving the button region which does not overlap
   with the plot region).
   
``rpgerasw``
............

::
   
   SUBROUTINE RPGERASW(X1,X2,Y1,Y2)

   Input: X1,X2,Y1,Y2
   
   Clear any rectangle defined by (X1,Y1) lower left corner
                                  (X2,Y2) upper right corner
   
   REAL X1 -> x-coordinate of the left hand edge of the rectangle to be
              cleared,in normalized device coordinates
   REAL X2 -> x-coordinate of the right hand edge of the rectangle to be
              cleared,in normalized device coordinates
   REAL Y1 -> y-coordinate of the bottom edge of the rectangle to be
              cleared,in normalized device coordinates
   REAL Y2 -> y-coordinate of the top edge of the rectangle to be
              cleared,in normalized device coordinates
   
   NOTE: this subroutine preserves the original viewport and window coordinate
         systems
