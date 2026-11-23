!------------------------------------------------------------------------------
! Version 24-November-2026                               file: fits2scratch.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This routine is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! SUBROUTINE FITS2SCRATCH(NF,FILENAME,LERR,LOK)
!
! Input: NF,FILENAME,LERR
! Output: LOK
!
! Reads a FITS image (the primary HDU, or the HDU given with the CFITSIO
! syntax FILENAME[n]) and writes it, header and data, in REDUCEME format into
! a scratch file connected to the logical unit NF, which is left rewound:
! INFILEX then reads it as any REDUCEME file (the scratch file is deleted
! when the unit is closed).
!
! Image: NAXIS=1 (spectrum: NSCAN=1), NAXIS=2, or NAXIS=3 with NAXIS3=1;
! NCHAN=NAXIS1, NSCAN=NAXIS2. Any BITPIX is accepted (CFITSIO applies BSCALE
! and BZERO and converts the data to REAL).
!
! REDUCEME header from the FITS keywords:
! - STWV=CRVAL1+(1-CRPIX1)*CDELT1, DISP=CDELT1 (CD1_1 if CDELT1 is missing;
!   CRPIX1=1 if missing); without CRVAL1 or CDELT1/CD1_1, STWV=1 and DISP=1
!   (pixels) with a warning; a logarithmic scale (CTYPE1 containing '-LOG')
!   is rejected;
! - AIRMASS; TIMEXPOS=EXPTIME (TIMEXPO, written by writefits, if EXPTIME is
!   missing); OBJECT; FITSFILE (the name of this file if the keyword is
!   missing); COMMENT=REDCOMM; missing numbers are set to 0;
! - ERRFRAME=T marks an error frame (' @ERROR@' is appended to the object
!   name in the scratch file, as in REDUCEME error frames).
!
! If this is the first FITS data file read (not an error frame: LERR=.FALSE.
! and ERRFRAME not true), its header cards selected by LKEEPCARD are kept to
! be copied into the FITS output files (SCRATCH2FITS).
!
! INTEGER       NF -> logical unit number for the scratch file
! CHARACTER*(*) FILENAME -> FITS file name (CFITSIO syntax allowed)
! LOGICAL       LERR -> .TRUE. if an error frame is expected
! LOGICAL       LOK -> .TRUE. if the file has been read; .FALSE. if it could
!               not be read (an error message has been printed)
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE FITS2SCRATCH(NF,FILENAME,LERR,LOK)
        IMPLICIT NONE
        INTEGER NF
        CHARACTER*(*) FILENAME
        LOGICAL LERR,LOK
!
! CFITSIO status for a keyword that does not exist
        INTEGER KEY_NO_EXIST
        PARAMETER (KEY_NO_EXIST=202)
!
        INTEGER IUNIT,ISTATUS,ISTATUS2,BLOCKSIZE,BITPIX,NAXIS,NAXES(3)
        INTEGER NSCAN0,NCHAN0,I,J,L,NKEYS,NMORE,NKEEP
        INTEGER FPIXELS(3),LPIXELS(3),INCS(3)
        DOUBLE PRECISION CRPIX1,CRVAL1,CDELT1
        REAL STWV0,DISP0,AIRMASS0,TIMEXPOS0
        REAL, ALLOCATABLE :: SCAN(:)
        CHARACTER*80 CTYPE1,CCOMMENT,CARD
        CHARACTER*80, ALLOCATABLE :: CARDS(:)
        CHARACTER*255 OBJECT0,FITSFILE0,COMMENT0
        LOGICAL LERRFRAME,ANYF,LCALIB,LDROPPED
!------------------------------------------------------------------------------
        LOK=.FALSE.
        L=TRUELEN(FILENAME)
        ISTATUS=0
        CALL FTGIOU(IUNIT,ISTATUS)
        CALL FTOPEN(IUNIT,FILENAME(1:L),0,BLOCKSIZE,ISTATUS)
        IF(ISTATUS.GT.0)THEN
          CALL FITSERR(ISTATUS,'cannot open the FITS file '//FILENAME(1:L))
          GOTO 900
        END IF
!------------------------------------------------------------------------------
! image dimensions
        CALL FTGIPR(IUNIT,3,BITPIX,NAXIS,NAXES,ISTATUS)
        IF(ISTATUS.GT.0)THEN
          CALL FITSERR(ISTATUS,'cannot read the image parameters of '//&
     &     FILENAME(1:L))
          GOTO 900
        END IF
        IF(NAXIS.EQ.1)THEN
          NCHAN0=NAXES(1)
          NSCAN0=1
        ELSEIF((NAXIS.EQ.2).OR.((NAXIS.EQ.3).AND.(NAXES(3).EQ.1)))THEN
          NCHAN0=NAXES(1)
          NSCAN0=NAXES(2)
        ELSEIF(NAXIS.EQ.0)THEN
          WRITE(*,101)'ERROR: this HDU of the FITS file contains no '//&
     &     'image data (an extension can be read with the syntax '//&
     &     'file.fits[n]).'
          GOTO 900
        ELSE
          WRITE(*,101)'ERROR: only FITS images with NAXIS=1 or 2 (or 3 '//&
     &     'with NAXIS3=1) can be read.'
          GOTO 900
        END IF
        IF((NCHAN0.LT.1).OR.(NSCAN0.LT.1))THEN
          WRITE(*,101)'ERROR: invalid image dimensions in the FITS file.'
          GOTO 900
        END IF
!------------------------------------------------------------------------------
! wavelength calibration
        CALL FTGKYS(IUNIT,'CTYPE1',CTYPE1,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          CTYPE1=' '
        END IF
        IF(INDEX(CTYPE1,'-LOG').GT.0)THEN
          WRITE(*,101)'ERROR: logarithmic wavelength scale (CTYPE1='''//&
     &     CTYPE1(1:TRUELEN(CTYPE1))//''') cannot be read by REDUCEME.'
          GOTO 900
        END IF
        LCALIB=.TRUE.
        CALL FTGKYD(IUNIT,'CRPIX1',CRPIX1,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          CRPIX1=1.D0
        END IF
        CALL FTGKYD(IUNIT,'CRVAL1',CRVAL1,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          LCALIB=.FALSE.
        END IF
        CALL FTGKYD(IUNIT,'CDELT1',CDELT1,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          CALL FTGKYD(IUNIT,'CD1_1',CDELT1,CCOMMENT,ISTATUS)
          IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
            ISTATUS=0
            LCALIB=.FALSE.
          END IF
        END IF
        IF(LCALIB)THEN
! CRVAL1 is the wavelength of pixel CRPIX1; STWV is that of pixel 1
          STWV0=REAL(CRVAL1+(1.D0-CRPIX1)*CDELT1)
          DISP0=REAL(CDELT1)
        ELSE
          WRITE(*,101)'WARNING: no wavelength calibration (CRVAL1, '//&
     &     'CDELT1) in the FITS header: STWV=1, DISP=1 (pixels).'
          STWV0=1.0
          DISP0=1.0
        END IF
!------------------------------------------------------------------------------
! other REDUCEME header values
        CALL FTGKYE(IUNIT,'AIRMASS',AIRMASS0,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          AIRMASS0=0.0
        END IF
        CALL FTGKYE(IUNIT,'EXPTIME',TIMEXPOS0,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          CALL FTGKYE(IUNIT,'TIMEXPO',TIMEXPOS0,CCOMMENT,ISTATUS)
          IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
            ISTATUS=0
            TIMEXPOS0=0.0
          END IF
        END IF
        CALL FTGKLS(IUNIT,'OBJECT',OBJECT0,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          OBJECT0=' '
        END IF
        CALL FTGKLS(IUNIT,'FITSFILE',FITSFILE0,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          FITSFILE0=FILENAME(1:L)
        END IF
        CALL FTGKLS(IUNIT,'REDCOMM',COMMENT0,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          COMMENT0=' '
        END IF
        CALL FTGKYL(IUNIT,'ERRFRAME',LERRFRAME,CCOMMENT,ISTATUS)
        IF(ISTATUS.EQ.KEY_NO_EXIST)THEN
          ISTATUS=0
          LERRFRAME=.FALSE.
        END IF
        IF(ISTATUS.GT.0)THEN
          CALL FITSERR(ISTATUS,'cannot read the header of '//FILENAME(1:L))
          GOTO 900
        END IF
!------------------------------------------------------------------------------
! header cards to be copied into the FITS output files (first FITS data file)
        IF((NCARDS_FITSHDR.LT.0).AND.(.NOT.LERR).AND.(.NOT.LERRFRAME))THEN
          CALL FTGHSP(IUNIT,NKEYS,NMORE,ISTATUS)
          ALLOCATE(CARDS(NKEYS))
          NKEEP=0
          LDROPPED=.FALSE.
          DO I=1,NKEYS
            CALL FTGREC(IUNIT,I,CARD,ISTATUS)
! a CONTINUE card belongs to the previous one
            IF(CARD(1:8).EQ.'CONTINUE')THEN
              IF(LDROPPED) CYCLE
            ELSE
              LDROPPED=.NOT.LKEEPCARD(CARD)
              IF(LDROPPED) CYCLE
            END IF
            NKEEP=NKEEP+1
            CARDS(NKEEP)=CARD
          END DO
          IF(ISTATUS.GT.0)THEN
            CALL FITSERR(ISTATUS,'cannot read the header of '//&
     &       FILENAME(1:L))
            GOTO 900
          END IF
          IF(ALLOCATED(CARDS_FITSHDR)) DEALLOCATE(CARDS_FITSHDR)
          ALLOCATE(CARDS_FITSHDR(MAX(NKEEP,1)))
          IF(NKEEP.GT.0) CARDS_FITSHDR(1:NKEEP)=CARDS(1:NKEEP)
          NCARDS_FITSHDR=NKEEP
        END IF
!------------------------------------------------------------------------------
! scratch file in REDUCEME format (see OUTFILEX)
        OPEN(NF,STATUS='SCRATCH',FORM='UNFORMATTED')
        WRITE(NF)'abcdefghijkl'
        WRITE(NF) NSCAN0,NCHAN0
        WRITE(NF) STWV0,DISP0
        WRITE(NF) AIRMASS0
        WRITE(NF) TIMEXPOS0
        J=TRUELEN(OBJECT0)
        IF(LERRFRAME)THEN
          IF(J.GT.0)THEN
            WRITE(NF) J+8
            WRITE(NF) OBJECT0(1:J)//' @ERROR@'
          ELSE
            WRITE(NF) 8
            WRITE(NF) ' @ERROR@'
          END IF
        ELSE
          WRITE(NF) J
          IF(J.GT.0) WRITE(NF) OBJECT0(1:J)
        END IF
        J=TRUELEN(FITSFILE0)
        WRITE(NF) J
        IF(J.GT.0) WRITE(NF) FITSFILE0(1:J)
        J=TRUELEN(COMMENT0)
        WRITE(NF) J
        IF(J.GT.0) WRITE(NF) COMMENT0(1:J)
! data, one scan per record: FITS pixel (j,i) is channel j of scan i; the
! image is read scan by scan through subsections (no 32-bit overflow of the
! linear pixel number in large images)
        ALLOCATE(SCAN(NCHAN0))
        FPIXELS=1
        LPIXELS=1
        INCS=1
        LPIXELS(1)=NCHAN0
        DO I=1,NSCAN0
          IF(NAXIS.GE.2)THEN
            FPIXELS(2)=I
            LPIXELS(2)=I
          END IF
! null value 0: no check of undefined pixels
          CALL FTGSVE(IUNIT,1,NAXIS,NAXES,FPIXELS,LPIXELS,INCS,0.0,SCAN,&
     &     ANYF,ISTATUS)
          IF(ISTATUS.GT.0)THEN
            CALL FITSERR(ISTATUS,'cannot read the data of '//FILENAME(1:L))
            CLOSE(NF)
            GOTO 900
          END IF
          WRITE(NF) (SCAN(J),J=1,NCHAN0)
        END DO
        REWIND(NF)
        LOK=.TRUE.
!------------------------------------------------------------------------------
! close the FITS file (also after an error: the status must be reset,
! since CFITSIO routines do nothing when it is not zero)
900     ISTATUS2=0
        CALL FTCLOS(IUNIT,ISTATUS2)
        ISTATUS2=0
        CALL FTFIOU(IUNIT,ISTATUS2)
101     FORMAT(A)
        END
