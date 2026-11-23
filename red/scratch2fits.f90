!------------------------------------------------------------------------------
! Version 24-November-2026                               file: scratch2fits.f90
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
! SUBROUTINE SCRATCH2FITS(NF,FILENAME)
!
! Input: NF,FILENAME
!
! Converts the REDUCEME file written by a program into the scratch file
! connected to the logical unit NF (opened by OUTFILEX) into the FITS file
! FILENAME, and closes (deletes) the scratch file.
!
! FITS file: BITPIX=-32; NAXIS1=NCHAN, NAXIS2=NSCAN (NAXIS=1 when NSCAN=1).
! Header:
! - the cards kept from the first FITS data file read (FITS2SCRATCH), if any;
! - CTYPE1='WAVE', CUNIT1='Angstrom', CRPIX1=1, CRVAL1=STWV, CDELT1=DISP;
! - OBJECT, EXPTIME (=TIMEXPOS), AIRMASS, ERRFRAME=T (error frames only),
!   FITSFILE and REDCOMM (=COMMENT) when they are not blank (long strings
!   are written with the CONTINUE convention);
! - DATE, and a HISTORY card with the REDUCEME version and the program.
!
! INTEGER       NF -> logical unit number of the scratch file
! CHARACTER*(*) FILENAME -> name of the FITS file to be created
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE SCRATCH2FITS(NF,FILENAME)
        USE REDLIB_M, ONLY: CREDUCEVERSION,THISPROGRAM
        IMPLICIT NONE
        INTEGER NF
        CHARACTER*(*) FILENAME
!
        INTEGER IUNIT,ISTATUS,ISTATUS2,NAXIS,NAXES(2)
        INTEGER NSCAN0,NCHAN0,I,J,L,LOBJ,LFITS,LCOMM,IOS
        INTEGER FPIXELS(2),LPIXELS(2)
        REAL STWV0,DISP0,AIRMASS0,TIMEXPOS0
        REAL, ALLOCATABLE :: SCAN(:)
        CHARACTER*12 CLAVE
        CHARACTER*255 OBJECT0,FITSFILE0,COMMENT0
        LOGICAL LERRFRAME
!------------------------------------------------------------------------------
        L=TRUELEN(FILENAME)
! REDUCEME header (see OUTFILEX)
        REWIND(NF)
        READ(NF) CLAVE
        READ(NF) NSCAN0,NCHAN0
        READ(NF) STWV0,DISP0
        READ(NF) AIRMASS0
        READ(NF) TIMEXPOS0
        OBJECT0=' '
        READ(NF) LOBJ
        IF(LOBJ.GT.0) READ(NF) OBJECT0(1:LOBJ)
        FITSFILE0=' '
        READ(NF) LFITS
        IF(LFITS.GT.0) READ(NF) FITSFILE0(1:LFITS)
        COMMENT0=' '
        READ(NF) LCOMM
        IF(LCOMM.GT.0) READ(NF) COMMENT0(1:LCOMM)
! error frames: ' @ERROR@' at the end of the object name
        LERRFRAME=.FALSE.
        IF(LOBJ.GE.8)THEN
          IF(OBJECT0(LOBJ-7:LOBJ).EQ.' @ERROR@')THEN
            LERRFRAME=.TRUE.
            OBJECT0(LOBJ-7:LOBJ)=' '
          END IF
        END IF
        LOBJ=TRUELEN(OBJECT0)
!------------------------------------------------------------------------------
        ISTATUS=0
        CALL FTGIOU(IUNIT,ISTATUS)
        CALL FTINIT(IUNIT,FILENAME(1:L),1,ISTATUS)
        IF(ISTATUS.GT.0)THEN
          CALL FITSERR(ISTATUS,'cannot create the FITS file '//&
     &     FILENAME(1:L))
          ISTATUS2=0
          CALL FTFIOU(IUNIT,ISTATUS2)
          CLOSE(NF)
          RETURN
        END IF
        NAXES(1)=NCHAN0
        NAXES(2)=NSCAN0
        IF(NSCAN0.EQ.1)THEN
          NAXIS=1
        ELSE
          NAXIS=2
        END IF
        CALL FTPHPS(IUNIT,-32,NAXIS,NAXES,ISTATUS)
! cards of the first FITS data file read
        DO I=1,NCARDS_FITSHDR
          CALL FTPREC(IUNIT,CARDS_FITSHDR(I),ISTATUS)
        END DO
! REDUCEME header; REAL values with 9 significant digits (exact round trip)
        CALL FTPKYS(IUNIT,'CTYPE1','WAVE','wavelength (linear scale)',&
     &   ISTATUS)
        CALL FTPKYS(IUNIT,'CUNIT1','Angstrom','wavelength unit',ISTATUS)
        CALL FTPKYE(IUNIT,'CRPIX1',1.0,-9,'reference pixel',ISTATUS)
        CALL FTPKYE(IUNIT,'CRVAL1',STWV0,-9,&
     &   'central wavelength of the reference pixel',ISTATUS)
        CALL FTPKYE(IUNIT,'CDELT1',DISP0,-9,'dispersion (Angstrom/pixel)',&
     &   ISTATUS)
        CALL FTPKLS(IUNIT,'OBJECT',OBJECT0(1:MAX(LOBJ,1)),'object name',&
     &   ISTATUS)
        CALL FTPKYE(IUNIT,'EXPTIME',TIMEXPOS0,-9,'exposure time (s)',&
     &   ISTATUS)
        CALL FTPKYE(IUNIT,'AIRMASS',AIRMASS0,-9,'airmass',ISTATUS)
        IF(LERRFRAME)THEN
          CALL FTPKYL(IUNIT,'ERRFRAME',.TRUE.,'REDUCEME error frame',&
     &     ISTATUS)
        END IF
        IF(LFITS.GT.0)THEN
          CALL FTPKLS(IUNIT,'FITSFILE',FITSFILE0(1:LFITS),&
     &     'original FITS file',ISTATUS)
        END IF
        IF(LCOMM.GT.0)THEN
          CALL FTPKLS(IUNIT,'REDCOMM',COMMENT0(1:LCOMM),&
     &     'REDUCEME comment',ISTATUS)
        END IF
! long string convention keyword, only if needed (68 characters per card)
        IF(MAX(LOBJ,LFITS,LCOMM).GT.68) CALL FTPLSW(IUNIT,ISTATUS)
        CALL FTPDAT(IUNIT,ISTATUS)
        CALL FTPHIS(IUNIT,CREDUCEVERSION(1:TRUELEN(CREDUCEVERSION))//&
     &   ': '//THISPROGRAM(1:TRUELEN(THISPROGRAM)),ISTATUS)
        IF(ISTATUS.GT.0)THEN
          CALL FITSERR(ISTATUS,'cannot write the header of '//FILENAME(1:L))
          GOTO 900
        END IF
!------------------------------------------------------------------------------
! data, one scan per record, written through subsections (no 32-bit
! overflow of the linear pixel number in large images)
        ALLOCATE(SCAN(NCHAN0))
        FPIXELS=1
        LPIXELS=1
        LPIXELS(1)=NCHAN0
        DO I=1,NSCAN0
          READ(NF,IOSTAT=IOS) (SCAN(J),J=1,NCHAN0)
          IF(IOS.NE.0)THEN
            WRITE(*,101)'ERROR: the program did not write all the scans '//&
     &       'of the FITS file '//FILENAME(1:L)
            GOTO 900
          END IF
          IF(NAXIS.EQ.2)THEN
            FPIXELS(2)=I
            LPIXELS(2)=I
          END IF
          CALL FTPSSE(IUNIT,1,NAXIS,NAXES,FPIXELS,LPIXELS,SCAN,ISTATUS)
        END DO
        CALL FTCLOS(IUNIT,ISTATUS)
        IF(ISTATUS.GT.0)THEN
          CALL FITSERR(ISTATUS,'cannot write the FITS file '//FILENAME(1:L))
        END IF
        ISTATUS=0
        CALL FTFIOU(IUNIT,ISTATUS)
        CLOSE(NF)
        RETURN
!------------------------------------------------------------------------------
! error: remove the incomplete FITS file
900     ISTATUS2=0
        CALL FTDELT(IUNIT,ISTATUS2)
        ISTATUS2=0
        CALL FTFIOU(IUNIT,ISTATUS2)
        CLOSE(NF)
101     FORMAT(A)
        END
