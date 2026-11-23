!------------------------------------------------------------------------------
! Version 24-November-2026                                 file: isfitsfile.f90
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
! LOGICAL FUNCTION ISFITSFILE(FILENAME)
!
! Input: FILENAME
! Output: ISFITSFILE
!
! Returns .TRUE. if FILENAME must be read as a FITS file: its first bytes are
! 'SIMPLE  =' (FITS), or the gzip signature (compressed files, which CFITSIO
! reads directly), or the name ends with an HDU specification in brackets
! (CFITSIO extended file name syntax, e.g. file.fits[1]); in the last case
! the existence of the file is not checked here.
!
! CHARACTER*(*) FILENAME -> file name
!
!Comment
!------------------------------------------------------------------------------
        LOGICAL FUNCTION ISFITSFILE(FILENAME)
        IMPLICIT NONE
        CHARACTER*(*) FILENAME
!
        INTEGER L,LU,IOS
        CHARACTER*9 CSTART
!------------------------------------------------------------------------------
        ISFITSFILE=.FALSE.
        L=TRUELEN(FILENAME)
        IF(L.EQ.0) RETURN
! HDU specification in brackets: only meaningful for FITS files
        IF((FILENAME(L:L).EQ.']').AND.(INDEX(FILENAME(1:L),'[').GT.1))THEN
          ISFITSFILE=.TRUE.
          RETURN
        END IF
! first bytes of the file (stream access, so that no record structure is
! assumed)
        OPEN(NEWUNIT=LU,FILE=FILENAME(1:L),STATUS='OLD',ACCESS='STREAM',&
     &   FORM='UNFORMATTED',ACTION='READ',IOSTAT=IOS)
        IF(IOS.NE.0) RETURN
        READ(LU,IOSTAT=IOS) CSTART
        CLOSE(LU)
        IF(IOS.NE.0) RETURN
        IF(CSTART.EQ.'SIMPLE  =')THEN
          ISFITSFILE=.TRUE.
! gzip signature: bytes 0x1F 0x8B
        ELSEIF((ICHAR(CSTART(1:1)).EQ.31).AND.&
     &   (ICHAR(CSTART(2:2)).EQ.139))THEN
          ISFITSFILE=.TRUE.
        END IF
        END
