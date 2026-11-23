!------------------------------------------------------------------------------
! Version 24-November-2026                                 file: isfitsname.f90
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
! LOGICAL FUNCTION ISFITSNAME(FILENAME)
!
! Input: FILENAME
! Output: ISFITSNAME
!
! Returns .TRUE. if the output file FILENAME must be written in FITS format,
! i.e. if its name ends with .fits, .fit or .fts (in upper or lower case).
!
! CHARACTER*(*) FILENAME -> file name
!
!Comment
!------------------------------------------------------------------------------
        LOGICAL FUNCTION ISFITSNAME(FILENAME)
        IMPLICIT NONE
        CHARACTER*(*) FILENAME
!
        INTEGER L,I,IC
        CHARACTER*5 CEXT
!------------------------------------------------------------------------------
        ISFITSNAME=.FALSE.
        L=TRUELEN(FILENAME)
        IF(L.LT.5) RETURN
! last 5 characters in lower case (ASCII)
        CEXT=FILENAME(L-4:L)
        DO I=1,5
          IC=ICHAR(CEXT(I:I))
          IF((IC.GE.65).AND.(IC.LE.90)) CEXT(I:I)=CHAR(IC+32)
        END DO
        IF(CEXT.EQ.'.fits')THEN
          ISFITSNAME=.TRUE.
        ELSEIF((CEXT(2:5).EQ.'.fit').OR.(CEXT(2:5).EQ.'.fts'))THEN
          ISFITSNAME=.TRUE.
        END IF
        END
