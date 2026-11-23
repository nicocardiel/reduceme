!------------------------------------------------------------------------------
! Version 24-November-2026                                 file: closefilex.f90
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
! SUBROUTINE CLOSEFILEX(NF)
!
! Input: NF
!
! Closes the file connected to the logical unit NF. If the unit is a FITS
! output file opened by OUTFILEX (a scratch file in REDUCEME format), the
! FITS file is written first (SCRATCH2FITS). For any other unit it is
! equivalent to CLOSE(NF), so it must be used to close every file opened
! with OUTFILEX (otherwise a FITS output file is never written).
!
! INTEGER NF -> logical unit number
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE CLOSEFILEX(NF)
        IMPLICIT NONE
        INTEGER NF
!
        INTEGER I,K
!------------------------------------------------------------------------------
        K=0
        DO I=1,NFITSOUT
          IF(IUNIT_FITSOUT(I).EQ.NF) K=I
        END DO
        IF(K.EQ.0)THEN
          CLOSE(NF)
          RETURN
        END IF
        CALL SCRATCH2FITS(NF,NAME_FITSOUT(K))
! remove the unit from the list of pending FITS output files
        DO I=K,NFITSOUT-1
          IUNIT_FITSOUT(I)=IUNIT_FITSOUT(I+1)
          NAME_FITSOUT(I)=NAME_FITSOUT(I+1)
        END DO
        NFITSOUT=NFITSOUT-1
        END
