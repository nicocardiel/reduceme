!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: guessef.f90
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
! SUBROUTINE GUESSEF(INFILE,OUTFILE)
!
! Input: INFILE
! Ouput: OUTFILE
!
! This subroutine determines the expected error file name OUTFILE from INFILE.
! An additional "e" character is located between the portion of the file name
! (INFILE) preceding the last period (if present) and the last period itself.
!
! CHARACTER*(*) INFILE
! CHARACTER*(*) OUTFILE
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE GUESSEF(INFILE,OUTFILE)
        IMPLICIT NONE
        CHARACTER*(*) INFILE
        CHARACTER*(*) OUTFILE
!
        INTEGER LIN,LOUT,I
!------------------------------------------------------------------------------
        LIN=TRUELEN(INFILE)
        LOUT=LEN(OUTFILE)
        IF(LIN.GT.LOUT-1) LIN=LOUT-1
        IF(LIN.GT.255) LIN=254
!
        DO I=1,LOUT
          OUTFILE(I:I)=' '
        END DO
!
        IF(INDEX(INFILE(1:LIN),'.').EQ.0)THEN
          OUTFILE(1:LIN)=INFILE(1:LIN)
          OUTFILE(LIN+1:LIN+1)='e'
          RETURN
        ENDIF
!
        DO I=LIN,1,-1
          IF(INFILE(I:I).EQ.'.')GOTO 10
          IF(INFILE(I:I).EQ.'/')THEN     !si hay barra de directorio, escapamos
            OUTFILE(1:LIN)=INFILE(1:LIN)
            OUTFILE(LIN+1:LIN+1)='e'
            RETURN
          END IF
        END DO
10      OUTFILE(1:I-1)=INFILE(1:I-1)
        OUTFILE(I:I)='e'
        OUTFILE(I+1:LIN+1)=INFILE(I:LIN)
!
        END
