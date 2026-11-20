!------------------------------------------------------------------------------
! Version 20-November-2026                                     File: showperc.f
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
! SUBROUTINE SHOWPERC(N1,N2,ISTEP,I,NEXTINFO)
!
! Input: N1,N2,ISTEP,I,NEXTINFO
! Output: NEXTINFO
!
! Display the percentage of work performed in a loop, which has been defined
! in the range from N1 to N2, with an incremental step ISTEP, being I the 
! current value of the loop control.
!
! INTEGER N1 -> first limit of the control variable of the loop
! INTEGER N2 -> last limit of the control variable of the loop
! INTEGER ISTEP -> the value by which the control variable is incremented
! INTEGER I -> current value of the control variable
! INTEGER NEXTINFO -> integer which stores the last fraction of 10 displayed;
!                     this variable is initialized in the first call of this
!                     routine
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE SHOWPERC(N1,N2,ISTEP,I,NEXTINFO)
        IMPLICIT NONE
        INTEGER N1,N2,ISTEP,I,NEXTINFO
!
        REAL FN1,FN2,FISTEP,FNEXTINFO
!------------------------------------------------------------------------------
        IF(I.LT.N1)THEN
          RETURN
        ELSEIF(I.EQ.N1)THEN
          WRITE(*,'(A)',ADVANCE='NO')'Working...'
          NEXTINFO=1
        END IF
!
        IF(NEXTINFO.GT.10) RETURN
!
        FN1=REAL(N1)
        FN2=REAL(N2)
        FISTEP=REAL(ISTEP)
        FNEXTINFO=REAL(NEXTINFO)
!
        IF((FN2-FN1)/FISTEP.GT.10.)THEN
          IF(I.GE.NINT(FN2/10.*FNEXTINFO))THEN
            IF(NEXTINFO.EQ.10)THEN
              WRITE(*,'(A8)')'100% OK!'
              NEXTINFO=NEXTINFO+1
            ELSE
              WRITE(*,'(I2,A3)',ADVANCE='NO')10*NEXTINFO,'%, '
              NEXTINFO=NEXTINFO+1
            END IF
          END IF
        ELSE
          IF(I.EQ.N2)THEN
            WRITE(*,'(I2,A4)') I,' OK!'
          ELSE
            WRITE(*,'(I2,A2)',ADVANCE='NO') I,', '
          END IF
        END IF
!
        END
