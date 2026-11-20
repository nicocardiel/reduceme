!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: welcome.f90
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
! SUBROUTINE WELCOME(CSTRING)
!
! Input: CSTRING
! Input (COMMON): THISPROGRAM,CREDUCEVERSION
!
! Write the welcome presentation of the programs. This subroutine also
! verifies whether the environment variable PGPLOT_DIR has been defined.
!
! CHARACTER*(*) CSTRING -> additional information to be shown as a centered
!               character string in the welcome presentation
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE WELCOME(CSTRING)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
!
        CHARACTER*(*) CSTRING
! local variables
        INTEGER I
        INTEGER LTAB
        INTEGER LP,LV,LR
        CHARACTER*1 CCONT
        CHARACTER*255 CENVIRONMENT
        CHARACTER*79 CBLANK
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
! comprobamos que estan definidas las variables de entorno
        CALL GET_ENVIRONMENT_VARIABLE('PGPLOT_DIR',CENVIRONMENT)
        IF(TRUELEN(CENVIRONMENT).EQ.0)THEN
          WRITE(*,101)'ERROR: PGPLOT_DIR is not defined.'
          WRITE(*,100,ADVANCE='NO')&
     &       'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') STOP
        END IF
!
!       CALL GETENV('reduceme_dir',CENVIRONMENT)
!       IF(TRUELEN(CENVIRONMENT).EQ.0)THEN
!         WRITE(*,101)'ERROR: reduceme_dir is not defined.'
!         WRITE(*,100)'Do you want to continue anyway (y/n) '
!         CCONT(1:1)=READC('n','yn')
!         IF(CCONT.EQ.'n') STOP
!       END IF
!------------------------------------------------------------------------------
        DO I=1,79
          CBLANK(I:I)=' '
        END DO
        LR=LEN(CREDUCEVERSION)
        LP=TRUELEN(THISPROGRAM)
        LV=TRUELEN(CSTRING)
        LTAB=(79-(11+LP))/2   !centramos la etiqueta con el nombre del programa
        WRITE(*,151)
        WRITE(*,100,ADVANCE='NO')CREDUCEVERSION
        WRITE(*,100,ADVANCE='NO')CBLANK(1:LTAB-LR)
        WRITE(*,100,ADVANCE='NO')'Welcome to '
        WRITE(*,100,ADVANCE='NO')THISPROGRAM(1:LP)
        WRITE(*,100,ADVANCE='NO')CBLANK(1:LTAB-LV-9)
        IF(2*LTAB+11+LP.EQ.78) WRITE(*,100,ADVANCE='NO')' '
        WRITE(*,100,ADVANCE='NO')'Version: '
        WRITE(*,101)CSTRING(1:LV)
        WRITE(*,152)
        WRITE(*,*)
!
100     FORMAT(A)
101     FORMAT(A)
151     FORMAT(79('*'))
152     FORMAT(79('-'))
        END
