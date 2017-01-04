C------------------------------------------------------------------------------
C Version 26-November-1996                                      File: welcome.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This routine is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C SUBROUTINE WELCOME(CSTRING)
C
C Input: CSTRING
C Input (COMMON): THISPROGRAM,CREDUCEVERSION
C
C Write the welcome presentation of the programs. This subroutine also
C verifies whether the environment variable PGPLOT_DIR has been defined.
C
C CHARACTER*(*) CSTRING -> additional information to be shown as a centered
C               character string in the welcome presentation
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE WELCOME(CSTRING)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        INTEGER TRUELEN
C
        CHARACTER*(*) CSTRING
C local variables
        INTEGER I
        INTEGER LTAB
        INTEGER LP,LV,LR
        CHARACTER*1 CCONT
        CHARACTER*255 CENVIRONMENT
        CHARACTER*79 CBLANK
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C comprobamos que estan definidas las variables de entorno
        CALL GETENV('PGPLOT_DIR',CENVIRONMENT)
        IF(TRUELEN(CENVIRONMENT).EQ.0)THEN
          WRITE(*,101)'ERROR: PGPLOT_DIR is not defined.'
          WRITE(*,100)'Do you want to continue anyway (y/n) '
          CCONT(1:1)=READC('n','yn')
          IF(CCONT.EQ.'n') STOP
        END IF
C
!       CALL GETENV('reduceme_dir',CENVIRONMENT)
!       IF(TRUELEN(CENVIRONMENT).EQ.0)THEN
!         WRITE(*,101)'ERROR: reduceme_dir is not defined.'
!         WRITE(*,100)'Do you want to continue anyway (y/n) '
!         CCONT(1:1)=READC('n','yn')
!         IF(CCONT.EQ.'n') STOP
!       END IF
C------------------------------------------------------------------------------
        DO I=1,79
          CBLANK(I:I)=' '
        END DO
        LR=LEN(CREDUCEVERSION)
        LP=TRUELEN(THISPROGRAM)
        LV=TRUELEN(CSTRING)
        LTAB=(79-(11+LP))/2   !centramos la etiqueta con el nombre del programa
        WRITE(*,151)
        WRITE(*,100)CREDUCEVERSION
        WRITE(*,100)CBLANK(1:LTAB-LR)
        WRITE(*,100)'Welcome to '
        WRITE(*,100)THISPROGRAM(1:LP)
        WRITE(*,100)CBLANK(1:LTAB-LV-9)
        IF(2*LTAB+11+LP.EQ.78) WRITE(*,100)' '
        WRITE(*,100)'Version: '
        WRITE(*,101)CSTRING(1:LV)
        WRITE(*,152)
        WRITE(*,*)
C
100     FORMAT(A,$)
101     FORMAT(A)
151     FORMAT(79('*'))
152     FORMAT(79('-'))
        END
