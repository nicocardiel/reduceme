!------------------------------------------------------------------------------
! Version 20-November-2026                                      File: rinterp.f
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
! SUBROUTINE CHEQUEA_FILEINDEX
!
! This subroutine verifies that there is an appropiate file containing
! the index definitions.
!
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE CHEQUEA_FILEINDEX
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
        INTEGER TRUEBEG,TRUELEN
!
        INTEGER I,II
        INTEGER L1,L2,N
        INTEGER ITI
        INTEGER NCONTI,NABSOR
        REAL WV(NWVMAX),FWV(NWVMAX/4)
        CHARACTER*8 CLABEL
        CHARACTER*60 CLDOS
        CHARACTER*255 REDUCEMEDIR
        LOGICAL LOGFILE1,LOGFILE2
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        INQUIRE(FILE='myindex.dat',EXIST=LOGFILE1)
        IF(LOGFILE1)THEN
          WRITE(*,101)'* Index definitions from file:'
          WRITE(*,101)'> ./myindex.dat'
          OPEN(20,FILE='./myindex.dat',STATUS='OLD',FORM='FORMATTED')
        ELSE
          CALL GET_ENVIRONMENT_VARIABLE('reduceme_dir',REDUCEMEDIR)
          L1=TRUEBEG(REDUCEMEDIR)
          L2=TRUELEN(REDUCEMEDIR)
          INQUIRE(FILE=REDUCEMEDIR(L1:L2)//'/files/index.dat',&
     &     EXIST=LOGFILE2)
          IF(LOGFILE2)THEN
            WRITE(*,101)'* Index definitions from file:'
            WRITE(*,100,ADVANCE='NO')'> '
            WRITE(*,100,ADVANCE='NO')REDUCEMEDIR(L1:L2)
            WRITE(*,101)'/files/index.dat'
            OPEN(20,FILE=REDUCEMEDIR(L1:L2)//'/files/index.dat',&
     &       STATUS='OLD',FORM='FORMATTED')
          ELSE
            WRITE(*,101)'FATAL ERROR: the following files do not exist:'
            WRITE(*,101)'-> myindex.dat'
            WRITE(*,100,ADVANCE='NO')'-> '
            WRITE(*,100,ADVANCE='NO')REDUCEMEDIR(L1:L2)
            WRITE(*,101)'/files/index.dat'
            STOP
          END IF
        END IF
!------------------------------------------------------------------------------
! saltamos las dos primeras lineas
        N=0
        READ(20,*,END=900)
        READ(20,*,END=900)
! leemos los indices
10      READ(20,133,END=12,ERR=900)CLABEL,ITI,CLDOS
        IF(ITI.LT.-1)THEN                               !es un indice pendiente
          NCONTI=-ITI
          DO I=1,NCONTI
            II=(2*I-1)
            READ(20,*,ERR=900) WV(II),WV(II+1)
          END DO
        ELSEIF(ITI.LE.100)THEN                        !no es un indice generico
          READ(CLDOS,'(6(1X,F9.3))') WV(1),WV(2),WV(3),WV(4),WV(5),WV(6)
        ELSEIF(ITI.GT.100)THEN                           !es un indice generico
          NCONTI=(ITI/100)                      !numero de regiones de continuo
          NABSOR=ITI-NCONTI*100             !numero de regiones con absorciones
          DO I=1,NCONTI
            II=(2*I-1)
            READ(20,*,ERR=900) WV(II),WV(II+1)
          END DO
          DO I=1,NABSOR
            II=(2*NCONTI)+(2*I-1)
            READ(20,*,ERR=900) WV(II),WV(II+1),FWV(I)
          END DO
        END IF
        N=N+1
        IF(MOD(N-1,8).NE.0)WRITE(*,100,ADVANCE='NO')', '
        WRITE(*,100,ADVANCE='NO')CLABEL
        IF(MOD(N,8).EQ.0) WRITE(*,*)
        GOTO 10
12      CLOSE(20)
        IF(MOD(N,8).NE.0) WRITE(*,*)
        WRITE(*,*)
133     FORMAT(A8,1X,I4,1X,A60)
!------------------------------------------------------------------------------
        IF(N.GT.NINDMAX)THEN
! en caso de fallar aqui, redimensionar NINDMAX en el fichero 'redlib.inc'
          WRITE(*,101)'FATAL ERROR: No. of defined indices > NINDMAX'
          STOP
        END IF
        IF(N.EQ.0)THEN
          WRITE(*,101)'FATAL ERROR: No. of defined indices = 0!'
          STOP
        END IF
        RETURN
!------------------------------------------------------------------------------
900     CONTINUE
        WRITE(*,101)'FATAL ERROR: while reading index data file.'
        STOP
!
100     FORMAT(A)
101     FORMAT(A)
        END
