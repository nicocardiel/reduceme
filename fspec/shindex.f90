!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: shindex.f90
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
! SUBROUTINE SHINDEX(LINDOK,MODE)
!
! Input: LINDOK,MODE
!
! Show a list of available indices (with an 80 character width format) in the
! subroutine SELINDEX. 
!
! LOGICAL LINDOK(NINDMAX) -> if .TRUE. the index is shown in the list
! INTEGER MODE -> if MODE=0 two extra options are displayed, namely
!                 -1:EXIT and 0:ALL.
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE SHINDEX(LINDOK,MODE)
        IMPLICIT NONE
!
        INCLUDE 'redlib.inc'
        LOGICAL LINDOK(NINDMAX)
        INTEGER MODE
!
        INTEGER I,N,L
        INTEGER K,NINDEXT,ITI
        REAL WV(NWVMAX),FWV(NWVMAX/4)
        CHARACTER*8 CLABEL
        CHARACTER*78 LINEA
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        DO I=1,LEN(LINEA)
          LINEA(I:I)=' '
        END DO
        WRITE(*,101)' *******************************************'//&
     &   '***********************************'
        CALL SELINDEX(0,WV,FWV,NINDEXT,CLABEL)!buscamos numero total de indices
!------------------>'12345678901231234567890123       !13 caracteres por indice
!------------------>'##:12345678  ##:12345678  '        !formato de cada indice
        IF(MODE.EQ.0)THEN
          LINEA(1:26)='-1:EXIT       0:ALL       '
          L=26                                !ultimo caracter en linea escrito
          N=2                                 !numero de indices en misma linea
        ELSE
          L=0                                 !ultimo caracter en linea escrito
          N=0                                 !numero de indices en misma linea
        END IF
        DO K=1,NINDEXT
          CALL SELINDEX(K,WV,FWV,ITI,CLABEL)
          IF(LINDOK(K))THEN
            N=N+1
            WRITE(LINEA(L+1:L+13),'(I2,A1,A8,2X)')K,':',CLABEL
            IF(K.EQ.NINDEXT)THEN
              WRITE(*,101)' '//LINEA
            ELSE
              IF(N.EQ.6)THEN
                WRITE(*,101)' '//LINEA
                DO I=1,LEN(LINEA)
                  LINEA(I:I)=' '
                END DO
                L=0
                N=0
              ELSE
                L=L+13
              END IF
            END IF
          ELSE
            IF(K.EQ.NINDEXT)THEN
              WRITE(*,101)' '//LINEA
            END IF
          END IF
        END DO
        WRITE(*,101)' *******************************************'//&
     &   '***********************************'
!
101     FORMAT(A)
        END
