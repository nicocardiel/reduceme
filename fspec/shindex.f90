C------------------------------------------------------------------------------
C Version 6-April-1999                                          File: shindex.f
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
C SUBROUTINE SHINDEX(LINDOK,MODE)
C
C Input: LINDOK,MODE
C
C Show a list of available indices (with an 80 character width format) in the
C subroutine SELINDEX. 
C
C LOGICAL LINDOK(NINDMAX) -> if .TRUE. the index is shown in the list
C INTEGER MODE -> if MODE=0 two extra options are displayed, namely
C                 -1:EXIT and 0:ALL.
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE SHINDEX(LINDOK,MODE)
        IMPLICIT NONE
C
        INCLUDE 'redlib.inc'
        LOGICAL LINDOK(NINDMAX)
        INTEGER MODE
C
        INTEGER I,N,L
        INTEGER K,NINDEXT,ITI
        REAL WV(NWVMAX),FWV(NWVMAX/4)
        CHARACTER*8 CLABEL
        CHARACTER*78 LINEA
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        DO I=1,LEN(LINEA)
          LINEA(I:I)=' '
        END DO
        WRITE(*,101)' *******************************************'//
     +   '***********************************'
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
        WRITE(*,101)' *******************************************'//
     +   '***********************************'
C
101     FORMAT(A)
        END
