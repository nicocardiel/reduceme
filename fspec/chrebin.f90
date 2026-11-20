!------------------------------------------------------------------------------
! Version 20-November-2026                                    File: chrebin.f90
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
! SUBROUTINE CHREBIN(CHANSHIFT,NCHAN,S,SS)
!
! Input: CHANSHIFT,NCHAN,S
! Output: SS
!
! This subroutine applies a constant channel shift to a spectrum.
!
! REAL CHANSHIFT -> channel shift to be applied
! INTEGER NCHAN  -> number of channels
! REAL S(NCHAN)  -> initial spectrum
! REAL SS(NCHAN) -> shifted spectrum
!
!Comment
!------------------------------------------------------------------------------
        SUBROUTINE CHREBIN(CHANSHIFT,NCHAN,S,SS)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        INTEGER J
        INTEGER NC1,NC2
        REAL S(NCMAX),SS(NCMAX)
        REAL CHANSHIFT
        REAL FNC1,FNC2
        REAL CA,CB,CC,X1,X2,AREA
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        IF(NCHAN.LT.3)THEN
          WRITE(*,100,ADVANCE='NO')&
     &       'FATAL ERROR: NCHAN.LT.3 (subroutine CHREBIN)'
          STOP
        END IF
!------------------------------------------------------------------------------
! caso trivial: CHANSHIFT=0
        IF(CHANSHIFT.EQ.0.0)THEN
          DO J=1,NCHAN
            SS(J)=S(J)
          END DO
          RETURN
        END IF
!------------------------------------------------------------------------------
        DO J=1,NCHAN                       !inicializamos el espectro de salida
          SS(J)=0.
        END DO
!------------------------------------------------------------------------------
! Comenzamos el bucle principal en numero de canales
        DO J=1,NCHAN           !recorremos todos los canales del nuevo espectro
! posiciones en espectro inicial que corresponden a al nuevo canal
          FNC1=REAL(J)-CHANSHIFT-1.0
          FNC2=REAL(J)-CHANSHIFT
! primer y ultimo canal a utilizar (numeros enteros)
          NC1=INT(FNC1)
          IF(FNC1.GE.0.0) NC1=NC1+1 !aqui usamos .GE.
          NC2=INT(FNC2)
          IF(FNC2.GT.0.0) NC2=NC2+1 !aqui usamos .GT.
          IF(NC2.NE.NC1+1)THEN
            WRITE(*,100,ADVANCE='NO')'NC1,NC2: '
            WRITE(*,*)NC1,NC2
            WRITE(*,100,ADVANCE='NO')'FATAL ERROR: NC2.NE.NC1+1'
            WRITE(*,101)' (subroutine CHREBIN)'
            STOP
          END IF
!---------
          IF(NC1.EQ.1)THEN
            CA=(S(NC1)+S(NC1+2))/2.-S(NC1+1)
            CB=(S(NC1+2)+3*S(NC1))/(-2.)+2.*S(NC1+1)
            CC=S(NC1)-CA/12.
            X1=FNC1-(AINT(FNC1)+0.5)
            X2=+0.5
            AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+&
     &           CB*(X2*X2-X1*X1)/2.+&
     &           CC*(X2-X1)
            SS(J)=AREA
          ELSEIF((NC1.GT.1).AND.(NC1.LT.NCHAN))THEN
            CA=(S(NC1-1)+S(NC1+1))/2.-S(NC1)
            CB=(S(NC1+1)-S(NC1-1))/2.
            CC=S(NC1)-CA/12.
            X1=FNC1-(AINT(FNC1)+0.5)
            X2=+0.5
            AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+&
     &           CB*(X2*X2-X1*X1)/2.+&
     &           CC*(X2-X1)
            SS(J)=AREA
          ELSEIF(NC1.EQ.NCHAN)THEN
            CA=(S(NC1-2)+S(NC1))/2.-S(NC1-1)
            CB=(3.*S(NC1)+S(NC1-2))/2.-2.*S(NC1-1)
            CC=S(NC1)-CA/12.
            X1=FNC1-(AINT(FNC1)+0.5)
            X2=+0.5
            AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+&
     &           CB*(X2*X2-X1*X1)/2.+&
     &           CC*(X2-X1)
            SS(J)=AREA
          END IF
!---------
          IF(NC2.EQ.1)THEN
            CA=(S(NC2)+S(NC2+2))/2.-S(NC2+1)
            CB=(S(NC2+2)+3*S(NC2))/(-2.)+2.*S(NC2+1)
            CC=S(NC2)-CA/12.
            X1=-0.5
            X2=FNC2-(AINT(FNC2)+0.5)
            AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+&
     &           CB*(X2*X2-X1*X1)/2.+&
     &           CC*(X2-X1)
            SS(J)=SS(J)+AREA
          ELSEIF((NC2.GT.1).AND.(NC2.LT.NCHAN))THEN
            CA=(S(NC2-1)+S(NC2+1))/2.-S(NC2)
            CB=(S(NC2+1)-S(NC2-1))/2.
            CC=S(NC2)-CA/12.
            X1=-0.5
            X2=FNC2-(AINT(FNC2)+0.5)
            AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+&
     &           CB*(X2*X2-X1*X1)/2.+&
     &           CC*(X2-X1)
            SS(J)=SS(J)+AREA
          ELSEIF(NC2.EQ.NCHAN)THEN
            CA=(S(NC2-2)+S(NC2))/2.-S(NC2-1)
            CB=(3.*S(NC2)+S(NC2-2))/2.-2.*S(NC2-1)
            CC=S(NC2)-CA/12.
            X1=-0.5
            X2=FNC2-(AINT(FNC2)+0.5)
            AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+&
     &           CB*(X2*X2-X1*X1)/2.+&
     &           CC*(X2-X1)
            SS(J)=SS(J)+AREA
          END IF
!---------
        END DO                            !final del bucle en numero de canales
!------------------------------------------------------------------------------
100     FORMAT(A)
101     FORMAT(A)
        END
