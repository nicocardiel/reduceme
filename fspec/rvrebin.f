C------------------------------------------------------------------------------
C Version 10-November-1997                                      File: rvrebin.f
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
C SUBROUTINE RVREBIN(RADVEL,NCHAN,S,SS,STWV,DISP)
C
C Input: RADVEL,NCHAN,S,STWV,DISP 
C Output: SS
C
C This subroutine applies a radial velocity shift to a spectrum.
C
C REAL RADVEL    -> radial velocity to be applied
C INTEGER NCHAN  -> number of channels
C REAL S(NCHAN)  -> initial spectrum
C REAL SS(NCHAN) -> shifted spectrum
C REAL STWV      -> central wavelength of the first pixel
C REAL DISP      -> dispersion (Angs./pixel) in the wavelength direction
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE RVREBIN(RADVEL,NCHAN,S,SS,STWV,DISP)
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
C
        REAL C
        PARAMETER (C=299792.46)
C
        INTEGER J,K
        INTEGER NC1,NC2
        REAL S(NCMAX),SS(NCMAX)
        REAL RADVEL
        REAL WLMIN,FACTOR
        REAL WL1,WL2
        REAL WL01,WL02
        REAL FNC1,FNC2
        REAL CA,CB,CC,X1,X2,AREA
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        IF(NCHAN.LT.3)THEN
          WRITE(*,101)'FATAL ERROR: NCHAN.LT.3 (subroutine RVREBIN)'
          STOP
        END IF
C------------------------------------------------------------------------------
C caso trivial: RADVEL=0
        IF(RADVEL.EQ.0.0)THEN
          DO J=1,NCHAN
            SS(J)=S(J)
          END DO
          RETURN
        END IF
C------------------------------------------------------------------------------
        DO J=1,NCHAN                       !inicializamos el espectro de salida
          SS(J)=0.
        END DO
C
        WLMIN=STWV-DISP/2.             !l.d.o. borde izquierdo del primer pixel
C------------------------------------------------------------------------------
        FACTOR=(1.+RADVEL/C)
        FACTOR=FACTOR/SQRT(1.-(RADVEL*RADVEL)/(C*C))    !correccion relativista
C------------------------------------------------------------------------------
C Comenzamos el bucle principal en numero de canales
        DO J=1,NCHAN           !recorremos todos los canales del nuevo espectro
C longitudes de onda del nuevo pixel (en el sistema del espectro de entrada,
C es decir, conservado STWV y DISP)
          WL1=WLMIN+REAL(J-1)*DISP
          WL2=WL1+DISP
C longitudes de onda del espectro inicial que tras variarlas en RADVEL dan
C lugar a WL1 y WL2
          WL01=WL1/FACTOR
          WL02=WL2/FACTOR
C posiciones en espectro inicial que corresponden a WL01 y WL02 (num. reales)
          FNC1=(WL01-WLMIN)/DISP
          FNC2=(WL02-WLMIN)/DISP
C primer y ultimo canal a utilizar (numeros enteros)
          NC1=INT(FNC1)
          IF(FNC1.GE.0.0) NC1=NC1+1 !aqui usamos .GE.
          NC2=INT(FNC2)
          IF(FNC2.GT.0.0) NC2=NC2+1 !aqui usamos .GT.
C sumamos todos los canales necesarios
c..............................................................................
          IF(NC1.EQ.NC2)THEN           !solo usa una fraccion de un mismo pixel
c-----------
            IF(NC1.EQ.1)THEN
              CA=(S(NC1)+S(NC1+2))/2.-S(NC1+1)
              CB=(S(NC1+2)+3*S(NC1))/(-2.)+2.*S(NC1+1)
              CC=S(NC1)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=FNC2-(AINT(FNC1)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            ELSEIF((NC1.GT.1).AND.(NC2.LT.NCHAN))THEN
              CA=(S(NC1-1)+S(NC1+1))/2.-S(NC1)
              CB=(S(NC1+1)-S(NC1-1))/2.
              CC=S(NC1)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=FNC2-(AINT(FNC1)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            ELSEIF(NC2.EQ.NCHAN)THEN
              CA=(S(NC2-2)+S(NC2))/2.-S(NC2-1)
              CB=(3.*S(NC2)+S(NC2-2))/2.-2.*S(NC2-1)
              CC=S(NC2)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=FNC2-(AINT(FNC1)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            END IF
c-----------
c..............................................................................
          ELSEIF(NC1+1.EQ.NC2)THEN                    !usa dos pixels contiguos
c-----------
            IF(NC1.EQ.1)THEN
              CA=(S(NC1)+S(NC1+2))/2.-S(NC1+1)
              CB=(S(NC1+2)+3*S(NC1))/(-2.)+2.*S(NC1+1)
              CC=S(NC1)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=+0.5
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            ELSEIF((NC1.GT.1).AND.(NC1.LT.NCHAN))THEN
              CA=(S(NC1-1)+S(NC1+1))/2.-S(NC1)
              CB=(S(NC1+1)-S(NC1-1))/2.
              CC=S(NC1)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=+0.5
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            ELSEIF(NC1.EQ.NCHAN)THEN
              CA=(S(NC1-2)+S(NC1))/2.-S(NC1-1)
              CB=(3.*S(NC1)+S(NC1-2))/2.-2.*S(NC1-1)
              CC=S(NC1)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=+0.5
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            END IF
c-----------
            IF(NC2.EQ.1)THEN
              CA=(S(NC2)+S(NC2+2))/2.-S(NC2+1)
              CB=(S(NC2+2)+3*S(NC2))/(-2.)+2.*S(NC2+1)
              CC=S(NC2)-CA/12.
              X1=-0.5
              X2=FNC2-(AINT(FNC2)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=SS(J)+AREA
            ELSEIF((NC2.GT.1).AND.(NC2.LT.NCHAN))THEN
              CA=(S(NC2-1)+S(NC2+1))/2.-S(NC2)
              CB=(S(NC2+1)-S(NC2-1))/2.
              CC=S(NC2)-CA/12.
              X1=-0.5
              X2=FNC2-(AINT(FNC2)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=SS(J)+AREA
            ELSEIF(NC2.EQ.NCHAN)THEN
              CA=(S(NC2-2)+S(NC2))/2.-S(NC2-1)
              CB=(3.*S(NC2)+S(NC2-2))/2.-2.*S(NC2-1)
              CC=S(NC2)-CA/12.
              X1=-0.5
              X2=FNC2-(AINT(FNC2)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=SS(J)+AREA
            END IF
c-----------
c..............................................................................
          ELSEIF(NC1+1.LT.NC2)THEN       !hay algun pixel completo entre medias
c-----------
            IF(NC1.EQ.1)THEN
              CA=(S(NC1)+S(NC1+2))/2.-S(NC1+1)
              CB=(S(NC1+2)+3*S(NC1))/(-2.)+2.*S(NC1+1)
              CC=S(NC1)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=+0.5
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            ELSEIF((NC1.GT.1).AND.(NC1.LT.NCHAN))THEN
              CA=(S(NC1-1)+S(NC1+1))/2.-S(NC1)
              CB=(S(NC1+1)-S(NC1-1))/2.
              CC=S(NC1)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=+0.5
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            ELSEIF(NC1.EQ.NCHAN)THEN
              CA=(S(NC1-2)+S(NC1))/2.-S(NC1-1)
              CB=(3.*S(NC1)+S(NC1-2))/2.-2.*S(NC1-1)
              CC=S(NC1)-CA/12.
              X1=FNC1-(AINT(FNC1)+0.5)
              X2=+0.5
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=AREA
            END IF
c-----------
            DO K=NC1+1,NC2-1
              IF((K.GE.1).AND.(K.LE.NCHAN))THEN
                SS(J)=SS(J)+S(K)      !sumamos los pixels intermedios completos
              END IF
            END DO
c-----------
            IF(NC2.EQ.1)THEN
              CA=(S(NC2)+S(NC2+2))/2.-S(NC2+1)
              CB=(S(NC2+2)+3*S(NC2))/(-2.)+2.*S(NC2+1)
              CC=S(NC2)-CA/12.
              X1=-0.5
              X2=FNC2-(AINT(FNC2)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=SS(J)+AREA
            ELSEIF((NC2.GT.1).AND.(NC2.LT.NCHAN))THEN
              CA=(S(NC2-1)+S(NC2+1))/2.-S(NC2)
              CB=(S(NC2+1)-S(NC2-1))/2.
              CC=S(NC2)-CA/12.
              X1=-0.5
              X2=FNC2-(AINT(FNC2)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=SS(J)+AREA
            ELSEIF(NC2.EQ.NCHAN)THEN
              CA=(S(NC2-2)+S(NC2))/2.-S(NC2-1)
              CB=(3.*S(NC2)+S(NC2-2))/2.-2.*S(NC2-1)
              CC=S(NC2)-CA/12.
              X1=-0.5
              X2=FNC2-(AINT(FNC2)+0.5)
              AREA=CA*(X2*X2*X2-X1*X1*X1)/3.+
     +             CB*(X2*X2-X1*X1)/2.+
     +             CC*(X2-X1)
              SS(J)=SS(J)+AREA
            END IF
c-----------
c..............................................................................
          ELSE                                 !si ocurre algo extranho paramos
            WRITE(*,101)'FATAL ERROR: in subroutine RVREBIN'
            STOP
          END IF
c..............................................................................
        END DO                            !final del bucle en numero de canales
C------------------------------------------------------------------------------
101     FORMAT(A,$)
        END
