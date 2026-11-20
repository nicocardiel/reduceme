C------------------------------------------------------------------------------
C Version 06-June-2003                                           File:ulogreb.f
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
C SUBROUTINE ULOGREB(COPC,S,N,CRVAL,CRPIX,CDELT,SS,M,STWV,DISP)
C
C Input: COPC,S,N,CRVAL,CRPIX,CDELT,STWV,DISP
C Output: SS,M
C
C Transform a spectrum S(N) in logarithmic wavelength scale into a linear
C wavelength scale.
C
C CHARACTER*1 COPC -> type of wavelength calibration of input spectrum
C                     COPC='1': CTYPE='WAVE'     (linear)
C                     COPC='2': CTYPE='WAVE-LOG' (log10)
C                     COPC='3': CTYPE='WAVE-LOG' (ln)
C                     COPC='4': CTYPE='wavenumber'
C REAL    S(N) -> input spectrum
C INTEGER N -> no. of points in input spectrum
C REAL    CRVAL,CRPIX,CDELT -> wavelength calibration of input spectrum
C                              if COPC='4', CRVAL and CRPIX are wavenumbers
C REAL    SS(M) -> output spectrum
C INTEGER M -> no. of points in output spectrum
C REAL    STWV, DISP -> linear wavelength calibration for output spectrum
C
Comment
C------------------------------------------------------------------------------
        SUBROUTINE ULOGREB(COPC,S,N,CRVAL,CRPIX,CDELT,SS,M,STWV,DISP)
        IMPLICIT NONE
        CHARACTER*1 COPC
        INTEGER N
        REAL S(N)
        REAL CRVAL,CRPIX,CDELT
        INTEGER M
        REAL SS(M)
        REAL STWV,DISP
C
        INTEGER K
        INTEGER J,J1,J2
        !OJO: hay que trabajar en doble precisión porque de lo contario
        !aparecen errores de redondeo al usar los logaritmos
        DOUBLE PRECISION W,W1,W2,W_CM
        DOUBLE PRECISION FJ1,FJ2
        DOUBLE PRECISION FFLUJO
C------------------------------------------------------------------------------
        DO K=1,M            !recorremos todos los pixels del espectro de salida
          W=DBLE(STWV)+DBLE(K-1)*DBLE(DISP)                 !l.d.o. del pixel K
          W1=W-0.5D0*DBLE(DISP)                 !l.d.o. borde izdo. del pixel K
          W2=W1+DBLE(DISP)                      !l.d.o. borde dcho. del pixel K
          IF(COPC.EQ.'1')THEN
            FJ1=(W1-DBLE(CRVAL))/DBLE(CDELT)+DBLE(CRPIX)   !pixel (número real)
            FJ2=(W2-DBLE(CRVAL))/DBLE(CDELT)+DBLE(CRPIX)   !pixel (número real)
            FFLUJO=1.0D0
          ELSEIF(COPC.EQ.'2')THEN
            FJ1=(DLOG10(W1)-DBLE(CRVAL))/DBLE(CDELT)+DBLE(CRPIX)  !pixel (real)
            FJ2=(DLOG10(W2)-DBLE(CRVAL))/DBLE(CDELT)+DBLE(CRPIX)  !pixel (real)
            FFLUJO=1.0D0
          ELSEIF(COPC.EQ.'3')THEN
            FJ1=(DLOG(W1)-DBLE(CRVAL))/DBLE(CDELT)+DBLE(CRPIX)    !pixel (real)
            FJ2=(DLOG(W2)-DBLE(CRVAL))/DBLE(CDELT)+DBLE(CRPIX)    !pixel (real)
            FFLUJO=1.0D0
          ELSEIF(COPC.EQ.'4')THEN
            FJ2=(1./W1-DBLE(CRVAL))/DBLE(CDELT)+DBLE(CRPIX)       !pixel (real)
            FJ1=(1./W2-DBLE(CRVAL))/DBLE(CDELT)+DBLE(CRPIX)       !pixel (real)
            W_CM=W/1.0D8 !l.d.o. en cm
            FFLUJO=1.0D0/(W_CM*W_CM)            !cambio de unidades en el flujo
          ELSE
            WRITE(*,101) 'FATAL ERROR: invalid COPC='//COPC
            STOP
          END IF
          J1=NINT(FJ1)                                   !pixel (número entero)
          J2=NINT(FJ2)                                   !pixel (número entero)
          SS(K)=0.0
          IF(J1.EQ.J2)THEN
            !el nuevo pixel está contenido por completo en un pixel del
            !espectro de entrada (eso sí, comprobamos que el número de pixel
            !está dentro del intervalo con datos)
            IF((J1.GE.1).AND.(J1.LE.N))THEN
              SS(K)=SS(K)+S(J1)*FFLUJO*REAL(FJ2-FJ1)
            END IF
          ELSE
            !en este caso el nuevo pixel se extiende, al menos, sobre dos
            !pixels del espectro de entrada; en esta situación vamos a calcular
            !la contribución de los pixels J1 y J2, y luego, si es necesario,
            !añadimos los pixels intermedios completos (si existen)
            IF((J1.GE.1).AND.(J1.LE.N))THEN
              SS(K)=SS(K)+S(J1)*FFLUJO*REAL(DBLE(J1)+0.5D0-FJ1)
            END IF
            IF((J2.GE.1).AND.(J2.LE.N))THEN
              SS(K)=SS(K)+S(J2)*FFLUJO*REAL(FJ2-(DBLE(J2)-0.5D0))
            END IF
            IF(J2-J1.GT.1)THEN
              DO J=J1+1,J2-1
                IF((J.GE.1).AND.(J.LE.N))THEN
                  SS(K)=SS(K)+S(J)*FFLUJO
                END IF
              END DO
            END IF
          END IF
        END DO
C------------------------------------------------------------------------------
101     FORMAT(A)
        END
