C------------------------------------------------------------------------------
C Version 6-December-1996                                       file: ecugal.f
C------------------------------------------------------------------------------
C Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
C Universidad Complutense de Madrid, 28040-Madrid, Spain
C E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
C------------------------------------------------------------------------------
C This program is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the Free
C Software Foundation; either version 2 of the License, or (at your option) any
C later version. See the file gnu-public-license.txt for details.
C------------------------------------------------------------------------------
Comment
C
C Program: ecugal
C Classification: miscellany
C Description: Transforms r.a. and dec. to galactic coordinates.
C
Comment
C
        PROGRAM ECUGAL
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'futils.inc'
        REAL READF
C
        REAL PI
        PARAMETER(PI=3.141593)
        REAL AR0,AR
        REAL DEC0,DEC
        REAL L,B
        REAL X,Y,Z
        REAL XX,YY,ZZ
        CHARACTER*1 CMORE
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='ecugal'
        CALL WELCOME('6-December-1996')
C
10      WRITE(*,100)'A.R. 1950.0 (HH.MMSS)'
        AR0=READF('@')
        CALL DECIMAL(AR0,AR)
        AR=AR*15*PI/180.
        WRITE(*,100)'DEC. 1950.0 (DD.MMSS)'
        DEC0=READF('@')
        CALL DECIMAL(DEC0,DEC)
        DEC=DEC*PI/180.
        X=COS(DEC)*COS(AR)
        Y=COS(DEC)*SIN(AR)
        Z=SIN(DEC)
        XX=-.06699*X-0.87276*Y-0.48354*Z
        YY=0.49273*X-0.45035*Y+0.74458*Z
        ZZ=-0.86760*X-0.18837*Y+0.46020*Z
        L=ATAN2(YY,XX)
        L=L*180./PI
        IF(L.LT.0.)L=L+360.
        B=ASIN(ZZ)
        B=B*180/PI
        WRITE(*,200)'l = ',L
        WRITE(*,200)'b = ',B
        WRITE(*,*)
        WRITE(*,100)'Continue (y/n) '
        CMORE(1:1)=READC('y','yn')
        IF(CMORE.EQ.'y') GOTO 10
        STOP
100     FORMAT(A,$)
200     FORMAT(A,F10.5)
        END
C **********************************************************************
C                                                     SUBROUTINE DECIMAL
C                                                     ******************
      SUBROUTINE DECIMAL(A,B)
C
C Paso de DD.MMSS a DD.dddd
C Utilizamos variables INTEGER*4 para evitar los errores de redondeo.
C Es posible hacer una llamada a la subrutina en la que los argumentos
C verdaderos correspondientes a las argumentos ficticios A y B
C coincidan.
C
      IMPLICIT NONE
C---> argumentos ficticios: ENTRADA
      REAL A
C---> argumentos ficticios: SALIDA
      REAL B
C---> variables locales
      INTEGER*4 AA,A1,A2,A3
      CHARACTER*1 SIGNO
C
      SIGNO='+'
      IF(A.LT.0)SIGNO='-'
      AA=ABS(NINT(A*10000))
      A1=AA/10000
      A2=AA-A1*10000
      A2=A2/100
      A3=AA-A1*10000-A2*100
      B=REAL(A1)+REAL(A2)/60.+REAL(A3)/3600.
      IF(SIGNO.EQ.'-')B=-B
      END
