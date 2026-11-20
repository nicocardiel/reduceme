!------------------------------------------------------------------------------
! Version 20-November-2026                                     file: ecugal.f90
!------------------------------------------------------------------------------
! Copyright N. Cardiel & J. Gorgas, Departamento de Astrofisica
! Universidad Complutense de Madrid, 28040-Madrid, Spain
! E-mail: ncl@astrax.fis.ucm.es or fjg@astrax.fis.ucm.es
!------------------------------------------------------------------------------
! This program is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the Free
! Software Foundation; either version 2 of the License, or (at your option) any
! later version. See the file gnu-public-license.txt for details.
!------------------------------------------------------------------------------
!Comment
!
! Program: ecugal
! Classification: miscellany
! Description: Transforms r.a. and dec. to galactic coordinates.
!
!Comment
!
        PROGRAM ECUGAL
        USE RED_M
        USE FUTILS_M
        USE REDLIB_M
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
!
        REAL PI
        PARAMETER(PI=3.141593)
        REAL AR0,AR
        REAL DEC0,DEC
        REAL L,B
        REAL X,Y,Z
        REAL XX,YY,ZZ
        CHARACTER*1 CMORE
!------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
!
        THISPROGRAM='ecugal'
        CALL WELCOME('20-November-2026')
!
10      WRITE(*,100,ADVANCE='NO')'A.R. 1950.0 (HH.MMSS)'
        AR0=READF('@')
        CALL DECIMAL(AR0,AR)
        AR=AR*15*PI/180.
        WRITE(*,100,ADVANCE='NO')'DEC. 1950.0 (DD.MMSS)'
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
        WRITE(*,100,ADVANCE='NO')'Continue (y/n) '
        CMORE(1:1)=READC('y','yn')
        IF(CMORE.EQ.'y') GOTO 10
        STOP
100     FORMAT(A)
200     FORMAT(A,F10.5)
        END
! **********************************************************************
!                                                     SUBROUTINE DECIMAL
!                                                     ******************
      SUBROUTINE DECIMAL(A,B)
!
! Paso de DD.MMSS a DD.dddd
! Utilizamos variables INTEGER*4 para evitar los errores de redondeo.
! Es posible hacer una llamada a la subrutina en la que los argumentos
! verdaderos correspondientes a las argumentos ficticios A y B
! coincidan.
!
      IMPLICIT NONE
!---> argumentos ficticios: ENTRADA
      REAL A
!---> argumentos ficticios: SALIDA
      REAL B
!---> variables locales
      INTEGER AA,A1,A2,A3
      CHARACTER*1 SIGNO
!
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
