!------------------------------------------------------------------------------
! Version 19-October-2007                                     File: fintgauss.f
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
! REAL FUNCTION FINTGAUSS(X1,X2,N,X0,FACTOR1)  -> integral of Gaussian
! REAL FUNCTION FINTGAUSSE(X1,X2,N,X0,FACTOR1) -> integral of Gaussian^2
!
! Input: X1,X2,N,X0,FACTOR1
! Output: FINTGAUSS (function) or FINTGAUSSE (function)
!
! Calculate the integral of a Gaussian of width SIGMA, center at wavelength X0,
! using Simpson's rule, between wavelengths X1 and X2, with N intervals 
! (N must be even).
!
! REAL X1,X2 -> integral limits (wavelengths)
! INTEGER N -> number of intervals to estimate integral
! REAL X0 -> Gaussian center (wavelength)
! REAL FACTOR1 -> = (c^2/(2*SIGMA^2))/(X0^2)
!
!Comment
!------------------------------------------------------------------------------
! Integral de una gaussiana de anchura SIGMA por el metodo de Simpson, entre 
! X1 y X2, con N intervalos (N debe ser par). Utilizamos la formula (1prim) del
! libro "Calculo Numerico Fundamental", Demidovich y Maron, pag. 656 (ojo, hay
! una errata en la formula de sigma2).
        REAL FUNCTION FINTGAUSS(X1,X2,N,X0,FACTOR1)
        IMPLICIT NONE
        REAL X1,X2
        INTEGER N
        REAL X0,FACTOR1
!
        DOUBLE PRECISION DX1,DX2
        DOUBLE PRECISION DX0,DFACTOR1
!
        INTEGER I
        DOUBLE PRECISION  X,Y1,Y2
        DOUBLE PRECISION H
        DOUBLE PRECISION SUM1,SUM2
!------------------------------------------------------------------------------
! chequeamos que N es par
        IF(MOD(N,2).NE.0)THEN
          STOP 'FATAL ERROR: N is odd in subroutine FINTGAUSS!'
        END IF
! pasamos a doble precision las variables REAL de entrada
        DX1=DBLE(X1)
        DX2=DBLE(X2)
        DX0=DBLE(X0)
        DFACTOR1=DBLE(FACTOR1)
!
        H=(DX2-DX1)/DBLE(N) !tamaño de cada intervalo
!
        SUM1=0.D0
        DO I=1,N-1,2
          X=DX1+H*DBLE(I)
          SUM1=SUM1+DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
!
        SUM2=0.D0          !sumamos terminos pares salvo el primero y el ultimo
        DO I=2,N-2,2
          X=DX1+H*DBLE(I)
          SUM2=SUM2+DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
!
        Y1=DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))       !funcion en el primer punto
        Y2=DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))       !funcion en el ultimo punto
!
        FINTGAUSS=REAL((Y1+Y2+4.D0*SUM1+2.D0*SUM2)*H/3.D0)            !solucion
!
        END
!
!******************************************************************************
! Integral de una gaussiana**2 de anchura SIGMA por el metodo de Simpson, entre 
! X1 y X2, con N intervalos (N debe ser par). Utilizamos la formula (1prim) del
! libro "Calculo Numerico Fundamental", Demidovich y Maron, pag. 656 (ojo, hay
! una errata en la formula de sigma2).
        REAL FUNCTION FINTGAUSSE(X1,X2,N,X0,FACTOR1)
        IMPLICIT NONE
        REAL X1,X2
        INTEGER N
        REAL X0,FACTOR1
!
        DOUBLE PRECISION DX1,DX2
        DOUBLE PRECISION DX0,DFACTOR1
!
        INTEGER I
        REAL X,Y1,Y2
        REAL SUM1,SUM2
!------------------------------------------------------------------------------
! chequeamos que N es par
        IF(MOD(N,2).NE.0)THEN
          STOP 'FATAL ERROR: N is odd in subroutine FINTGAUSSE!'
        END IF
! pasamos a doble precision las variables REAL de entrada
        DX1=DBLE(X1)
        DX2=DBLE(X2)
        DX0=DBLE(X0)
        DFACTOR1=DBLE(FACTOR1)
!
        SUM1=0.D0                                     !sumamos terminos impares
        DO I=1,N-1,2
          X=DX1+(DX2-DX1)*DBLE(I)/DBLE(N)
          SUM1=SUM1+&
     &     DEXP(DFACTOR1*(X-DX0)*(X-DX0))*DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
!
        SUM2=0.D0          !sumamos terminos pares salvo el primero y el ultimo
        DO I=2,N-2,2
          X=DX1+(DX2-DX1)*DBLE(I)/DBLE(N)
          SUM2=SUM2+&
     &     DEXP(DFACTOR1*(X-DX0)*(X-DX0))*DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
!
        Y1=DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))*&
     &   DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))
        Y2=DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))*&
     &   DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))
!
        FINTGAUSSE=&
     &   DBLE((Y1+Y2+4.D0*SUM1+2.D0*SUM2)*(DX2-DX1)/DBLE(3*N))        !solucion
!
        END
