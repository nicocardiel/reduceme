C------------------------------------------------------------------------------
C Version 19-October-2007                                     File: fintgauss.f
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
C REAL FUNCTION FINTGAUSS(X1,X2,N,X0,FACTOR1)  -> integral of Gaussian
C REAL FUNCTION FINTGAUSSE(X1,X2,N,X0,FACTOR1) -> integral of Gaussian^2
C
C Input: X1,X2,N,X0,FACTOR1
C Output: FINTGAUSS (function) or FINTGAUSSE (function)
C
C Calculate the integral of a Gaussian of width SIGMA, center at wavelength X0,
C using Simpson's rule, between wavelengths X1 and X2, with N intervals 
C (N must be even).
C
C REAL X1,X2 -> integral limits (wavelengths)
C INTEGER N -> number of intervals to estimate integral
C REAL X0 -> Gaussian center (wavelength)
C REAL FACTOR1 -> = (c^2/(2*SIGMA^2))/(X0^2)
C
Comment
C------------------------------------------------------------------------------
C Integral de una gaussiana de anchura SIGMA por el metodo de Simpson, entre 
C X1 y X2, con N intervalos (N debe ser par). Utilizamos la formula (1prim) del
C libro "Calculo Numerico Fundamental", Demidovich y Maron, pag. 656 (ojo, hay
C una errata en la formula de sigma2).
        REAL FUNCTION FINTGAUSS(X1,X2,N,X0,FACTOR1)
        IMPLICIT NONE
        REAL X1,X2
        INTEGER N
        REAL X0,FACTOR1
C
        DOUBLE PRECISION DX1,DX2
        DOUBLE PRECISION DX0,DFACTOR1
C
        INTEGER I
        DOUBLE PRECISION  X,Y1,Y2
        DOUBLE PRECISION H
        DOUBLE PRECISION SUM1,SUM2
C------------------------------------------------------------------------------
C chequeamos que N es par
        IF(MOD(N,2).NE.0)THEN
          STOP 'FATAL ERROR: N is odd in subroutine FINTGAUSS!'
        END IF
C pasamos a doble precision las variables REAL de entrada
        DX1=DBLE(X1)
        DX2=DBLE(X2)
        DX0=DBLE(X0)
        DFACTOR1=DBLE(FACTOR1)
C
        H=(DX2-DX1)/DBLE(N) !tamaño de cada intervalo
C
        SUM1=0.D0
        DO I=1,N-1,2
          X=DX1+H*DBLE(I)
          SUM1=SUM1+DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
C
        SUM2=0.D0          !sumamos terminos pares salvo el primero y el ultimo
        DO I=2,N-2,2
          X=DX1+H*DBLE(I)
          SUM2=SUM2+DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
C
        Y1=DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))       !funcion en el primer punto
        Y2=DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))       !funcion en el ultimo punto
C
        FINTGAUSS=REAL((Y1+Y2+4.D0*SUM1+2.D0*SUM2)*H/3.D0)            !solucion
C
        END
C
C******************************************************************************
C Integral de una gaussiana**2 de anchura SIGMA por el metodo de Simpson, entre 
C X1 y X2, con N intervalos (N debe ser par). Utilizamos la formula (1prim) del
C libro "Calculo Numerico Fundamental", Demidovich y Maron, pag. 656 (ojo, hay
C una errata en la formula de sigma2).
        REAL FUNCTION FINTGAUSSE(X1,X2,N,X0,FACTOR1)
        IMPLICIT NONE
        REAL X1,X2
        INTEGER N
        REAL X0,FACTOR1
C
        DOUBLE PRECISION DX1,DX2
        DOUBLE PRECISION DX0,DFACTOR1
C
        INTEGER I
        REAL X,Y1,Y2
        REAL SUM1,SUM2
C------------------------------------------------------------------------------
C chequeamos que N es par
        IF(MOD(N,2).NE.0)THEN
          STOP 'FATAL ERROR: N is odd in subroutine FINTGAUSSE!'
        END IF
C pasamos a doble precision las variables REAL de entrada
        DX1=DBLE(X1)
        DX2=DBLE(X2)
        DX0=DBLE(X0)
        DFACTOR1=DBLE(FACTOR1)
C
        SUM1=0.D0                                     !sumamos terminos impares
        DO I=1,N-1,2
          X=DX1+(DX2-DX1)*DBLE(I)/DBLE(N)
          SUM1=SUM1+
     +     DEXP(DFACTOR1*(X-DX0)*(X-DX0))*DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
C
        SUM2=0.D0          !sumamos terminos pares salvo el primero y el ultimo
        DO I=2,N-2,2
          X=DX1+(DX2-DX1)*DBLE(I)/DBLE(N)
          SUM2=SUM2+
     +     DEXP(DFACTOR1*(X-DX0)*(X-DX0))*DEXP(DFACTOR1*(X-DX0)*(X-DX0))
        END DO
C
        Y1=DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))*
     +   DEXP(DFACTOR1*(DX1-DX0)*(DX1-DX0))
        Y2=DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))*
     +   DEXP(DFACTOR1*(DX2-DX0)*(DX2-DX0))
C
        FINTGAUSSE=
     +   DBLE((Y1+Y2+4.D0*SUM1+2.D0*SUM2)*(DX2-DX1)/DBLE(3*N))        !solucion
C
        END
