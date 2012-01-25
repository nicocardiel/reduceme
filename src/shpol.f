C------------------------------------------------------------------------------
C Version 8-December-1996                                         file: shpol.f
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
C Program: shpol
C Classification: wavelengths
C Description: Determines the new coefficients of a polynomial after changing 
C the x-axis origin.
C
Comment
C
C Calcula los coeficientes de un polinomio tras realizar un cambio de origen
C
        PROGRAM SHPOL
C
        IMPLICIT NONE
        INCLUDE 'futils.inc'
        INCLUDE 'iofile.inc'
        INCLUDE 'redlib.inc'
        REAL READF
C
        INTEGER I,J,K,KK
        INTEGER NDEG
        REAL A(0:19),B(0:19)
        REAL COMB,XOFF
        CHARACTER*1 CSAVE
        CHARACTER*75 INFILE,OUTFILE
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
C
        THISPROGRAM='shpol'
        CALL WELCOME('8-December-1996')
C------------------------------------------------------------------------------
        WRITE(*,100)'File name (file with polynomial coefficients)'
        INFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
C
        WRITE(*,101)'----------------------------'
        K=-1
10      CONTINUE
        READ(20,*,END=20) KK,A(K+1)
        K=K+1
        WRITE(*,*)K,A(K)
        GOTO 10
20      CONTINUE
        CLOSE(20)
        WRITE(*,101)'----------------------------'
C
        NDEG=K
C
30      WRITE(*,100)'X offset (999=EXIT) '
        XOFF=READF('999')
        IF(XOFF.EQ.999.) STOP
        DO I=0,NDEG
          B(I)=0.
        END DO
C
        DO I=0,NDEG
          DO J=0,I
            B(I-J)=B(I-J)+A(I)*COMB(I,J)*(XOFF**REAL(J))
          END DO
        END DO
C
        WRITE(*,100)'Generate output file name (y/n) '
        CSAVE(1:1)=READC('n','yn')
C
        WRITE(*,100)'Output file name'
        IF(CSAVE.EQ.'y') OUTFILE=OUTFILEX(30,'@',0,0,.0,.0,3,.FALSE.)
        WRITE(*,101)'----------------------------'
        DO K=0,NDEG
          WRITE(*,*)K,B(K)
          IF(CSAVE.EQ.'y') WRITE(30,*)K,B(K)
        END DO
        IF(CSAVE.EQ.'y') CLOSE(30)
        WRITE(*,101)'----------------------------'
        GOTO 30
100     FORMAT(A,$)
101     FORMAT(A)
        END
C
C******************************************************************************
C
        REAL FUNCTION COMB(N,M)
        IMPLICIT NONE
        INTEGER N,M
        REAL FACT
C
        COMB=FACT(N)/(FACT(M)*FACT(N-M))
        END
C
C******************************************************************************
C
        REAL FUNCTION FACT(K)
        IMPLICIT NONE
        INTEGER K
        INTEGER J
C
        FACT=1.
        IF(K.LT.0)THEN
          STOP 'Factorial of negative number.'
        ELSEIF(K.GT.0)THEN
          DO J=1,K
            FACT=FACT*REAL(J)
          END DO
        END IF
        END
