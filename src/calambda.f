C------------------------------------------------------------------------------
C Version 28-November-1996                                     file: calambda.f
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
C Program: calambda
C Classification: wavelengths
C Description: determine the wavelength as a function of the channel number,
C using the wavelength calibration polynomial.
C
Comment
C
C Calcula la longitud de onda de un determinado pixel a partir del polinomio
C de calibracion en l.d.o.
C
        PROGRAM CALAMBDA
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READI
        INTEGER READILIM
C
        INTEGER K
        INTEGER NTERMS
        INTEGER NP
        REAL POL1,POL2
        REAL A(20)
        CHARACTER*75 POLFILE
C------------------------------------------------------------------------------
        CALL AVOID_WARNINGS(STWV,DISP,NSCAN,NCHAN)
        OUTFILEX=OUTFILEX
C
        THISPROGRAM='calambda'
        CALL WELCOME('28-November-1996')
        CALL SHOWHLP('explanation')
C
        WRITE(*,100)'Polynomial file name'
        POLFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
C
        NTERMS=0
10      READ(20,*,END=20)K,A(NTERMS+1)
        NTERMS=NTERMS+1
        GOTO 10
20      CLOSE(20)
C
        WRITE(*,100)'NCHAN'
        NCHAN=READI('@')
C
        WRITE(*,100)'Pol. degree: '
        WRITE(*,*) NTERMS-1
C
        NP=1
        POL1=A(NTERMS)
        DO K=NTERMS-1,1,-1
          POL1=POL1*REAL(NP)+A(K)
        END DO
        WRITE(*,'(A,I5,A,$)')'Pixel: ',NP,'  Lambda: '
        WRITE(*,*) POL1
        NP=NCHAN
        POL2=A(NTERMS)
        DO K=NTERMS-1,1,-1
          POL2=POL2*REAL(NP)+A(K)
        END DO
        WRITE(*,'(A,I5,A,$)')'Pixel: ',NP,'  Lambda: '
        WRITE(*,*) POL2
        WRITE(*,100)'Dispersion: '
        WRITE(*,*)(POL2-POL1)/REAL(NCHAN-1)
C
30      WRITE(*,100)'Pixel (0=EXIT) '
        NP=READILIM('0',0,NCHAN)
        IF(NP.EQ.0)STOP
        POL1=A(NTERMS)
        DO K=NTERMS-1,1,-1
          POL1=POL1*REAL(NP)+A(K)
        END DO
        WRITE(*,'(A,I5,A,$)')'Pixel: ',NP,'  Lambda: '
        WRITE(*,*) POL1
        GOTO 30
C
        STOP
100     FORMAT(A,$)
        END
