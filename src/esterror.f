C------------------------------------------------------------------------------
C Version 29-September-2005                                    file: esterror.f
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
C Program: esterror
C Classification: error handling
C Description: Creates an error file from an initial image, taking into
C account the r.m.s. measured in each spectrum in a given wavelength range.
C
Comment
C 
        PROGRAM ESTERROR
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NREGMAX
        PARAMETER (NREGMAX=NCMAX)
        REAL C                                      !velocidad de la luz (km/s)
        PARAMETER (C=2.9979246E+5)
C
        REAL FMEAN0,FMEAN1
C
        INTEGER I,J,K,L,JJ
        INTEGER NREG,NPIXELS
        INTEGER NPOLDEG_REGION(NREGMAX)
        INTEGER NC1(NREGMAX),NC2(NREGMAX)
        INTEGER NFIT
        REAL A(NCMAX,NSMAX),B(NCMAX,NSMAX),ERR(NCMAX,NSMAX)
        REAL RADVEL,RCVEL,RCVEL1
        REAL LDO1(NREGMAX),LDO2(NREGMAX)
        REAL XF(NCMAX),YF(NCMAX),YFF(NCMAX)
        REAL CFIT(10),CHISQR
        REAL PIXEL(NCMAX)
        REAL FMEAN,FSIGMA
        REAL SNRATIO
        CHARACTER*75 INFILE,OUTFILE
        LOGICAL IFCHAN(NCMAX)
C------------------------------------------------------------------------------
        THISPROGRAM='esterror'
        CALL WELCOME('29-September-2005')
C leemos la imagen inicial
        WRITE(*,100)'Input file name'
        INFILE=INFILEX(20,'@',NSCAN,NCHAN,STWV,DISP,1,.FALSE.)
        DO I=1,NSCAN
          READ(20) (A(J,I),J=1,NCHAN)
        END DO
        CLOSE(20)
C------------------------------------------------------------------------------
C definimos region en la que vamos a medir el r.m.s. en los espectros
        WRITE(*,100)'Radial velocity (km/sec) '
        RADVEL=READF('0.0')
        RCVEL=RADVEL/C
        RCVEL1=1.0+RCVEL
        RCVEL1=RCVEL1/SQRT(1.-RCVEL*RCVEL)
C
        WRITE(*,100)'Number of regions to measure r.m.s.'
        NREG=READILIM('@',1,NCMAX)
C
        DO J=1,NCHAN
          IFCHAN(J)=.FALSE. !no pixel has yet been defined to be used
        END DO
C
        DO K=1,NREG
          WRITE(*,100) '>>> Region #'
          WRITE(*,*) K
          WRITE(*,100) 'Wavelength 1'
          LDO1(K)=READF('@')
          WRITE(*,100) 'Wavelength 2'
          LDO2(K)=READF('@')
          IF(LDO1(K).GT.LDO2(K)) STOP 'WV1.GT.WV2'
          NC1(K)=1+NINT((LDO1(K)*RCVEL1-STWV)/DISP)
          NC2(K)=1+NINT((LDO2(K)*RCVEL1-STWV)/DISP)
          WRITE(*,100) 'NC1,NC2: '
          WRITE(*,*) NC1(K),NC2(K)
          IF(NC1(K).LT.1) STOP 'NC1.LT.1'
          IF(NC2(K).GT.NCHAN) STOP 'NC2.GT.NCHAN'
          DO J=NC1(K),NC2(K)
            IF(IFCHAN(J)) STOP 'Channel already used'
            IFCHAN(J)=.TRUE.
          END DO
          WRITE(*,100) 'Polynomial degree to be applied '
          NPOLDEG_REGION(K)=READILIM('0',0,9)
          IF(NPOLDEG_REGION(K)+1.GT.NC2(K)-NC1(K)+1) 
     +     STOP 'Not enough points for fit'
          WRITE(*,*)
        END DO
C------------------------------------------------------------------------------
C ajustamos regiones con polinomios y sustraemos dichos polinomios
        DO I=1,NSCAN
          DO K=1,NREG
            NFIT=NC2(K)-NC1(K)+1
            DO J=NC1(K),NC2(K)
              JJ=J-NC1(K)+1
              XF(JJ)=REAL(J)
              YF(JJ)=A(J,I)
            END DO
            CALL POLFIT(XF,YF,YF,NFIT,NPOLDEG_REGION(K)+1,0,CFIT,CHISQR)
            DO J=NC1(K),NC2(K)
              JJ=J-NC1(K)+1
              YFF(JJ)=CFIT(NPOLDEG_REGION(K)+1)
              DO L=NPOLDEG_REGION(K),1,-1
                YFF(JJ)=YFF(JJ)*XF(JJ)+CFIT(L)
              END DO
              B(J,I)=A(J,I)-YFF(JJ)
            END DO
          END DO
        END DO
C medimos r.m.s. en las regiones definidas
        DO I=1,NSCAN
          NPIXELS=0
          DO K=1,NREG
            DO J=NC1(K),NC2(K)
              NPIXELS=NPIXELS+1
              PIXEL(NPIXELS)=B(J,I)
            END DO
          END DO
          FMEAN=FMEAN0(NPIXELS,PIXEL,FSIGMA) !FSIGMA es el r.m.s.
          NPIXELS=0
          DO K=1,NREG
            DO J=NC1(K),NC2(K)
              NPIXELS=NPIXELS+1
              PIXEL(NPIXELS)=A(J,I)
            END DO
          END DO
          FMEAN=FMEAN1(NPIXELS,PIXEL) !FMEAN es la señal promedio
          SNRATIO=FMEAN/FSIGMA
          DO J=1,NCHAN
            IF(SNRATIO.GT.0.0)THEN
              ERR(J,I)=A(J,I)/SNRATIO
            ELSE
              ERR(J,I)=A(J,I)/0.000001
            END IF
          END DO
        END DO
C------------------------------------------------------------------------------
C salvamos el fichero de errores
        WRITE(*,100)'Output file name '
        CALL GUESSEF(INFILE,OUTFILE)
        OUTFILE=OUTFILEX(30,OUTFILE,NSCAN,NCHAN,STWV,DISP,1,.TRUE.)
        DO I=1,NSCAN
          WRITE(30) (ERR(J,I),J=1,NCHAN)
        END DO
        CLOSE(30)
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
        END
