C------------------------------------------------------------------------------
C Version 23-May-2006                                          file: predchan.f
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
C Program: predchan
C Classification: wavelengths
C Description: determine the final position (channel and wavelength) of a pixel
C in the wavelength axis corresponding to a given channel position (known 
C before the C-distortion correction and the wavelength calibration).
C
Comment
C
        PROGRAM PREDCHAN
        IMPLICIT NONE
        INCLUDE 'redlib.inc'
        INCLUDE 'iofile.inc'
        INTEGER READILIM
        REAL READF
C
        INTEGER NCDISTMAX
        PARAMETER (NCDISTMAX=9) !maximum number of C-distortion iterations
        INTEGER NDEGMAX
        PARAMETER (NDEGMAX=19) !maximum polynomial degree
C
        INTEGER I,J,K,L
        INTEGER I0
        INTEGER IBEG,IEND
        INTEGER NCFILES,NDEGC,NDEGW
        INTEGER IDUM
        REAL COEFF(0:NDEGMAX)
        REAL COEFF_C(0:NDEGMAX,NSMAX)
        REAL COEFF_W(0:NDEGMAX)
        REAL FDUM
        REAL XJ0,WPOL,XJ0_CCOR,XJ0_NEW,DXJ0
!       REAL LTWV !longitud de onda del centro del ultimo pixel
        REAL LDO_INI,LDO_FIN
        REAL XJ_INI,XJ_FIN
        CHARACTER*50 CDUMMY
        CHARACTER*255 CFILE
        LOGICAL LOOP
C------------------------------------------------------------------------------
        OUTFILEX=OUTFILEX
        NDEGC=0      !avoid compilation warning
        XJ0_CCOR=0.0 !avoid compilation warning
C
        THISPROGRAM='predchan'
        CALL WELCOME('23-May-2006')
C
        WRITE(*,100) 'Original NSCAN '
        NSCAN=READILIM('@',1,NSMAX)
        WRITE(*,100) 'Original NCHAN '
        NCHAN=READILIM('@',1,NCMAX)
        WRITE(*,100) 'STWV '
        STWV=READF('@')
        WRITE(*,100) 'DISP '
        DISP=READF('@')
C------------------------------------------------------------------------------
        WRITE(*,100) 'Number of iterations in C-distortion correction'//
     +   ' (0=none) ' 
        NCFILES=READILIM('0',0,NCDISTMAX)
C
        IF(NCFILES.GT.0)THEN
          WRITE(*,100)'Spectra for which the corr. was derived: '//
     +       'from... to... '
          WRITE(CDUMMY,'(I1,A1,I8)') 1,',',NSCAN
          CALL RMBLANK(CDUMMY,CDUMMY,L)
          CALL READ2I(CDUMMY(1:L),IBEG,IEND)
          WRITE(*,100) 'Order of polynomial '
          NDEGC=READILIM('@',0,NDEGMAX)
          DO I=IBEG,IEND
            DO K=0,NDEGC
              COEFF_C(K,I)=0.0
            END DO
          END DO
          DO J=1,NCFILES
            WRITE(*,'(A,I1,$)') 'C-distortion file name #',J
            CFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
            DO I=IBEG,IEND
              READ(20,*) FDUM,(COEFF(K),K=0,NDEGC)
              IF(FDUM.NE.REAL(I)) STOP 'FATAL ERROR: FDUM.NE.scan_no.'
              DO K=0,NDEGC
                COEFF_C(K,I)=COEFF_C(K,I)+COEFF(K)
              END DO
            END DO
            CLOSE(20)
          END DO
        END IF
C------------------------------------------------------------------------------
        WRITE(*,100) 'Wavelength calibration polynomial file name '
        CFILE=INFILEX(20,'@',0,0,0.,0.,3,.FALSE.)
        K=0
10      READ(20,*,END=12) IDUM,COEFF_W(K)
        K=K+1
        GOTO 10
12      CLOSE(20)
        NDEGW=K-1
        WRITE(*,100) '>>> Polynomial degree (wavelength calibration): '
        WRITE(*,*) NDEGW
C------------------------------------------------------------------------------
        LDO_INI=COEFF_W(NDEGW)
        LDO_FIN=COEFF_W(NDEGW)
        DO K=NDEGW-1,0,-1
          LDO_INI=LDO_INI*0.5+COEFF_W(K)
          LDO_FIN=LDO_FIN*(REAL(NCHAN)+0.5)+COEFF_W(K)
        END DO
        XJ_INI=(LDO_INI-STWV)/DISP+1.0
        XJ_FIN=(LDO_FIN-STWV)/DISP+1.0
        WRITE(*,100) '>>> First pixel with data:'
        IF(XJ_INI.GT.0.5)THEN
          WRITE(*,*) NINT(XJ_INI)
        ELSE
          WRITE(*,*) 1
        END IF
        WRITE(*,100) '>>> Last  pixel with data:'
        IF(XJ_FIN.LT.REAL(NCHAN)-0.5)THEN
          WRITE(*,*) NINT(XJ_FIN)
        ELSE
          WRITE(*,*) NCHAN
        END IF
C------------------------------------------------------------------------------
        I0=0
        XJ0=0.0
        LOOP=.TRUE.
        DO WHILE(LOOP)
          WRITE(CDUMMY,*) I0
          WRITE(*,100) 'Scan (INTEGER; 0=EXIT) '
          I0=READILIM(CDUMMY,0,NSCAN)
          IF(I0.GT.0)THEN
            WRITE(CDUMMY,*) XJ0
            WRITE(*,100) 'Channel (REAL) '
            XJ0=READF(CDUMMY)
            IF(NCFILES.GT.0)THEN
              DXJ0=COEFF_C(NDEGC,I0)
              DO K=NDEGC-1,0,-1
                DXJ0=DXJ0*XJ0+COEFF_C(K,I0)
              END DO
              XJ0_CCOR=XJ0-DXJ0
            END IF
            WPOL=COEFF_W(NDEGW)
            DO K=NDEGW-1,0,-1
              WPOL=WPOL*XJ0_CCOR+COEFF_W(K)
            END DO
            XJ0_NEW=(WPOL-STWV)/DISP+1.0
            WRITE(*,100) '>>> New wavelength.......: '
            WRITE(*,*) WPOL
            WRITE(*,100) '>>> New channel (real)...: '
            WRITE(*,*) XJ0_NEW
            WRITE(*,100) '>>> New channel (integer): '
            WRITE(*,*) NINT(XJ0_NEW)
          ELSE
            LOOP=.FALSE.
          END IF
        END DO
C------------------------------------------------------------------------------
        STOP
100     FORMAT(A,$)
        END
